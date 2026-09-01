#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""External-validation concordance: AMRFinderPlus + ResFinder vs phenotype (M13).

Head-to-head genotype-vs-phenotype validation of the KB antibiotics. Two
established genotypic AMR callers are run on the genome assemblies (on the HPC,
in ``amr-tools.sif``); this script (a) ``prep`` writes the genome list + paths
the SLURM job needs, and (b) ``post`` parses every per-genome tool output into a
per-antibiotic resistant/susceptible call, then scores each caller against the
EUCAST/CLSI phenotype (ground truth) with the clinical metrics in
``lib/concordance`` — balanced accuracy, sensitivity, specificity, Cohen's κ and
the FDA major/very-major error bands — and cross-compares the two callers
(κ + McNemar).

Two-container flow (like step 14):
    amr.sif       16_external_concordance.py --mode prep   # genome list + paths.sh
    amr-tools.sif amrfinder … ; python -m resfinder …      # per genome (SLURM loop)
    amr.sif       16_external_concordance.py --mode post    # parse + metrics

Genotypic mapping
-----------------
* **ResFinder** already emits a per-antibiotic phenotype prediction
  (``pheno_table_<species>.txt``) via its curated database — parsed directly.
* **AMRFinderPlus** emits a determinant table with ``Class``/``Subclass`` drug
  labels; a genome is genotypic-R for an antibiotic if any AMR-type row's
  Class/Subclass matches that antibiotic's keyword set (``AFP_KEYWORDS``). This
  uses NCBI's own curation faithfully — how well genotype predicts phenotype
  then falls out of the concordance metrics (e.g. AMPICILLIN via any
  ``BETA-LACTAM`` β-lactamase, but CEFOTAXIME only via ``CEPHALOSPORIN`` ESBL/
  AmpC, so a narrow TEM-1 correctly does NOT imply cefotaxime resistance).

Output (results/{organism}/external_validation/):
    16_concordance_{organism}.csv       — one row per (antibiotic, caller) vs phenotype
    16_concordance_summary_{organism}.json
"""

import argparse
import csv
import datetime
import json
import os
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import concordance as C  # noqa: E402
from lib.config import load_config, resolve_path, get_target  # noqa: E402
from lib.logging_utils import get_logger  # noqa: E402
from lib.registry import load_amrfinder_keywords  # noqa: E402

# Antibiotic -> AMRFinderPlus Class/Subclass keyword set (upper-case tokens).
# Registry-driven (config/registry/antibiotics.yaml `amrfinder_keywords`; audit
# Issue 9) so adding an antibiotic needs no code change; a built-in fallback keeps
# the script runnable if the registry section is absent. A narrow β-lactamase
# (Subclass BETA-LACTAM) implies ampicillin-R but NOT cefotaxime-R (that needs an
# ESBL/AmpC -> Subclass CEPHALOSPORIN).
_AFP_FALLBACK = {
    "ampicillin":    {"AMPICILLIN", "BETA-LACTAM"},
    "cefotaxime":    {"CEFOTAXIME", "CEPHALOSPORIN"},
    "ciprofloxacin": {"CIPROFLOXACIN", "FLUOROQUINOLONE", "QUINOLONE"},
}
AFP_KEYWORDS = load_amrfinder_keywords() or _AFP_FALLBACK
DEFAULT_ANTIBIOTICS = list(AFP_KEYWORDS)

# Tool versions for provenance (audit Issue 8). These are the M13 BASELINE — the
# whole headline is "our model beats AMRFinderPlus" (K. pneu ciprofloxacin: bACC
# 0.926 vs 0.538), which means nothing unless the KB says WHICH AMRFinderPlus.
#
# SOFTWARE and DATABASE versions are different facts and both matter: the same
# AMRFinderPlus binary calls different genes off a newer DB. The old default here
# was "AMRFinderPlus 2026-05-15.1" — that is the DATABASE version stamped as if
# it were the software (the software is 4.2.7), so the KB recorded a date where a
# version belongs.
#
# Asked of the tools themselves rather than hardcoded: a literal default rots
# silently, and this one had. Env overrides stay for the case where the binary
# is not on PATH (e.g. running from another container).
def _probe(cmd, args, env_var, label):
    """Ask the tool for its own version; fall back to the env var, then None.

    The exit code is checked, not just the output: `python -m resfinder --version`
    with resfinder absent exits 1 and prints "No module named resfinder" to
    stderr — reading that as a version would stamp an error message into the KB
    as provenance, which is worse than admitting it is unknown.
    """
    v = os.environ.get(env_var)
    if v:
        return v
    try:
        import subprocess
        r = subprocess.run([cmd] + args, capture_output=True, text=True,
                           check=False, timeout=30)
        if r.returncode != 0:
            return None
        out = (r.stdout or r.stderr or "").strip().splitlines()
        if out:
            return f"{label} {out[0].strip()}"
    except Exception:
        pass
    return None


def _amrfinder_db_version():
    """AMRFinderPlus DB version — a separate fact from the software version.
    The DB lives on scratch (not baked into the image) and is named by date."""
    v = os.environ.get("AMR_AFP_DB_VERSION")
    if v:
        return v
    db = PROJECT_ROOT / "data" / "external" / "amrfinder_db"
    try:
        dated = sorted(p.name for p in db.iterdir()
                       if p.is_dir() and p.name[:4].isdigit())
        return dated[-1] if dated else None
    except Exception:
        return None


AFP_VERSION = _probe("amrfinder", ["--version"], "AMR_AFP_VERSION", "AMRFinderPlus")
AFP_DB_VERSION = _amrfinder_db_version()
RF_VERSION = _probe("python", ["-m", "resfinder", "--version"], "AMR_RF_VERSION", "ResFinder")

# evidence_source strings: software + DB together, so a KB row states exactly what
# produced it (e.g. "AMRFinderPlus 4.2.7 (DB 2026-05-15.1)").
AFP_SOURCE = " ".join(filter(None, [
    AFP_VERSION or "AMRFinderPlus (version unknown)",
    f"(DB {AFP_DB_VERSION})" if AFP_DB_VERSION else None]))
RF_SOURCE = RF_VERSION or "ResFinder (version unknown)"


def _tokens(field):
    """Split an AMRFinderPlus Class/Subclass cell into upper-case tokens."""
    if not field or field.upper() in ("NA", ""):
        return set()
    return {t.strip().upper() for t in field.replace(",", "/").split("/") if t.strip()}


# Agents AMRFinderPlus cannot resolve, so it is not scored on them at all.
# Without this they come out as an all-negative predictor at balanced accuracy
# 0.500, which reads as "the tool performed at chance" when the truth is "the
# tool was not asked", and that is unfair to the tool in a head-to-head table.
AFP_NOT_ASSESSABLE = {
    # AMRFinderPlus labels every van determinant `GLYCOPEPTIDE | VANCOMYCIN` and
    # publishes no TEICOPLANIN subclass, so it cannot separate vanA (vancomycin
    # and teicoplanin) from vanB (vancomycin only).
    "teicoplanin",
    # AMRFinderPlus resolves beta-lactamases but not inhibitor combinations: it
    # publishes no SULBACTAM subclass, so it cannot say whether sulbactam
    # restores activity against the enzymes it did find.
    "ampicillin_sulbactam",
    # No MONOBACTAM or AZTREONAM subclass appears anywhere in the delivered
    # output, so aztreonam draws no call at all: 0% called against an 88.7%
    # phenotype rate in K. pneumoniae. Scoring that as an all-negative predictor
    # would report the tool as failing a question it was never able to answer.
    "aztreonam",
}


def parse_amrfinder(tsv_path, antibiotics=DEFAULT_ANTIBIOTICS, keywords=AFP_KEYWORDS):
    """One AMRFinderPlus TSV -> {antibiotic: 0/1}. R if any AMR-type row's
    Class/Subclass tokens intersect the antibiotic's keyword set. Agents in
    AFP_NOT_ASSESSABLE are omitted, so downstream code reports no call rather
    than a spurious susceptible one."""
    calls = {ab: 0 for ab in antibiotics if ab not in AFP_NOT_ASSESSABLE}
    with open(tsv_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("Type") or row.get("Element type") or "").strip().upper() != "AMR":
                continue
            # `--plus` adds elements NCBI flags as of interest -- stress response,
            # efflux, virulence -- which are not by themselves resistance
            # determinants. Counting them as calls inflates the tool's false
            # positives: mepA, the chromosomal MATE pump carried by 2501 of 2505
            # S. aureus genomes here, is a `plus` element, and taking it at face
            # value called tetracycline in 100% of them against a 23.2%
            # phenotype rate. Only the curated `core` set is scored.
            if (row.get("Scope") or "core").strip().lower() != "core":
                continue
            toks = _tokens(row.get("Class")) | _tokens(row.get("Subclass"))
            for ab in calls:
                if toks & keywords.get(ab, set()):
                    calls[ab] = 1
    return calls


# ResFinder writes antibiotic names as free text with '+' and spaces
# ("amoxicillin+clavulanic acid"), while this project keys them with underscores
# ("amoxicillin_clavulanic_acid"). A literal lower-case comparison therefore
# matched neither, and every combination agent silently produced no ResFinder
# call at all. Matching is done on the token SET so separator style stops
# mattering, and never on a substring: "ampicillin_sulbactam" must NOT match
# ResFinder's plain "ampicillin", because the inhibitor changes the phenotype.
def _ab_tokens(name):
    return frozenset(t for t in re.split(r"[^a-z0-9]+", str(name).strip().lower()) if t)


# Combinations ResFinder reports one component at a time. It publishes no row
# for the combination itself, so the call has to be assembled from the parts.
# Rule: resistant if EITHER component is called resistant — a genotypic
# convention, since a single acquired sul or dfr determinant is enough to
# abolish the synergy the combination depends on. This is a decision taken here,
# not something ResFinder reports, and it is stated in the thesis as such.
RF_COMPONENTS = {
    "trimethoprim_sulfamethoxazole": ("trimethoprim", "sulfamethoxazole"),
}


def parse_resfinder(pheno_table_path, antibiotics=DEFAULT_ANTIBIOTICS):
    """One ResFinder pheno_table -> {antibiotic: 0/1}. Reads the '# Antimicrobial
    <TAB> Class <TAB> WGS-predicted phenotype …' rows directly.

    Antibiotics ResFinder does not report at all (e.g. ampicillin_sulbactam,
    oxacillin) are simply absent from the returned dict, so downstream code
    scores them as 'no ResFinder call' rather than as a susceptible call.
    """
    wanted = {_ab_tokens(ab): ab for ab in antibiotics}
    comp_of = {}
    for ab, parts in RF_COMPONENTS.items():
        if ab in antibiotics:
            for part in parts:
                comp_of.setdefault(_ab_tokens(part), []).append(ab)

    calls, component_hits = {}, {}
    with open(pheno_table_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            toks = _ab_tokens(parts[0])
            is_r = 1 if parts[2].strip().lower().startswith("resistant") else 0
            if toks in wanted:
                calls[wanted[toks]] = is_r
            for ab in comp_of.get(toks, ()):
                component_hits.setdefault(ab, []).append(is_r)

    # Assemble component-built combinations, but never overwrite a direct row.
    for ab, hits in component_hits.items():
        if ab not in calls:
            calls[ab] = 1 if any(hits) else 0
    return calls


def load_phenotype(metadata_file, antibiotics):
    """{genome_id: {antibiotic: 0/1/None}} from amr_phenotypes.csv (blank -> None)."""
    import pandas as pd
    # Genome ID as str at read time (audit Issue 5): "562.10" parsed as float ->
    # 562.1 would mismatch the {gid}.fna filename and the model_preds join.
    gid_col = pd.read_csv(metadata_file, nrows=0).columns[0]
    df = pd.read_csv(metadata_file, encoding="utf-8", dtype={gid_col: str})
    # BV-BRC writes combination agents with a slash ("trimethoprim/
    # sulfamethoxazole"); this project keys them with underscores. A literal
    # r.get(ab) therefore returned None for every combination, so those models
    # were scored against no phenotype at all and reported n=0. Same token-set
    # match as parse_resfinder, for the same reason.
    col_of = {_ab_tokens(c): c for c in df.columns}
    cols = {ab: col_of.get(_ab_tokens(ab)) for ab in antibiotics}
    missing = sorted(ab for ab, c in cols.items() if c is None)
    if missing:
        print(f"  note: no phenotype column for {', '.join(missing)}")
    pheno = {}
    for _, r in df.iterrows():
        gid = str(r[gid_col])
        row = {}
        for ab in antibiotics:
            col = cols.get(ab)
            v = r.get(col) if col else None
            row[ab] = None if (v is None or (isinstance(v, float) and v != v)) else int(v)
        pheno[gid] = row
    return pheno


def do_prep(organism, antibiotics, config, out_dir, logger):
    metadata_file = resolve_path("metadata_file", organism=organism, config=config)
    raw_genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    pheno = load_phenotype(metadata_file, antibiotics)
    # genomes with a label for ANY target antibiotic AND a present assembly
    genomes = [g for g, row in pheno.items()
               if any(row[ab] is not None for ab in antibiotics)
               and (raw_genomes_dir / f"{g}.fna").exists()]
    genomes.sort()
    (out_dir / "16_genomes.txt").write_text("\n".join(genomes) + "\n", encoding="utf-8")
    paths_sh = out_dir / "16_paths.sh"
    paths_sh.write_text("\n".join([
        f'GENOMES_DIR="{raw_genomes_dir}"',
        f'GENOME_LIST="{out_dir / "16_genomes.txt"}"',
        f'AFP_DIR="{out_dir / "amrfinder"}"',
        f'RF_DIR="{out_dir / "resfinder"}"',
        f'OUT_DIR="{out_dir}"']) + "\n", encoding="utf-8")
    logger.info("prep: %d genomes with a target-antibiotic label + assembly", len(genomes))
    logger.info("  ✓ %s", out_dir / "16_genomes.txt")
    logger.info("  ✓ %s (source paths for the SLURM job)", paths_sh)


def _resfinder_pheno_file(rf_genome_dir):
    """The species-specific pheno_table (has cefotaxime/ciprofloxacin); fall back
    to the generic one."""
    hits = sorted(rf_genome_dir.glob("pheno_table_*.txt"))
    if hits:
        return hits[0]
    generic = rf_genome_dir / "pheno_table.txt"
    return generic if generic.exists() else None


def head_to_head(genomes, pheno, afp, rf, model_calls, antibiotics):
    """3-way concordance on the model's held-out TEST genomes (leakage-free):
    model vs AMRFinderPlus vs ResFinder vs phenotype on the identical genome set,
    plus model-vs-tool κ/McNemar. Only antibiotics with model predictions (from
    06's saved test split) are included; the tools are re-scored on exactly those
    genomes so all three predictors share one sample."""
    out = {}
    for ab in antibiotics:
        mcall = model_calls.get(ab)
        if not mcall:
            continue
        common = [g for g in genomes if g in mcall and g in afp and g in rf
                  and ab in afp[g] and ab in rf[g]
                  and pheno.get(g, {}).get(ab) is not None]
        if not common:
            continue
        yt = [pheno[g][ab] for g in common]
        ym = [mcall[g] for g in common]
        ya = [afp[g][ab] for g in common]
        yr = [rf[g][ab] for g in common]
        out[ab] = {
            "n_common_test_genomes": len(common),
            "model": C.score_pair(yt, ym),
            "amrfinderplus": C.score_pair(yt, ya),
            "resfinder": C.score_pair(yt, yr),
            "model_vs_amrfinderplus": {"cohen_kappa": C.cohen_kappa(ym, ya),
                                       "mcnemar": C.mcnemar(ym, ya)},
            "model_vs_resfinder": {"cohen_kappa": C.cohen_kappa(ym, yr),
                                   "mcnemar": C.mcnemar(ym, yr)},
        }
    return out


def do_post(organism, antibiotics, out_dir, config, logger):
    metadata_file = resolve_path("metadata_file", organism=organism, config=config)
    pheno = load_phenotype(metadata_file, antibiotics)
    afp_dir, rf_dir = out_dir / "amrfinder", out_dir / "resfinder"

    genomes = [g.strip() for g in (out_dir / "16_genomes.txt").read_text().splitlines() if g.strip()] \
        if (out_dir / "16_genomes.txt").exists() else sorted(pheno)

    # gather per-genome caller calls (None if that genome's output is missing)
    afp, rf = {}, {}
    n_afp = n_rf = 0
    for g in genomes:
        t = afp_dir / f"afp_{g}.tsv"
        if t.exists():
            afp[g] = parse_amrfinder(t, antibiotics); n_afp += 1
        p = _resfinder_pheno_file(rf_dir / f"rf_{g}")
        if p:
            rf[g] = parse_resfinder(p, antibiotics); n_rf += 1
    logger.info("post: parsed AMRFinderPlus %d / ResFinder %d of %d genomes",
                n_afp, n_rf, len(genomes))

    rows, summary = [], {"organism": organism,
                         "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                         "n_genomes": len(genomes), "n_amrfinder": n_afp, "n_resfinder": n_rf,
                         "antibiotics": {}}
    for ab in antibiotics:
        # align on genomes with a phenotype for this antibiotic
        evald = [g for g in genomes if pheno.get(g, {}).get(ab) is not None]
        y_true = [pheno[g][ab] for g in evald]
        y_afp = [afp.get(g, {}).get(ab) for g in evald]
        y_rf = [rf.get(g, {}).get(ab) for g in evald]
        ab_doc = {"n_evaluable": len(evald), "n_resistant_phenotype": sum(y_true)}
        for caller, y in (("amrfinderplus", y_afp), ("resfinder", y_rf)):
            s = C.score_pair(y_true, y)
            ab_doc[caller] = s
            rows.append({"antibiotic": ab, "caller": caller, **{k: s[k] for k in
                         ("n", "sensitivity", "specificity", "balanced_accuracy",
                          "cohen_kappa", "major_error_rate", "very_major_error_rate")}})
        # caller-vs-caller agreement (paired, same genomes)
        ab_doc["amrfinder_vs_resfinder"] = {
            "cohen_kappa": C.cohen_kappa(y_afp, y_rf), "mcnemar": C.mcnemar(y_afp, y_rf)}
        summary["antibiotics"][ab] = ab_doc
        logger.info("  %s: n=%d  AFP bACC=%s κ=%s  RF bACC=%s κ=%s", ab, len(evald),
                    _r(ab_doc["amrfinderplus"]["balanced_accuracy"]),
                    _r(ab_doc["amrfinderplus"]["cohen_kappa"]),
                    _r(ab_doc["resfinder"]["balanced_accuracy"]),
                    _r(ab_doc["resfinder"]["cohen_kappa"]))

    # ---- model-vs-tool head-to-head on the model's held-out test genomes ----
    import pandas as pd
    model_calls = {}
    for ab in antibiotics:
        f = out_dir / f"16_model_preds_{ab}.csv"
        if f.exists():
            mp = pd.read_csv(f)
            model_calls[ab] = dict(zip(mp["Genome ID"].astype(str),
                                       mp["model_pred"].astype(int)))
    if model_calls:
        summary["head_to_head_model_test_genomes"] = head_to_head(
            genomes, pheno, afp, rf, model_calls, antibiotics)
        for ab, h in summary["head_to_head_model_test_genomes"].items():
            logger.info("  H2H %s (n=%d): model bACC=%s κ=%s | AFP bACC=%s | RF bACC=%s",
                        ab, h["n_common_test_genomes"],
                        _r(h["model"]["balanced_accuracy"]), _r(h["model"]["cohen_kappa"]),
                        _r(h["amrfinderplus"]["balanced_accuracy"]),
                        _r(h["resfinder"]["balanced_accuracy"]))
    else:
        logger.info("  (no 16_model_preds_*.csv yet — run 06 to enable model head-to-head)")

    csv_path = out_dir / f"16_concordance_{organism}.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["antibiotic", "caller", "n", "sensitivity",
                           "specificity", "balanced_accuracy", "cohen_kappa",
                           "major_error_rate", "very_major_error_rate"])
        w.writeheader()
        w.writerows(rows)
    summary_path = out_dir / f"16_concordance_summary_{organism}.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    logger.info("  ✓ %s", csv_path)
    logger.info("  ✓ %s", summary_path)
    return summary


def _r(x):
    return "NA" if x is None else f"{x:.3f}"


def write_kb_evidence(db_path, summary, logger, organism):
    """Persist the concordance result into amrk.db `validation_evidence` (M11):
    one row per (antibiotic, caller) vs phenotype + the model head-to-head, linked
    to that antibiotic's model run. Idempotent (clears prior concordance rows)."""
    import sqlite3
    conn = sqlite3.connect(str(db_path))
    try:
        # Key the run lookup by (organism, antibiotic), not by antibiotic alone:
        # the same agent is modelled in several organisms, so an antibiotic-only
        # map collapses six runs into one and attributes every organism's
        # concordance to whichever run happened to be last.
        runs, models = {}, {}
        for mid, ab, org in conn.execute(
                "SELECT m.model_id, m.antibiotic, r.organism FROM models m "
                "JOIN pipeline_runs r ON r.run_id = m.run_id"):
            models[(org, ab)] = mid
        for rid, ab, org in conn.execute(
                "SELECT m.run_id, m.antibiotic, r.organism FROM models m "
                "JOIN pipeline_runs r ON r.run_id = m.run_id"):
            runs[(org, ab)] = rid
        etypes = ("concordance_amrfinderplus", "concordance_resfinder", "head_to_head_model")
        # Clear only THIS organism's rows. Deleting every concordance row on each
        # call meant a six-organism sweep left only the last organism behind.
        own_runs = [r for (o, _), r in runs.items() if o == organism]
        if own_runs:
            conn.execute(
                f"DELETE FROM validation_evidence WHERE evidence_type IN "
                f"({','.join('?' * len(etypes))}) AND pipeline_run_id IN "
                f"({','.join('?' * len(own_runs))})", (*etypes, *own_runs))
        own_models = [m for (o, _), m in models.items() if o == organism]
        if own_models:
            conn.execute("DELETE FROM external_concordance WHERE model_id IN "
                         f"({','.join('?' * len(own_models))})", own_models)
        n = 0
        for ab, doc in summary["antibiotics"].items():
            rid = runs.get((organism, ab))
            mid = models.get((organism, ab))
            # An antibiotic with no model in this organism has nothing to attach
            # a concordance row to. Writing one anyway left a NULL run_id that no
            # organism-scoped delete could clear, so re-running the sweep grew the
            # table instead of replacing it.
            if rid is None:
                continue
            for caller, et, src in (
                    ("amrfinderplus", "concordance_amrfinderplus", AFP_SOURCE),
                    ("resfinder", "concordance_resfinder", RF_SOURCE)):
                s = doc[caller]
                conn.execute(
                    "INSERT INTO validation_evidence(unitig_id, evidence_type, "
                    "evidence_source, evidence_score, pipeline_run_id) VALUES (NULL,?,?,?,?)",
                    (et, f"{src} vs EUCAST/CLSI (bACC={_r(s['balanced_accuracy'])}, "
                     f"kappa={_r(s['cohen_kappa'])}, n={s['n']})", s["cohen_kappa"], rid))
                n += 1
                # The purpose-built table: one row per (model, caller) with the
                # full clinical metric set, which validation_evidence cannot hold.
                if mid is not None and s.get("n"):
                    conn.execute(
                        "INSERT INTO external_concordance(model_id, caller, reference, "
                        "n_test, sensitivity, specificity, balanced_accuracy, cohen_kappa, "
                        "major_error_rate, very_major_error_rate) VALUES (?,?,?,?,?,?,?,?,?,?)",
                        (mid, caller, "EUCAST/CLSI phenotype (BV-BRC)", s["n"],
                         s.get("sensitivity"), s.get("specificity"),
                         s.get("balanced_accuracy"), s.get("cohen_kappa"),
                         s.get("major_error_rate"), s.get("very_major_error_rate")))
        for ab, h in summary.get("head_to_head_model_test_genomes", {}).items():
            if (organism, ab) not in runs:
                continue
            m = h["model"]
            conn.execute(
                "INSERT INTO validation_evidence(unitig_id, evidence_type, "
                "evidence_source, evidence_score, pipeline_run_id) VALUES (NULL,?,?,?,?)",
                ("head_to_head_model",
                 f"unitig model vs AMRFinderPlus/ResFinder on held-out test "
                 f"(bACC={_r(m['balanced_accuracy'])}, kappa={_r(m['cohen_kappa'])}, "
                 f"n={h['n_common_test_genomes']})", m["cohen_kappa"], runs[(organism, ab)]))
            n += 1
        conn.commit()
        logger.info("  ✓ wrote %d concordance evidence rows to KB (%s)", n, db_path)
    finally:
        conn.close()


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="External-validation concordance (M13).")
    ap.add_argument("--mode", choices=["prep", "post"], required=True)
    ap.add_argument("--organism", default=get_target(config=config)[0])
    ap.add_argument("--antibiotics", default=",".join(DEFAULT_ANTIBIOTICS),
                    help="comma-separated (default: ampicillin,cefotaxime,ciprofloxacin)")
    ap.add_argument("--write-kb", action="store_true",
                    help="(post) also write concordance to amrk.db validation_evidence (M11)")
    ap.add_argument("--db", default=None, help="KB path (default: results/{org}/kb/amrk.db)")
    args = ap.parse_args()
    organism = args.organism
    antibiotics = [a.strip() for a in args.antibiotics.split(",") if a.strip()]
    out_dir = PROJECT_ROOT / "results" / organism / "external_validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = get_logger("m13-concordance")

    if args.mode == "prep":
        do_prep(organism, antibiotics, config, out_dir, logger)
    else:
        summary = do_post(organism, antibiotics, out_dir, config, logger)
        if args.write_kb:
            db_path = Path(args.db) if args.db else (
                PROJECT_ROOT / "results" / organism / "kb" / "amrk.db")
            if db_path.exists():
                write_kb_evidence(db_path, summary, logger, organism)
            else:
                logger.warning("--write-kb: KB not found at %s (skipped)", db_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
