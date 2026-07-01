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
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import concordance as C  # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402
from lib.logging_utils import get_logger  # noqa: E402

# Antibiotic -> AMRFinderPlus Class/Subclass keyword set (upper-case tokens).
# A narrow β-lactamase (Subclass BETA-LACTAM) implies ampicillin-R but NOT
# cefotaxime-R (which needs an ESBL/AmpC -> Subclass CEPHALOSPORIN).
AFP_KEYWORDS = {
    "ampicillin":    {"AMPICILLIN", "BETA-LACTAM"},
    "cefotaxime":    {"CEFOTAXIME", "CEPHALOSPORIN"},
    "ciprofloxacin": {"CIPROFLOXACIN", "FLUOROQUINOLONE", "QUINOLONE"},
}
DEFAULT_ANTIBIOTICS = ["ampicillin", "cefotaxime", "ciprofloxacin"]


def _tokens(field):
    """Split an AMRFinderPlus Class/Subclass cell into upper-case tokens."""
    if not field or field.upper() in ("NA", ""):
        return set()
    return {t.strip().upper() for t in field.replace(",", "/").split("/") if t.strip()}


def parse_amrfinder(tsv_path, antibiotics=DEFAULT_ANTIBIOTICS, keywords=AFP_KEYWORDS):
    """One AMRFinderPlus TSV -> {antibiotic: 0/1}. R if any AMR-type row's
    Class/Subclass tokens intersect the antibiotic's keyword set."""
    calls = {ab: 0 for ab in antibiotics}
    with open(tsv_path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if (row.get("Type") or row.get("Element type") or "").strip().upper() != "AMR":
                continue
            toks = _tokens(row.get("Class")) | _tokens(row.get("Subclass"))
            for ab in antibiotics:
                if toks & keywords.get(ab, set()):
                    calls[ab] = 1
    return calls


def parse_resfinder(pheno_table_path, antibiotics=DEFAULT_ANTIBIOTICS):
    """One ResFinder pheno_table -> {antibiotic: 0/1}. Reads the '# Antimicrobial
    <TAB> Class <TAB> WGS-predicted phenotype …' rows directly."""
    wanted = {ab.lower(): ab for ab in antibiotics}
    calls = {}
    with open(pheno_table_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip().lower()
            if name in wanted:
                calls[wanted[name]] = 1 if parts[2].strip().lower().startswith("resistant") else 0
    return calls


def load_phenotype(metadata_file, antibiotics):
    """{genome_id: {antibiotic: 0/1/None}} from amr_phenotypes.csv (blank -> None)."""
    import pandas as pd
    df = pd.read_csv(metadata_file, encoding="utf-8")
    gid_col = df.columns[0]
    pheno = {}
    for _, r in df.iterrows():
        gid = str(r[gid_col])
        row = {}
        for ab in antibiotics:
            v = r.get(ab)
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


def _r(x):
    return "NA" if x is None else f"{x:.3f}"


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="External-validation concordance (M13).")
    ap.add_argument("--mode", choices=["prep", "post"], required=True)
    ap.add_argument("--organism", default=config.get("project", {}).get("organism", "ecoli"))
    ap.add_argument("--antibiotics", default=",".join(DEFAULT_ANTIBIOTICS),
                    help="comma-separated (default: ampicillin,cefotaxime,ciprofloxacin)")
    args = ap.parse_args()
    organism = args.organism
    antibiotics = [a.strip() for a in args.antibiotics.split(",") if a.strip()]
    out_dir = PROJECT_ROOT / "results" / organism / "external_validation"
    out_dir.mkdir(parents=True, exist_ok=True)
    logger = get_logger("m13-concordance")

    if args.mode == "prep":
        do_prep(organism, antibiotics, config, out_dir, logger)
    else:
        do_post(organism, antibiotics, out_dir, config, logger)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
