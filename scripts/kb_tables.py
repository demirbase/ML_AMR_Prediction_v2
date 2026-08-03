#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Export tidy, machine-readable CSV summary tables from the unified AMR-KB
(+ a few per-run result JSONs the KB doesn't store: PFER, pyseer counts, H2).

Feeds kb_figures.py and any external plotting / thesis tables. Read-only.

Usage:
    python scripts/kb_tables.py --db results/kb/amrk.db \
        --results results --out results/tables

Outputs (results/tables/):
    models_summary.csv    one row/model: perf (lineage-CV, ROC, MCC…) + provenance
    kb_overview.csv       one row/model: n_stable, n_confirmed, PFER, pyseer-sig, H2…
    biomarkers.csv        one row/(model,unitig): gain/SHAP/stable + gene/tier/ARO +
                          MDA + pyseer-p + prevalence  (the full evidence table)
    mechanisms.csv        one row/(model, on-target confirmed gene) — class-filtered
"""
import argparse
import glob
import json
import os
import re
import sqlite3
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

# CARD aro_drug_class keyword(s) that count as ON-TARGET for a registry class.
CLASS_TO_ARO_KEYWORD = {
    "penicillins": ("penam", "penicillin"),
    "cephalosporins": ("cephalosporin", "cephamycin"),
    "beta_lactams_carbapenems_others": ("carbapenem", "monobactam"),
    "quinolones": ("fluoroquinolone", "quinolone"),
    "aminoglycosides": ("aminoglycoside",),
    "tetracyclines": ("tetracycline",),
    "folate_pathway_inhibitors": ("sulfonamide", "diaminopyrimidine"),
}


def _org(run_id):
    return run_id.split("__")[0]


def _write(rows, cols, path):
    import csv
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})
    print(f"  ✓ {path}  ({len(rows)} rows)")


def _load_json(pattern):
    hits = glob.glob(pattern)
    return json.load(open(hits[0])) if hits else {}


def main():
    ap = argparse.ArgumentParser(description="Export tidy CSV tables from the AMR-KB.")
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--results", default="results", help="root holding {org}/{ab}/05_explainability")
    ap.add_argument("--out", default="results/tables")
    args = ap.parse_args()

    c = sqlite3.connect(args.db)
    c.row_factory = sqlite3.Row
    out = Path(args.out)

    cls_of = dict(c.execute("SELECT antibiotic, drug_class FROM antibiotics").fetchall())
    models = c.execute("SELECT * FROM models ORDER BY model_id").fetchall()

    # ---- models_summary + kb_overview -------------------------------------
    msum, overview = [], []
    for m in models:
        mid, ab, rid = m["model_id"], m["antibiotic"], m["run_id"]
        pr = c.execute("SELECT * FROM pipeline_runs WHERE run_id=?", (rid,)).fetchone()
        org = pr["organism"] if pr else _org(rid)
        rdir = f"{args.results}/{org}/{ab}/05_explainability"
        st = _load_json(f"{rdir}/13_stability_summary_{ab}.json")
        ps = _load_json(f"{rdir}/14_pyseer_summary_{ab}.json")
        vm = _load_json(f"{rdir}/08_validation_metrics_{ab}.json")

        msum.append(dict(
            model_id=mid, organism=org, antibiotic=ab, drug_class=cls_of.get(ab),
            n_genomes=pr["n_genomes"] if pr else None,
            lineage_cv_auc=m["auc_mean_seeds"], lineage_cv_std=m["auc_std_seeds"],
            roc_auc_singlesplit=m["roc_auc"], pr_auc=m["pr_auc"], mcc=m["mcc"],
            balanced_accuracy=m["balanced_accuracy"], accuracy=m["accuracy"],
            n_trees=m["n_trees"], operating_threshold=m["operating_threshold"],
            git_commit=pr["git_commit"] if pr else None,
            random_seed=pr["random_seed"] if pr else None,
            card_version=pr["card_version"] if pr else None, run_id=rid,
        ))
        n_stable = c.execute("SELECT COUNT(*) FROM unitig_model_scores WHERE model_id=? AND stable=1", (mid,)).fetchone()[0]
        n_conf = c.execute("SELECT COUNT(*) FROM blast_annotations WHERE model_id=? AND tier='confirmed'", (mid,)).fetchone()[0]
        n_pysig = c.execute("SELECT COUNT(*) FROM validation_evidence WHERE evidence_type='pyseer_lmm' AND evidence_score<=0.05 AND pipeline_run_id=?", (rid,)).fetchone()[0]
        n_rsnp = c.execute("SELECT COUNT(*) FROM variant_snp_check WHERE model_id=? AND allele_class='resistant_allele'", (mid,)).fetchone()[0]
        overview.append(dict(
            model_id=mid, organism=org, antibiotic=ab, drug_class=cls_of.get(ab),
            lineage_cv_auc=m["auc_mean_seeds"], mcc=m["mcc"],
            cpss_n_stable=n_stable, pfer_bound=st.get("pfer_bound"),
            avg_selected_per_fit=st.get("avg_selected_per_fit"),
            n_confirmed_card=n_conf, resistant_allele_snp=n_rsnp,
            pyseer_significant=ps.get("n_significant"),
            cpss_stable_and_pyseer_sig=ps.get("n_cpss_stable_significant"),
            recovery_rate=vm.get("known_mechanism_recovery_rate"),
            H2_pass=vm.get("H2_pass"), novel_fraction=vm.get("novel_candidate_fraction"),
        ))
    _write(msum, list(msum[0].keys()), out / "models_summary.csv")
    _write(overview, list(overview[0].keys()), out / "kb_overview.csv")

    # ---- biomarkers (full per-unitig evidence table) -----------------------
    bio = []
    q = """
      SELECT s.model_id, s.unitig_id, u.sequence, s.gain, s.mean_abs_shap,
             s.selection_frequency, s.stable, s.composite_score, s.selection_method,
             b.gene_symbol, b.tier, b.source_db, b.identity_pct, b.coverage, b.aro_accession,
             b.aro_gene_family, b.aro_drug_class, b.aro_resistance_mechanism,
             f.prevalence_resistant, f.prevalence_susceptible,
             -- Derive the gap when the KB has none: step 10 did not emit
             -- delta_prevalence before 2026-08, so every KB written earlier stores NULL
             -- even though both prevalences are there.
             COALESCE(f.delta_prevalence,
                      f.prevalence_resistant - f.prevalence_susceptible) AS delta_prevalence,
             f.odds_ratio, f.fisher_p, f.discriminative,
             et.evidence_tier, et.n_evidence_layers, et.evidence_layers, et.is_novel_candidate
      FROM (
          -- One row per (unitig, model). unitig_model_scores is keyed on
          -- (unitig_id, model_id, selection_method), so a biomarker found by BOTH the
          -- gain path (07/07b) and CPSS (13) legitimately has two rows — but a table
          -- of biomarkers must list it once, or it is counted twice everywhere. The
          -- two rows populate disjoint columns (gain/in_gain_topn vs
          -- selection_frequency/stable/mean_abs_shap/composite_score), so MAX() merges
          -- them without losing anything, and selection_method records both.
          SELECT unitig_id, model_id,
                 MAX(gain)                AS gain,
                 MAX(in_gain_topn)        AS in_gain_topn,
                 MAX(selection_frequency) AS selection_frequency,
                 MAX(stable)              AS stable,
                 MAX(composite_score)     AS composite_score,
                 MAX(mean_abs_shap)       AS mean_abs_shap,
                 GROUP_CONCAT(DISTINCT selection_method) AS selection_method
          FROM unitig_model_scores GROUP BY unitig_id, model_id
      ) s
      LEFT JOIN (
          -- ONE BLAST row per (unitig, model). A unitig can carry both a CARD and an
          -- NCBI hit, and joining the raw table multiplied those biomarkers into
          -- several rows: biomarkers.csv came out at 3693 rows for 3571 biomarkers and
          -- every tier count was inflated (strong_novel read 88 instead of the KB's 23).
          -- Best hit = CARD first (gene_symbol/tier/ARO come from there), then higher
          -- identity, then lower e-value.
          SELECT * FROM (
              SELECT *, ROW_NUMBER() OVER (
                         PARTITION BY unitig_id, model_id
                         ORDER BY (source_db='card') DESC,
                                  COALESCE(identity_pct, -1) DESC,
                                  COALESCE(evalue, 1e9) ASC) AS _rn
              FROM blast_annotations
          ) WHERE _rn = 1
      ) b ON b.unitig_id=s.unitig_id AND b.model_id=s.model_id
      LEFT JOIN unitig_background_frequency f ON f.unitig_id=s.unitig_id AND f.model_id=s.model_id
      LEFT JOIN unitig_evidence_tier et ON et.unitig_id=s.unitig_id AND et.model_id=s.model_id
      LEFT JOIN unitigs u ON u.unitig_id=s.unitig_id
    """
    id2ab = {m["model_id"]: (m["antibiotic"], _org(m["run_id"])) for m in models}
    mda = {}
    for r in c.execute("SELECT unitig_id, evidence_score, pipeline_run_id FROM validation_evidence WHERE evidence_type='permutation_mda'"):
        mda[(r["unitig_id"], r["pipeline_run_id"])] = r["evidence_score"]
    pysc = {}
    for r in c.execute("SELECT unitig_id, evidence_score, pipeline_run_id FROM validation_evidence WHERE evidence_type='pyseer_lmm'"):
        pysc[(r["unitig_id"], r["pipeline_run_id"])] = r["evidence_score"]
    run_of = {m["model_id"]: m["run_id"] for m in models}
    for r in c.execute(q):
        ab, org = id2ab[r["model_id"]]
        rid = run_of[r["model_id"]]
        d = dict(r)
        d["antibiotic"] = ab
        d["organism"] = org
        d["drug_class"] = cls_of.get(ab)
        d["mda_auc_drop"] = mda.get((r["unitig_id"], rid))
        d["pyseer_lrt_p"] = pysc.get((r["unitig_id"], rid))
        bio.append(d)
    bcols = ["model_id", "organism", "antibiotic", "drug_class", "unitig_id", "sequence",
             "gain", "mean_abs_shap", "selection_frequency", "stable", "composite_score",
                 "selection_method",
             "mda_auc_drop", "pyseer_lrt_p", "prevalence_resistant", "prevalence_susceptible",
             "delta_prevalence", "odds_ratio", "fisher_p", "discriminative",
             "gene_symbol", "tier", "source_db", "identity_pct", "coverage", "aro_accession",
             "aro_gene_family", "aro_drug_class", "aro_resistance_mechanism",
             "evidence_tier", "n_evidence_layers", "evidence_layers", "is_novel_candidate"]
    # Guard: biomarkers.csv is one row per (model, unitig). A LEFT JOIN that fans out
    # silently double-counts biomarkers in every downstream table, figure and thesis
    # number, and nothing else would notice — so fail loudly instead.
    _pairs = {(d["model_id"], d["unitig_id"]) for d in bio}
    if len(_pairs) != len(bio):
        raise SystemExit(
            f"ERROR: biomarkers rows ({len(bio)}) != distinct (model, unitig) pairs "
            f"({len(_pairs)}) — a join is fanning out; tier counts would be inflated."
        )
    _write(bio, bcols, out / "biomarkers.csv")

    # ---- mechanisms: on-target confirmed genes (class-filtered) ------------
    mech = []
    for m in models:
        mid, ab = m["model_id"], m["antibiotic"]
        org = _org(m["run_id"])
        kws = CLASS_TO_ARO_KEYWORD.get(cls_of.get(ab, ""), ())
        for r in c.execute(
            "SELECT gene_symbol, aro_gene_family, aro_drug_class, aro_resistance_mechanism, "
            "MAX(identity_pct) idp, MIN(evalue) ev, COUNT(*) n "
            "FROM blast_annotations WHERE model_id=? AND tier IN ('confirmed','candidate') "
            "AND gene_symbol IS NOT NULL GROUP BY gene_symbol ORDER BY n DESC", (mid,)):
            adc = (r["aro_drug_class"] or "").lower()
            on_target = any(k in adc for k in kws) if kws else None
            mech.append(dict(
                model_id=mid, organism=org, antibiotic=ab, drug_class=cls_of.get(ab),
                gene_symbol=r["gene_symbol"], aro_gene_family=r["aro_gene_family"],
                aro_drug_class=r["aro_drug_class"], mechanism=r["aro_resistance_mechanism"],
                on_target=on_target, max_identity=r["idp"], min_evalue=r["ev"], n_unitigs=r["n"]))
    _write(mech, list(mech[0].keys()), out / "mechanisms.csv")

    print(f"\nDONE — 4 tidy tables in {out}/ ({len(models)} models).")


if __name__ == "__main__":
    main()
