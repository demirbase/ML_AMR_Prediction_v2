#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Six thesis tables that the delivered artefacts imply but do not contain.

    lineage_summary        population structure per organism, next to the CV inflation
                           it explains -- lineage counts and singletons exist nowhere else
    evidence_accounting    the seven-produced / six-counted / four-firing ledger, as rows
    headline_biomarkers    biomarkers.csv (3,571 rows) reduced to what a results chapter
                           can print, under the mandatory CARD tier filter
    provenance_tools       what actually ran: nine tool versions, seed, config hashes
    hyperparameters        the Optuna choice per model, joined to the delivered model
    limitations            METHODOLOGY 5.3's nine limitations with the numbers re-measured

Two rules this script follows deliberately.

FIRST: every number is recomputed from the delivered artefacts. Where a limitation
cannot be recomputed locally -- BV-BRC collection dates, country shares, the
max_target_seqs flag -- the row says so in `evidence_source` instead of quietly
repeating prose. That distinction is the point of the limitations table: a reader
should be able to tell which caveats were verified and which are testimony.

Re-measuring found one discrepancy already: METHODOLOGY 5.3 states PFER max ~14,
while the delivered kb_overview.csv maxes at 12.9 across 21 of 45 models above 1.
The table carries the measured value and flags the difference.

SECOND: a table is never written from missing inputs. lineage_summary needs the
PopPUNK cluster CSVs; if they are absent the table is skipped with a message rather
than written empty over a good one -- the failure mode that cost figures 15 and 16
a review round.

    python scripts/kb_tables_thesis.py --db results/kb/amrk.db \
        --tables results/tables --data data/processed --runs runs \
        [--only lineage_summary,limitations]
"""
import argparse
import glob
import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

# layer token in unitig_evidence_tier.evidence_layers -> the validation_evidence type
# that produces it, and the step that writes it. label_permutation has no token: it
# grades a model, not a biomarker, which is exactly why 7 are produced and 6 counted.
LAYER_MAP = [
    ("blast",      "blast",                "08_blast_annotation",          "yes"),
    ("prevalence", "background_frequency",  "10_kmer_background_frequency", "yes"),
    ("snp",        "snp",                   "11_variant_snp_check",         "yes"),
    ("mda",        "permutation_mda",       "12_permutation_test",          "yes"),
    ("cpss",       "stability_selection",   "13_stability_selection",       "yes"),
    ("pyseer",     "pyseer_lmm",            "14_pyseer_lmm",                "yes"),
    (None,         "label_permutation",     "12b_label_permutation_test",   "no"),
]


# ------------------------------------------------------------------ B1 lineage
def t_lineage_summary(ctx):
    files = sorted(glob.glob(str(Path(ctx["data"]) / "*" / "lineage" / "poppunk_clusters.csv")))
    if not files:
        print("  skip lineage_summary — no poppunk_clusters.csv under "
              f"{ctx['data']} (nothing written)")
        return None
    ms, cv = ctx["models_summary"], ctx["cv_comparison"]
    rows = []
    for f in files:
        org = Path(f).parts[-3]
        d = pd.read_csv(f)
        col = [c for c in d.columns if "luster" in c][-1]
        vc = d[col].value_counts()
        p = vc / vc.sum()
        mm, mc = ms[ms.organism == org], cv[cv.organism == org]
        rows.append({
            "organism": org,
            "n_genomes_clustered": len(d),
            "n_lineages": int(d[col].nunique()),
            "largest_lineage_n": int(vc.iloc[0]),
            "clonality_pct": round(100 * vc.iloc[0] / len(d), 1),
            "n_singleton_lineages": int((vc == 1).sum()),
            "singleton_pct_of_lineages": round(100 * (vc == 1).sum() / d[col].nunique(), 1),
            "simpson_diversity": round(1 - float((p ** 2).sum()), 4),
            "shannon_diversity": round(float(-(p * np.log(p)).sum()), 4),
            "n_models": int(len(mm)),
            "mean_lineage_cv_auc": round(float(mm.lineage_cv_auc.mean()), 4) if len(mm) else None,
            "mean_random_cv_auc": round(float(mc.random_cv_auc.mean()), 4) if len(mc) else None,
            "mean_inflation": round(float(mc.inflation.mean()), 4) if len(mc) else None,
        })
    out = pd.DataFrame(rows).sort_values("clonality_pct", ascending=False)
    ctx["_lineage_stats"] = _lineage_inflation_stats(out)
    return out


def _lineage_inflation_stats(d):
    """Correlate five structure measures against CV inflation -- and report all five.

    Reporting only the strongest would be selection: with n=6 organisms, five candidate
    measures and no pre-registration, the best p-value is not evidence. Simpson diversity
    happens to hold up under Spearman where largest-lineage share does not, and that is
    interesting, but it is only defensible if the reader can see the whole set.
    """
    try:
        from scipy import stats
    except ImportError:
        return None
    res = {}
    for col in ["clonality_pct", "simpson_diversity", "shannon_diversity",
                "n_lineages", "n_singleton_lineages"]:
        r, pr = stats.pearsonr(d[col], d.mean_inflation)
        rho, ps = stats.spearmanr(d[col], d.mean_inflation)
        # p-values keep six decimals: rounded to four, Simpson's 0.004543 becomes
        # 0.0045, and a downstream "%.3f" then renders it 0.004 (round-half-even) while
        # the same p computed from raw renders 0.005. The figure and the text would
        # disagree on a number the thesis quotes.
        res[col] = {"pearson_r": round(float(r), 4), "pearson_p": round(float(pr), 6),
                    "spearman_rho": round(float(rho), 4), "spearman_p": round(float(ps), 6)}
    return {
        "n_organisms": int(len(d)),
        "target": "mean_inflation (random-CV AUC minus lineage-CV AUC, averaged per organism)",
        "correlations": res,
        "caveat_n": "n=6 organisms. No coefficient here is confirmatory; the thesis should "
                    "call this a trend.",
        "caveat_selection": "Five measures were tested without pre-registration. Simpson "
                            "diversity gives the strongest association and is the only one "
                            "significant under Spearman, but reporting it alone would be "
                            "selection on the outcome -- all five are published together.",
        "note_consistency": "clonality_pct reproduces the value recorded during production "
                            "(Pearson r 0.914 here vs 0.912 recorded; rounding of the "
                            "clonality percentages).",
    }


# --------------------------------------------------------------- B2 evidence
def t_evidence_accounting(ctx):
    c = ctx["conn"]
    produced = dict(c.execute(
        "select evidence_type, count(*) from validation_evidence group by 1"))
    n_graded = c.execute("select count(*) from unitig_evidence_tier").fetchone()[0]
    rows = []
    for token, etype, step, counted in LAYER_MAP:
        if token is None:
            graded, note = 0, ("model-level: validates a model against its label-permutation "
                               "null, grades no biomarker — this is why 7 are produced but 6 "
                               "are counted")
            verdict = "not applicable to biomarkers"
        else:
            graded = c.execute(
                "select count(*) from unitig_evidence_tier "
                "where ','||evidence_layers||',' like ?", (f"%,{token},%",)).fetchone()[0]
            if graded == 0 and token == "snp":
                verdict, note = "produced, never fires", (
                    "wiring fault: step 11 scans step 07's FASTA while the graded universe "
                    "comes from steps 10+13, so the two sets do not intersect. 21 "
                    "resistant_allele findings are attached to no tier. Fixable, but the fix "
                    "changes the KB and therefore needs a new Zenodo version")
            elif graded == 0 and token == "mda":
                verdict, note = "produced, never fires", (
                    "real negative, not a fault: the candidate sets overlap fully (2409/2409) "
                    "but no unitig passes BH-FDR q<0.05. R=100 permutations is underpowered "
                    "for ~2400 candidates, so this must NOT be written as 'no feature matters'")
            else:
                verdict, note = "fires", "counted into the evidence tier"
        rows.append({
            "layer_token": token or "(none)",
            "kb_evidence_type": etype,
            "producing_step": step,
            "produced_evidence_rows": produced.get(etype, 0),
            "counted_into_tier": counted,
            "graded_pairs": graded,
            "share_of_graded_pct": round(100 * graded / n_graded, 1) if n_graded else None,
            "verdict": verdict,
            "note": note,
        })
    return pd.DataFrame(rows)


# -------------------------------------------------------------- B3 headline
def t_headline_biomarkers(ctx, top_n=3):
    b, ms = ctx["biomarkers"], ctx["models_summary"]
    # The CARD tier filter is not optional: 3,007 of 3,611 hits are tier='none', and
    # without it an A. baumannii model returns mecA.
    f = b[b.tier.isin(["confirmed", "candidate"])].copy()
    f = f.sort_values(["model_id", "n_evidence_layers", "composite_score",
                       "selection_frequency", "unitig_id"],
                      ascending=[True, False, False, False, True])
    # Rank DISTINCT genes: several unitigs routinely map to the same determinant, and a
    # top-3 of unitigs printed ErmB three times where the reader wanted three genes.
    # Keep the best-supported unitig per gene and record how many backed it.
    f["n_unitigs_same_gene"] = f.groupby(["model_id", "gene_symbol"]).unitig_id.transform("size")
    f = f.drop_duplicates(["model_id", "gene_symbol"], keep="first")
    f["rank_in_model"] = f.groupby("model_id").cumcount() + 1
    keep = f[f.rank_in_model <= top_n].merge(
        ms[["model_id", "lineage_cv_auc", "n_genomes"]], on="model_id", how="left")
    cols = ["model_id", "organism", "antibiotic", "drug_class", "n_genomes",
            "lineage_cv_auc", "rank_in_model", "gene_symbol", "n_unitigs_same_gene",
            "aro_gene_family",
            "aro_resistance_mechanism", "identity_pct", "coverage", "tier",
            "evidence_tier", "n_evidence_layers", "evidence_layers",
            "delta_prevalence", "odds_ratio", "pyseer_lrt_p", "selection_frequency",
            "composite_score", "unitig_id"]
    return keep[[c for c in cols if c in keep.columns]].sort_values(
        ["organism", "antibiotic", "rank_in_model"])


# ------------------------------------------------------------ B4 provenance
def t_provenance_tools(ctx):
    c = ctx["conn"]
    n_runs = c.execute("select count(*) from pipeline_runs").fetchone()[0]
    tools = [("CARD", "card_version"), ("KMC", "kmc_version"),
             ("XGBoost", "xgboost_version"), ("unitig-caller", "unitig_caller_version"),
             ("bcalm", "bcalm_version"), ("PopPUNK", "poppunk_version"),
             ("graph-tool", "graph_tool_version"), ("BLAST", "blast_version"),
             ("pyseer", "pyseer_version")]
    rows = []
    for label, col in tools:
        vals = [r[0] for r in c.execute(f"select distinct {col} from pipeline_runs")]
        val = vals[0] if len(vals) == 1 else " | ".join(str(v) for v in vals)
        note = ""
        if str(val) == "None":
            note = ("honest NULL: bcalm exposes no version CLI" if label == "bcalm"
                    else "not recorded")
            val = "not reported"      # else read_csv turns the text 'None' into NaN
        elif len(vals) > 1:
            note = f"NOT uniform across runs: {len(vals)} distinct values"
        rows.append({"item": label, "category": "tool", "value": val,
                     "n_runs_recording": n_runs, "note": note})
    for label, col, note in [
        ("random seed", "random_seed", "single seed across the whole panel"),
        ("min_support", "min_support", "unitig prevalence floor, fixed"),
    ]:
        vals = [r[0] for r in c.execute(f"select distinct {col} from pipeline_runs")]
        rows.append({"item": label, "category": "parameter",
                     "value": vals[0] if len(vals) == 1 else " | ".join(map(str, vals)),
                     "n_runs_recording": n_runs,
                     "note": note if len(vals) == 1 else f"{len(vals)} distinct values"})
    n_cfg = c.execute("select count(distinct config_hash) from pipeline_runs").fetchone()[0]
    commits = [r[0][:7] for r in c.execute("select distinct git_commit from pipeline_runs")]
    dirty = c.execute("select count(*) from pipeline_runs where git_dirty = 1").fetchone()[0]
    rows += [
        {"item": "config hash", "category": "provenance", "value": f"{n_cfg} distinct",
         "n_runs_recording": n_runs, "note": "one per run — the run's resolved configuration"},
        {"item": "git commit", "category": "provenance", "value": ", ".join(sorted(commits)),
         "n_runs_recording": n_runs,
         "note": f"the panel spans {len(commits)} commits; state this rather than implying one"},
        {"item": "git_dirty", "category": "provenance", "value": f"1 on {dirty}/{n_runs} runs",
         "n_runs_recording": n_runs,
         "note": "the working tree was patched beyond the recorded commit, so the hash alone "
                 "does not restore the code (METHODOLOGY 5.3, item 5)"},
        {"item": "KB schema", "category": "provenance",
         "value": c.execute("select kb_schema_version from kb_metadata").fetchone()[0],
         "n_runs_recording": n_runs, "note": "version-locked by test_version_alignment"},
    ]
    return pd.DataFrame(rows)


# ---------------------------------------------------------- B5 hyperparameters
def t_hyperparameters(ctx):
    files = glob.glob(str(Path(ctx["runs"]) / "*" / "*" / "*" / "run_metadata.json"))
    if not files:
        print(f"  skip hyperparameters — no run_metadata.json under {ctx['runs']}")
        return None
    c = ctx["conn"]
    known = {r[0] for r in c.execute("select run_id from pipeline_runs")}
    kb = pd.read_sql_query(
        "select model_id, run_id, n_trees, operating_threshold, auc_mean_seeds from models", c)
    rows = []
    for f in files:
        d = json.load(open(f))
        # Filtering on the KB's run_ids drops the two June leftovers in this directory
        # without having to special-case them by date or path.
        if d.get("run_id") not in known or not d.get("params"):
            continue
        r = {"run_id": d["run_id"], "organism": d.get("organism"),
             "antibiotic": d.get("antibiotic"), "n_trials": d.get("n_trials"),
             "best_score_optuna": d.get("best_score"), "n_genomes": d.get("n_genomes")}
        r.update({k: v for k, v in d["params"].items()})
        rows.append(r)
    if not rows:
        print("  skip hyperparameters — no run_metadata.json matched a KB run_id")
        return None
    out = pd.DataFrame(rows).drop_duplicates("run_id").merge(kb, on="run_id", how="left")
    front = ["model_id", "organism", "antibiotic", "n_genomes", "n_trials",
             "best_score_optuna", "auc_mean_seeds", "n_trees", "operating_threshold"]
    rest = [c_ for c_ in out.columns if c_ not in front + ["run_id"]]
    return out[front + rest + ["run_id"]].sort_values(["organism", "antibiotic"])


# ------------------------------------------------------------- B6 limitations
def t_limitations(ctx):
    c, k, ms = ctx["conn"], ctx["kb_overview"], ctx["models_summary"]
    q = lambda s: c.execute(s).fetchone()[0]
    pfer_max = round(float(k.pfer_bound.max()), 1)
    pfer_over = int((k.pfer_bound > 1).sum())
    n_conc = q("select count(*) from external_concordance")
    n_conc_models = q("select count(distinct model_id) from external_concordance")
    n_overlap = q("select count(*) from unitig_antibiotic_overlap")
    dirty = q("select count(*) from pipeline_runs where git_dirty = 1")
    n_runs = q("select count(*) from pipeline_runs")
    n_org_reg = q("select count(*) from organisms")
    n_org_panel = q("select count(distinct organism) from pipeline_runs")
    efm = round(float(ms[ms.organism == "enterococcus_faecium"].lineage_cv_auc.mean()), 3)
    snp_rows = q("select count(*) from variant_snp_check")
    snp_graded = q("select count(*) from unitig_evidence_tier"
                   " where evidence_layers like '%snp%'")
    nov = ctx["novel_ctx"]
    plasmid = int((nov.replicon_call == "plasmid").sum()) if nov is not None else None
    mixed = int((nov.replicon_call == "mixed").sum()) if nov is not None else None

    L = [
        (1, "No temporal or geographic hold-out",
         "Concordance against AMRFinderPlus and ResFinder ran on 2026-09-01 and covers the panel, "
         "but it scores the tools on the model's OWN held-out split, which is a chunk split rather "
         "than a lineage-aware one — a design that favours the model by the margin section 4.3 "
         "measures. Temporal validation remains impossible (BV-BRC AMR phenotypes end in 2021; "
         "≥2023 isolates: E. coli 28, K. pneumoniae 13, A. baumannii 11, S. aureus 0). Geographic "
         "validation was not performed and collections are country-dominated (E. coli 58% Norway, "
         "A. baumannii 63% USA).",
         f"recomputed: external_concordance holds {n_conc} rows over "
         f"{n_conc_models} of 45 models. Collection dates and country shares are NOT recomputable "
         f"from the delivered artefacts — they come from METHODOLOGY 5.3.",
         "partly recomputed", "3.x methods + 5.1 discussion"),
        (2, "Labels are BV-BRC as published",
         "No MIC re-interpretation was attempted: raw MIC completeness falls to 9% (S. aureus) "
         "and units are mixed.",
         "not recomputable locally — raw MIC fields are not in the KB.",
         "from METHODOLOGY 5.3", "3.x methods"),
        (3, f"PFER exceeds 1 in {pfer_over} of the {len(k)} models",
         "In those stable sets the expected number of false positives is above one.",
         f"recomputed: pfer_bound max = {pfer_max} over {pfer_over} of {len(k)} models above 1. "
         f"METHODOLOGY 5.3 read 'max ~14' until this table was built and is now corrected; "
         f"figure 02 had rendered {pfer_max} correctly all along, so the prose was the only "
         f"place the wrong number survived.",
         "recomputed", "4.4 stability results"),
        (4, "Co-carriage is linkage, not causation",
         "sul/qacEdelta1 co-occurrence is a class-1 integron; a novel K. pneumoniae gentamicin "
         "unitig maps to an MCR-1-carrying plasmid without being mcr-1. All KB claims are "
         "associational.",
         f"recomputed: of 23 novel biomarkers, {plasmid} sit on plasmids and {mixed} give a "
         f"mixed replicon signal (mobile-element signature) — novel_ncbi_context.csv, figure 36.",
         "recomputed", "4.9 novel candidates + 5.x"),
        (5, "git_dirty = 1 on every run",
         "Each run records its commit, but the working tree was modified at execution time, so "
         "the commit hash alone does not reconstruct the code. Seed, config hash, CARD version "
         "and tool versions are recorded and consistent; bcalm reports an honest NULL.",
         f"recomputed: git_dirty = 1 on {dirty}/{n_runs} runs.",
         "recomputed", "3.x reproducibility + FAIR R1.2"),
        (6, "nt replicon proportions are bounded samples",
         "The remote pass used max_target_seqs = 50, so replicon proportions are proportions of "
         "the alignments BLAST retained, not a census of nt.",
         "not recomputable locally — max_target_seqs is a runtime flag, not a stored field.",
         "from METHODOLOGY 5.3", "4.9 novel candidates"),
        (7, "The organisms table holds one more organism than the panel",
         "Enterobacter cloacae is registered in the reference table but has no trained model.",
         f"recomputed: organisms = {n_org_reg} rows, pipeline_runs covers {n_org_panel} "
         f"organisms. Needs a footnote wherever '7 organisms' could be inferred.",
         "recomputed", "3.1 dataset + 4.1 KB overview"),
        (8, "Genome QC enforces two of the four criteria it measures",
         "CheckM2/QUAST ran on all 17,742 assemblies and 98.7% (17,516) pass, but only "
         "completeness ≥95% and contamination ≤5% are enforced. N50 ≥50 kb and contigs ≤500 are "
         "computed then deliberately not applied: an N50 gate would remove 1,305 of 2,078 "
         "E. faecium genomes (63%), selecting on assembly provenance rather than genome quality. "
         "Assembly fragmentation is therefore an uncontrolled covariate.",
         f"recomputed: E. faecium — the organism most affected — carries the highest panel mean "
         f"AUC at {efm}. The 17,742/17,516 and 1,305/2,078 counts come from METHODOLOGY 5.3; "
         f"figure 12 plots all four criteria and marks which two are enforced.",
         "partly recomputed", "3.2 QC + 5.x discussion"),
        (9, "pyseer p-values are strongly inflated against a uniform null",
         "Genomic inflation λ = 0.7–79.3, median 2.6. Some is real — thousands of unitigs in "
         "tight LD tag the same locus — but the QQ plot cannot separate that from stratification "
         "the kinship term did not absorb. Figure 29 is a diagnostic of the scan, not evidence "
         "that population structure was controlled.",
         f"not recomputed here (λ is computed inside figure 29 from the per-model pyseer "
         f"outputs). Related and recomputed: variant_snp_check holds {snp_rows} rows, of which "
         f"{snp_graded} reach a tier since the 2026-09-01 loader fix, and "
         f"unitig_antibiotic_overlap holds {n_overlap} rows.",
         "from METHODOLOGY 5.3 + figure 29", "4.4 significance + 5.x"),
    ]
    return pd.DataFrame(L, columns=["n", "limitation", "detail", "evidence",
                                    "evidence_source", "affects"])


BUILDERS = {
    "lineage_summary": t_lineage_summary,
    "evidence_accounting": t_evidence_accounting,
    "headline_biomarkers": t_headline_biomarkers,
    "provenance_tools": t_provenance_tools,
    "hyperparameters": t_hyperparameters,
    "limitations": t_limitations,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True)
    ap.add_argument("--tables", required=True, help="tidy tables dir (read and write)")
    ap.add_argument("--data", default="data/processed", help="for the PopPUNK cluster CSVs")
    ap.add_argument("--runs", default="runs", help="for run_metadata.json hyperparameters")
    ap.add_argument("--out", help="output dir (default: --tables)")
    ap.add_argument("--only", help="comma list: " + ",".join(BUILDERS))
    a = ap.parse_args()

    tdir, out = Path(a.tables), Path(a.out or a.tables)
    out.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(a.db)
    # float_precision="round_trip" is not a nicety. pandas' default CSV float parser
    # is accurate to within an ULP, not exact, and which ULP it lands on varies by
    # version: on this laptop it read composite_score 23.890940407379954 back as
    # 23.89094040737995, and in the HPC container it did not. One value out of 480 was
    # enough to make headline_biomarkers.csv differ between the two machines while
    # every other table matched byte for byte. round_trip parses to the float whose
    # repr IS the source text, so a column copied from a tidy table is reproduced
    # exactly, on any machine, by any pandas.
    read = (lambda n: pd.read_csv(tdir / n, float_precision="round_trip")
            if (tdir / n).exists() else None)
    ctx = {"conn": conn, "data": a.data, "runs": a.runs,
           "biomarkers": read("biomarkers.csv"),
           "models_summary": read("models_summary.csv"),
           "cv_comparison": read("cv_comparison.csv"),
           "kb_overview": read("kb_overview.csv"),
           "novel_ctx": read("novel_ncbi_context.csv")}

    want = [s.strip() for s in a.only.split(",")] if a.only else list(BUILDERS)
    for name in want:
        if name not in BUILDERS:
            raise SystemExit(f"unknown table '{name}'; choose from {', '.join(BUILDERS)}")
        df = BUILDERS[name](ctx)
        if df is None:
            continue
        dest = out / f"{name}.csv"
        df.to_csv(dest, index=False)
        print(f"  wrote {dest.name}: {len(df)} rows x {len(df.columns)} cols")
        if name == "lineage_summary" and ctx.get("_lineage_stats"):
            js = out / "lineage_summary_stats.json"
            js.write_text(json.dumps(ctx["_lineage_stats"], indent=2), encoding="utf-8")
            print(f"  wrote {js.name}: 5 structure measures vs CV inflation")
    conn.close()


if __name__ == "__main__":
    main()
