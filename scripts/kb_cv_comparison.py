#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Random-vs-lineage-aware CV comparison (the table reviewers ask for).

Both numbers come from step 07b with EVERYTHING held fixed except the grouping:
same 5 folds, same stratification, same step-04 hyperparameters, same models.
The canonical run uses StratifiedGroupKFold over PopPUNK lineages; the comparison
run (AMR_CV_MODE=random) uses plain StratifiedKFold and writes '*_randomcv.csv'.
So the gap between them is attributable to lineage leakage and nothing else.

    lineage : results/{org}/{ab}/04_evaluation/10_repeated_holdout_summary_{ab}.csv
    random  : same directory, ..._{ab}_randomcv.csv

Usage:
    python scripts/kb_cv_comparison.py --results results \
        --tables results/tables --out results/figures

Outputs:
    results/tables/cv_comparison.csv   one row per model + the panel summary in the log
    results/figures/07_cv_random_vs_lineage.png (+pdf)
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

LINEAGE_METHOD = "lineage_group_kfold"
RANDOM_METHOD = "random_stratified_kfold"


def _read_summary(path, expect_method):
    """(mean_auc, std_auc, cv_method) from a 07b summary CSV, or None if absent.

    The cv_method check is the point of the exercise: if a file that is supposed to
    hold lineage-CV actually holds a fallback holdout (or vice versa), the comparison
    would silently measure the wrong thing.
    """
    if not path.exists() or path.stat().st_size == 0:
        return None
    df = pd.read_csv(path)
    if "seed" not in df.columns or "roc_auc" not in df.columns:
        return None
    method = str(df["cv_method"].dropna().iloc[0]) if "cv_method" in df.columns else ""
    if not method.startswith(expect_method):
        print(f"  ⚠ SKIP {path.name}: cv_method='{method}' (expected {expect_method}*)")
        return None
    rows = df.set_index("seed")["roc_auc"]
    mean = float(rows.get("MEAN", np.nan))
    std = float(rows.get("STD", np.nan))
    return mean, std, method


def main():
    ap = argparse.ArgumentParser(description="Random vs lineage-aware CV comparison.")
    ap.add_argument("--results", default="results")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--out", default="results/figures")
    args = ap.parse_args()

    tables = Path(args.tables)
    ms = pd.read_csv(tables / "models_summary.csv")   # organism/antibiotic/drug_class

    rows = []
    for _, m in ms.iterrows():
        org, ab = m["organism"], m["antibiotic"]
        d = Path(args.results) / org / ab / "04_evaluation"
        lin = _read_summary(d / f"10_repeated_holdout_summary_{ab}.csv", LINEAGE_METHOD)
        rnd = _read_summary(d / f"10_repeated_holdout_summary_{ab}_randomcv.csv", RANDOM_METHOD)
        if lin is None or rnd is None:
            print(f"  ⚠ missing pair for {org}/{ab} — excluded")
            continue
        rows.append({
            "organism": org, "antibiotic": ab, "drug_class": m.get("drug_class"),
            "n_genomes": m.get("n_genomes"),
            "random_cv_auc": round(rnd[0], 4), "random_cv_std": round(rnd[1], 4),
            "lineage_cv_auc": round(lin[0], 4), "lineage_cv_std": round(lin[1], 4),
            "inflation": round(rnd[0] - lin[0], 4),
            "inflation_pct": round(100.0 * (rnd[0] - lin[0]) / lin[0], 2) if lin[0] else None,
        })

    if not rows:
        sys.exit("ERROR: no model has BOTH a lineage and a random CV summary.")

    # Tie-break on (organism, antibiotic). Sorting on inflation alone is stable, so a
    # tie inherits the input order -- which comes from directory enumeration and differs
    # between filesystems: K. pneumoniae ciprofloxacin and gentamicin both sit at
    # inflation 0.0384 and swapped rows between the laptop and the HPC. Same numbers,
    # different file, so the artefact stopped being byte-reproducible across machines.
    df = (pd.DataFrame(rows)
          .sort_values(["inflation", "organism", "antibiotic"],
                       ascending=[False, True, True])
          .reset_index(drop=True))
    tables.mkdir(parents=True, exist_ok=True)
    out_csv = tables / "cv_comparison.csv"
    df.to_csv(out_csv, index=False, quoting=csv.QUOTE_MINIMAL)
    print(f"  ✓ {out_csv}  ({len(df)} models)")

    # ---- panel summary ----------------------------------------------------
    infl = df["inflation"].to_numpy()
    print("\n" + "=" * 70)
    print("RANDOM vs LINEAGE-AWARE CV")
    print("=" * 70)
    print(f"  models compared      : {len(df)}")
    print(f"  random-CV   mean AUC : {df['random_cv_auc'].mean():.3f}")
    print(f"  lineage-CV  mean AUC : {df['lineage_cv_auc'].mean():.3f}")
    print(f"  mean inflation       : {infl.mean():+.3f} "
          f"(median {np.median(infl):+.3f}, max {infl.max():+.3f})")
    print(f"  random > lineage in  : {(infl > 0).sum()}/{len(df)} models")
    wilcox_note = ""
    try:
        from scipy.stats import wilcoxon
        stat, p_w = wilcoxon(df["random_cv_auc"], df["lineage_cv_auc"])
        print(f"  Wilcoxon signed-rank : W={stat:.1f}, p={p_w:.2e}")
        # Carried onto the figure: this p-value IS the claim the figure makes, and a
        # reader should not have to find it in a log or the running text.
        wilcox_note = f" · Wilcoxon signed-rank p = {p_w:.1e}"
    except Exception as e:                                  # pragma: no cover
        print(f"  (Wilcoxon unavailable: {e})")
    worst = df.iloc[0]
    print(f"  largest gap          : {worst['antibiotic']} ({worst['organism']}) "
          f"{worst['random_cv_auc']:.3f} -> {worst['lineage_cv_auc']:.3f}")
    print("=" * 70)

    # ---- figure -----------------------------------------------------------
    try:
        from kb_figures import _abbr, _colour, _short
    except Exception:                                       # pragma: no cover
        def _abbr(o): return o[:2].title()
        def _colour(o): return "#2c7fb8"
        def _short(a): return a
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(9.5, 0.42 * len(df) + 1.6))
    for yi, (_, r) in zip(y, df.iterrows()):
        ax.plot([r["lineage_cv_auc"], r["random_cv_auc"]], [yi, yi],
                color="lightgrey", lw=2, zorder=1)
        ax.scatter(r["random_cv_auc"], yi, color="#999999", s=30, zorder=2)
        ax.scatter(r["lineage_cv_auc"], yi, color=_colour(r["organism"]), s=48, zorder=3)
    ax.axvline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{_short(a)} ({_abbr(o)})"
                        for a, o in zip(df.antibiotic, df.organism)], fontsize=7.5)
    ax.set_xlabel("ROC-AUC")
    ax.set_xlim(min(0.40, float(df["lineage_cv_auc"].min()) - 0.04), 1.0)
    ax.set_title("Removing the lineage grouping inflates the AUC\n"
                 f"same folds, same models — mean inflation {infl.mean():+.3f} "
                 f"in {(infl > 0).sum()}/{len(df)} models{wilcox_note}\n"
                 "(sorted by gap)", fontsize=11)
    # One legend entry per organism, so the colour actually resolves to a name. The
    # legend used to say "dot coloured by organism" while offering no key for it.
    org_handles = [Line2D([], [], marker="o", ls="", color=_colour(o), label=_abbr(o))
                   for o in sorted(df["organism"].unique())]
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color="#999999", label="random 5-fold CV (lineage-blind)"),
        Line2D([], [], marker="o", ls="", color="#444444",
               label="lineage-aware CV (reported)"),
    ] + org_handles, fontsize=8, loc="lower left", frameon=False, ncol=2)
    ax.invert_yaxis()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "07_cv_random_vs_lineage.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / "07_cv_random_vs_lineage.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out}/07_cv_random_vs_lineage.png (+pdf)")


if __name__ == "__main__":
    main()
