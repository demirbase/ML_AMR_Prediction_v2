#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thesis figures, part 2 — model, tuning and statistical validation.

Reads the per-model artefacts every step already writes (Optuna history, ROC/PR/
calibration curves, bootstrap CIs, permutation nulls, CPSS frequencies, pyseer
associations) and turns them into panel-level figures.

⚠️ Curve figures (ROC/PR/calibration/probability) come from step 06, which scores a
SINGLE random split — they describe the fitted model, not its generalisation. The
generalisation metric is the lineage-CV AUC (kb_figures.py fig 01). Titles say so.

Usage:
    python scripts/kb_figures_model.py --results results --tables results/tables \
        --out results/figures [--only hpo,params,roc,pr,calib,opchar,forest,probdist,mda,qq,manhattan,cpss]

Figures (results/figures/):
    20_optuna_convergence     HPO search behaviour across models
    21_hyperparameters        where the 45 tunings landed
    22_roc_curves             per-organism ROC panels (single split)
    23_pr_curves              per-organism PR panels (single split)
    24_calibration            reliability curves
    25_operating_characteristics  sensitivity vs specificity at the operating point
    26_auc_forest             single-split AUC with bootstrap CI, vs lineage-CV
    27_probability_separation R/S predicted-probability separation
    28_mda_permutation        MDA importance vs its permutation null
    29_pyseer_qq              LMM p-value calibration (QQ)
    30_pyseer_manhattan       association strength across unitigs
    31_cpss_frequencies       CPSS selection-frequency distributions
"""
import argparse
import glob
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from kb_figures import _abbr, _colour, _display, _save, _short, _sortkey  # noqa: E402

ORG_ORDER = ["ecoli", "kpneumoniae", "staphylococcus_aureus",
             "acinetobacter_baumannii", "pseudomonas_aeruginosa", "enterococcus_faecium"]


def _organisms(ms):
    return [o for o in ORG_ORDER if o in set(ms.organism)] + \
           sorted(set(ms.organism) - set(ORG_ORDER))


def _p(results, org, ab, sub, name):
    return Path(results) / org / ab / sub / name


def _read(path):
    return pd.read_csv(path) if Path(path).exists() else None


def _grid(n, ncols=3, w=4.3, h=3.5):
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(w * ncols, h * nrows))
    axes = np.atleast_1d(axes).ravel()
    for a in axes[n:]:
        a.axis("off")
    return fig, axes


# ------------------------------------------------------------------ HPO
def fig_optuna(results, ms, out):
    """Search convergence: best-so-far objective per trial, one line per model."""
    fig, ax = plt.subplots(figsize=(9, 5))
    n = 0
    for r in ms.itertuples():
        df = _read(_p(results, r.organism, r.antibiotic, "03_model_optimization",
                      f"01_optuna_history_{r.antibiotic}.csv"))
        if df is None or "Objective_Value" not in df:
            continue
        v = df.loc[df.get("State", "COMPLETE") == "COMPLETE", "Objective_Value"].to_numpy()
        if v.size == 0:
            continue
        ax.plot(np.arange(1, v.size + 1), np.maximum.accumulate(v), lw=1,
                alpha=0.55, color=_colour(r.organism))
        n += 1
    ax.set_xlabel("Optuna trial"); ax.set_ylabel("best objective so far (validation AUC)")
    # Describe what the curves show, not what the config allows: patience=15 exists, but
    # most searches keep finding small improvements and run the full 30 trials, so the
    # earlier caption ("early stopping is why runs end before trial 30") contradicted
    # the very lines it labelled.
    ax.set_title(f"Hyperparameter search plateaus within a few trials ({n} models)\n"
                 "most of the achievable objective is reached by trial ~5; "
                 "the remainder are marginal gains", fontsize=10.5)
    fig.tight_layout()
    _save(fig, out, "20_optuna_convergence")


def fig_params(results, ms, out):
    """Which hyperparameters mattered, aggregated over the panel."""
    imp = {}
    for r in ms.itertuples():
        df = _read(_p(results, r.organism, r.antibiotic, "03_model_optimization",
                      f"02_optuna_importance_{r.antibiotic}.csv"))
        if df is None or "Hyperparameter" not in df:
            continue
        for h, v in zip(df["Hyperparameter"], df["Importance"]):
            imp.setdefault(h, []).append(float(v))
    if not imp:
        print("  (params: no importance files — skipped)"); return
    order = sorted(imp, key=lambda k: -np.median(imp[k]))
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    # `vert` is deprecated in matplotlib 3.11 and removed in 3.13; vertical is the
    # default, so simply not passing it works on every version.
    ax.boxplot([imp[k] for k in order], widths=0.6, showfliers=False)
    for i, k in enumerate(order, start=1):
        ax.scatter(np.random.default_rng(0).normal(i, 0.05, len(imp[k])), imp[k],
                   s=12, alpha=0.5, color="#2c7fb8")
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel("Optuna importance")
    ax.set_title("Which hyperparameters drive the objective (across models)", fontsize=11)
    fig.tight_layout()
    _save(fig, out, "21_hyperparameters")


# --------------------------------------------------------------- curves
def _curve_panel(results, ms, orgs, out, fname, sub, pattern, xcol, ycol,
                 xlabel, ylabel, title, diag=False, name=""):
    fig, axes = _grid(len(orgs))
    for ax, org in zip(axes, orgs):
        sub_ms = ms[ms.organism == org]
        drawn = 0
        for r in sub_ms.itertuples():
            df = _read(_p(results, org, r.antibiotic, sub, pattern.format(ab=r.antibiotic)))
            if df is None or xcol not in df or ycol not in df:
                continue
            ax.plot(df[xcol], df[ycol], lw=1.1, alpha=0.85, label=_short(r.antibiotic))
            drawn += 1
        if diag:
            ax.plot([0, 1], [0, 1], ls="--", c="grey", lw=0.8)
        ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
        ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
        ax.set_title(f"{_display(org)} ({drawn})", fontsize=9.5, style="italic")
        if drawn:
            ax.legend(fontsize=6, frameon=False, loc="lower right" if diag else "best")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    _save(fig, out, fname)
    del name


def fig_roc(results, ms, orgs, out):
    _curve_panel(results, ms, orgs, out, "22_roc_curves", "04_evaluation",
                 "02_roc_curve_{ab}.csv", "False_Positive_Rate", "True_Positive_Rate",
                 "false positive rate", "true positive rate",
                 "ROC curves — SINGLE-SPLIT fit (generalisation is the lineage-CV AUC)",
                 diag=True)


def fig_pr(results, ms, orgs, out):
    _curve_panel(results, ms, orgs, out, "23_pr_curves", "04_evaluation",
                 "03_pr_curve_{ab}.csv", "Recall", "Precision",
                 "recall", "precision",
                 "Precision-recall curves — SINGLE-SPLIT fit")


def fig_calibration(results, ms, orgs, out):
    _curve_panel(results, ms, orgs, out, "24_calibration", "04_evaluation",
                 "05_calibration_curve_{ab}.csv", "Mean_Predicted_Probability",
                 "Fraction_Positives", "mean predicted probability",
                 "observed fraction positive",
                 "Calibration (reliability) curves — SINGLE-SPLIT fit", diag=True)


def fig_opchar(results, ms, out):
    """Sensitivity vs specificity at each model's operating threshold."""
    rows = []
    for r in ms.itertuples():
        df = _read(_p(results, r.organism, r.antibiotic, "04_evaluation",
                      f"06_comprehensive_metrics_{r.antibiotic}.csv"))
        if df is None or "Sensitivity_Recall_TPR" not in df:
            continue
        rows.append({"organism": r.organism, "antibiotic": r.antibiotic,
                     "sens": float(df["Sensitivity_Recall_TPR"].iloc[0]),
                     "spec": float(df["Specificity_TNR"].iloc[0]),
                     "mcc": float(df["MCC"].iloc[0]) if "MCC" in df else np.nan})
    if not rows:
        print("  (opchar: no metrics files — skipped)"); return
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7.4, 6))
    for org, g in df.groupby("organism"):
        ax.scatter(g.spec, g.sens, s=48, color=_colour(org), edgecolor="k", lw=0.35,
                   label=_display(org))
    for t in df.itertuples():
        if t.sens < 0.75 or t.spec < 0.75:                # label only the weak corner
            ax.annotate(f"{_short(t.antibiotic)} ({_abbr(t.organism)})", (t.spec, t.sens),
                        textcoords="offset points", xytext=(6, 2), fontsize=6.5, color="#555")
    ax.axhline(0.9, ls=":", c="grey", lw=0.8); ax.axvline(0.9, ls=":", c="grey", lw=0.8)
    ax.margins(0.12)            # annotations near the corners were being clipped
    ax.set_xlabel("specificity (TNR)"); ax.set_ylabel("sensitivity (TPR)")
    ax.set_title("Operating characteristics at each model's threshold\n(single split; dotted = 0.9)",
                 fontsize=10.5)
    ax.legend(fontsize=7.5, frameon=False, loc="lower left")
    fig.tight_layout()
    _save(fig, out, "25_operating_characteristics")


def fig_forest(results, ms, out):
    """Single-split AUC with bootstrap CI, next to the lineage-CV value: the CI
    quantifies sampling noise WITHIN a split — it says nothing about lineage leakage,
    which is why the two can disagree far beyond the interval."""
    rows = []
    for r in _sortkey(ms).itertuples():
        df = _read(_p(results, r.organism, r.antibiotic, "04_evaluation",
                      f"08_bootstrap_ci_{r.antibiotic}.csv"))
        if df is None or "ROC_AUC" not in df:
            continue
        rows.append({"organism": r.organism, "antibiotic": r.antibiotic,
                     "auc": float(df["ROC_AUC"].iloc[0]),
                     "lo": float(df["ROC_AUC_CI_low"].iloc[0]),
                     "hi": float(df["ROC_AUC_CI_high"].iloc[0]),
                     "lineage": r.lineage_cv_auc})
    if not rows:
        print("  (forest: no bootstrap files — skipped)"); return
    df = pd.DataFrame(rows)
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(9, 0.32 * len(df) + 1.6))
    ax.hlines(y, df.lo, df.hi, color="#bbbbbb", lw=2)
    ax.scatter(df.auc, y, s=34, color=[_colour(o) for o in df.organism], zorder=3,
               edgecolor="k", lw=0.3)
    ax.scatter(df.lineage, y, marker="|", s=90, color="black", linewidths=1.4, zorder=4)
    ax.axvline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(df.antibiotic, df.organism)],
                       fontsize=7)
    ax.invert_yaxis(); ax.set_xlabel("ROC-AUC")
    ax.set_title("Single-split AUC with 95 % bootstrap CI (dot+bar)\n"
                 "black tick = lineage-CV AUC — the gap is leakage, not sampling noise",
                 fontsize=10.5)
    fig.tight_layout()
    _save(fig, out, "26_auc_forest")


def fig_probdist(results, ms, orgs, out):
    """How well the predicted probabilities separate R from S."""
    fig, axes = _grid(len(orgs))
    for ax, org in zip(axes, orgs):
        sub = ms[ms.organism == org]
        pooled_r, pooled_s = [], []
        for r in sub.itertuples():
            df = _read(_p(results, org, r.antibiotic, "04_evaluation",
                          f"04_probability_distribution_{r.antibiotic}.csv"))
            if df is None or "True_Label" not in df:
                continue
            pooled_r += df.loc[df.True_Label == 1, "Predicted_Probability"].tolist()
            pooled_s += df.loc[df.True_Label == 0, "Predicted_Probability"].tolist()
        if not pooled_r:
            ax.axis("off"); continue
        ax.hist(pooled_s, bins=40, range=(0, 1), alpha=0.6, color="#7fbf7f",
                density=True, label=f"susceptible (n={len(pooled_s):,})")
        ax.hist(pooled_r, bins=40, range=(0, 1), alpha=0.6, color="#d62728",
                density=True, label=f"resistant (n={len(pooled_r):,})")
        ax.set_xlabel("predicted probability"); ax.set_ylabel("density")
        ax.set_title(_display(org), fontsize=9.5, style="italic")
        ax.legend(fontsize=6.5, frameon=False)
    fig.suptitle("Predicted-probability separation, pooled over each organism's models", fontsize=12)
    fig.tight_layout()
    _save(fig, out, "27_probability_separation")


# ---------------------------------------------------------- validation
def fig_mda(results, ms, out):
    """MDA permutation importance: how many candidates survive their own null.
    Few do — correlated unitigs cover for each other when one is permuted, which is
    exactly why CPSS + pyseer carry the selection argument instead."""
    rows, pooled = [], []
    for r in ms.itertuples():
        df = _read(_p(results, r.organism, r.antibiotic, "05_explainability",
                      f"12_permutation_test_{r.antibiotic}.csv"))
        if df is None or "mda_auc_drop" not in df:
            continue
        sig = int(df.get("permutation_significant", pd.Series(dtype=int)).sum())
        pooled += pd.to_numeric(df["mda_auc_drop"], errors="coerce").dropna().tolist()
        rows.append({"organism": r.organism, "antibiotic": r.antibiotic,
                     "tested": len(df), "significant": sig,
                     "median_drop": float(df["mda_auc_drop"].median()),
                     "max_drop": float(df["mda_auc_drop"].max())})
    if not rows:
        print("  (mda: no permutation files — skipped)"); return
    df = pd.DataFrame(rows)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    x = np.arange(len(df))
    a1.bar(x, df.max_drop, color=[_colour(o) for o in df.organism], edgecolor="k", lw=0.35)
    a1.axhline(0, c="grey", lw=0.8)
    a1.set_xticks(x); a1.set_xticklabels([f"{_short(a)} ({_abbr(o)})"
                                          for a, o in zip(df.antibiotic, df.organism)],
                                         rotation=90, fontsize=6)
    a1.set_ylabel("largest single-feature AUC drop")
    a1.set_title("MDA permutation: strongest individual effect", fontsize=10)
    # A histogram of the per-model significant COUNT is a single bar at zero — true but
    # empty. Show the distribution of the effects themselves, which is what "no single
    # unitig is indispensable" actually means.
    if pooled:
        a2.hist(pooled, bins=60, color="#756bb1", edgecolor="k")
        a2.set_yscale("log")
        a2.axvline(0, c="grey", lw=0.9)
        n_sig = int(df.significant.sum())
        a2.set_xlabel("AUC drop when a single feature is permuted")
        a2.set_ylabel("features (log)")
        # The power caveat is not optional. Read without it, "0 pass FDR" invites the
        # conclusion that no feature matters, when the honest reading is that this test
        # cannot resolve one: R=100 permutations against ~2,400 candidates puts a coarse
        # floor under the attainable q-value.
        a2.set_title("Almost no single unitig is individually indispensable\n"
                     f"{len(pooled):,} candidates tested · {n_sig} pass FDR in "
                     f"{int((df.significant > 0).sum())}/{len(df)} models — linked features "
                     "substitute for each other (hence CPSS/pyseer)\n"
                     "underpowered by design: 100 permutations vs "
                     f"{len(pooled):,} candidates bounds how small q can get, so this is "
                     "'not resolvable here', not 'no effect'", fontsize=8.5)
    fig.tight_layout()
    _save(fig, out, "28_mda_permutation")


def _pyseer_assoc(results, org, ab):
    """All tested variants (not just the significant ones) from pyseer's raw output."""
    hits = glob.glob(str(_p(results, org, ab, "05_explainability", f"14_pyseer_assoc_{ab}.txt")))
    if not hits:
        return None
    try:
        df = pd.read_csv(hits[0], sep="\t")
    except Exception:
        return None
    col = "lrt-pvalue" if "lrt-pvalue" in df.columns else (
        "filter-pvalue" if "filter-pvalue" in df.columns else None)
    if col is None:
        return None
    p = pd.to_numeric(df[col], errors="coerce").dropna()
    return p[p > 0]


def _lambda_gc(p):
    """Genomic inflation factor: median observed chi-square over its null median.

    lambda ~ 1 is calibrated; lambda >> 1 means the bulk of the distribution is
    shifted, not just the tail.
    """
    p = np.asarray(p, dtype=float)
    p = p[(p > 0) & (p <= 1)]
    if p.size == 0:
        return float("nan")
    try:
        from scipy.stats import chi2, norm
        obs = float(np.median(chi2.isf(p, 1)))
        return obs / float(chi2.ppf(0.5, 1))
    except Exception:                                        # pragma: no cover
        from math import sqrt
        from statistics import NormalDist
        nd = NormalDist()
        z = np.array([nd.inv_cdf(max(pi / 2, 1e-300)) for pi in p])
        return float(np.median(z ** 2)) / 0.4549364231195724


def fig_pyseer_qq(results, ms, orgs, out):
    """QQ of the LMM p-values, with the genomic inflation factor per model.

    A well-calibrated mixed model tracks the diagonal in the bulk and lifts only in
    the tail. These curves lift throughout, so the figure must NOT be captioned as
    evidence that the random effect absorbed population structure — it previously
    was, which asserted the opposite of what it shows. Two things drive the lift and
    the plot alone cannot separate them: genuine signal carried by thousands of
    unitigs in tight LD across the same locus, and any stratification the kinship
    term did not remove. Lambda is printed so the magnitude is at least quantified.
    """
    fig, axes = _grid(len(orgs))
    lambdas = []
    for ax, org in zip(axes, orgs):
        drawn, xmax = 0, 1.0
        for r in ms[ms.organism == org].itertuples():
            p = _pyseer_assoc(results, org, r.antibiotic)
            if p is None or len(p) < 50:
                continue
            obs = -np.log10(np.sort(p.to_numpy()))
            exp = -np.log10(np.linspace(1 / len(obs), 1, len(obs)))
            lam = _lambda_gc(p.to_numpy())
            if np.isfinite(lam):
                lambdas.append(lam)
            ax.plot(exp, obs, lw=0.9, alpha=0.8,
                    label=f"{_short(r.antibiotic)}  λ={lam:.1f}")
            xmax = max(xmax, float(exp.max()))
            drawn += 1
        # The expected axis can only reach log10(n_tested) ≈ 3.7 for 5 000 variants while
        # the observed tail reaches 250+. Forcing a square view (the old max-of-both
        # limit) squashed the expected axis to nothing and every curve looked like a
        # vertical line. Scale each axis to its own range and draw y=x over x only.
        ax.plot([0, xmax], [0, xmax], ls="--", c="grey", lw=0.8)
        ax.set_xlim(0, xmax * 1.05)
        ax.set_xlabel("expected  −log10(p)"); ax.set_ylabel("observed  −log10(p)")
        ax.set_title(f"{_display(org)} ({drawn})", fontsize=9.5, style="italic")
        if drawn:
            ax.legend(fontsize=6, frameon=False)
    rng_note = (f"  ·  genomic inflation λ = {min(lambdas):.1f}–{max(lambdas):.1f} "
                f"(median {np.median(lambdas):.1f})" if lambdas else "")
    fig.suptitle("pyseer LMM p-values, QQ vs the uniform null" + rng_note + "\n"
                 "the lift is not confined to the tail: unitigs in tight LD tag the same "
                 "locus thousands of times, and the plot cannot separate that from "
                 "residual stratification —\nso this is a diagnostic of the association "
                 "scan, not a certificate that the kinship term absorbed population "
                 "structure", fontsize=10)
    fig.tight_layout()
    _save(fig, out, "29_pyseer_qq")


def fig_pyseer_manhattan(results, ms, orgs, out):
    """Association strength across the tested unitigs, with each model's Bonferroni
    threshold (unitigs have no coordinates, so the x-axis is rank, not position)."""
    fig, axes = _grid(len(orgs))
    for ax, org in zip(axes, orgs):
        drawn = 0
        for r in ms[ms.organism == org].itertuples():
            p = _pyseer_assoc(results, org, r.antibiotic)
            if p is None or len(p) < 50:
                continue
            v = -np.log10(p.to_numpy())
            ax.scatter(np.arange(len(v)), v, s=2, alpha=0.35, label=_short(r.antibiotic))
            js = _p(results, org, r.antibiotic, "05_explainability",
                    f"14_pyseer_summary_{r.antibiotic}.json")
            if js.exists():
                thr = json.loads(js.read_text()).get("bonferroni_threshold")
                if thr:
                    ax.axhline(-np.log10(thr), ls="--", lw=0.8, c="#d62728")
            drawn += 1
        ax.set_xlabel("unitig (tested order)"); ax.set_ylabel("−log10(p)")
        ax.set_title(f"{_display(org)} ({drawn})", fontsize=9.5, style="italic")
        if drawn:
            ax.legend(fontsize=6, frameon=False, markerscale=4)
    fig.suptitle("pyseer LMM associations — red line = Bonferroni threshold "
                 "(x is rank: unitigs have no genome coordinate)", fontsize=12)
    fig.tight_layout()
    _save(fig, out, "30_pyseer_manhattan")


def fig_cpss(results, ms, out):
    """CPSS selection frequencies vs the π≥0.6 stability threshold."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.6))
    pooled, per_model = [], []
    for r in ms.itertuples():
        df = _read(_p(results, r.organism, r.antibiotic, "05_explainability",
                      f"13_stability_selection_{r.antibiotic}.csv"))
        if df is None or "selection_frequency" not in df:
            continue
        pooled += df["selection_frequency"].tolist()
        per_model.append({"organism": r.organism, "antibiotic": r.antibiotic,
                          "n_stable": int((df["selection_frequency"] >= 0.6).sum()),
                          "n_candidates": len(df)})
    if not pooled:
        print("  (cpss: no stability files — skipped)"); return
    a1.hist(pooled, bins=40, range=(0, 1), color="#31a354", edgecolor="k")
    a1.axvline(0.6, ls="--", c="#d62728", lw=1.2)
    a1.text(0.62, 0.86, "π = 0.6\n(stable)", transform=a1.transAxes,
            fontsize=8, color="#d62728")
    a1.set_yscale("log")
    a1.set_xlabel("CPSS selection frequency"); a1.set_ylabel("features (log)")
    a1.set_title(f"Selection frequency across all candidates ({len(pooled):,})", fontsize=10)
    pm = pd.DataFrame(per_model)
    x = np.arange(len(pm))
    a2.bar(x, pm.n_stable, color=[_colour(o) for o in pm.organism], edgecolor="k", lw=0.35)
    a2.set_xticks(x); a2.set_xticklabels([f"{_short(a)} ({_abbr(o)})"
                                          for a, o in zip(pm.antibiotic, pm.organism)],
                                         rotation=90, fontsize=6)
    a2.set_ylabel("stable features (π≥0.6)")
    a2.set_title("Stable set size per model", fontsize=10)
    fig.tight_layout()
    _save(fig, out, "31_cpss_frequencies")


def main():
    ap = argparse.ArgumentParser(description="Thesis figures — model, tuning, validation.")
    ap.add_argument("--results", default="results")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--out", default="results/figures")
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    ms = pd.read_csv(Path(args.tables) / "models_summary.csv")
    orgs = _organisms(ms)
    out = Path(args.out)
    only = {s.strip() for s in args.only.split(",") if s.strip()}

    todo = [
        ("hpo",       lambda: fig_optuna(args.results, ms, out)),
        ("params",    lambda: fig_params(args.results, ms, out)),
        ("roc",       lambda: fig_roc(args.results, ms, orgs, out)),
        ("pr",        lambda: fig_pr(args.results, ms, orgs, out)),
        ("calib",     lambda: fig_calibration(args.results, ms, orgs, out)),
        ("opchar",    lambda: fig_opchar(args.results, ms, out)),
        ("forest",    lambda: fig_forest(args.results, ms, out)),
        ("probdist",  lambda: fig_probdist(args.results, ms, orgs, out)),
        ("mda",       lambda: fig_mda(args.results, ms, out)),
        ("qq",        lambda: fig_pyseer_qq(args.results, ms, orgs, out)),
        ("manhattan", lambda: fig_pyseer_manhattan(args.results, ms, orgs, out)),
        ("cpss",      lambda: fig_cpss(args.results, ms, out)),
    ]
    for name, fn in todo:
        if only and name not in only:
            continue
        try:
            fn()
        except Exception as e:
            print(f"  ✗ {name} failed: {type(e).__name__}: {e}")
    print("DONE.")


if __name__ == "__main__":
    main()
