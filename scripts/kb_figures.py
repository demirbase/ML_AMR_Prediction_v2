#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thesis figures for the unified AMR-KB, built from the tidy tables that
kb_tables.py exports (+ raw 12b null CSVs for the permutation histogram).

Run kb_tables.py FIRST, then:
    python scripts/kb_figures.py --tables results/tables --results results \
        --out figures [--only performance,cpss_pfer,cross_org,mechanism,null_hist]

Each figure is saved as PNG (200 dpi) + PDF. Colours: E. coli blue, K. pneumoniae
red (+ a 3rd colour auto-assigned for any further organism). Edit CLASS_ORDER /
palette below to taste — this is a scaffold you own, not a black box.
"""
import argparse
import glob
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

CLASS_ORDER = ["penicillins", "cephalosporins", "beta_lactams_carbapenems_others",
               "quinolones", "aminoglycosides", "tetracyclines",
               "folate_pathway_inhibitors"]
PALETTE = {"ecoli": "#2c7fb8", "kpneumoniae": "#de2d26",
           "paeruginosa": "#31a354", "saureus": "#756bb1"}
_EXTRA = ["#e6ab02", "#a6761d", "#666666"]


def _colour(org, _cache={}):
    if org in PALETTE:
        return PALETTE[org]
    return _cache.setdefault(org, _EXTRA[len(_cache) % len(_EXTRA)])


def _short(ab):
    return ab.replace("_", "/")[:18]


def _abbr(org):
    return {"ecoli": "Ec", "kpneumoniae": "Kp", "paeruginosa": "Pa", "saureus": "Sa"}.get(org, org[:2].title())


def _sortkey(df):
    df = df.copy()
    df["_c"] = df["drug_class"].map({c: i for i, c in enumerate(CLASS_ORDER)}).fillna(99)
    return df.sort_values(["_c", "organism", "antibiotic"])


def _legend(ax, orgs):
    ax.legend(handles=[Patch(color=_colour(o), label={"ecoli": "E. coli", "kpneumoniae": "K. pneumoniae"}.get(o, o))
                       for o in orgs], fontsize=9, loc="lower left")


def _save(fig, out, name):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {out/name}.png (+pdf)")


def fig_performance(tables, out):
    df = _sortkey(pd.read_csv(tables / "models_summary.csv"))
    x = np.arange(len(df))
    col = [_colour(o) for o in df["organism"]]
    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.bar(x, df["lineage_cv_auc"], yerr=df["lineage_cv_std"], color=col,
           capsize=3, edgecolor="black", linewidth=0.4, alpha=0.9)
    ax.axhline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_ylim(0.4, 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_short(a)}\n({_abbr(o)})" for a, o in zip(df.antibiotic, df.organism)], rotation=90, fontsize=7.5)
    ax.set_ylabel("Lineage-aware CV ROC-AUC (mean ± SD)")
    ax.set_title("Per-antibiotic generalisation performance")
    _legend(ax, df["organism"].unique())
    _save(fig, out, "01_performance_lineageCV")


def fig_cpss_pfer(tables, out):
    df = _sortkey(pd.read_csv(tables / "kb_overview.csv"))
    x = np.arange(len(df))
    col = [_colour(o) for o in df["organism"]]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True, gridspec_kw={"hspace": 0.08})
    a1.bar(x, df["cpss_n_stable"], color=col, edgecolor="k", lw=0.4, alpha=0.9)
    a1.set_ylabel("CPSS stable unitigs (π≥0.6)")
    a1.set_title("CPSS stability selection — stable biomarker count & PFER bound")
    a2.bar(x, df["pfer_bound"], color=col, edgecolor="k", lw=0.4, alpha=0.9)
    a2.set_yscale("log")
    a2.axhline(1, ls="--", c="grey", lw=0.8)
    a2.set_ylabel("PFER bound (E[false positives], log)")
    a2.set_xticks(x)
    a2.set_xticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(df.antibiotic, df.organism)], rotation=90, fontsize=7.5)
    _legend(a1, df["organism"].unique())
    _save(fig, out, "02_cpss_pfer")


def fig_cross_org(tables, out):
    """Drugs assayed in ≥2 organisms → do the confirmed on-target genes agree?"""
    mech = pd.read_csv(tables / "mechanisms.csv")
    ot = mech[mech["on_target"] == True]  # noqa: E712
    shared = [ab for ab, g in ot.groupby("antibiotic") if g["organism"].nunique() >= 2]
    if not shared:
        print("  (cross_org: no drug shared across organisms yet — skipped)")
        return
    fig, ax = plt.subplots(figsize=(10, 0.5 + 0.5 * len(shared)))
    for i, ab in enumerate(sorted(shared)):
        for org, g in ot[ot.antibiotic == ab].groupby("organism"):
            fams = sorted(set(g["aro_gene_family"].dropna()))[:4]
            ax.text(0.02 if org == "ecoli" else 0.52, i, f"{_abbr(org)}: {', '.join(fams) or '—'}",
                    fontsize=8, va="center", color=_colour(org))
        ax.text(-0.01, i, _short(ab), fontsize=8, va="center", ha="right", fontweight="bold")
    ax.set_ylim(-0.5, len(shared) - 0.5)
    ax.axis("off")
    ax.set_title("Cross-organism mechanism concordance (on-target confirmed gene families)")
    _save(fig, out, "03_cross_organism")


def fig_mechanism(tables, out):
    mech = pd.read_csv(tables / "mechanisms.csv")
    ot = mech[mech["on_target"] == True]  # noqa: E712
    top = (ot.sort_values("n_unitigs", ascending=False)
             .groupby(["model_id", "organism", "antibiotic", "drug_class"], as_index=False)
             .agg(genes=("gene_symbol", lambda s: ", ".join(sorted(set(s))[:3]))))
    top = _sortkey(top)
    fig, ax = plt.subplots(figsize=(11, 0.4 + 0.42 * len(top)))
    for i, r in enumerate(top.itertuples()):
        ax.text(0.0, i, f"{_short(r.antibiotic)} ({_abbr(r.organism)})", fontsize=8.5, va="center", fontweight="bold", color=_colour(r.organism))
        ax.text(0.42, i, r.genes or "—", fontsize=8.5, va="center")
    ax.set_ylim(-0.5, len(top) - 0.5)
    ax.axis("off")
    ax.set_title("On-target confirmed resistance mechanism per model (class-filtered)")
    _save(fig, out, "04_mechanism_panel")


def fig_null_hist(tables, results, out):
    files = sorted(glob.glob(f"{results}/*/*/05_explainability/12b_label_permutation_nulls_*.csv"))
    if not files:
        print("  (null_hist: no 12b null CSVs — skipped)")
        return
    n = len(files)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3 * ncol, 2.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, f in zip(axes, files):
        ab = os.path.basename(f).replace("12b_label_permutation_nulls_", "").replace(".csv", "")
        d = pd.read_csv(f)
        col = d.columns[0]
        nulls = pd.to_numeric(d[col], errors="coerce").dropna()
        summ = json.load(open(f.replace("_nulls_", "_summary_").replace(".csv", ".json"))) if os.path.exists(f.replace("_nulls_", "_summary_").replace(".csv", ".json")) else {}
        real = summ.get("real_roc_auc") or summ.get("real_test_roc_auc")
        ax.hist(nulls, bins=15, color="#999999", edgecolor="white")
        if real:
            ax.axvline(real, color="#d7301f", lw=2)
        ax.set_title(_short(ab), fontsize=8)
        ax.set_xlim(0.4, 1.0)
        ax.tick_params(labelsize=6)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Label-permutation null vs REAL ROC-AUC (red line) — model-level significance", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, out, "05_label_permutation_nulls")


FIGS = {"performance": lambda t, r, o: fig_performance(t, o),
        "cpss_pfer": lambda t, r, o: fig_cpss_pfer(t, o),
        "cross_org": lambda t, r, o: fig_cross_org(t, o),
        "mechanism": lambda t, r, o: fig_mechanism(t, o),
        "null_hist": lambda t, r, o: fig_null_hist(t, r, o)}


def main():
    ap = argparse.ArgumentParser(description="Thesis figures from the AMR-KB tidy tables.")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--only", default=None, help="comma list: " + ",".join(FIGS))
    args = ap.parse_args()
    tables, out = Path(args.tables), Path(args.out)
    want = args.only.split(",") if args.only else list(FIGS)
    for name in want:
        if name not in FIGS:
            print(f"  ! unknown figure '{name}'")
            continue
        FIGS[name](tables, args.results, out)
    print("DONE.")


if __name__ == "__main__":
    main()
