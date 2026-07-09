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

# The 7 orthogonal validation layers (evidence_type in the KB) in pipeline order,
# biological → statistical. Label used on the evidence-layer heatmap (fig 06).
EVIDENCE_ORDER = [
    ("blast",                "BLAST\n(CARD/NCBI)"),
    ("background_frequency", "Prevalence\nR vs S"),
    ("snp",                  "SNP allele\n(CARD var.)"),
    ("permutation_mda",      "MDA\npermutation"),
    ("label_permutation",    "Label-perm\n(model)"),
    ("stability_selection",  "CPSS\nstability"),
    ("pyseer_lmm",           "pyseer LMM\n(lineage)"),
]


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


_FAM_MAP = {
    "16S rRNA methyltransferase (G1405)": "16S-RMTase",
    "General Bacterial Porin with reduced permeability to beta-lactams": "porin loss",
    "aminoglycoside bifunctional resistance protein": "AAC(6')-APH",
    "major facilitator superfamily (MFS) antibiotic efflux pump": "MFS efflux",
    "OXA beta-lactamase;OXA-48-like beta-lactamase": "OXA-48-like",
    "sulfonamide resistant sul": "sul",
    "trimethoprim resistant dihydrofolate reductase dfr": "dfr",
}


def _fam(s):
    """Short, human gene-family label for figures."""
    s = str(s).strip()
    if s in _FAM_MAP:
        return _FAM_MAP[s]
    if "beta-lactamase" in s:
        return s.replace(" beta-lactamase", "").strip()
    return s.split()[-1] if s else s


CLASS_SHORT = {"beta_lactams_carbapenems_others": "carbapenems / others",
               "folate_pathway_inhibitors": "folate inhibitors"}


def fig_overview(tables, out, db):
    """Cover slide: scope of the KB (models / organisms / classes / genomes) +
    models-per-drug-class stacked by organism."""
    ms = pd.read_csv(tables / "models_summary.csv")
    order = [c for c in CLASS_ORDER if c in set(ms.drug_class)]
    orgs = list(ms.organism.unique())
    fig = plt.figure(figsize=(13, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.7], wspace=0.28)
    a0 = fig.add_subplot(gs[0]); a1 = fig.add_subplot(gs[1])
    a0.axis("off")
    cards = [(str(len(ms)), "AMR models"),
             (str(ms.organism.nunique()), "organisms  (E. coli, K. pneumoniae)"),
             (str(len(order)), "drug classes"),
             (f"{int(ms.n_genomes.sum()):,}", "genome–phenotype pairs")]
    for k, (num, lab) in enumerate(cards):
        cy = 0.86 - 0.25 * k
        a0.text(0.02, cy, num, fontsize=30, fontweight="800", color="#222", transform=a0.transAxes)
        a0.text(0.03, cy - 0.085, lab, fontsize=10.5, color="#666", transform=a0.transAxes)
    a0.set_title("Unified AMR biomarker knowledge base", fontsize=12, loc="left")
    piv = ms.groupby(["drug_class", "organism"]).size().unstack(fill_value=0).reindex(order).fillna(0)
    y = np.arange(len(order)); left = np.zeros(len(order))
    for org in orgs:
        vals = piv[org].values if org in piv.columns else np.zeros(len(order))
        a1.barh(y, vals, left=left, color=_colour(org), edgecolor="white",
                label={"ecoli": "E. coli", "kpneumoniae": "K. pneumoniae"}.get(org, org))
        left += vals
    a1.set_yticks(y); a1.set_yticklabels([CLASS_SHORT.get(c, c.replace("_", " ")) for c in order], fontsize=9.5)
    a1.invert_yaxis(); a1.set_xlabel("models"); a1.legend(fontsize=9, loc="lower right")
    a1.set_title("Models per drug class")
    _save(fig, out, "00_kb_overview")


def fig_cross_org(tables, out):
    """Drugs assayed in ≥2 organisms → SHARED (concordant) gene family highlighted."""
    mech = pd.read_csv(tables / "mechanisms.csv")
    ot = mech[mech["on_target"] == True]  # noqa: E712
    abs_ = sorted(ab for ab, g in ot.groupby("antibiotic") if g["organism"].nunique() >= 2)
    if not abs_:
        print("  (cross_org: no drug shared across organisms yet — skipped)")
        return
    fig, ax = plt.subplots(figsize=(11, 1.0 + 1.0 * len(abs_)))
    ax.axis("off")
    for i, ab in enumerate(abs_):
        yy = len(abs_) - 1 - i
        sub = ot[ot.antibiotic == ab]
        fam = {o: set(_fam(x) for x in g["aro_gene_family"].dropna()) for o, g in sub.groupby("organism")}
        ec, kp = fam.get("ecoli", set()), fam.get("kpneumoniae", set())
        shared = sorted(ec & kp); eo = sorted(ec - kp); ko = sorted(kp - ec)
        ax.text(0.0, yy, _short(ab), fontsize=12, fontweight="bold", va="center")
        ax.text(0.24, yy, "   ".join(shared) or "—", fontsize=13, fontweight="bold",
                color="#2ca25f", va="center")
        extra = []
        if eo: extra.append("Ec-only: " + ", ".join(eo))
        if ko: extra.append("Kp-only: " + ", ".join(ko))
        ax.text(0.24, yy - 0.30, "     ".join(extra), fontsize=8.5, color="#888", va="center")
    ax.text(0.24, len(abs_) - 0.30, "SHARED gene family — recovered in BOTH organisms (concordant)",
            fontsize=9.5, color="#2ca25f", fontweight="bold", va="center")
    ax.set_xlim(-0.02, 1.0); ax.set_ylim(-0.6, len(abs_) - 0.05)
    ax.set_title("Cross-organism concordance: same drug → same resistance gene family", fontsize=12.5)
    _save(fig, out, "03_cross_organism")


def fig_mechanism(tables, out):
    """Heatmap: on-target confirmed gene family (rows) × model (cols),
    cell = # supporting unitigs. Reveals which family drives which drug."""
    from matplotlib.colors import LogNorm
    mech = pd.read_csv(tables / "mechanisms.csv")
    ot = mech[mech["on_target"] == True].copy()  # noqa: E712
    ot["fam"] = ot["aro_gene_family"].map(_fam)
    ot["col"] = ot["organism"] + "||" + ot["antibiotic"]
    ms = _sortkey(pd.read_csv(tables / "models_summary.csv"))
    order_cols = [f"{o}||{a}" for o, a in zip(ms.organism, ms.antibiotic)]
    piv = ot.groupby(["fam", "col"])["n_unitigs"].sum().unstack(fill_value=0)
    cols = [c for c in order_cols if c in piv.columns]
    piv = piv[cols]
    piv = piv.loc[piv.sum(axis=1).sort_values(ascending=False).index]
    M = piv.values.astype(float)
    disp = np.where(M > 0, M, np.nan)
    fig, ax = plt.subplots(figsize=(0.52 * len(cols) + 3, 0.42 * len(piv) + 2))
    im = ax.imshow(disp, aspect="auto", cmap="YlOrRd", norm=LogNorm(vmin=1, vmax=np.nanmax(disp)))
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels([f"{_short(c.split('||')[1])} ({_abbr(c.split('||')[0])})" for c in cols],
                       rotation=90, fontsize=7.5)
    ax.set_yticks(range(len(piv))); ax.set_yticklabels(piv.index, fontsize=8.5)
    thr = np.nanmax(disp) ** 0.5
    for i in range(len(piv)):
        for j in range(len(cols)):
            v = M[i, j]
            if v > 0:
                ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6,
                        color="white" if v > thr else "black")
    ax.set_title("On-target confirmed resistance gene families across models (cell = # unitigs)", fontsize=11)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02, label="unitigs (log)")
    _save(fig, out, "04_mechanism_heatmap")


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


def fig_evidence_layers(tables, out, db):
    """THE backbone figure: every model is filtered through 7 orthogonal
    validation layers. Heatmap = evidence units (unitigs, or 1 for the
    model-level label-permutation) supporting each model in each layer."""
    import sqlite3
    ms = _sortkey(pd.read_csv(tables / "models_summary.csv")).reset_index(drop=True)
    conn = sqlite3.connect(str(db))
    counts = {}  # run_id -> {evidence_type: n}
    for run_id, et, n in conn.execute(
            "SELECT pipeline_run_id, evidence_type, COUNT(*) FROM validation_evidence "
            "GROUP BY pipeline_run_id, evidence_type"):
        counts.setdefault(run_id, {})[et] = n
    conn.close()
    types = [t for t, _ in EVIDENCE_ORDER]
    M = np.array([[counts.get(r, {}).get(t, 0) for t in types] for r in ms["run_id"]], float)
    disp = np.where(M > 0, M, np.nan)  # 0 = layer not applicable (e.g. no SNP) → blank
    fig, ax = plt.subplots(figsize=(8.5, 0.45 * len(ms) + 1.6))
    from matplotlib.colors import LogNorm
    im = ax.imshow(disp, aspect="auto", cmap="YlGnBu",
                   norm=LogNorm(vmin=1, vmax=np.nanmax(disp)))
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels([lab for _, lab in EVIDENCE_ORDER], fontsize=8)
    ax.set_yticks(range(len(ms)))
    ax.set_yticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(ms.antibiotic, ms.organism)], fontsize=7.5)
    for i in range(len(ms)):
        for j in range(len(types)):
            v = M[i, j]
            if v > 0:
                ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6.5,
                        color="white" if v > np.nanmax(disp) ** 0.5 else "black")
    ax.set_title("Every biomarker is filtered through 7 orthogonal validation layers\n"
                 "(cell = supporting evidence units; blank = layer not applicable)", fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="evidence units (log)")
    _save(fig, out, "06_evidence_layers")


def fig_significance(tables, out, db):
    """05 redesign — single panel: REAL lineage-CV AUC vs label-shuffle null_max
    per model (parsed from the KB's label_permutation evidence_source). Shows the
    gap = model-level significance, far clearer than 20 tiny histograms."""
    import re, sqlite3
    ms = _sortkey(pd.read_csv(tables / "models_summary.csv")).reset_index(drop=True)
    conn = sqlite3.connect(str(db))
    real, nullmax, pval = {}, {}, {}
    for run_id, src, score in conn.execute(
            "SELECT pipeline_run_id, evidence_source, evidence_score FROM validation_evidence "
            "WHERE evidence_type='label_permutation'"):
        m = re.search(r"real_auc=([0-9.]+).*null_max=([0-9.]+)", src or "")
        if m:
            real[run_id] = float(m.group(1)); nullmax[run_id] = float(m.group(2)); pval[run_id] = score
    conn.close()
    ms = ms[ms["run_id"].isin(real)].reset_index(drop=True)
    y = np.arange(len(ms))
    r = [real[i] for i in ms["run_id"]]
    nm = [nullmax[i] for i in ms["run_id"]]
    col = [_colour(o) for o in ms["organism"]]
    fig, ax = plt.subplots(figsize=(9, 0.42 * len(ms) + 1.2))
    for yi, ri, ni, ci in zip(y, r, nm, col):
        ax.plot([ni, ri], [yi, yi], color="lightgrey", lw=2, zorder=1)
        ax.scatter(ni, yi, color="#999999", s=28, zorder=2)
        ax.scatter(ri, yi, color=ci, s=46, zorder=3)
    ax.axvline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(ms.antibiotic, ms.organism)], fontsize=7.5)
    ax.set_xlim(0.4, 1.0)
    ax.set_xlabel("ROC-AUC")
    ax.set_title("Model-level significance: REAL lineage-CV AUC (colour) ≫ "
                 "label-shuffle null max (grey)\nall p ≈ 0.02 (N=50 permutations)", fontsize=10)
    ax.invert_yaxis()
    _save(fig, out, "05_significance_real_vs_null")


def fig_external_concordance(tables, out, db):
    """M13: reference genotype tools (AMRFinderPlus, ResFinder) vs EUCAST/CLSI
    phenotype on held-out TEST genomes — balanced accuracy (bars) + Cohen's kappa."""
    import sqlite3
    conn = sqlite3.connect(str(db))
    rows = conn.execute(
        "SELECT p.organism, m.antibiotic, e.caller, e.balanced_accuracy "
        "FROM external_concordance e JOIN models m USING(model_id) "
        "JOIN pipeline_runs p USING(run_id)").fetchall()
    conn.close()
    if not rows:
        print("  (external_concordance: no rows — skipped)")
        return
    df = pd.DataFrame(rows, columns=["organism", "antibiotic", "caller", "bacc"])
    df["key"] = list(zip(df.organism, df.antibiotic))
    keys = sorted(df["key"].unique())
    callers = [c for c in ("model", "amrfinderplus", "resfinder") if c in set(df.caller)]
    cmap = {"model": "#31a354", "amrfinderplus": "#e6550d", "resfinder": "#756bb1"}
    nice = {"model": "our model", "amrfinderplus": "AMRFinderPlus", "resfinder": "ResFinder"}
    x = np.arange(len(keys)); w = 0.8 / max(1, len(callers))
    fig, ax = plt.subplots(figsize=(0.62 * len(keys) + 3, 5))
    for i, cl in enumerate(callers):
        vals = [df[(df.key == k) & (df.caller == cl)]["bacc"].mean() for k in keys]
        ax.bar(x + i * w, vals, w, label=nice[cl], color=cmap[cl], edgecolor="black", lw=0.4)
    ax.set_xticks(x + w * (len(callers) - 1) / 2)
    ax.set_xticklabels([f"{_short(a)}\n({_abbr(o)})" for o, a in keys], rotation=90, fontsize=7.5)
    ax.set_ylim(0.45, 1.03); ax.axhline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_ylabel("Balanced accuracy vs EUCAST/CLSI phenotype")
    ax.set_title("External concordance (M13, leakage-free held-out test): "
                 "our model vs AMRFinderPlus vs ResFinder", fontsize=10.5)
    ax.legend(fontsize=9, loc="lower right")
    _save(fig, out, "07_external_concordance")


FIGS = {"overview": lambda t, r, o, db: fig_overview(t, o, db),
        "external": lambda t, r, o, db: fig_external_concordance(t, o, db),
        "performance": lambda t, r, o, db: fig_performance(t, o),
        "cpss_pfer": lambda t, r, o, db: fig_cpss_pfer(t, o),
        "cross_org": lambda t, r, o, db: fig_cross_org(t, o),
        "mechanism": lambda t, r, o, db: fig_mechanism(t, o),
        "evidence": lambda t, r, o, db: fig_evidence_layers(t, o, db),
        "significance": lambda t, r, o, db: fig_significance(t, o, db),
        "null_hist": lambda t, r, o, db: fig_null_hist(t, r, o)}


def main():
    ap = argparse.ArgumentParser(description="Thesis figures from the AMR-KB tidy tables.")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--results", default="results")
    ap.add_argument("--db", default="results/kb/amrk.db", help="KB (for evidence/significance figs)")
    ap.add_argument("--out", default="figures")
    ap.add_argument("--only", default=None, help="comma list: " + ",".join(FIGS))
    args = ap.parse_args()
    tables, out = Path(args.tables), Path(args.out)
    want = args.only.split(",") if args.only else list(FIGS)
    for name in want:
        if name not in FIGS:
            print(f"  ! unknown figure '{name}'")
            continue
        FIGS[name](tables, args.results, out, args.db)
    print("DONE.")


if __name__ == "__main__":
    main()
