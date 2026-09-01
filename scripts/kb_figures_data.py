#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thesis figures, part 1 — data → QC → population structure → features.

Everything upstream of the model: what went in, what quality control removed, how
each organism's population is structured, and what the feature space looks like.
kb_figures.py covers the model/KB end; this covers the beginning of the pipeline.

Usage:
    python scripts/kb_figures_data.py --results results --data data/processed \
        --tables results/tables --db results/kb/amrk.db --out results/figures \
        [--only qc,contiguity,passrate,composition,balance,lineage,clonality,structure,resistance,features,lengths]

Figures (results/figures/):
    10_genome_qc_scatter          CheckM2 completeness x contamination, per organism
    11_assembly_contiguity        QUAST N50 x contig count
    12_qc_pass_rates              pass rate + failure reasons
    13_dataset_composition        per-model R/S balance
    14_balance_vs_auc             minority fraction vs lineage-CV AUC
    15_lineage_size_distribution  PopPUNK lineage rank-size curves
    16_lineage_resistance         resistance rate inside the largest lineages
    17_clonality_vs_inflation     clonal dominance vs random-CV inflation
    18_feature_counts             unitigs retained per model
    19_unitig_lengths             unitig length distribution
"""
import argparse
import json
import sqlite3
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


def _qc_dir(results, org):
    return Path(results) / org / "global_exploration" / "genome_qc"


def _grid(n, ncols=3, w=4.3, h=3.4):
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(w * ncols, h * nrows))
    axes = np.atleast_1d(axes).ravel()
    for a in axes[n:]:
        a.axis("off")
    return fig, axes


# ---------------------------------------------------------------- A1 / A2 / A3
def fig_qc_scatter(results, orgs, out):
    """CheckM2 completeness x contamination — the gate that actually excludes."""
    fig, axes = _grid(len(orgs))
    for ax, org in zip(axes, orgs):
        f = _qc_dir(results, org) / "checkm2" / "quality_report.tsv"
        if not f.exists():
            ax.set_title(f"{_abbr(org)} — no CheckM2"); ax.axis("off"); continue
        df = pd.read_csv(f, sep="\t")
        ok = (df["Completeness"] >= 95) & (df["Contamination"] <= 5)
        ax.scatter(df.loc[ok, "Completeness"], df.loc[ok, "Contamination"], s=6,
                   alpha=0.45, color=_colour(org), label=f"pass ({int(ok.sum())})")
        # Black, not red: K. pneumoniae's own palette colour IS red, so red crosses were
        # indistinguishable from its passing genomes in that panel.
        ax.scatter(df.loc[~ok, "Completeness"], df.loc[~ok, "Contamination"], s=16,
                   alpha=0.9, color="black", marker="x", linewidths=0.9,
                   label=f"fail ({int((~ok).sum())})")
        ax.axvline(95, ls="--", c="grey", lw=0.8)
        ax.axhline(5, ls="--", c="grey", lw=0.8)
        ax.set_yscale("symlog", linthresh=1)
        # Contamination is a percentage >= 0; symlog's default view drew a negative
        # decade under the axis, which is not a value the metric can take.
        ax.set_ylim(bottom=0)
        ax.set_xlim(0, 102)
        ax.set_xlabel("completeness (%)"); ax.set_ylabel("contamination (%)")
        ax.set_title(_display(org), fontsize=10, style="italic")
        ax.legend(fontsize=7, loc="upper left", frameon=False)
    fig.suptitle("Genome QC — CheckM2 gate (completeness ≥95 %, contamination ≤5 %)", fontsize=12)
    fig.tight_layout()
    _save(fig, out, "10_genome_qc_scatter")


def fig_contiguity(results, orgs, out):
    """QUAST N50 x contigs — reported, but ADVISORY: fragmented assemblies keep
    their AMR content, so gating on contiguity would have discarded ~63 % of
    E. faecium for no biological reason."""
    fig, axes = _grid(len(orgs))
    for ax, org in zip(axes, orgs):
        f = _qc_dir(results, org) / "quast" / "transposed_report.tsv"
        if not f.exists():
            ax.set_title(f"{_abbr(org)} — no QUAST"); ax.axis("off"); continue
        df = pd.read_csv(f, sep="\t")
        n50 = next((c for c in df.columns if c.strip() == "N50"), None)
        nct = next((c for c in df.columns if c.strip() == "# contigs"), None)
        if not n50 or not nct:
            ax.set_title(f"{_abbr(org)} — columns?"); ax.axis("off"); continue
        ax.scatter(df[nct], df[n50], s=7, alpha=0.45, color=_colour(org))
        ax.axhline(50000, ls="--", c="grey", lw=0.8)
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("# contigs"); ax.set_ylabel("N50 (bp)")
        ax.set_title(f"{_display(org)}  (n={len(df)})", fontsize=10, style="italic")
    fig.suptitle("Assembly contiguity (QUAST) — advisory, NOT an exclusion gate", fontsize=12)
    fig.tight_layout()
    _save(fig, out, "11_assembly_contiguity")


def fig_pass_rates(results, orgs, out):
    """Pass rate + why genomes failed."""
    rows = []
    for org in orgs:
        f = _qc_dir(results, org) / f"02d_genome_qc_summary_{org}.json"
        if not f.exists():
            continue
        s = json.loads(f.read_text(encoding="utf-8"))
        rows.append({"organism": org, "n": s.get("n_genomes"), "pass": s.get("n_pass"),
                     "fail": s.get("n_fail"), "rate": s.get("pass_rate"),
                     "completeness": s.get("n_fail_completeness"),
                     "contamination": s.get("n_fail_contamination"),
                     "n50": s.get("n_fail_n50"),
                     "contigs": s.get("n_fail_contigs")})
    if not rows:
        print("  (pass-rate: no 02d summaries — skipped)"); return
    df = pd.DataFrame(rows)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.2))
    x = np.arange(len(df))
    a1.bar(x, 100 * df["rate"], color=[_colour(o) for o in df.organism],
           edgecolor="k", lw=0.4)
    for xi, r, n in zip(x, df["rate"], df["n"]):
        a1.text(xi, 100 * r + 0.4, f"{100*r:.1f}%\n(n={n})", ha="center", fontsize=8)
    a1.set_xticks(x); a1.set_xticklabels([_abbr(o) for o in df.organism])
    a1.set_ylim(80, 103); a1.set_ylabel("genomes passing QC (%)")
    tot_p, tot_n = int(df["pass"].sum()), int(df["pass"].sum() + df["fail"].sum())
    a1.set_title("CheckM2 pass rate — the ENFORCED gate\n"
                 f"{tot_p:,}/{tot_n:,} = {100*tot_p/tot_n:.1f}% pass "
                 "(note the truncated y-axis)", fontsize=9.5)
    # Four criteria are COMPUTED; only two are ENFORCED. The assembly-contiguity pair
    # was evaluated and then deliberately not applied, because an N50>=50 kb gate
    # removed 63% of the E. faecium collection (1,305 of 2,078) — short-contig draft
    # assemblies are the norm for that species in BV-BRC, so the gate was selecting on
    # assembly provenance rather than on genome quality. Plotting only the two enforced
    # reasons, with no trace of the other two, let the figure imply that contiguity was
    # never a question. It has to be visible, and visibly not applied.
    w = 0.2
    a2.bar(x - 1.5 * w, df["completeness"], width=w, color="#2c7fb8",
           edgecolor="k", lw=0.4, label="completeness <95 % (enforced)")
    a2.bar(x - 0.5 * w, df["contamination"], width=w, color="#d62728",
           edgecolor="k", lw=0.4, label="contamination >5 % (enforced)")
    a2.bar(x + 0.5 * w, df["n50"], width=w, color="#bdbdbd", edgecolor="k", lw=0.4,
           hatch="//", label="N50 <50 kb (computed, NOT enforced)")
    a2.bar(x + 1.5 * w, df["contigs"], width=w, color="#f0f0f0", edgecolor="k", lw=0.4,
           hatch="//", label="contigs >500 (computed, NOT enforced)")
    a2.set_yscale("symlog", linthresh=10)
    a2.set_xticks(x); a2.set_xticklabels([_abbr(o) for o in df.organism])
    a2.set_ylabel("genomes failing the criterion (symlog)")
    a2.legend(fontsize=7, frameon=False, ncol=2)
    worst = df.loc[df["n50"].idxmax()]
    a2.set_title("Why genomes were excluded — and what was measured but not applied\n"
                 f"an N50 gate would have removed {int(worst['n50']):,} of "
                 f"{int(worst['n'])} {_abbr(worst['organism'])} genomes "
                 f"({100*worst['n50']/worst['n']:.0f}%)", fontsize=9)
    fig.tight_layout()
    _save(fig, out, "12_qc_pass_rates")


# ------------------------------------------------------------------- A4 / A5
def _label_counts(data, org, ab):
    f = Path(data) / org / ab / "matrix_unitig" / f"y_{ab}.csv"
    if not f.exists():
        return None
    y = pd.read_csv(f)["label"].to_numpy()
    return int((y == 1).sum()), int((y == 0).sum())


def fig_composition(data, ms, out):
    """Per-model resistant/susceptible composition — the class balance every AUC
    has to be read against."""
    df = _sortkey(ms).reset_index(drop=True)
    counts = [_label_counts(data, r.organism, r.antibiotic) for r in df.itertuples()]
    keep = [i for i, c in enumerate(counts) if c]
    if not keep:
        print("  (composition: no y_*.csv found — skipped)"); return
    df = df.iloc[keep].reset_index(drop=True)
    R = np.array([counts[i][0] for i in keep]); S = np.array([counts[i][1] for i in keep])
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.bar(x, R, color="#d62728", edgecolor="k", lw=0.3, label="resistant")
    ax.bar(x, S, bottom=R, color="#7fbf7f", edgecolor="k", lw=0.3, label="susceptible")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(df.antibiotic, df.organism)],
                       rotation=90, fontsize=7)
    ax.set_ylabel("genomes with a phenotype")
    ax.set_title("Dataset composition per model (QC-passed, lineage-covered genomes)")
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    _save(fig, out, "13_dataset_composition")


def fig_balance_vs_auc(data, ms, out):
    """Does class imbalance explain performance? (It mostly does not — which is
    why the low scores need the population-structure explanation instead.)"""
    rows = []
    for r in ms.itertuples():
        c = _label_counts(data, r.organism, r.antibiotic)
        if not c:
            continue
        R, S = c
        rows.append({"organism": r.organism, "antibiotic": r.antibiotic,
                     "minority_frac": min(R, S) / (R + S),
                     "n": R + S, "auc": r.lineage_cv_auc})
    if not rows:
        print("  (balance: no labels — skipped)"); return
    df = pd.DataFrame(rows)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    for org, g in df.groupby("organism"):
        a1.scatter(g.minority_frac, g.auc, s=42, color=_colour(org),
                   edgecolor="k", lw=0.3, label=_display(org))
        a2.scatter(g.n, g.auc, s=42, color=_colour(org), edgecolor="k", lw=0.3)
    r1 = np.corrcoef(df.minority_frac, df.auc)[0, 1]
    r2 = np.corrcoef(df.n, df.auc)[0, 1]
    a1.set_xlabel("minority-class fraction"); a1.set_ylabel("lineage-CV ROC-AUC")
    a1.set_title(f"Class balance vs performance (r = {r1:.2f})", fontsize=10)
    a2.set_xscale("log"); a2.set_xlabel("genomes in the model")
    a2.set_title(f"Sample size vs performance (r = {r2:.2f})", fontsize=10)
    a1.legend(fontsize=7.5, frameon=False, loc="lower right")
    fig.tight_layout()
    _save(fig, out, "14_balance_vs_auc")


# -------------------------------------------------------------- B1 / B3 / B2
def _clusters(data, org):
    f = Path(data) / org / "lineage" / "poppunk_clusters.csv"
    if not f.exists():
        return None
    return pd.read_csv(f, dtype={"Genome ID": str})


def fig_lineage_sizes(data, orgs, out):
    """Rank-size curves: how much of each organism sits in its biggest lineage.
    This is the property lineage-aware CV exists to respect."""
    # Bail out BEFORE creating the figure when no organism has cluster data. Skipping
    # each organism individually and saving anyway wrote a fully empty pair of axes over
    # a good figure — this file lives in data/processed/, which is routinely pruned from
    # the laptop, so "regenerate everything" silently destroyed the artefact instead of
    # leaving it alone. Every other generator here skips; this one has to as well.
    have = [(org, _clusters(data, org)) for org in orgs]
    have = [(org, cl) for org, cl in have if cl is not None]
    if not have:
        print("  (lineage sizes: no PopPUNK cluster files under "
              f"{data} — skipped, existing figure left untouched)")
        return
    fig, ax = plt.subplots(figsize=(8.5, 5))
    txt = []
    for org, cl in have:
        sizes = cl["Cluster"].value_counts().to_numpy()
        frac = 100 * sizes / sizes.sum()
        ax.plot(np.arange(1, len(frac) + 1), frac, marker="o", ms=2.5, lw=1.2,
                color=_colour(org), label=f"{_display(org)}  ({len(sizes)} lineages)")
        txt.append((org, frac[0], len(sizes), int((sizes == 1).sum())))
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("lineage rank"); ax.set_ylabel("% of the organism's genomes")
    ax.set_title("PopPUNK lineage size distribution (log-log)")
    ax.legend(fontsize=8, frameon=False)
    note = "  ·  ".join(f"{_abbr(o)} largest {f:.1f}%" for o, f, _, _ in txt)
    ax.text(0.5, -0.16, note, transform=ax.transAxes, ha="center", fontsize=8, color="#555")
    fig.tight_layout()
    _save(fig, out, "15_lineage_size_distribution")


def fig_lineage_resistance(data, ms, orgs, out):
    """Resistance rate INSIDE the largest lineages — clonal confounding, seen directly.
    A lineage that is ~all-resistant lets a model score by recognising the clone;
    holding that lineage out is precisely what lineage-aware CV does."""
    # Same guard as fig_lineage_sizes: an all-blank grid is worse than no output,
    # because _save overwrites whatever was there.
    # The per-panel branch below blanks a panel when the label/genome CSVs are
    # absent, and _save then overwrites a good figure with an empty grid — the
    # exact failure this guard exists to prevent, so it has to test the label
    # files too, not just the cluster assignments.
    def _panel_has_data(org):
        if _clusters(data, org) is None:
            return False
        sub = ms[ms.organism == org]
        if sub.empty:
            return False
        ab = sub.sort_values(["n_genomes", "antibiotic"],
                             ascending=[False, True]).iloc[0]["antibiotic"]
        gdir = Path(data) / org / ab / "matrix_unitig"
        return (gdir / f"genomes_{ab}.csv").exists() and (gdir / f"y_{ab}.csv").exists()

    if not any(_panel_has_data(org) for org in orgs):
        print("  (lineage resistance: no genomes_*/y_* CSVs — skipped, existing figure kept)")
        return
    if all(_clusters(data, org) is None for org in orgs):
        print("  (lineage resistance: no PopPUNK cluster files under "
              f"{data} — skipped, existing figure left untouched)")
        return
    fig, axes = _grid(len(orgs), ncols=3, w=4.6, h=3.4)
    for ax, org in zip(axes, orgs):
        cl = _clusters(data, org)
        sub = ms[ms.organism == org]
        if cl is None or sub.empty:
            ax.axis("off"); continue
        # the organism's largest model = the most representative phenotype
        # Deterministic pick when two antibiotics have the same genome count.
        ab = sub.sort_values(["n_genomes", "antibiotic"],
                             ascending=[False, True]).iloc[0]["antibiotic"]
        gdir = Path(data) / org / ab / "matrix_unitig"
        gf, yf = gdir / f"genomes_{ab}.csv", gdir / f"y_{ab}.csv"
        if not (gf.exists() and yf.exists()):
            ax.axis("off"); continue
        g = pd.read_csv(gf, dtype=str)["Genome ID"].astype(str)
        y = pd.read_csv(yf)["label"].to_numpy()
        d = pd.DataFrame({"Genome ID": g, "y": y}).merge(
            cl.astype({"Genome ID": str}), on="Genome ID", how="inner")
        top = d["Cluster"].value_counts().head(10)
        rate = [100 * d.loc[d.Cluster == c, "y"].mean() for c in top.index]
        xs = np.arange(len(top))
        ax.bar(xs, rate, color=_colour(org), edgecolor="k", lw=0.4)
        ax.axhline(100 * d["y"].mean(), ls="--", c="grey", lw=1,
                   label=f"organism-wide {100*d['y'].mean():.0f}%")
        for xi, c in zip(xs, top.values):
            ax.text(xi, 2, f"n={c}", ha="center", fontsize=6.5, rotation=90, color="white")
        ax.set_xticks(xs); ax.set_xticklabels([str(c) for c in top.index], fontsize=7)
        ax.set_ylim(0, 105); ax.set_ylabel("% resistant")
        ax.set_xlabel("lineage (10 largest)")
        ax.set_title(f"{_display(org)} — {_short(ab)}", fontsize=9.5, style="italic")
        ax.legend(fontsize=7, frameon=False)
    fig.suptitle("Resistance is clone-structured: R-rate inside the largest lineages", fontsize=12)
    fig.tight_layout()
    _save(fig, out, "16_lineage_resistance")


def fig_clonality_vs_inflation(data, tables, orgs, out):
    """Clonal dominance vs how much a lineage-blind CV inflates the AUC."""
    cvf = Path(tables) / "cv_comparison.csv"
    if not cvf.exists():
        print("  (clonality: run kb_cv_comparison.py first — skipped)"); return
    cv = pd.read_csv(cvf).groupby("organism")["inflation"].mean()
    rows = []
    for org in orgs:
        cl = _clusters(data, org)
        if cl is None or org not in cv.index:
            continue
        sizes = cl["Cluster"].value_counts()
        rows.append({"organism": org, "largest_pct": 100 * sizes.iloc[0] / sizes.sum(),
                     "inflation": cv[org]})
    if len(rows) < 3:
        print("  (clonality: too few organisms — skipped)"); return
    df = pd.DataFrame(rows)
    from scipy.stats import pearsonr, spearmanr
    r, p = pearsonr(df.largest_pct, df.inflation)
    rho, ps = spearmanr(df.largest_pct, df.inflation)
    fig, ax = plt.subplots(figsize=(7.2, 5))
    for t in df.itertuples():
        ax.scatter(t.largest_pct, t.inflation, s=110, color=_colour(t.organism),
                   edgecolor="k", lw=0.5, zorder=3)
        ax.annotate(_display(t.organism), (t.largest_pct, t.inflation),
                    textcoords="offset points", xytext=(9, -3), fontsize=8.5, style="italic")
    z = np.polyfit(df.largest_pct, df.inflation, 1)
    xs = np.linspace(df.largest_pct.min() - 2, df.largest_pct.max() + 4, 50)
    ax.plot(xs, np.polyval(z, xs), ls="--", c="grey", lw=1, zorder=1)
    ax.margins(x=0.16)          # the right-most label ran into the axis edge
    ax.set_xlabel("largest lineage (% of the organism's genomes)")
    ax.set_ylabel("mean AUC inflation when the lineage grouping is removed")
    ax.set_title("The more clonal the organism, the more a random split flatters it\n"
                 f"Pearson r = {r:.3f} (p = {p:.3f}) · Spearman ρ = {rho:.3f} (p = {ps:.3f}) · n = {len(df)}",
                 fontsize=10.5)
    ax.text(0.02, 0.02, "n = 6 organisms: treat as a trend, not an estimate",
            transform=ax.transAxes, fontsize=8, color="#777")
    fig.tight_layout()
    _save(fig, out, "17_clonality_vs_inflation")


# ----------------------------------------------------------------- C1 / C2
MEASURE_LABELS = [
    ("clonality_pct",        "largest lineage (% of genomes)"),
    ("simpson_diversity",    "Simpson diversity of lineages"),
    ("shannon_diversity",    "Shannon diversity of lineages"),
    ("n_lineages",           "number of PopPUNK lineages"),
    ("n_singleton_lineages", "singleton lineages"),
]


def fig_structure_vs_inflation(tables, out):
    """All five structure measures against CV inflation -- deliberately all five.

    Figure 17 plots one measure, largest-lineage share, and reports r and rho beside it.
    That is the measure the production run happened to pick, and it is significant under
    Pearson (r +0.914, p 0.011) but NOT under Spearman (rho +0.771, p 0.072). Simpson
    diversity is significant under both (r -0.944 p 0.005, rho -0.829 p 0.042).

    Publishing Simpson alone would be selection on the outcome: five measures, six
    organisms, no pre-registration. So this figure shows the whole set, marks which
    survive a rank test, and says the count out loud in the title. A reader can then see
    that the direction is consistent across every measure -- more clonal, more inflation
    -- while no single p-value here carries weight on its own.

    Reads lineage_summary.csv / lineage_summary_stats.json so the figure and the table
    can never disagree; kb_tables_thesis.py computes both.
    """
    tf, sf = Path(tables) / "lineage_summary.csv", Path(tables) / "lineage_summary_stats.json"
    if not tf.exists():
        print("  (structure: lineage_summary.csv missing — run kb_tables_thesis.py — skipped)")
        return
    d = pd.read_csv(tf)
    stats = json.loads(sf.read_text()) if sf.exists() else {"correlations": {}}
    cor = stats.get("correlations", {})

    fig, axes = _grid(6, ncols=3, w=4.5, h=3.6)
    n_sig = 0
    for ax, (col, lab) in zip(axes, MEASURE_LABELS):
        if col not in d.columns:
            ax.axis("off"); continue
        x, y = d[col].to_numpy(float), d["mean_inflation"].to_numpy(float)
        for _, t in d.iterrows():
            ax.scatter(t[col], t["mean_inflation"], s=95, zorder=3,
                       color=_colour(t["organism"]), edgecolor="k", lw=0.5)
            ax.annotate(_abbr(t["organism"]), (t[col], t["mean_inflation"]),
                        textcoords="offset points", xytext=(6, 4), fontsize=8)
        if len(x) > 2:
            ax.plot(np.sort(x), np.polyval(np.polyfit(x, y, 1), np.sort(x)),
                    ls="--", c="grey", lw=1, zorder=1)
        c = cor.get(col)
        if c:
            sig = c["spearman_p"] < 0.05
            n_sig += bool(sig)
            ax.set_title(
                f"r = {c['pearson_r']:+.3f} (p {c['pearson_p']:.3f})\n"
                f"rho = {c['spearman_rho']:+.3f} (p {c['spearman_p']:.3f})"
                + ("  \u2713 rank test" if sig else ""),
                fontsize=9, color="#238b45" if sig else "#444444")
        ax.set_xlabel(lab, fontsize=9)
        ax.set_ylabel("mean AUC inflation", fontsize=9)
        ax.margins(x=0.18, y=0.22)
        ax.tick_params(labelsize=8)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)

    ax = axes[5]
    ax.axis("off")
    ax.text(0.0, 0.98,
            "Why all five are shown\n\n"
            f"n = {len(d)} organisms. Five measures were tested against the same target "
            "without pre-registration, so the best p-value among them is not evidence.\n\n"
            "Only Simpson diversity survives a rank test; largest-lineage share, the "
            "measure figure 17 plots, does not (p 0.072).\n\n"
            "What is robust is the direction, which every measure agrees on: the less "
            "diverse the population, the more a lineage-blind split flatters the model.\n\n"
            "Report this as a trend, not an estimate.",
            transform=ax.transAxes, fontsize=8.6, va="top", ha="left", wrap=True,
            bbox=dict(boxstyle="round,pad=0.6", fc="#f7f7f7", ec="#cccccc", lw=0.8))

    fig.suptitle(
        "Population structure vs cross-validation inflation — five measures, all reported\n"
        f"{n_sig} of {len(MEASURE_LABELS)} reach p < 0.05 under Spearman; with n = {len(d)} "
        "and five measures tried, treat the direction as the finding, not the coefficients",
        fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.90])
    _save(fig, out, "40_structure_vs_inflation")


def fig_feature_counts(db, ms, out):
    """Unitigs retained per model (post support-filter) — the feature space size."""
    if not Path(db).exists():
        print("  (features: no KB — skipped)"); return
    conn = sqlite3.connect(db)
    nf = dict(conn.execute(
        "SELECT m.antibiotic || '|' || p.organism, m.n_features FROM models m "
        "JOIN pipeline_runs p ON p.run_id = m.run_id").fetchall())
    conn.close()
    df = _sortkey(ms).reset_index(drop=True)
    vals = [nf.get(f"{r.antibiotic}|{r.organism}", np.nan) for r in df.itertuples()]
    if all(v != v for v in vals):
        print("  (features: n_features empty — skipped)"); return
    x = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(13, 4.6))
    ax.bar(x, np.array(vals) / 1e6, color=[_colour(o) for o in df.organism],
           edgecolor="k", lw=0.35)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(df.antibiotic, df.organism)],
                       rotation=90, fontsize=7)
    ax.set_ylabel("unitig features (millions)")
    ax.set_title("Feature-space size per model — unitigs surviving the support filter")
    fig.tight_layout()
    _save(fig, out, "18_feature_counts")


def fig_unitig_lengths(data, ms, out, sample=40000, stride=25):
    """Unitig length distribution — why 'blastn-short' is the right BLAST task.

    Prefers the full features.txt. Where only ``features_sample.txt`` is present
    — the systematic 1-in-25 extract taken on the HPC because the full feature
    files total ~3.8 GB — that file is read with stride 1, since it has already
    been strided. The sampling is stated on the figure so the two cases are not
    confused.
    """
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    drawn, all_lens, from_sample = 0, [], False
    for org, g in ms.groupby("organism"):
        lens = []
        for r in g.itertuples():
            mdir = Path(data) / r.organism / r.antibiotic / "matrix_unitig"
            f, step = mdir / "features.txt", stride
            if not f.exists():
                f, step = mdir / "features_sample.txt", 1   # already 1-in-25
                if not f.exists():
                    continue
                from_sample = True
            # Stride through the file instead of reading its head: features.txt is
            # written in matrix-column order, so the first N lines are one contiguous
            # slice of the feature space, not a sample of it.
            with open(f, encoding="utf-8", errors="replace") as fh:
                for i, line in enumerate(fh):
                    if i % step:
                        continue
                    lens.append(len(line.split("\t")[0]))
                    if len(lens) >= sample:
                        break
            break                                   # one model per organism is enough
        if not lens:
            continue
        all_lens += lens
        ax.hist(lens, bins=60, range=(20, 120), histtype="step", lw=1.5,
                density=True, color=_colour(org), label=f"{_display(org)} (n={len(lens):,})")
        drawn += 1
    if not drawn:
        plt.close(fig); print("  (lengths: no features.txt/features_sample.txt — skipped)"); return
    hi = float(np.percentile(all_lens, 99.5)) if all_lens else 120
    ax.set_xlim(min(all_lens) - 2, max(35, hi) + 5)
    ax.set_xlabel("unitig length (bp)"); ax.set_ylabel("density")
    src = "systematic 1-in-25 sample" if from_sample else "1-in-%d sample" % stride
    ax.set_title(f"Unitig length distribution ({src})\n"
                 "short unitigs are why BLAST runs in 'blastn-short' mode", fontsize=10.5)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    _save(fig, out, "19_unitig_lengths")


def main():
    ap = argparse.ArgumentParser(description="Thesis figures — data, QC, population, features.")
    ap.add_argument("--results", default="results")
    ap.add_argument("--data", default="data/processed")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--out", default="results/figures")
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    ms = pd.read_csv(Path(args.tables) / "models_summary.csv")
    orgs = _organisms(ms)
    out = Path(args.out)
    only = {s.strip() for s in args.only.split(",") if s.strip()}

    todo = [
        ("qc",          lambda: fig_qc_scatter(args.results, orgs, out)),
        ("contiguity",  lambda: fig_contiguity(args.results, orgs, out)),
        ("passrate",    lambda: fig_pass_rates(args.results, orgs, out)),
        ("composition", lambda: fig_composition(args.data, ms, out)),
        ("balance",     lambda: fig_balance_vs_auc(args.data, ms, out)),
        ("lineage",     lambda: fig_lineage_sizes(args.data, orgs, out)),
        ("resistance",  lambda: fig_lineage_resistance(args.data, ms, orgs, out)),
        ("clonality",   lambda: fig_clonality_vs_inflation(args.data, args.tables, orgs, out)),
        ("structure",   lambda: fig_structure_vs_inflation(args.tables, out)),
        ("features",    lambda: fig_feature_counts(args.db, ms, out)),
        ("lengths",     lambda: fig_unitig_lengths(args.data, ms, out)),
    ]
    for name, fn in todo:
        if only and name not in only:
            continue
        try:
            fn()
        except Exception as e:                       # one bad figure must not kill the set
            print(f"  ✗ {name} failed: {type(e).__name__}: {e}")
    print("DONE.")


if __name__ == "__main__":
    main()
