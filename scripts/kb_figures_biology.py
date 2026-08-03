#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thesis figures, part 3 — biology and the knowledge base itself.

Everything here reads results/tables/biomarkers.csv (one row per (model, unitig),
already de-duplicated by kb_tables.py) plus the KB for the tier counts.

Usage:
    python scripts/kb_figures_biology.py --tables results/tables \
        --db results/kb/amrk.db --out results/figures [--only novel,blast,prevalence,funnel]

Figures (results/figures/):
    32_novel_candidates    the strong_novel biomarkers, the KB's own contribution
    33_blast_identity      identity x coverage, with the confidence-tier cutoffs
    34_prevalence          resistant vs susceptible prevalence per biomarker
    35_evidence_funnel     45 models -> biomarkers -> tiers -> strong_novel
"""
import argparse
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

from kb_figures import _abbr, _colour, _save, _short  # noqa: E402

TIER_COLOURS = {"confirmed": "#238b45", "strong_novel": "#d62728",
                "candidate": "#fdae61", "weak": "#9ecae1", "none": "#cccccc"}
TIER_ORDER = ["confirmed", "strong_novel", "candidate", "weak", "none"]


def fig_novel(bio, out):
    """The strong_novel set: CPSS-stable AND pyseer-significant AND no CARD hit.
    These are the KB's actual contribution — biomarkers a database lookup misses."""
    nov = bio[bio.evidence_tier == "strong_novel"].copy()
    if nov.empty:
        print("  (novel: no strong_novel rows — skipped)"); return
    nov["neglogp"] = -np.log10(pd.to_numeric(nov.get("pyseer_lrt_p"), errors="coerce")
                               .replace(0, np.nan))
    # Effect size, in order of preference. delta_prevalence is the meaningful one and
    # IS available for all 23 (they sit in both the gain and CPSS paths) — it only
    # looked empty because step 10 never emitted the column, so the KB stored NULL;
    # that is fixed at the source and derived in kb_tables. The fallbacks stay for KBs
    # built before the fix, and a value that is genuinely missing is skipped, never
    # drawn as a zero-length bar.
    strength, s_label = None, ""
    for col, lab in (("delta_prevalence", "prevalence(R) − prevalence(S)"),
                     ("mean_abs_shap", "mean |TreeSHAP|"),
                     ("composite_score", "composite score"),
                     ("gain", "XGBoost gain")):
        v = pd.to_numeric(nov.get(col), errors="coerce")
        if v is not None and v.notna().any():
            strength, s_label = v, lab          # keep NaN as NaN — see below
            break
    if strength is None:
        strength, s_label = pd.Series(np.nan, index=nov.index), "(no effect-size column)"
    n_missing = int(strength.isna().sum())
    # Size on MAGNITUDE: delta_prevalence is signed (two biomarkers here are enriched in
    # susceptible genomes, i.e. negative), and scaling a marker size by a negative number
    # yields a negative size — silently invalid.
    smax = float(strength.abs().max()) if strength.notna().any() else 1.0
    smax = smax or 1.0

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.4),
                                 gridspec_kw={"width_ratios": [1.15, 1]})
    sizes = 40 + 260 * (strength.abs() / smax).fillna(0)   # NaN -> smallest dot, not zero-value
    for org, g in nov.groupby("organism"):
        a1.scatter(g["selection_frequency"], g["neglogp"], s=sizes.loc[g.index],
                   color=_colour(org), edgecolor="k", lw=0.4, alpha=0.85,
                   label=f"{_abbr(org)} ({len(g)})")
    a1.axvline(0.6, ls="--", c="grey", lw=0.8)
    a1.set_xlabel("CPSS selection frequency"); a1.set_ylabel("pyseer LMM  −log10(p)")
    a1.set_title(f"{len(nov)} strong_novel biomarkers\n"
                 f"stable + LMM-significant + no CARD hit (dot size = |{s_label}|)",
                 fontsize=10.5)
    a1.legend(fontsize=8, frameon=False, title="organism", title_fontsize=8)

    # Bars only for biomarkers whose effect size is actually defined. Plotting a NaN as
    # a zero-length bar would again show "not measured" as "no effect" — the same trap
    # delta_prevalence set earlier in this figure.
    lab_all = [f"{_short(r.antibiotic)} ({_abbr(r.organism)})" for r in nov.itertuples()]
    have = np.where(strength.notna().to_numpy())[0]
    order = have[np.argsort(strength.to_numpy()[have])]
    y = np.arange(len(order))
    a2.barh(y, strength.to_numpy()[order],
            color=[_colour(nov.iloc[i]["organism"]) for i in order], edgecolor="k", lw=0.35)
    a2.set_yticks(y); a2.set_yticklabels([lab_all[i] for i in order], fontsize=7)
    a2.set_xlabel(s_label)
    note = f" · {n_missing} of {len(nov)} without a value (not shown)" if n_missing else ""
    a2.set_title(f"Discriminative gap of each novel biomarker{note}\n"
                 f"{s_label}", fontsize=9.5)
    fig.tight_layout()
    _save(fig, out, "32_novel_candidates")


def fig_blast(bio, out):
    """Where the confidence tiers actually fall in identity/coverage space."""
    d = bio.dropna(subset=["identity_pct", "coverage"]).copy()
    if d.empty:
        print("  (blast: no identity/coverage — skipped)"); return
    d["cov_pct"] = 100 * pd.to_numeric(d["coverage"], errors="coerce")
    src = d["source_db"].astype(str).str.lower() if "source_db" in d.columns \
        else pd.Series("card", index=d.index)
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    # CARD and NCBI hits share this space but only the CARD pass feeds the tier cutoffs
    # — a strong_novel biomarker has NO CARD hit, so the point plotted for it is its NCBI
    # hit. Drawing both with one marker implied all of them were graded against CARD.
    for tier in TIER_ORDER:
        for source, marker, alpha in (("card", "o", 0.55), ("ncbi", "x", 0.7)):
            g = d[(d.evidence_tier == tier) & (src == source)]
            if g.empty:
                continue
            ax.scatter(g["identity_pct"], g["cov_pct"], s=18 if marker == "o" else 22,
                       alpha=alpha, marker=marker, linewidths=0.7,
                       color=TIER_COLOURS.get(tier, "#888888"),
                       label=f"{tier} · {source.upper()} ({len(g)})",
                       edgecolor="none" if marker == "o" else None)
    ax.axvline(95, ls="--", c="grey", lw=0.8); ax.axhline(95, ls="--", c="grey", lw=0.8)
    ax.axvline(90, ls=":", c="grey", lw=0.7); ax.axhline(80, ls=":", c="grey", lw=0.7)
    ax.set_xlabel("BLAST identity (%)"); ax.set_ylabel("query coverage (%)")
    lo = float(np.nanmin(d["identity_pct"]))
    ax.set_xlim(max(60, lo - 2), 101)
    ax.set_title(f"Biomarker BLAST hits and the confidence-tier cutoffs "
                 f"({len(d):,} of {len(bio):,} biomarkers have a hit)\n"
                 "dashed = confirmed 95/95, dotted = candidate 90/80 — "
                 "the cutoffs grade the CARD pass (circles); NCBI hits (×) give context\n"
                 "strong_novel points are CARD ALIGNMENTS that did not qualify as gene "
                 "hits — that is what makes them novel", fontsize=9)
    ax.legend(fontsize=7, frameon=False, loc="lower left", ncol=2)
    fig.tight_layout()
    _save(fig, out, "33_blast_identity")


def fig_prevalence(bio, out):
    """Prevalence in resistant vs susceptible genomes. Distance from the diagonal is
    the discriminative signal; tier colour shows whether CARD already knew about it."""
    d = bio.dropna(subset=["prevalence_resistant", "prevalence_susceptible"]).copy()
    if d.empty:
        print("  (prevalence: no background-frequency rows — skipped)"); return
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    for tier in TIER_ORDER:
        g = d[d.evidence_tier == tier]
        if g.empty:
            continue
        ax.scatter(g["prevalence_susceptible"], g["prevalence_resistant"], s=16,
                   alpha=0.5, color=TIER_COLOURS.get(tier, "#888888"),
                   label=f"{tier} ({len(g)})", edgecolor="none")
    nov = d[d.evidence_tier == "strong_novel"]
    ax.scatter(nov["prevalence_susceptible"], nov["prevalence_resistant"], s=70,
               facecolor="none", edgecolor="#d62728", lw=1.4, zorder=5)
    ax.plot([0, 1], [0, 1], ls="--", c="grey", lw=0.9)
    ax.set_xlabel("prevalence in susceptible genomes")
    ax.set_ylabel("prevalence in resistant genomes")
    ax.set_title("Biomarker prevalence, resistant vs susceptible\n"
                 "distance from the diagonal = discriminative power", fontsize=10.5)
    # Step 10 scores step 09's candidate list, so biomarkers found ONLY by CPSS have no
    # prevalence row (1,162 of 3,571). State the coverage rather than letting a reader
    # read their absence as an absence of signal.
    ax.text(0.02, 0.96, f"{len(d):,} of {len(bio):,} biomarkers have prevalence stats\n"
                        "(step 10 scores the step-09 candidate list; biomarkers found ONLY by CPSS are not in it)",
            transform=ax.transAxes, fontsize=7.5, color="#666", va="top")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    _save(fig, out, "34_prevalence")


def fig_funnel(bio, db, out):
    """The whole KB as one funnel: models -> graded biomarkers -> tiers -> novel."""
    n_models = int(bio["model_id"].nunique())
    n_bio = len(bio)
    counts = bio["evidence_tier"].value_counts()
    n_card = int(bio["gene_symbol"].notna().sum())
    n_stable = int(pd.to_numeric(bio.get("stable"), errors="coerce").fillna(0).sum())
    if Path(db).exists():
        with sqlite3.connect(db) as c:
            n_models = c.execute("SELECT COUNT(*) FROM models").fetchone()[0] or n_models

    # NOT a nested funnel: 'CARD hit' and 'CPSS-stable' are cross-cutting attributes
    # (2,045 stable > 1,567 with a CARD hit), so tapering the bars by position would
    # assert a subset relation that does not hold. Bar width encodes the COUNT.
    stages = [
        (n_bio, "biomarkers graded", "#41ab5d"),
        (n_stable, "CPSS-stable", "#756bb1"),
        (n_card, "with a CARD hit", "#fdae61"),
        (int(counts.get("confirmed", 0)), "confirmed", "#238b45"),
        (int(counts.get("strong_novel", 0)), "strong_novel", "#d62728"),
    ]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8),
                                 gridspec_kw={"width_ratios": [1.25, 1]})
    a1.axis("off")
    top = max(v for v, _, _ in stages) or 1
    a1.text(0.5, 0.97, f"{n_models} AMR models  →", ha="center", va="center",
            fontsize=12, fontweight="bold", color="#2c7fb8")
    for i, (val, lab, col) in enumerate(stages):
        w = max(0.16, val / top)
        yy = 0.86 - (i + 0.5) / (len(stages) + 0.4)
        a1.barh(yy, w, height=0.11, left=0, color=col, edgecolor="white")
        inside = w > 0.42
        a1.text(w - 0.02 if inside else w + 0.02, yy, f"{val:,}  {lab}",
                ha="right" if inside else "left", va="center", fontsize=11,
                color="white" if inside else "#333", fontweight="bold")
    a1.set_xlim(0, 1.05); a1.set_ylim(0, 1.02)
    a1.set_title("The knowledge base in numbers\n"
                 "(bar length = count; these are overlapping attributes, not nested subsets)",
                 fontsize=11)

    tiers = [t for t in TIER_ORDER if t in counts.index]
    vals = [int(counts[t]) for t in tiers]
    a2.bar(np.arange(len(tiers)), vals,
           color=[TIER_COLOURS[t] for t in tiers], edgecolor="k", lw=0.4)
    for i, v in enumerate(vals):
        a2.text(i, v + max(vals) * 0.015, f"{v:,}", ha="center", fontsize=9)
    a2.set_xticks(np.arange(len(tiers)))
    a2.set_xticklabels(tiers, rotation=20, ha="right", fontsize=9)
    a2.set_ylabel("biomarkers")
    a2.set_title("Evidence tiers (7 orthogonal layers)", fontsize=11)
    fig.tight_layout()
    _save(fig, out, "35_evidence_funnel")


def main():
    ap = argparse.ArgumentParser(description="Thesis figures — biology and the KB.")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--out", default="results/figures")
    ap.add_argument("--only", default="")
    args = ap.parse_args()

    bio = pd.read_csv(Path(args.tables) / "biomarkers.csv")
    out = Path(args.out)
    only = {s.strip() for s in args.only.split(",") if s.strip()}

    todo = [
        ("novel",      lambda: fig_novel(bio, out)),
        ("blast",      lambda: fig_blast(bio, out)),
        ("prevalence", lambda: fig_prevalence(bio, out)),
        ("funnel",     lambda: fig_funnel(bio, args.db, out)),
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
