#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Thesis figures for the unified AMR-KB, built from the tidy tables that
kb_tables.py exports (+ raw 12b null CSVs for the permutation histogram).

Run kb_tables.py FIRST, then:
    python scripts/kb_figures.py --tables results/tables --results results \
        --out figures [--only performance,cpss_pfer,cross_org,mechanism,null_hist,combos]

Each figure is saved as PNG (200 dpi) + PDF. Colours come from PALETTE (registry slugs); any organism not listed there gets an
auto-assigned colour, and display names come from the registry. Edit CLASS_ORDER /
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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Preferred display order for the classes that carry the thesis's headline results;
# _class_order() appends everything else alphabetically, so this list never filters.
CLASS_ORDER = ["penicillins", "cephalosporins", "carbapenems", "monobactams",
               "quinolones", "aminoglycosides", "tetracyclines", "glycylcyclines",
               "glycopeptides", "macrolides", "lincosamides", "phenicols",
               "polymyxins", "folate_pathway_inhibitors"]
# Keys are REGISTRY SLUGS (organisms.yaml / pipeline_runs.organism). They used to be
# short forms ("saureus", "paeruginosa") that never matched the real slugs, so four of
# the six organisms silently fell through to the auto-assigned _EXTRA colours.
PALETTE = {"ecoli": "#2c7fb8", "kpneumoniae": "#de2d26",
           "staphylococcus_aureus": "#756bb1", "acinetobacter_baumannii": "#e6ab02",
           "pseudomonas_aeruginosa": "#31a354", "enterococcus_faecium": "#a6761d"}
_EXTRA = ["#666666", "#1b9e77", "#d95f02"]

# The six layers `classify_evidence_tier()` folds into a grade, in its own order, keyed
# by the token it writes into unitig_evidence_tier.evidence_layers. NOT the seven
# evidence_type values in validation_evidence: `label_permutation` is model-level and
# grades no biomarker, so it belongs to figure 05, not to this per-biomarker grid.
# Same palette and order as kb_figures_biology, so a tier is one colour across the set.
TIER_COLOURS = {"confirmed": "#238b45", "strong_novel": "#d62728",
                "candidate": "#74c476", "weak": "#c6dbef", "none": "#eeeeee"}
TIER_ORDER = ["confirmed", "strong_novel", "candidate", "weak", "none"]

TIER_LAYER_ORDER = [
    ("blast",      "BLAST\n(CARD)"),
    ("prevalence", "Prevalence\nR vs S"),
    ("snp",        "SNP allele\n(CARD var.)"),
    ("mda",        "MDA\npermutation"),
    ("cpss",       "CPSS\nstability"),
    ("pyseer",     "pyseer LMM\n(lineage)"),
]


def _colour(org, _cache={}):
    if org in PALETTE:
        return PALETTE[org]
    return _cache.setdefault(org, _EXTRA[len(_cache) % len(_EXTRA)])


def _short(ab):
    return ab.replace("_", "/")[:18]


def _display(org, _cache={}):
    """'Escherichia coli' for a slug — from the registry, never hardcoded, so the
    figures follow the panel instead of naming two organisms forever."""
    if org not in _cache:
        name = org
        try:
            from lib.registry import get_organism
            name = (get_organism(org) or {}).get("display_name") or org
        except Exception:
            name = {"ecoli": "Escherichia coli",
                    "kpneumoniae": "Klebsiella pneumoniae"}.get(org, org.replace("_", " ").title())
        _cache[org] = name
    return _cache[org]


def _abbr(org):
    """'Ec' from 'Escherichia coli' — genus+species initials of the display name."""
    parts = _display(org).split()
    if len(parts) >= 2:
        return parts[0][0].upper() + parts[1][0].lower()
    return org[:2].title()


def _class_order(series):
    """Drug classes present in the data: the curated CLASS_ORDER first, then anything
    else alphabetically. CLASS_ORDER lists 7 classes while the panel now spans 14, and
    filtering *to* it silently dropped half the classes from the overview figure —
    a figure must never quietly narrow the KB it claims to summarise. CLASS_ORDER now
    lists all 14, but the append-the-rest behaviour is what keeps that safe."""
    present = set(str(c) for c in series.dropna())
    known = [c for c in CLASS_ORDER if c in present]
    return known + sorted(present - set(known))


def _sortkey(df):
    df = df.copy()
    df["_c"] = df["drug_class"].map({c: i for i, c in enumerate(_class_order(df["drug_class"]))}).fillna(99)
    return df.sort_values(["_c", "organism", "antibiotic"])


def _legend(ax, orgs, outside=False):
    """`outside=True` parks the legend right of the axes: with 45 bars a legend drawn
    inside covers real data (it sat on top of the low-AUC bars, which are exactly the
    ones a reader needs to see)."""
    handles = [Patch(color=_colour(o), label=_display(o)) for o in orgs]
    if outside:
        ax.legend(handles=handles, fontsize=9, loc="upper left",
                  bbox_to_anchor=(1.005, 1.0), frameon=False)
    else:
        ax.legend(handles=handles, fontsize=9, loc="lower left")


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
    # Clip the upper whisker at 1.0. ROC-AUC cannot exceed 1, but mean+SD can, and the
    # unclipped whisker ran past the axis and was cut off mid-line — which reads as a
    # rendering fault rather than as a wide interval.
    auc = df["lineage_cv_auc"].to_numpy(float)
    sd = df["lineage_cv_std"].to_numpy(float)
    yerr = np.vstack([sd, np.minimum(sd, 1.0 - auc)])
    ax.bar(x, auc, yerr=yerr, color=col, capsize=3, edgecolor="black",
           linewidth=0.4, alpha=0.9)
    ax.axhline(0.5, ls="--", c="grey", lw=0.8)
    # Park the chance label in clear air on the right; at x=0.3 it sat behind the first
    # bar and was unreadable.
    ax.text(len(df) - 1.2, 0.505, "chance", fontsize=8, color="grey",
            va="bottom", ha="right")
    # Floor below the weakest model (0.429 for A. baumannii ceftazidime): a 0.4 floor
    # clipped that bar to an invisible sliver, hiding the panel's most informative
    # result — the clonally-confounded model lineage-CV is supposed to expose.
    ax.set_ylim(min(0.40, float(auc.min()) - 0.06), 1.02)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{_short(a)}\n({_abbr(o)})" for a, o in zip(df.antibiotic, df.organism)], rotation=90, fontsize=7.5)
    ax.set_ylabel("Lineage-aware CV ROC-AUC (mean ± SD)")

    # Call out any model at or below chance. This is the thesis's clonal-confounding
    # exhibit, and an unlabelled short bar in a row of 45 is easy to read past.
    below = np.where(auc <= 0.5)[0]
    for i in below:
        ax.annotate(f"{auc[i]:.3f}\nbelow chance", xy=(x[i], auc[i]),
                    xytext=(x[i] + 2.6, 0.455), fontsize=7.5, color="#b2182b",
                    fontweight="bold", ha="center", va="center",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec="#b2182b", lw=0.6, alpha=0.95),
                    arrowprops=dict(arrowstyle="->", color="#b2182b", lw=1.0))
    n_below = len(below)
    sub = (f" — {n_below} model{'s' if n_below != 1 else ''} at or below chance"
           if n_below else "")
    ax.set_title(f"Per-antibiotic generalisation performance, lineage-aware CV "
                 f"(n={len(df)} models, mean {auc.mean():.3f}){sub}", fontsize=11)
    _legend(ax, df["organism"].unique(), outside=True)
    _save(fig, out, "01_performance_lineageCV")


def fig_cpss_pfer(tables, out):
    df = _sortkey(pd.read_csv(tables / "kb_overview.csv"))
    x = np.arange(len(df))
    col = [_colour(o) for o in df["organism"]]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True, gridspec_kw={"hspace": 0.08})
    a1.bar(x, df["cpss_n_stable"], color=col, edgecolor="k", lw=0.4, alpha=0.9)
    a1.set_ylabel("CPSS stable unitigs (π≥0.6)")
    # Quantify the PFER>1 models. HANDOFF/METHODOLOGY carry this as a stated
    # limitation ("some models exceed one expected false positive"), and the figure is
    # where a reader should be able to check it rather than take it on trust.
    pf = pd.to_numeric(df["pfer_bound"], errors="coerce")
    n_over = int((pf > 1).sum())
    a1.set_title("CPSS stability selection — stable biomarker count & PFER bound\n"
                 f"{n_over} of {len(df)} models exceed PFER = 1 "
                 f"(max {pf.max():.1f}): in those stable sets more than one selected "
                 "unitig is expected to be a false positive", fontsize=10.5)
    a2.bar(x, df["pfer_bound"], color=col, edgecolor="k", lw=0.4, alpha=0.9)
    a2.set_yscale("log")
    a2.axhline(1, ls="--", c="grey", lw=0.8)
    a2.text(len(df) - 0.5, 1.06, "PFER = 1", fontsize=7.5, color="grey",
            ha="right", va="bottom")
    a2.set_ylabel("PFER bound (E[false positives], log)")
    a2.set_xticks(x)
    a2.set_xticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(df.antibiotic, df.organism)], rotation=90, fontsize=7.5)
    _legend(a1, df["organism"].unique(), outside=True)
    _save(fig, out, "02_cpss_pfer")


_FAM_MAP = {
    "16S rRNA methyltransferase (G1405)": "16S-RMTase",
    "General Bacterial Porin with reduced permeability to beta-lactams": "porin loss",
    "aminoglycoside bifunctional resistance protein": "AAC(6')-APH",
    "major facilitator superfamily (MFS) antibiotic efflux pump": "MFS efflux",
    "resistance-nodulation-cell division (RND) antibiotic efflux pump": "RND efflux",
    "ATP-binding cassette (ABC) antibiotic efflux pump": "ABC efflux",
    "small multidrug resistance (SMR) antibiotic efflux pump": "SMR efflux",
    "OXA beta-lactamase;OXA-48-like beta-lactamase": "OXA-48-like",
    "sulfonamide resistant sul": "sul",
    "trimethoprim resistant dihydrofolate reductase dfr": "dfr",
    "Erm 23S ribosomal RNA methyltransferase": "Erm (23S rRNA MTase)",
    "macrolide phosphotransferase (MPH)": "MPH",
    "chloramphenicol acetyltransferase (CAT)": "CAT",
    "streptothricin acetyltransferase (SAT)": "SAT",
    "msr-type ABC-F protein": "Msr (ABC-F)",
    "methicillin resistant PBP2": "PBP2 (mecA)",
    "tetracycline-resistant ribosomal protection protein": "ribosomal protection",
    "ADC beta-lactamases pending classification for carbapenemase activity": "ADC",
    "Van ligase;glycopeptide resistance gene cluster": "vanA ligase",
    "AAC(6');AAC(6')-Ib-cr": "AAC(6')",
}


# Last words that carry no information on their own. Taking the final token of an ARO
# family name is a decent shortener ("...APH(3')" -> "APH(3')"), but for families that
# end in a generic noun it produced labels like "protein" and "pump" on the figures —
# unreadable, and indistinguishable between families.
_GENERIC_TAIL = {"protein", "proteins", "pump", "pumps", "enzyme", "gene", "genes",
                 "family", "system", "transporter", "determinant", "cluster"}


def _fam(s):
    """Short, human gene-family label for figures."""
    s = str(s).strip()
    if s in _FAM_MAP:
        return _FAM_MAP[s]
    if not s:
        return s
    # ARO composes some families as semicolon-joined qualifiers. Splitting on
    # whitespace alone left labels like "cluster;vanS" and "AAC(6');AAC(6')-Ib-cr"
    # on the axes, so resolve the two shapes that occur here before anything else:
    #   "glycopeptide resistance gene cluster;vanS" -> "vanS (van cluster)"
    #   "OXA beta-lactamase;OXA-23-like beta-lactamase" -> "OXA-23-like"
    if ";" in s:
        head, _, tail = s.partition(";")
        head, tail = head.strip(), tail.strip()
        if "cluster" in head.lower() and tail:
            return f"{tail} (van cluster)" if tail.lower().startswith("van") else tail
        # Same family qualified twice: keep the more specific (longer) side.
        pick = tail if len(tail) >= len(head) else head
        return _fam(pick) if pick != s else pick
    if "beta-lactamase" in s:
        return s.replace(" beta-lactamase", "").strip()
    parts = s.split()
    if parts[-1].lower() in _GENERIC_TAIL:
        # Drop the trailing generic words, then keep the two that identify the family:
        #   "tetracycline-resistant ribosomal protection protein" -> "ribosomal protection"
        #   "glycopeptide resistance gene cluster"                -> "glycopeptide resistance"
        while parts and parts[-1].lower() in _GENERIC_TAIL:
            parts.pop()
        return " ".join(parts[-2:]) if parts else s
    return parts[-1]


CLASS_SHORT = {"folate_pathway_inhibitors": "folate inhibitors"}


def fig_overview(tables, out, db):
    """Cover slide: scope of the KB (models / organisms / classes / genomes) +
    models-per-drug-class stacked by organism."""
    ms = pd.read_csv(tables / "models_summary.csv")
    order = _class_order(ms.drug_class)
    orgs = list(ms.organism.unique())
    fig = plt.figure(figsize=(13, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.85, 1.7], wspace=0.28)
    a0 = fig.add_subplot(gs[0]); a1 = fig.add_subplot(gs[1])
    a0.axis("off")
    cards = [(str(len(ms)), "AMR models"),
             (str(ms.organism.nunique()), "ESKAPEE organisms"),
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
                label=_display(org))
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
    fig, ax = plt.subplots(figsize=(13, 1.0 + 1.0 * len(abs_)))
    ax.axis("off")
    for i, ab in enumerate(abs_):
        yy = len(abs_) - 1 - i
        sub = ot[ot.antibiotic == ab]
        # {gene family -> organisms that recovered it}. This used to intersect two
        # hardcoded organisms (ecoli/kpneumoniae) and label everything else "-only",
        # which silently ignored the other four organisms of the panel: a drug like
        # ciprofloxacin is assayed in five.
        fam_orgs = {}
        for o, g in sub.groupby("organism"):
            for f in {_fam(x) for x in g["aro_gene_family"].dropna()}:
                fam_orgs.setdefault(f, set()).add(o)
        shared = sorted((f for f, o in fam_orgs.items() if len(o) >= 2),
                        key=lambda f: (-len(fam_orgs[f]), f))
        single = sorted((f for f, o in fam_orgs.items() if len(o) == 1), key=str)
        ax.text(0.0, yy, _short(ab), fontsize=12, fontweight="bold", va="center")
        txt = "   ".join(f"{f} ({','.join(sorted(_abbr(o) for o in fam_orgs[f]))})"
                         for f in shared) or "—"
        ax.text(0.26, yy, txt, fontsize=12, fontweight="bold", color="#2ca25f", va="center")
        if single:
            per = ", ".join(f"{_abbr(next(iter(fam_orgs[f])))}: {f}" for f in single[:6])
            if len(single) > 6:
                per += f", +{len(single) - 6} more"
            ax.text(0.26, yy - 0.30, "single-organism — " + per,
                    fontsize=8.5, color="#888", va="center")
    n_org = ot.groupby("antibiotic")["organism"].nunique().reindex(abs_)
    ax.text(0.26, len(abs_) - 0.30,
            f"SHARED gene family — recovered in ≥2 organisms (concordant); "
            f"drugs span up to {int(n_org.max())} organisms",
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
    # Say how many models are drawn. Columns exist only for models with at least one
    # on-target hit, and a reader who is not told that reads the missing ones as
    # "no resistance genes found" rather than "off-target hits only".
    n_models_total = len(ms)
    ax.set_title("On-target confirmed resistance gene families across models "
                 "(cell = # unitigs)\n"
                 f"{len(cols)} of {n_models_total} models have an on-target hit; "
                 "the rest recovered only off-target or no CARD genes",
                 fontsize=10.5)
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
    n_perm = 0
    for ax, f in zip(axes, files):
        ab = os.path.basename(f).replace("12b_label_permutation_nulls_", "").replace(".csv", "")
        # The organism is in the path, not the filename. Without it four panels were
        # titled "tetracycline" and four "ciprofloxacin", with no way to tell them apart.
        org = Path(f).parts[-4] if len(Path(f).parts) >= 4 else ""
        d = pd.read_csv(f)
        # 12b writes `permutation,null_roc_auc`. Taking columns[0] histogrammed the
        # permutation INDEX (1..N), not the null AUC — and the hardcoded 0.4-1.0 xlim
        # then pushed that 1..50 range entirely off-axis, so all 45 panels rendered as
        # empty boxes with a lone red line and nobody could tell the data was wrong.
        # Name the column; fall back to the last one only if the schema changes.
        col = "null_roc_auc" if "null_roc_auc" in d.columns else d.columns[-1]
        nulls = pd.to_numeric(d[col], errors="coerce").dropna()
        if nulls.empty:
            ax.axis("off")
            continue
        summ = json.load(open(f.replace("_nulls_", "_summary_").replace(".csv", ".json"))) if os.path.exists(f.replace("_nulls_", "_summary_").replace(".csv", ".json")) else {}
        real = summ.get("real_roc_auc") or summ.get("real_test_roc_auc")
        n_perm = max(n_perm, len(nulls))
        # Two reasons this panel used to render as an empty box with one red line:
        # a white edgecolor on bars only ~2 px wide at this panel size painted over the
        # fill entirely, and a hardcoded xlim of 0.4-1.0 squeezed a null that lives in
        # ~0.50-0.65 into a quarter of the axis. Fill without an edge, and let each
        # panel frame its own data.
        ax.hist(nulls, bins=12, color="#9e9ac8", edgecolor="none")
        if real:
            ax.axvline(real, color="#d7301f", lw=1.6)
        lo = float(min(nulls.min(), real if real else nulls.min()))
        hi = float(max(nulls.max(), real if real else nulls.max()))
        pad = max(0.02, 0.08 * (hi - lo))
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_title(f"{_short(ab)} ({_abbr(org)})" if org else _short(ab), fontsize=7.5)
        ax.tick_params(labelsize=6)
        ax.set_yticks([])
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("Label-permutation null (purple) vs the model's REAL ROC-AUC (red line)\n"
                 f"{n} models · N={n_perm} shuffles each — every real AUC sits outside its "
                 "own null, so p is at the 1/(N+1) floor for all of them\n"
                 "each panel is scaled to its own null: the nulls differ in location and "
                 "width, which a shared axis hid", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _save(fig, out, "05_label_permutation_nulls")


def fig_evidence_layers(tables, out, db):
    """THE backbone figure: how many of each model's biomarkers PASSED each tier layer.

    Rewritten. It used to count rows in ``validation_evidence``, which is the number of
    candidates *tested* — not the number that passed. Because BLAST, prevalence and MDA
    all test the same candidate set, those three columns were identical in every row
    (51/51/51, 59/59/59 …), and the figure asserted 40–70 MDA hits per model while
    figure 28 correctly reported that **zero** unitigs clear the MDA FDR in any model.
    The SNP column had a second fault: it counted step-11 unitigs that are not
    biomarkers of that model at all (they never enter ``unitig_model_scores``).

    The tier ladder keeps its own record of what actually passed —
    ``unitig_evidence_tier.evidence_layers`` — so read that. Zeros are drawn as an
    explicit grey ``0`` rather than left blank: two layers are zero everywhere, and a
    blank cell reads as "not applicable" when the honest statement is "measured, none
    passed". See METHODOLOGY.md §5.2.
    """
    import sqlite3
    ms = _sortkey(pd.read_csv(tables / "models_summary.csv")).reset_index(drop=True)
    conn = sqlite3.connect(str(db))
    counts, totals = {}, {}
    for run_id, layers, n in conn.execute(
            "SELECT p.run_id, e.evidence_layers, COUNT(*) "
            "FROM unitig_evidence_tier e "
            "JOIN models m ON m.model_id = e.model_id "
            "JOIN pipeline_runs p ON p.run_id = m.run_id "
            "GROUP BY p.run_id, e.evidence_layers"):
        totals[run_id] = totals.get(run_id, 0) + n
        for tok in (layers or "").split(","):
            tok = tok.strip()
            if tok:
                counts.setdefault(run_id, {})[tok] = counts.get(run_id, {}).get(tok, 0) + n
    conn.close()

    keys = [k for k, _ in TIER_LAYER_ORDER]
    M = np.array([[counts.get(r, {}).get(k, 0) for k in keys] for r in ms["run_id"]], float)
    disp = np.where(M > 0, M, np.nan)
    fig, ax = plt.subplots(figsize=(8.8, 0.45 * len(ms) + 2.0))
    from matplotlib.colors import LogNorm
    im = ax.imshow(disp, aspect="auto", cmap="YlGnBu",
                   norm=LogNorm(vmin=1, vmax=np.nanmax(disp)))
    ax.set_xticks(range(len(keys)))
    ax.set_xticklabels([lab for _, lab in TIER_LAYER_ORDER], fontsize=8)
    ax.set_yticks(range(len(ms)))
    ax.set_yticklabels([f"{_short(a)} ({_abbr(o)})  n={totals.get(r, 0)}"
                        for a, o, r in zip(ms.antibiotic, ms.organism, ms.run_id)],
                       fontsize=7.5)
    thr = np.nanmax(disp) ** 0.5
    for i in range(len(ms)):
        for j in range(len(keys)):
            v = M[i, j]
            ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6.5,
                    color=("white" if v > thr else "black") if v > 0 else "#b0b0b0")
    dead = [lab.replace("\n", " ") for k, lab in TIER_LAYER_ORDER if M[:, keys.index(k)].sum() == 0]
    note = (f"  ·  {' and '.join(dead)} pass nowhere (0/45 models)" if dead else "")
    ax.set_title("Biomarkers passing each tier layer, per model\n"
                 "cell = biomarkers of that model passing that layer; "
                 "row label n = biomarkers graded" + note, fontsize=10)
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="biomarkers passing (log)")
    _save(fig, out, "06_evidence_layers")


def fig_evidence_combinations(tables, out, db):
    """Which evidence layers fire TOGETHER, and what grade each combination earns.

    Figure 06 counts each layer separately, per model. That marginal view cannot show
    co-occurrence: it says prevalence fired 2,055 times and pyseer 1,172 times without
    saying that 558 biomarkers rest on those two and nothing else. This is the UpSet of
    the same ladder, over all 3,571 graded (unitig, model) pairs.

    Worth drawing because the mapping turns out to be deterministic: every one of the 14
    observed combinations earns exactly one tier, so colouring the bars by tier makes
    `classify_evidence_tier()` readable without opening the code. Two things become
    visible that prose keeps having to assert:

      * `blast` alone earns `candidate` while a lone statistical layer earns `weak` --
        a CARD hit outweighs one statistical signal.
      * `strong_novel` is not a separate test. It is exactly the cell where three
        statistical layers fire and BLAST does not: 23 biomarkers with real evidence and
        no known CARD gene behind them.

    The `snp` and `mda` rows are empty across every column, and are drawn greyed with
    the reason rather than omitted -- a missing row reads as "not measured" when the
    honest statement is "measured, never fires". See METHODOLOGY.md 5.2.
    """
    import sqlite3
    from matplotlib.patches import Rectangle
    conn = sqlite3.connect(str(db))
    raw = conn.execute("SELECT COALESCE(evidence_layers,''), evidence_tier, COUNT(*) "
                       "FROM unitig_evidence_tier GROUP BY 1,2").fetchall()
    conn.close()

    combos = {}
    for layers, tier, n in raw:
        toks = tuple(t for t in (layers or "").split(",") if t.strip())
        combos.setdefault(toks, {})[tier] = combos.get(toks, {}).get(tier, 0) + n
    # If a combination ever earned two grades the colour would be a lie, so say so.
    ambiguous = {k: v for k, v in combos.items() if len(v) > 1}
    items = sorted(((k, max(v, key=v.get), sum(v.values())) for k, v in combos.items()),
                   key=lambda x: -x[2])
    total = sum(n for _, _, n in items)

    keys = [k for k, _ in TIER_LAYER_ORDER]
    marg = {k: sum(n for toks, _, n in items if k in toks) for k in keys}
    row_order = ([kv for kv in TIER_LAYER_ORDER if marg[kv[0]]]
                 + [kv for kv in TIER_LAYER_ORDER if not marg[kv[0]]])
    n_live = sum(1 for kv in TIER_LAYER_ORDER if marg[kv[0]])

    fig, (ax, axm) = plt.subplots(
        2, 1, figsize=(12.4, 7.4), sharex=True,
        gridspec_kw={"height_ratios": [3, 2], "hspace": 0.06})

    xs = range(len(items))
    for x, (toks, tier, n) in zip(xs, items):
        ax.bar(x, n, width=0.62, color=TIER_COLOURS.get(tier, "#888888"),
               edgecolor="k", lw=0.4, zorder=3)
        ax.annotate(f"{n:,}", (x, n), textcoords="offset points", xytext=(0, 3),
                    ha="center", fontsize=8, zorder=4)
    ax.set_yscale("log")
    ax.set_ylim(0.7, max(n for _, _, n in items) * 2.4)
    ax.set_ylabel("biomarkers (log)")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.5, zorder=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)

    handles = [plt.Rectangle((0, 0), 1, 1, fc=TIER_COLOURS.get(t, "#888"), ec="k", lw=0.4)
               for t in TIER_ORDER]
    ax.legend(handles, TIER_ORDER, ncol=5, fontsize=8, frameon=False,
              loc="upper right", title="evidence tier earned", title_fontsize=8)

    # The novel cell is the whole point of the KB; point at it.
    for x, (toks, tier, n) in zip(xs, items):
        if tier == "strong_novel":
            ax.add_patch(Rectangle((x - 0.46, 0.7), 0.92, n / 0.7, fill=False,
                                   ec=TIER_COLOURS["strong_novel"], lw=1.4, ls="--",
                                   zorder=5, clip_on=False))
            ax.annotate("three statistical layers,\nno CARD hit = novel",
                        (x, n), textcoords="offset points", xytext=(0, 52),
                        ha="center", fontsize=8, color=TIER_COLOURS["strong_novel"],
                        arrowprops=dict(arrowstyle="->", lw=1.0, shrinkB=16,
                                        color=TIER_COLOURS["strong_novel"]))

    dead = []
    for y, (k, lab) in enumerate(row_order):
        alive = marg[k] > 0
        if not alive:
            dead.append(lab.replace("\n", " "))
        axm.axhline(y, color="#eeeeee", lw=8, zorder=0)
        for x, (toks, _, _) in zip(xs, items):
            on = k in toks
            axm.scatter(x, y, s=54, zorder=3,
                        color="#333333" if on else ("#dddddd" if alive else "#f2f2f2"))
    for x, (toks, _, _) in zip(xs, items):
        ys = [i for i, (kk, _) in enumerate(row_order) if kk in toks]
        if len(ys) > 1:
            axm.plot([x, x], [min(ys), max(ys)], color="#333333", lw=1.6, zorder=2)

    axm.set_yticks(range(len(row_order)))
    axm.set_yticklabels(
        [f"{lab.replace(chr(10), ' ')}  ({marg[k]:,})" if marg[k]
         else f"{lab.replace(chr(10), ' ')}  (0 — never fires)"
         for k, lab in row_order], fontsize=8.5)
    for tick, (k, _) in zip(axm.get_yticklabels(), row_order):
        if not marg[k]:
            tick.set_color("#b00020")
    if n_live < len(row_order):        # rule off the layers that never fire
        axm.axhline(n_live - 0.5, color="#b00020", lw=0.8, ls=":", zorder=4)
    axm.set_ylim(len(row_order) - 0.5, -0.5)
    axm.set_xlim(-0.7, len(items) - 0.3)
    axm.set_xticks(list(xs))
    axm.set_xticklabels([f"{len(t)}" if t else "0" for t, _, _ in items], fontsize=8)
    axm.set_xlabel("layers firing in the combination  ·  one column per observed combination")
    for side in ("top", "right", "left"):
        axm.spines[side].set_visible(False)
    axm.tick_params(axis="y", length=0)

    sub = (f"{total:,} graded (unitig, model) pairs · {len(items)} observed combinations · "
           f"each earns exactly one tier, so the colour IS the grading rule")
    if ambiguous:
        sub = (f"{total:,} pairs · WARNING: {len(ambiguous)} combination(s) span more than "
               f"one tier; bar colour shows the majority grade")
    if dead:
        sub += "\nmeasured but never fires: " + ", ".join(dead)
    fig.suptitle("How the evidence layers combine, and what each combination is worth\n"
                 + sub, fontsize=11.5)
    fig.subplots_adjust(top=0.84)
    _save(fig, out, "39_evidence_combinations")


def fig_significance(tables, out, db):
    """05 — model-level significance: the observed AUC of step 12b's label-permutation
    test vs its shuffled-label null, per model.

    The observed value is 12b's OWN split AUC, NOT the lineage-CV score: the two differ
    sharply exactly where it matters (A. baumannii ceftazidime is 0.91 here and 0.429
    under lineage-CV). This figure used to call it 'REAL lineage-CV AUC', which invited
    the reader to conclude that a clonally-confounded model generalises. The lineage-CV
    value is now drawn as a separate black tick so both are visible and distinct."""
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
    lcv = list(ms["lineage_cv_auc"])
    fig, ax = plt.subplots(figsize=(9.5, 0.42 * len(ms) + 1.4))
    for yi, ri, ni, ci, li in zip(y, r, nm, col, lcv):
        ax.plot([ni, ri], [yi, yi], color="lightgrey", lw=2, zorder=1)
        ax.scatter(ni, yi, color="#999999", s=28, zorder=2)
        ax.scatter(ri, yi, color=ci, s=46, zorder=3)
        ax.scatter(li, yi, marker="|", color="black", s=90, linewidths=1.4, zorder=4)
    ax.axvline(0.5, ls="--", c="grey", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{_short(a)} ({_abbr(o)})" for a, o in zip(ms.antibiotic, ms.organism)], fontsize=7.5)
    ax.set_xlim(min(0.4, float(min(lcv)) - 0.03), 1.0)
    ax.set_xlabel("ROC-AUC")
    ax.set_title("Model-level significance: observed AUC (colour) ≫ label-shuffle null max (grey)\n"
                 "black tick = lineage-aware CV AUC (the reported, generalisation metric)\n"
                 "every model sits at p = 1/(N+1) \u2248 0.02, the SMALLEST value "
                 "N=50 permutations can return — a floor, not a measured difference "
                 "between models", fontsize=9.5)
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color="#444444", label="observed AUC (12b split)"),
        Line2D([], [], marker="o", ls="", color="#999999", label="label-shuffle null max"),
        Line2D([], [], marker="|", ls="", color="black", markeredgewidth=1.4, label="lineage-CV AUC"),
    ], fontsize=8, loc="lower left", frameon=False)
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
        "combos": lambda t, r, o, db: fig_evidence_combinations(t, o, db),
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
