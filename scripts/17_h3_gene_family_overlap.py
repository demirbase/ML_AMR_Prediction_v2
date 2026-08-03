#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H3 — cross-antibiotic biomarker overlap, tested the way E2 says it must be.

H3: two antibiotics of the SAME class share resistance determinants more than two
antibiotics of different classes.

Why this is not a Fisher test on raw unitigs (docs/literature/E2.md): unitigs are in
tight linkage — one plasmid contributes thousands of correlated unitigs — so the
i.i.d. assumption collapses and the p-values come out spuriously microscopic. E2's
prescription is followed here:

  1. COMPONENTISE first. The unit of analysis is the ARO **gene family** a biomarker
     maps to, not the unitig. Two drugs can share a mechanism yet no exact unitig
     (segmentation is genome-set dependent), and a single family absorbs its own
     linked unitigs.
  2. NULL UNIVERSE = what both analyses could actually have recovered, i.e. the gene
     families observed across that ORGANISM's models — NOT the union of the two sets
     being compared (which inflates enrichment by construction), and not the pan-genome
     or all of CARD (which deflates it). Both models of a pair see the same genomes and
     the same CARD snapshot, so their reachable spaces coincide.
  3. Fisher exact (one-sided, enrichment) on the components.
  4. Monte-Carlo check: draw same-sized family sets from the universe to confirm the
     analytic p is not an artefact of the small-N regime.
  5. **Benjamini-Yekutieli**, not BH: cross-resistance makes the pair tests negatively
     as well as positively dependent, which violates BH's PRDS assumption.
  6. Report **Overlap Coefficient** (k/min) and **Fold Enrichment**, not Jaccard —
     our sets differ in size by an order of magnitude and Jaccard just tracks the
     size imbalance.
  7. H3 itself is then a within- vs cross-class contrast of the overlap coefficients
     (Mann-Whitney U), because that is the claim, not any single pair's p-value.

Usage:
    python scripts/17_h3_gene_family_overlap.py --db results/kb/amrk.db \
        --tables results/tables --figures results/figures [--permutations 1000]

Outputs:
    results/tables/h3_gene_family_overlap.csv   one row per (organism, ab1, ab2)
    results/tables/h3_summary.json              the within-vs-cross contrast
    results/figures/08_h3_overlap.png (+pdf)
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

try:
    # Same organism abbreviations as every other figure (Ab/Ec/Ef/Kp/Pa/Sa, derived from
    # the registry). A local slug slice produced "Ac"/"En"/"Ps"/"St" here and nowhere else.
    from kb_figures import _abbr, _short
except Exception:                                            # pragma: no cover
    def _abbr(org):
        return org[:2].title()

    def _short(ab):
        return str(ab).replace("_", "/")[:18]

TIERS = ("confirmed", "candidate")


def family_sets(conn, tiers=TIERS):
    """{(organism, antibiotic): {aro_gene_family}} — the components of H3."""
    ph = ",".join("?" * len(tiers))
    q = (f"SELECT p.organism, m.antibiotic, b.aro_gene_family "
         f"FROM blast_annotations b "
         f"JOIN models m ON m.model_id = b.model_id "
         f"JOIN pipeline_runs p ON p.run_id = m.run_id "
         f"WHERE b.tier IN ({ph}) AND b.aro_gene_family IS NOT NULL "
         f"AND TRIM(b.aro_gene_family) != ''")
    sets = {}
    for org, ab, fam in conn.execute(q, tiers):
        sets.setdefault((org, ab), set()).add(fam.strip())
    return sets


def fisher_greater(k, K, n, N):
    """One-sided Fisher exact (enrichment) = hypergeometric survival P(X >= k)."""
    try:
        from scipy.stats import hypergeom
        return float(hypergeom.sf(k - 1, N, K, n))
    except Exception:                                        # pragma: no cover
        from math import comb
        if N <= 0 or K > N or n > N:
            return float("nan")
        total = comb(N, n)
        if total == 0:
            return float("nan")
        return sum(comb(K, i) * comb(N - K, n - i)
                   for i in range(k, min(K, n) + 1)) / total


def mc_pvalue(k, K, n, N, B, rng):
    """Empirical P(overlap >= k) when two sets of sizes K and n are drawn at random
    from a universe of N components. Guards the analytic p in the small-N regime."""
    if N <= 0 or K <= 0 or n <= 0:
        return float("nan")
    univ = np.arange(N)
    hits = 0
    for _ in range(B):
        a = rng.choice(univ, size=min(K, N), replace=False)
        b = rng.choice(univ, size=min(n, N), replace=False)
        if len(np.intersect1d(a, b, assume_unique=True)) >= k:
            hits += 1
    return (hits + 1) / (B + 1)          # add-one: never reports p = 0


def benjamini_yekutieli(pvals):
    """BY-adjusted p-values (BH scaled by the harmonic number c(m)).

    BH is not valid here: cross-resistance induces both positive and negative
    dependence between pair tests, so the PRDS condition BH needs does not hold.
    BY is valid under ARBITRARY dependence — the conservative, defensible choice.
    """
    p = np.asarray(pvals, dtype=float)
    ok = ~np.isnan(p)
    out = np.full(p.shape, np.nan)
    m = int(ok.sum())
    if m == 0:
        return out
    c_m = np.sum(1.0 / np.arange(1, m + 1))
    idx = np.argsort(p[ok])
    ranked = p[ok][idx]
    adj = ranked * m * c_m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]     # enforce monotonicity
    res = np.empty(m)
    res[idx] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


def main():
    ap = argparse.ArgumentParser(description="H3 gene-family overlap test (E2 framework).")
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--figures", default="results/figures")
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    conn = sqlite3.connect(args.db)
    cls_of = dict(conn.execute("SELECT antibiotic, drug_class FROM antibiotics").fetchall())
    sets = family_sets(conn)
    conn.close()
    if not sets:
        sys.exit("ERROR: no ARO gene families in the KB — run 08/09 + populate first.")

    rng = np.random.default_rng(args.seed)
    organisms = sorted({o for o, _ in sets})
    rows = []
    for org in organisms:
        abs_ = sorted(ab for (o, ab) in sets if o == org)
        # Universe: components reachable by this organism's analyses (see docstring).
        universe = set().union(*(sets[(org, ab)] for ab in abs_))
        N = len(universe)
        for i, a1 in enumerate(abs_):
            for a2 in abs_[i + 1:]:
                A, B_ = sets[(org, a1)], sets[(org, a2)]
                k = len(A & B_)
                K, n = len(A), len(B_)
                oc = k / min(K, n) if min(K, n) else float("nan")
                fe = ((k / n) / (K / N)) if (n and K and N) else float("nan")
                rows.append({
                    "organism": org, "ab1": a1, "ab2": a2,
                    "class1": cls_of.get(a1), "class2": cls_of.get(a2),
                    "same_class": bool(cls_of.get(a1) and cls_of.get(a1) == cls_of.get(a2)),
                    "n_families_ab1": K, "n_families_ab2": n,
                    "n_shared": k, "universe": N,
                    "overlap_coefficient": round(oc, 4) if oc == oc else None,
                    "fold_enrichment": round(fe, 3) if fe == fe else None,
                    "fisher_p": fisher_greater(k, K, n, N),
                    "mc_p": mc_pvalue(k, K, n, N, args.permutations, rng),
                    "shared_families": "; ".join(sorted(A & B_)),
                })

    df = pd.DataFrame(rows)
    df["fisher_p_BY"] = benjamini_yekutieli(df["fisher_p"].to_numpy())
    df = df.sort_values(["organism", "fisher_p"]).reset_index(drop=True)
    tables = Path(args.tables); tables.mkdir(parents=True, exist_ok=True)
    df.to_csv(tables / "h3_gene_family_overlap.csv", index=False)
    print(f"  ✓ {tables/'h3_gene_family_overlap.csv'}  ({len(df)} pairs)")

    # ---- H3: within-class vs cross-class overlap ---------------------------
    within = df.loc[df.same_class, "overlap_coefficient"].dropna()
    cross = df.loc[~df.same_class, "overlap_coefficient"].dropna()
    summary = {
        "n_pairs": int(len(df)),
        "n_within_class": int(len(within)), "n_cross_class": int(len(cross)),
        "mean_overlap_within": round(float(within.mean()), 4) if len(within) else None,
        "mean_overlap_cross": round(float(cross.mean()), 4) if len(cross) else None,
        "n_significant_BY_0.05": int((df["fisher_p_BY"] <= 0.05).sum()),
        "permutations": args.permutations,
        "universe_definition": "ARO gene families observed across the organism's models",
        "correction": "Benjamini-Yekutieli (arbitrary dependence)",
        # Two caveats that must travel WITH the result, not be discovered by a reviewer:
        # (1) the panel was curated to maximise CLASS coverage, which by construction
        #     leaves very few same-class pairs — H3's within group is small;
        # (2) each pair compares a handful of gene families drawn from a small universe,
        #     so a single pair has almost no power and BY (correctly conservative under
        #     dependence) leaves none individually significant. H3 rests on the CONTRAST
        #     between the two groups, not on any pair's p-value.
        "caveat_within_group_size": (
            "few same-class pairs by construction: the panel was curated for class "
            "coverage, trimming redundant same-class drugs"),
        "caveat_per_pair_power": (
            "per-pair tests are underpowered at these set sizes; the H3 claim is the "
            "within- vs cross-class contrast, not individual pair significance"),
    }
    if len(within) and len(cross):
        try:
            from scipy.stats import mannwhitneyu
            u, p = mannwhitneyu(within, cross, alternative="greater")
            summary["mannwhitney_U"] = float(u)
            summary["mannwhitney_p_within_gt_cross"] = float(p)
            summary["H3_supported"] = bool(p <= 0.05)
        except Exception as e:                               # pragma: no cover
            summary["mannwhitney_error"] = str(e)
    (tables / "h3_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"  ✓ {tables/'h3_summary.json'}")
    print("\n" + "=" * 64)
    for k_, v in summary.items():
        print(f"  {k_:32s} {v}")
    print("=" * 64)

    # ---- figure: within vs cross + the significant pairs --------------------
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 5.2),
                                 gridspec_kw={"width_ratios": [0.8, 1.6]})
    if len(within) and len(cross):
        # matplotlib >=3.9 renamed boxplot's `labels` to `tick_labels`; support both so
        # the figure does not depend on which matplotlib the container resolved.
        try:
            a1.boxplot([cross.values, within.values],
                       tick_labels=["cross-class", "within-class"],
                       widths=0.55, showfliers=False)
        except TypeError:                                    # pragma: no cover
            a1.boxplot([cross.values, within.values],
                       labels=["cross-class", "within-class"],
                       widths=0.55, showfliers=False)
        for xi, vals in ((1, cross.values), (2, within.values)):
            a1.scatter(np.random.default_rng(1).normal(xi, 0.06, len(vals)), vals,
                       s=18, alpha=0.55, color="#2c7fb8" if xi == 1 else "#de2d26")
    a1.set_ylabel("Overlap coefficient  (shared / min set size)")
    ttl = "H3: within-class pairs share more determinants"
    if "mannwhitney_p_within_gt_cross" in summary:
        ttl += f"\nMann-Whitney p = {summary['mannwhitney_p_within_gt_cross']:.3g}"
    a1.set_title(ttl, fontsize=10)

    # Rank by how MANY families are shared, not by the coefficient: with sets of one or
    # two families k/min(K,n) is 1.0 for any overlap at all, so ranking on it produced 18
    # bars of identical length that separated nothing. Pairs where the smaller set has a
    # single family are dropped for the same reason — their coefficient is uninformative.
    cand = df.dropna(subset=["overlap_coefficient"]).copy()
    cand = cand[(cand[["n_families_ab1", "n_families_ab2"]].min(axis=1) >= 2)
                & (cand["n_shared"] >= 1)]
    top = cand.sort_values(["n_shared", "fold_enrichment"], ascending=False).head(18)
    if top.empty:
        a2.axis("off")
        a2.text(0.5, 0.5, "no pair has ≥2 families on both sides", ha="center", fontsize=10)
    else:
        y = np.arange(len(top))
        col = ["#de2d26" if s else "#2c7fb8" for s in top.same_class]
        a2.barh(y, top["n_shared"], color=col, edgecolor="k", lw=0.4)
        a2.set_yticks(y)
        a2.set_yticklabels([f"{_short(r.ab1)}–{_short(r.ab2)} ({_abbr(r.organism)})"
                            for r in top.itertuples()], fontsize=7)
        for yi, r in zip(y, top.itertuples()):
            a2.text(r.n_shared + 0.06, yi, f"OC {r.overlap_coefficient:.2f} · FE {r.fold_enrichment:g}",
                    va="center", fontsize=6.5, color="#555")
        a2.invert_yaxis()
        a2.set_xlim(0, top["n_shared"].max() * 1.55)
        a2.set_xlabel("shared gene families (k)")
        a2.set_title("Pairs sharing the most gene families "
                     "(red = same class, blue = cross-class)\n"
                     "OC = overlap coefficient, FE = fold enrichment", fontsize=9.5)
    figs = Path(args.figures); figs.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(figs / "08_h3_overlap.png", dpi=200, bbox_inches="tight")
    fig.savefig(figs / "08_h3_overlap.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  ✓ {figs}/08_h3_overlap.png (+pdf)")


if __name__ == "__main__":
    main()
