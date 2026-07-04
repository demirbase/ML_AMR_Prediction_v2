#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 10 — k-mer background frequency / discriminativeness (ROADMAP §1.1).

A k-mer can be stable, high-Gain AND map to a known ARG, yet still be a poor
marker if it is present in (nearly) every genome — resistant and susceptible
alike. BLAST tells us *what gene* a k-mer belongs to; it does NOT tell us
whether the k-mer *discriminates* resistance. This step closes that gap.

For each candidate k-mer (07_kb_candidates: gain top-N ∪ stable set) it streams
the out-of-core matrix once and computes:
    prevalence_resistant   = present in R genomes / n_R
    prevalence_susceptible = present in S genomes / n_S
    prevalence_overall
    odds_ratio + Fisher's exact p (2×2 present/absent × R/S)
    discriminative flag      (|Δprevalence| >= min_delta AND fisher_p < alpha)
    fisher_q + discriminative_fdr   (Benjamini-Hochberg FDR over the candidate
                                     set; q < alpha — ROADMAP §0.2)

This distinguishes a genuine resistance marker from a ubiquitous / lineage
(conserved) sequence and feeds the KB `kmer_background_frequency` record. A
k-mer that is "confirmed" by BLAST but NOT discriminative is flagged — it is
likely a wildtype gene region or a clonal-lineage signal, not a resistance
determinant (cf. the gyrA caveat in step 09).

Output:
    results/{org}/{ab}/05_explainability/10_kmer_background_frequency_{ab}.csv
        = 07_kb_candidates columns + the prevalence/discriminativeness columns
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import load_npz
from scipy.stats import fisher_exact

from lib.config import load_config, resolve_path, get_target

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Discriminativeness thresholds (citeable in Methods)
MIN_PREVALENCE_DELTA = 0.10   # |prev_R - prev_S| must be at least this
FISHER_ALPHA = 0.05           # Fisher's exact significance (raw p AND BH-FDR q)


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg FDR-adjusted q-values (ROADMAP §0.2 — multiple-testing
    correction over the candidate set). NaN p-values pass through as NaN and are
    excluded from the ranking. Returns an array aligned to ``pvals``.
    """
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    finite = np.where(np.isfinite(p))[0]
    if finite.size == 0:
        return q
    pf = p[finite]
    order = np.argsort(pf)
    ranked = pf[order]
    m = pf.size
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]   # enforce monotonicity
    qf = np.empty(m)
    qf[order] = np.clip(adj, 0, 1)
    q[finite] = qf
    return q


def compute_kmer_stats(present_r, n_r, present_s, n_s,
                       min_delta=MIN_PREVALENCE_DELTA, alpha=FISHER_ALPHA):
    """
    Pure function: prevalence + discriminativeness for one k-mer.

    Args are presence counts and group sizes for resistant (r) / susceptible (s)
    genomes. Returns a dict of prevalences, odds ratio, Fisher p and a
    discriminative flag. Unit-testable without any matrix I/O.
    """
    present_r, n_r = int(present_r), int(n_r)
    present_s, n_s = int(present_s), int(n_s)
    absent_r, absent_s = n_r - present_r, n_s - present_s

    prev_r = present_r / n_r if n_r else float('nan')
    prev_s = present_s / n_s if n_s else float('nan')
    prev_all = ((present_r + present_s) / (n_r + n_s)) if (n_r + n_s) else float('nan')

    # 2×2 contingency: rows = present/absent, cols = R/S
    try:
        odds_ratio, fisher_p = fisher_exact([[present_r, present_s],
                                             [absent_r, absent_s]])
    except ValueError:
        odds_ratio, fisher_p = float('nan'), float('nan')

    delta = abs(prev_r - prev_s) if (prev_r == prev_r and prev_s == prev_s) else float('nan')
    discriminative = bool(delta == delta and delta >= min_delta
                          and fisher_p == fisher_p and fisher_p < alpha)

    return {
        'present_resistant': present_r, 'n_resistant': n_r,
        'prevalence_resistant': prev_r,
        'present_susceptible': present_s, 'n_susceptible': n_s,
        'prevalence_susceptible': prev_s,
        'prevalence_overall': prev_all,
        'odds_ratio': float(odds_ratio),
        'fisher_p': float(fisher_p),
        'enriched_in': ('resistant' if (prev_r == prev_r and prev_s == prev_s and prev_r > prev_s)
                        else 'susceptible' if (prev_r == prev_r and prev_s == prev_s and prev_s > prev_r)
                        else 'equal'),
        'discriminative': discriminative,
    }


def count_presence_by_label(indices, chunk_files, y_all):
    """
    Stream the matrix chunks once and count, for each feature index, how many
    resistant / susceptible genomes carry that k-mer (presence = value > 0).

    Returns {feature_index: (present_R, present_S)}. One chunk in RAM at a time.
    """
    idx_list = sorted(set(int(i) for i in indices))
    pres_r = {i: 0 for i in idx_list}
    pres_s = {i: 0 for i in idx_list}
    row = 0
    for f in chunk_files:
        Xc = load_npz(f).tocsc()
        n = Xc.shape[0]
        y = y_all[row:row + n]
        # Slice only the candidate columns -> dense (n_rows × n_candidates)
        sub = (Xc[:, idx_list].toarray() > 0)
        r_mask = (y == 1)
        s_mask = (y == 0)
        col_r = sub[r_mask].sum(axis=0)
        col_s = sub[s_mask].sum(axis=0)
        for j, i in enumerate(idx_list):
            pres_r[i] += int(col_r[j])
            pres_s[i] += int(col_s[j])
        row += n
        del Xc, sub
    return {i: (pres_r[i], pres_s[i]) for i in idx_list}


def main():
    config = load_config()
    antibiotic = get_target(config=config)[1]
    organism = get_target(config=config)[0]

    matrix_dir = resolve_path('matrix_dir', organism=organism, antibiotic=antibiotic, config=config)
    explain_dir = resolve_path('dir_05_explainability', organism=organism,
                               antibiotic=antibiotic, config=config)

    kb_path = explain_dir / f"07_kb_candidates_{antibiotic}.csv"
    if not kb_path.exists():
        print(f"ERROR: KB candidates not found: {kb_path}\n  Run 09_biological_summary.py first.")
        sys.exit(1)
    kb = pd.read_csv(kb_path)

    print("=" * 80)
    print(f"K-MER BACKGROUND FREQUENCY / DISCRIMINATIVENESS: {antibiotic.upper()} ({organism})")
    print("=" * 80)

    if kb.empty:
        out_path = explain_dir / f"10_kmer_background_frequency_{antibiotic}.csv"
        kb.to_csv(out_path, index=False)
        print(f"  No candidate k-mers; wrote empty {out_path.name}.")
        return

    # feature index from 'feature_id' (e.g. "f19862101")
    # drop the 'f' prefix (feature ids are 'f<index>'); str[1:] not lstrip('f')
    # which would also strip any further leading 'f' chars (audit Issue 25).
    kb['feature_index'] = kb['feature_id'].astype(str).str[1:].astype(int)

    y_path = matrix_dir / f"y_{antibiotic}.csv"
    chunk_files = sorted(matrix_dir.glob(f"X_{antibiotic}_part_*.npz"),
                         key=lambda x: int(x.stem.split('_')[-1]))
    if not y_path.exists() or not chunk_files:
        print(f"ERROR: matrix/labels missing in {matrix_dir}")
        sys.exit(1)
    y_all = pd.read_csv(y_path)['label'].values.astype(int)
    n_r, n_s = int((y_all == 1).sum()), int((y_all == 0).sum())
    print(f"  Genomes: {len(y_all)} (resistant={n_r}, susceptible={n_s}) | "
          f"candidate k-mers: {len(kb)}")

    presence = count_presence_by_label(kb['feature_index'].tolist(), chunk_files, y_all)

    stat_rows = []
    for _, r in kb.iterrows():
        pr, ps = presence.get(int(r['feature_index']), (0, 0))
        stat_rows.append(compute_kmer_stats(pr, n_r, ps, n_s))
    stats = pd.DataFrame(stat_rows)
    # BH-FDR across the candidate set (ROADMAP §0.2): adds fisher_q + an
    # FDR-corrected discriminative flag (|Δprev| >= delta AND q < alpha). The raw
    # 'discriminative' (per-test p) column is kept for backward compatibility.
    stats['fisher_q'] = benjamini_hochberg(stats['fisher_p'].values)
    _delta = (stats['prevalence_resistant'] - stats['prevalence_susceptible']).abs()
    stats['discriminative_fdr'] = (
        (_delta >= MIN_PREVALENCE_DELTA) & (stats['fisher_q'] < FISHER_ALPHA)
    ).fillna(False).astype(bool)
    out = pd.concat([kb.reset_index(drop=True), stats], axis=1)

    out_path = explain_dir / f"10_kmer_background_frequency_{antibiotic}.csv"
    out.to_csv(out_path, index=False)

    # ---- summary -----------------------------------------------------------
    stable = out[out['stable']] if 'stable' in out else out
    n_disc = int(out['discriminative'].sum())
    n_disc_fdr = int(out['discriminative_fdr'].sum())
    n_stable_disc = int(stable['discriminative'].sum()) if len(stable) else 0
    # confirmed by BLAST but NOT discriminative -> likely wildtype/lineage
    if 'confidence_tier' in out:
        conf_not_disc = out[(out['confidence_tier'] == 'confirmed') & (~out['discriminative'])]
    else:
        conf_not_disc = out.iloc[0:0]

    print(f"  Discriminative (|Δprev|≥{MIN_PREVALENCE_DELTA:g}, Fisher p<{FISHER_ALPHA:g}): "
          f"{n_disc}/{len(out)}  (stable: {n_stable_disc}/{len(stable)})")
    print(f"  Discriminative after BH-FDR (q<{FISHER_ALPHA:g}): {n_disc_fdr}/{len(out)}")
    if len(conf_not_disc):
        print(f"  ⚠ {len(conf_not_disc)} CONFIRMED-by-BLAST k-mer(s) are NOT discriminative "
              f"(ubiquitous / likely wildtype or lineage signal):")
        for _, r in conf_not_disc.iterrows():
            print(f"      {r['kmer']} [{r.get('card_gene','')}] "
                  f"prev_R={r['prevalence_resistant']:.2f} prev_S={r['prevalence_susceptible']:.2f}")
    print(f"  Saved: {out_path.name}")
    print("=" * 80)


if __name__ == "__main__":
    main()
