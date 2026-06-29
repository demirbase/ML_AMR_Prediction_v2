#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stability selection (CPSS) + SHAP — Step 13  (ROADMAP §0.1)

Why
---
The per-feature MDA (step 12) is flat because unitigs are *redundant*: a single
resistance gene spawns many overlapping unitigs, so no one of them is individually
indispensable. **Complementary Pairs Stability Selection (CPSS)** answers the right
question instead — *which unitigs are selected reproducibly across many random
subsamples?* — and comes with a finite-sample false-discovery (PFER) bound
(Meinshausen & Bühlmann 2010; Shah & Samworth 2013). This is the thesis-grade
replacement for ad-hoc seed stability.

Pipeline (staged, so CPSS runs on a tractable feature set)
----------------------------------------------------------
1. **Chi² prefilter** — one streamed pass over the full ~4.9M-unitig matrix scores
   each unitig's association with the R/S label (2x2 chi-square); keep the top
   ``--n-candidates`` (default 5000). Zero-variance / ubiquitous unitigs drop out.
2. **CPSS** — ``B`` complementary pairs (default 100): each draws a random 50%
   subsample and its complement (2B subsamples total); fit the frozen-HP base
   learner (8 trees, from the experiment config) on each and record which
   candidates it *selects* (appears in a tree split, gain>0). Selection frequency
   = (# subsamples selecting it) / 2B.
3. **Stable set** — selection frequency >= ``--pi`` (default 0.6); report the
   Meinshausen-Bühlmann **PFER bound**.
4. **SHAP** — XGBoost's built-in TreeSHAP (``pred_contribs`` — no `shap` package)
   on a final model over the candidates → mean |SHAP| per unitig (directional
   importance, replaces raw Gain).

Output (results/{org}/{ab}/05_explainability/)
    13_stability_selection_{ab}.csv  — per-candidate chi2, selection_freq, stable, mean_abs_shap, kmer
    13_stability_summary_{ab}.json   — params, n_stable, PFER bound, top unitigs
"""

import argparse
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import xgboost as xgb
from scipy.sparse import load_npz, vstack

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

# Reuse 07b's frozen-HP loader, matrix paths and features.txt mapper.
_s = importlib.import_module("07b_feature_stability")


def chi2_prefilter(chunk_files, y_all, chunk_size, top_k):
    """One streamed pass: per-unitig 2x2 chi-square vs the R/S label.

    Returns (top_idx, chi2_top) — the global column indices of the top_k unitigs
    by chi-square and their scores. Streams one chunk at a time (never densifies
    the full matrix)."""
    n_feat = None
    present_R = present_S = None
    nR = nS = 0
    for f in chunk_files:
        cid = int(f.stem.split("_")[-1])
        X = load_npz(f).tocsr()
        start = cid * chunk_size
        yc = y_all[start:start + X.shape[0]]
        if n_feat is None:
            n_feat = X.shape[1]
            present_R = np.zeros(n_feat, dtype=np.float64)
            present_S = np.zeros(n_feat, dtype=np.float64)
        rmask = yc == 1
        present_R += np.asarray(X[rmask].sum(axis=0)).ravel()
        present_S += np.asarray(X[~rmask].sum(axis=0)).ravel()
        nR += int(rmask.sum()); nS += int((~rmask).sum())
        del X
    N = nR + nS
    a, b = present_R, present_S          # present & R, present & S
    c, d = nR - a, nS - b                # absent & R, absent & S
    # chi-square with a tiny denom guard; zero-variance columns -> 0.
    denom = (a + b) * (c + d) * (a + c) * (b + d)
    chi2 = np.where(denom > 0, N * (a * d - b * c) ** 2 / denom, 0.0)
    k = min(top_k, n_feat)
    top_idx = np.argpartition(chi2, -k)[-k:]
    top_idx = top_idx[np.argsort(chi2[top_idx])[::-1]]   # high -> low
    return top_idx, chi2[top_idx]


def load_candidate_matrix(chunk_files, chunk_size, top_idx):
    """Stream the chunks and keep only the candidate columns -> X_k (csr)."""
    cols = np.asarray(top_idx)
    parts = []
    for f in sorted(chunk_files, key=lambda x: int(x.stem.split("_")[-1])):
        parts.append(load_npz(f).tocsr()[:, cols])
    return vstack(parts, format="csr")


def cpss(X_k, y, params, total_trees, B, rng):
    """Complementary-pairs stability selection over the candidate columns.

    Returns (selection_freq[K], avg_selected). For each of B pairs we fit the base
    learner on a random 50% subsample AND its complement (2B fits); a candidate is
    'selected' by a fit if it appears in a tree split (gain>0)."""
    n, K = X_k.shape
    counts = np.zeros(K, dtype=np.int64)
    n_selected = []
    half = n // 2
    for b in range(B):
        perm = rng.permutation(n)
        for rows in (perm[:half], perm[half:2 * half]):  # complementary pairs
            yr = y[rows]
            if len(np.unique(yr)) < 2:
                continue
            pos = int(yr.sum()); pw = (len(yr) - pos) / pos if pos else 1.0
            d = xgb.DMatrix(X_k[rows], label=yr,
                            weight=np.where(yr == 1, pw, 1.0))
            model = xgb.train(params, d, num_boost_round=total_trees)
            sel = {int(k[1:]) for k in model.get_score(importance_type="gain")}
            n_selected.append(len(sel))
            for j in sel:
                counts[j] += 1
    n_fits = 2 * B
    return counts / n_fits, (float(np.mean(n_selected)) if n_selected else 0.0)


def shap_importance(X_k, y, params, total_trees):
    """Mean |TreeSHAP| per candidate from a final model (XGBoost built-in)."""
    pos = int(y.sum()); pw = (len(y) - pos) / pos if pos else 1.0
    dtrain = xgb.DMatrix(X_k, label=y, weight=np.where(y == 1, pw, 1.0))
    model = xgb.train(params, dtrain, num_boost_round=total_trees)
    contribs = model.predict(xgb.DMatrix(X_k), pred_contribs=True)  # (n, K+1)
    return np.abs(contribs[:, :-1]).mean(axis=0)                    # drop bias col


def pfer_bound(avg_selected, n_candidates, pi):
    """Meinshausen-Bühlmann upper bound on the expected number of false positives:
    PFER <= q^2 / ((2*pi - 1) * p). Valid for pi > 0.5."""
    if pi <= 0.5:
        return float("nan")
    return (avg_selected ** 2) / ((2 * pi - 1) * n_candidates)


def main():
    ap = argparse.ArgumentParser(description="CPSS stability selection + SHAP.")
    ap.add_argument("--n-candidates", type=int, default=5000, help="Chi² prefilter top-K")
    ap.add_argument("--B", type=int, default=100, help="complementary pairs (2B fits)")
    ap.add_argument("--pi", type=float, default=0.6, help="stability threshold")
    ap.add_argument("--base-trees", type=int, default=10,
                    help="trees in the CPSS base selector — kept small/sparse and "
                         "decoupled from the final model so few features are picked "
                         "per fit, which is what makes the PFER bound tight")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    organism, antibiotic = _s.ORGANISM, _s.TARGET_ANTIBIOTIC
    matrix_dir, config = _s.MATRIX_DIR, _s.config
    chunk_size = _s.CHUNK_SIZE
    out_dir = _s.resolve_path("dir_05_explainability", organism=organism,
                              antibiotic=antibiotic, config=config)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(f"STABILITY SELECTION (CPSS) + SHAP  —  {organism} / {antibiotic}")
    print(f"  K={args.n_candidates}  B={args.B} (2B fits)  pi={args.pi}  base_trees={args.base_trees}")
    print("=" * 78)

    y_all = np.asarray(
        _s.pd.read_csv(matrix_dir / f"y_{antibiotic}.csv", encoding="utf-8")["label"].values
    ).astype(int)
    chunk_files = sorted(matrix_dir.glob(f"X_{antibiotic}_part_*.npz"),
                         key=lambda x: int(x.stem.split("_")[-1]))
    if not chunk_files:
        print(f"ERROR: no matrix chunks in {matrix_dir}"); sys.exit(1)

    params, total_trees = _s.load_fixed_params()
    params = dict(params); params.setdefault("tree_method", "hist"); params["max_bin"] = 2

    print("  [1/4] Chi² prefilter (streamed)...", flush=True)
    top_idx, chi2_top = chi2_prefilter(chunk_files, y_all, chunk_size, args.n_candidates)
    print(f"        kept {len(top_idx)} candidates (chi2 max={chi2_top[0]:.1f}, "
          f"min={chi2_top[-1]:.1f})")

    print("  [2/4] loading candidate matrix...", flush=True)
    X_k = load_candidate_matrix(chunk_files, chunk_size, top_idx)

    print(f"  [3/4] CPSS — {2*args.B} fits (base selector: {args.base_trees} trees)...", flush=True)
    sel_freq, avg_sel = cpss(X_k, y_all, params, args.base_trees, args.B, rng)

    print("  [4/4] SHAP (TreeSHAP, built-in)...", flush=True)
    shap_imp = shap_importance(X_k, y_all, params, total_trees)

    kmers = _s.map_indices_to_kmers(set(int(i) for i in top_idx))
    stable = sel_freq >= args.pi
    pfer = pfer_bound(avg_sel, len(top_idx), args.pi)

    df = _s.pd.DataFrame({
        "feature_index": top_idx,
        "kmer": [kmers.get(int(i), "") for i in top_idx],
        "chi2": chi2_top,
        "selection_frequency": sel_freq,
        "stable": stable.astype(int),
        "mean_abs_shap": shap_imp,
    }).sort_values(["stable", "selection_frequency", "mean_abs_shap"],
                   ascending=[False, False, False])
    csv_path = out_dir / f"13_stability_selection_{antibiotic}.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8")

    # Emit the stable set as FASTA so it can be BLAST-validated against CARD
    # (biological closure of the statistically-stable selection).
    stable_df = df[(df["stable"] == 1) & (df["kmer"].astype(str).str.len() > 0)]
    fasta_path = out_dir / f"13_stable_features_{antibiotic}.fasta"
    with open(fasta_path, "w", encoding="utf-8") as fh:
        for rank, (_, r) in enumerate(stable_df.iterrows(), 1):
            fh.write(f">stable_{rank}|freq_{r['selection_frequency']:.2f}|"
                     f"shap_{r['mean_abs_shap']:.4g}|fidx_{int(r['feature_index'])}\n")
            fh.write(f"{r['kmer']}\n")

    n_stable = int(stable.sum())
    summary = {
        "antibiotic": antibiotic, "organism": organism,
        "n_candidates": int(len(top_idx)), "B": args.B, "n_fits": 2 * args.B,
        "pi": args.pi, "seed": args.seed,
        "n_stable": n_stable,
        "avg_selected_per_fit": avg_sel,
        "pfer_bound": pfer,
        "base_trees": args.base_trees,
        "n_trees": int(total_trees),
        "top_stable": df[df["stable"] == 1].head(15)[
            ["kmer", "selection_frequency", "mean_abs_shap"]].to_dict("records"),
    }
    (out_dir / f"13_stability_summary_{antibiotic}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 78)
    print(f"  stable unitigs (freq>={args.pi}): {n_stable}/{len(top_idx)} "
          f"| avg selected/fit={avg_sel:.1f} | PFER bound={pfer:.3f}")
    print(f"  ✓ {csv_path.name}\n  ✓ 13_stability_summary_{antibiotic}.json"
          f"\n  ✓ {fasta_path.name} ({n_stable} stable unitigs for BLAST)")
    print("=" * 78)


if __name__ == "__main__":
    main()
