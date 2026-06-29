#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Label-permutation significance test — Step 12b  (ROADMAP §1.7 / must-have M9)

Why (complements 12's MDA)
--------------------------
12's per-feature MDA is conservative under unitig redundancy (many unitigs tag
the same ARG, so permuting one barely moves AUC). This test instead asks the
**model-level** question — *is the whole model's skill real, or could a model
this good arise by chance?* — via the classical permutation test (Ojala &
Garriga 2010): shuffle the labels, retrain with the SAME frozen hyper-parameters
(no retuning), evaluate on the same held-out split, and build the **null ROC-AUC
distribution**. The real AUC's empirical p-value is ``(1 + #{null>=real})/(N+1)``.

Key efficiency point
--------------------
The feature matrix is IDENTICAL across permutations — only the labels change. So
we build ONE in-core ``DMatrix`` for train and one for test, then each
permutation just resets ``set_label``/``set_weight`` and refits 8 trees. XGBoost
caches the (label-independent) histogram index on the DMatrix after the first
fit, so permutations 2..N are seconds each — no per-perm matrix rebuild (which
made the naive version I/O-bound, ~17 min/perm, and tripped the HPC low-eff
killer). max_bin=2 on binary features → identical binning to 05's QuantileDMatrix
→ the real AUC reproduces the ~0.953 headline.

Output (results/{org}/{ab}/05_explainability/)
    12b_label_permutation_summary_{ab}.json  — real_auc, null mean/std/max, empirical_p, significant
    12b_label_permutation_nulls_{ab}.csv     — per-permutation null AUC (for the thesis histogram)
"""

import argparse
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import xgboost as xgb
from scipy.sparse import load_npz, vstack
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

# Reuse 07b's frozen-HP loader / split helpers and 06's exact held-out split.
_s = importlib.import_module("07b_feature_stability")
_ev = importlib.import_module("06_evaluation")


def build_chunk_split(chunk_files, test_filenames):
    """Sample-level (train_mask, test_mask) from the experiment config's chunk
    split — the same held-out test set 06 reports the headline AUC on."""
    offsets, n_total = _s.chunk_offsets(chunk_files)
    test_names = set(test_filenames)
    test_mask = np.zeros(n_total, dtype=bool)
    for f, start, end in offsets:
        if f.name in test_names:
            test_mask[start:end] = True
    if not test_mask.any():
        raise RuntimeError("No chunk matched the config test_files — split empty.")
    return n_total, ~test_mask, test_mask


def main():
    ap = argparse.ArgumentParser(description="Label-permutation null test (M9).")
    ap.add_argument("--n-permutations", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    N = args.n_permutations
    rng = np.random.default_rng(args.seed)

    organism, antibiotic = _s.ORGANISM, _s.TARGET_ANTIBIOTIC
    matrix_dir, config = _s.MATRIX_DIR, _s.config
    out_dir = _s.resolve_path("dir_05_explainability", organism=organism,
                              antibiotic=antibiotic, config=config)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(f"LABEL-PERMUTATION NULL TEST  —  {organism} / {antibiotic}  |  N={N}")
    print("=" * 78)

    y_all = np.asarray(
        _s.pd.read_csv(matrix_dir / f"y_{antibiotic}.csv", encoding="utf-8")["label"].values
    ).astype(int)
    chunk_files = sorted(matrix_dir.glob(f"X_{antibiotic}_part_*.npz"),
                         key=lambda x: int(x.stem.split("_")[-1]))
    if not chunk_files:
        print(f"ERROR: no matrix chunks in {matrix_dir}"); sys.exit(1)

    test_filenames, _thr = _ev.load_test_files_from_config()
    n_total, train_mask, test_mask = build_chunk_split(chunk_files, test_filenames)

    params, total_trees = _s.load_fixed_params()
    max_bin = int(params.get("max_bin", 2))
    print(f"  rows: total={n_total} train={int(train_mask.sum())} test={int(test_mask.sum())}")
    print(f"  trees={total_trees} max_bin={max_bin} tree_method={params.get('tree_method')}")

    # Build the TRAIN matrix ONCE as a compact in-core QuantileDMatrix (max_bin=2
    # -> ~1 byte/nnz), streamed chunk-by-chunk so the raw ~4.9M-wide matrix is
    # never fully materialised in RAM (a plain CSV/CSR copy OOM-killed a 120G
    # node). QuantileDMatrix supports set_label/set_weight, so each permutation
    # just swaps labels and refits — no rebuild. pos_weight is applied per-perm
    # via set_weight, so build unweighted here.
    print("  building train QuantileDMatrix (once, streamed)...", flush=True)
    dtrain = _s.build_quantile_dmatrix(chunk_files, y_all, _s.CHUNK_SIZE,
                                       max_bin=max_bin, row_mask=train_mask,
                                       pos_weight=None)
    # Build the test matrix in ASCENDING chunk order so it aligns with
    # y[test_mask]. The experiment config's test_files can be listed in a
    # non-ascending order; using that order (load_test_data) would misalign the
    # per-row labels and collapse the AUC to ~0.5. Stream + select test rows.
    test_parts = []
    for f in chunk_files:
        cid = int(f.stem.split("_")[-1])
        X = load_npz(f).tocsr()
        local = test_mask[cid * _s.CHUNK_SIZE: cid * _s.CHUNK_SIZE + X.shape[0]]
        if local.any():
            test_parts.append(X[local])
    dtest = xgb.DMatrix(vstack(test_parts, format="csr"))
    assert dtest.num_row() == int(test_mask.sum()), "test row count mismatch"

    def fit_eval(y_vec):
        ytr, yte = y_vec[train_mask], y_vec[test_mask]
        pos = int(ytr.sum()); pw = (len(ytr) - pos) / pos if pos else 1.0
        dtrain.set_label(ytr)
        dtrain.set_weight(np.where(ytr == 1, pw, 1.0))
        model = xgb.train(params, dtrain, num_boost_round=total_trees)
        if len(np.unique(yte)) < 2:
            return np.nan
        return roc_auc_score(yte, model.predict(dtest))

    real_auc = fit_eval(y_all)
    print(f"\n  REAL test ROC-AUC = {real_auc:.4f}\n  running {N} label permutations...",
          flush=True)

    null = np.empty(N)
    for r in range(N):
        null[r] = fit_eval(rng.permutation(y_all))     # shuffle ALL labels
        if (r + 1) % 10 == 0 or r == 0:
            print(f"    perm {r+1:>3}/{N}  null_auc={null[r]:.4f}  "
                  f"(running max={null[:r+1].max():.4f})", flush=True)

    n_ge = int(np.sum(null >= real_auc))
    p_emp = (1 + n_ge) / (N + 1)
    summary = {
        "antibiotic": antibiotic, "organism": organism,
        "n_permutations": N, "seed": args.seed,
        "real_roc_auc": float(real_auc),
        "null_auc_mean": float(null.mean()), "null_auc_std": float(null.std()),
        "null_auc_min": float(null.min()), "null_auc_max": float(null.max()),
        "n_null_ge_real": n_ge, "empirical_p": p_emp,
        "significant": bool(p_emp < 0.05),
        "split_method": "experiment_config_chunk_split", "n_trees": int(total_trees),
    }
    _s.pd.DataFrame({"permutation": np.arange(1, N + 1), "null_roc_auc": null}).to_csv(
        out_dir / f"12b_label_permutation_nulls_{antibiotic}.csv", index=False, encoding="utf-8")
    (out_dir / f"12b_label_permutation_summary_{antibiotic}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 78)
    print(f"  REAL AUC={real_auc:.4f} | null mean={null.mean():.4f} "
          f"max={null.max():.4f} | p={p_emp:.4g} "
          f"({'SIGNIFICANT' if p_emp < 0.05 else 'n.s.'})")
    print(f"  ✓ 12b_label_permutation_summary_{antibiotic}.json")
    print(f"  ✓ 12b_label_permutation_nulls_{antibiotic}.csv")
    print("=" * 78)


if __name__ == "__main__":
    main()
