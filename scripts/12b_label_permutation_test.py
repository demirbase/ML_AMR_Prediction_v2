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
(no retuning — ROADMAP's frozen-HP variant), evaluate on the same held-out
split, and build the **null ROC-AUC distribution**. The real AUC's empirical
p-value is ``(1 + #{null >= real}) / (N + 1)``.

Method
------
* Single fixed train/test split = the experiment config's chunk split (the exact
  split 06 reports the headline AUC on), so the real AUC reproduces ~0.953.
* Frozen params + tree count from ``config_{ab}.yaml`` (n_estimators = 8); reuse
  07b's ``train_one_seed`` / ``eval_one_seed`` (same full-data boosting as 05).
* Each permutation reshuffles the ENTIRE label vector (so train+test stay
  consistently permuted) and retrains from scratch — this is the expensive,
  fully-honest null (HPC job; see the SLURM snippet in HANDOFF §0.2).

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

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import env_bool  # noqa: E402

# Reuse 07b's out-of-core training/eval (single source of truth for the boosting
# regime) and 06's exact held-out split. Digit-leading names need importlib.
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
    return offsets, n_total, ~test_mask, test_mask


def main():
    ap = argparse.ArgumentParser(description="Label-permutation null test (M9).")
    ap.add_argument("--n-permutations", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    N = args.n_permutations
    rng = np.random.default_rng(args.seed)

    organism, antibiotic = _s.ORGANISM, _s.TARGET_ANTIBIOTIC
    matrix_dir, models_dir = _s.MATRIX_DIR, _s.MODELS_DIR
    config = _s.config
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
    offsets, n_total, train_mask, test_mask = build_chunk_split(chunk_files, test_filenames)
    params, total_trees = _s.load_fixed_params()
    max_bin = int(params.get("max_bin", 2))
    use_extmem = env_bool("AMR_EXTERNAL_MEMORY", config["training"].get("external_memory", True))
    cache_dir = models_dir / "_xgb_cache_perm"

    print(f"  rows: total={n_total} train={int(train_mask.sum())} test={int(test_mask.sum())}")
    print(f"  trees={total_trees} | max_bin={max_bin} | external_memory={use_extmem}\n")

    def fit_eval(y_vec):
        model = _s.train_one_seed(params, total_trees, chunk_files, y_vec,
                                  train_mask, max_bin, cache_dir, use_extmem)
        auc, _, _ = _s.eval_one_seed(model, offsets, y_vec, test_mask)
        return auc

    # Real model (same pipeline) — should reproduce the ~0.953 headline.
    real_auc = fit_eval(y_all)
    print(f"  REAL test ROC-AUC = {real_auc:.4f}\n  running {N} label permutations...")

    null = np.empty(N)
    for r in range(N):
        y_perm = rng.permutation(y_all)            # shuffle ALL labels
        null[r] = fit_eval(y_perm)
        if (r + 1) % 10 == 0 or r == 0:
            print(f"    perm {r+1:>3}/{N}  null_auc={null[r]:.4f}  "
                  f"(running max={null[:r+1].max():.4f})")

    n_ge = int(np.sum(null >= real_auc))
    p_emp = (1 + n_ge) / (N + 1)
    summary = {
        "antibiotic": antibiotic, "organism": organism,
        "n_permutations": N, "seed": args.seed,
        "real_roc_auc": float(real_auc),
        "null_auc_mean": float(null.mean()), "null_auc_std": float(null.std()),
        "null_auc_min": float(null.min()), "null_auc_max": float(null.max()),
        "n_null_ge_real": n_ge,
        "empirical_p": p_emp,
        "significant": bool(p_emp < 0.05),
        "split_method": "experiment_config_chunk_split",
        "n_trees": int(total_trees),
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
