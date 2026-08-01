#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 07b — Feature stability + generalisation CV (out-of-core).

Quantifies how reproducible the top k-mer/unitig features are across resampling
AND gives an honest generalisation ROC-AUC mean±std (ROADMAP §0.1 M2 / §1.5).

DESIGN:
    Splits come from build_cv_splits():
      • PREFERRED — lineage-aware StratifiedGroupKFold (PopPUNK labels from 02c):
        an entire lineage stays on one side, so the AUC is NOT inflated by
        lineage leakage and stability is measured across lineage-disjoint folds.
      • FALLBACK — legacy 5-seed repeated holdout when no lineage file exists.
    For each split:
        - sample-level train/test masks (lineage-grouped or stratified holdout)
        - train XGBoost with the FIXED hyperparameters from step 04
          (config/experiments/{organism}/config_{antibiotic}.yaml) — HPO is done
          ONCE and held fixed across seeds; it never sees a seed's test split.
        - evaluate ROC-AUC on the held-out 20%
        - extract the top-N (analysis.top_n_features) Gain k-mers
    Aggregate:
        - ROC-AUC mean ± std across the 5 seeds
        - per-k-mer selection_frequency = (#seeds in top-N) / len(SEEDS)
          ("stable" >= 0.6) and mean Gain
        - mean pairwise Jaccard similarity of the 5 top-N k-mer sets

MEMORY MODEL: streamed, never materialised. A sample-level boolean mask selects
the train/test rows that fall inside each chunk's offset range, and training
reuses the same regime as 05 — standard full-data boosting over a streaming
(Ext)QuantileDMatrix (lib.xgb_data), one chunk read at a time.

Outputs:
    models/{organism}/{antibiotic}/seed{S}/xgboost_{antibiotic}_seed{S}.json
    results/{organism}/{antibiotic}/04_evaluation/10_repeated_holdout_summary_{antibiotic}.csv
    results/{organism}/{antibiotic}/05_explainability/06_feature_stability_{antibiotic}.csv
"""

import gc
import os
import shutil
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
import yaml
from scipy.sparse import load_npz
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit

from lib.config import load_config, resolve_path, env_bool, get_target
from lib.xgb_data import build_quantile_dmatrix, global_pos_weight

SEEDS = [42, 123, 777, 1024, 2025]
TEST_SIZE = 0.20

# AMR_CV_MODE=random runs the SAME 5-fold machinery with the lineage grouping removed
# (plain StratifiedKFold), for the random-vs-lineage-aware comparison reviewers expect.
# It is a measurement mode, not a pipeline mode: every artefact it writes is suffixed
# so it can never overwrite the canonical lineage-CV outputs the KB is populated from.
CV_MODE = os.environ.get("AMR_CV_MODE", "lineage").strip().lower()
RANDOM_CV = CV_MODE == "random"
OUT_SUFFIX = "_randomcv" if RANDOM_CV else ""

PROJECT_ROOT = Path(__file__).resolve().parent.parent
config = load_config()
TARGET_ANTIBIOTIC = get_target(config=config)[1]
ORGANISM = get_target(config=config)[0]
TOP_N = config['analysis']['top_n_features']
CHUNK_SIZE = config['preprocessing']['chunk_size']

MATRIX_DIR = resolve_path('matrix_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
MODELS_DIR = resolve_path('models_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
EVAL_DIR = resolve_path('dir_04_evaluation', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
EXPLAIN_DIR = resolve_path('dir_05_explainability', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)


def load_fixed_params():
    """Load the fixed hyperparameters tuned once in step 04 (no per-seed retune)."""
    cfg_path = resolve_path('experiment_config', organism=ORGANISM,
                            antibiotic=TARGET_ANTIBIOTIC, config=config)
    if not cfg_path.exists():
        print(f"ERROR: tuned config not found: {cfg_path}\n  Run 04_optimization.py first.")
        sys.exit(1)
    with open(cfg_path, 'r', encoding='utf-8') as f:
        ab_cfg = yaml.safe_load(f)

    base = {
        'objective': config['xgboost_params']['objective'],
        'eval_metric': config['xgboost_params'].get('eval_metric', 'auc'),
        'tree_method': config['xgboost_params']['tree_method'],
        'device': config['xgboost_params']['device'],
        'verbosity': 0,
        'max_bin': 2,
    }
    base.update(ab_cfg.get('xgboost_params', {}))
    best = dict(ab_cfg.get('best_params', {}))
    total_trees = max(1, int(best.pop('n_estimators', 100)))
    base.update(best)
    base.pop('scale_pos_weight', None)
    base.setdefault('base_score', 0.5)   # avoid pure-chunk base_score error (see 05)
    return base, total_trees


def chunk_offsets(chunk_files):
    """Return [(file, start, end), ...] reading ONLY each chunk's shape (no densify)."""
    offsets, cur = [], 0
    for f in chunk_files:
        with np.load(f) as d:               # CSR npz stores 'shape' without loading data
            rows = int(d['shape'][0])
        offsets.append((f, cur, cur + rows))
        cur += rows
    return offsets, cur


def train_one_seed(params, total_trees, chunk_files, y_all, train_mask, max_bin,
                   cache_dir, use_extmem):
    """Full-data boosting on this seed's train rows.

    Matches 05's regime: one DMatrix over all train-split rows (selected via the
    sample-level mask), global neg/pos class weighting, standard boosting. Built
    in-core when ``use_extmem`` is False (config training.external_memory), else
    as an ExtMemQuantileDMatrix spilling to ``cache_dir`` (removed after).
    """
    pos_weight = global_pos_weight(chunk_files, y_all, CHUNK_SIZE, row_mask=train_mask)
    cache_prefix = None
    if use_extmem:
        shutil.rmtree(cache_dir, ignore_errors=True)
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_prefix = str(cache_dir / "cache")
    dtrain = None
    try:
        dtrain = build_quantile_dmatrix(chunk_files, y_all, CHUNK_SIZE,
                                        max_bin=max_bin, row_mask=train_mask,
                                        pos_weight=pos_weight,
                                        cache_prefix=cache_prefix)
        model = xgb.train(params, dtrain, num_boost_round=total_trees)
    finally:
        dtrain = None              # release so the cache files unlock
        gc.collect()
        shutil.rmtree(cache_dir, ignore_errors=True)
    return model


def eval_one_seed(model, offsets, y_all, test_mask):
    """Stream the test-rows of each chunk and return ROC-AUC."""
    y_true, y_prob = [], []
    for f, start, end in offsets:
        local = test_mask[start:end]
        if not local.any():
            continue
        X = load_npz(f)[local]
        y_prob.extend(model.predict(xgb.DMatrix(X)))
        y_true.extend(y_all[start:end][local])
        del X
        gc.collect()
    y_true, y_prob = np.array(y_true), np.array(y_prob)
    if len(np.unique(y_true)) < 2:
        return np.nan, y_true, y_prob
    return roc_auc_score(y_true, y_prob), y_true, y_prob


def top_feature_indices(model):
    """Return (set of top-N feature indices, {index: gain}) from the model."""
    imp = model.get_score(importance_type='gain')
    if not imp:
        return set(), {}
    top = sorted(imp.items(), key=lambda kv: kv[1], reverse=True)[:TOP_N]
    idx_gain = {int(name[1:]): float(g) for name, g in top}
    return set(idx_gain), idx_gain


def map_indices_to_kmers(indices):
    """Map feature indices -> k-mer sequences via features.txt (single pass)."""
    features_file = MATRIX_DIR / "features.txt"
    mapping = {}
    if not indices or not features_file.exists():
        return mapping
    needed = set(indices)
    with open(features_file, 'r', encoding='utf-8') as f:
        for line_idx, line in enumerate(f):
            if line_idx in needed:
                mapping[line_idx] = line.split()[0]
                if len(mapping) == len(needed):
                    break
    return mapping


def mean_pairwise_jaccard(sets):
    """Mean Jaccard over all pairs; empty union -> 0.0 for that pair."""
    vals = []
    for a, b in combinations(sets, 2):
        union = a | b
        vals.append(len(a & b) / len(union) if union else 0.0)
    return float(np.mean(vals)) if vals else float('nan')


def build_cv_splits(y_all, n_total, genomes_csv, lineage_csv, n_splits, seed=42):
    """Build the resampling splits as sample-level (train_mask, test_mask) pairs.

    Lineage-aware by default (ROADMAP §0.1 M2): when PopPUNK labels exist
    (02c -> poppunk_clusters.csv), use **StratifiedGroupKFold** so an entire
    lineage stays on one side — the reported AUC mean±std is then an HONEST
    generalisation estimate (no lineage leakage) and the per-feature selection
    frequency is measured across lineage-disjoint folds. Falls back to the
    legacy 5-seed repeated holdout when no lineage file is present (e.g. the
    synthetic integration test), preserving prior behaviour.

    Returns (splits, method_label, split_labels).
    """
    if RANDOM_CV:
        # Deliberately lineage-BLIND: same n_splits, same stratification, same models —
        # the ONLY difference from the canonical path is that lineages are ignored, so
        # the AUC gap between the two is attributable to lineage leakage and nothing
        # else. This is the comparison a reviewer asks for; it is never the KB metric.
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        masks = []
        for tr_idx, te_idx in skf.split(np.zeros(n_total), y_all):
            tr = np.zeros(n_total, dtype=bool); tr[tr_idx] = True
            te = np.zeros(n_total, dtype=bool); te[te_idx] = True
            masks.append((tr, te))
        return masks, f"random_stratified_kfold_{n_splits}fold", list(range(len(masks)))

    if lineage_csv.exists() and genomes_csv.exists():
        try:
            from lib.lineage import load_lineage, group_kfold_masks, no_group_leakage
            groups = load_lineage(genomes_csv, lineage_csv)
            n_clusters = len(set(groups.tolist()))
            if len(groups) == n_total and n_clusters >= n_splits:
                masks = group_kfold_masks(y_all, groups, n_splits=n_splits,
                                          stratified=True, seed=seed)
                # watertight guard: no lineage may span train+test in any fold
                if not all(no_group_leakage(tr, te, groups) for tr, te in masks):
                    raise RuntimeError("lineage leakage detected in group-kfold masks")
                return masks, f"lineage_group_kfold_{n_splits}fold", list(range(len(masks)))
            print(f"  ⚠ lineage labels unusable (aligned {len(groups)} vs {n_total} rows, "
                  f"{n_clusters} clusters < {n_splits} folds); using 5-seed holdout.")
        except Exception as e:
            print(f"  ⚠ lineage CV unavailable ({e}); using 5-seed holdout.")

    masks = []
    for s in SEEDS:
        sss = StratifiedShuffleSplit(n_splits=1, test_size=TEST_SIZE, random_state=s)
        tr_idx, te_idx = next(sss.split(np.zeros(n_total), y_all))
        tr = np.zeros(n_total, dtype=bool); tr[tr_idx] = True
        te = np.zeros(n_total, dtype=bool); te[te_idx] = True
        masks.append((tr, te))
    return masks, "repeated_holdout_5seed", list(SEEDS)


def main():
    print("=" * 80)
    print(f"FEATURE STABILITY & GENERALISATION CV: {TARGET_ANTIBIOTIC.upper()} ({ORGANISM})")
    print("=" * 80)

    params, total_trees = load_fixed_params()

    y_path = MATRIX_DIR / f"y_{TARGET_ANTIBIOTIC}.csv"
    if not y_path.exists():
        print(f"ERROR: label file not found: {y_path}")
        sys.exit(1)
    y_all = pd.read_csv(y_path, encoding='utf-8')['label'].values.astype(int)

    chunk_files = sorted(MATRIX_DIR.glob(f"X_{TARGET_ANTIBIOTIC}_part_*.npz"),
                         key=lambda x: int(x.stem.split('_')[-1]))
    if not chunk_files:
        print(f"ERROR: no matrix chunks in {MATRIX_DIR}")
        sys.exit(1)

    offsets, n_total = chunk_offsets(chunk_files)
    if n_total != len(y_all):
        print(f"WARNING: matrix rows ({n_total}) != labels ({len(y_all)}); using min.")
    print(f"  Samples: {n_total} | chunks: {len(offsets)} | trees/seed: {total_trees} | top-N: {TOP_N}")

    max_bin = int(params.get('max_bin', 2))
    use_extmem = env_bool('AMR_EXTERNAL_MEMORY', config['training'].get('external_memory', True))
    print(f"  external_memory={use_extmem}")

    # Resampling scheme: lineage-aware StratifiedGroupKFold when PopPUNK labels
    # exist (honest, no lineage leakage), else the legacy 5-seed holdout.
    genomes_csv = MATRIX_DIR / f"genomes_{TARGET_ANTIBIOTIC}.csv"
    try:
        lineage_dir = resolve_path('lineage_dir', organism=ORGANISM, config=config)
    except KeyError:
        lineage_dir = resolve_path('data_dir', config=config) / "processed" / ORGANISM / "lineage"
    n_splits = int(config.get('lineage', {}).get('n_splits', 5))
    splits, cv_method, split_labels = build_cv_splits(
        y_all, n_total, genomes_csv, lineage_dir / "poppunk_clusters.csv", n_splits)
    print(f"  CV scheme: {cv_method} ({len(splits)} splits)")
    if RANDOM_CV:
        print("  " + "=" * 74)
        print("  MEASUREMENT MODE (AMR_CV_MODE=random): lineage-BLIND CV on purpose.")
        print("  The AUC produced here is the INFLATED comparison baseline; it is not")
        print(f"  the KB metric and is written to *{OUT_SUFFIX}.csv so the canonical")
        print("  lineage-CV outputs stay untouched. Do NOT populate from these files.")
        print("  " + "=" * 74)
    elif not cv_method.startswith("lineage_group_kfold"):
        print("  " + "!" * 74)
        print("  ⚠ WARNING: NO lineage-aware CV (PopPUNK clusters absent/unusable).")
        print("  ⚠ The reported AUC is a 5-seed holdout and may be lineage-INFLATED —")
        print("  ⚠ it is NOT a lineage-corrected metric. Run 02c PopPUNK for this")
        print(f"  ⚠ organism ({ORGANISM}) before trusting auc_mean_seeds as lineage-CV.")
        print("  " + "!" * 74)

    # Splits run sequentially; each split's full-data boosting already saturates
    # the allocated cores (one DMatrix over all train rows). The 'seed' column in
    # the summary is kept for backward compatibility (it holds the fold id here).
    seed_rows, seed_sets, gain_accum = [], [], {}
    for label, (train_mask, test_mask) in zip(split_labels, splits):
        print(f"\n--- SPLIT {label} ---")
        try:
            cache_dir = MODELS_DIR / f"_xgb_cache_split{label}"
            model = train_one_seed(params, total_trees, chunk_files, y_all,
                                   train_mask, max_bin, cache_dir, use_extmem)
            auc, _, _ = eval_one_seed(model, offsets, y_all, test_mask)
            idx_set, idx_gain = top_feature_indices(model)

            seed_rows.append({'seed': label, 'roc_auc': auc, 'n_top_features': len(idx_set)})
            seed_sets.append(idx_set)
            for idx, g in idx_gain.items():
                gain_accum.setdefault(idx, []).append(g)

            split_dir = MODELS_DIR / f"split{label}{OUT_SUFFIX}"
            split_dir.mkdir(parents=True, exist_ok=True)
            model.save_model(str(split_dir / f"xgboost_{TARGET_ANTIBIOTIC}_split{label}.json"))
            print(f"  ROC-AUC: {auc:.4f} | top features: {len(idx_set)}")
        except Exception as e:
            print(f"  ✗ split {label} failed: {e}")

    if not seed_rows:
        print("\nERROR: no seed completed; nothing to aggregate.")
        sys.exit(1)

    # ---- aggregate ----------------------------------------------------------
    aucs = np.array([r['roc_auc'] for r in seed_rows], dtype=float)
    auc_mean, auc_std = float(np.nanmean(aucs)), float(np.nanstd(aucs))
    jaccard = mean_pairwise_jaccard(seed_sets) if len(seed_sets) > 1 else float('nan')

    EVAL_DIR.mkdir(parents=True, exist_ok=True)
    EXPLAIN_DIR.mkdir(parents=True, exist_ok=True)

    summary = pd.DataFrame(seed_rows)
    summary = pd.concat([summary, pd.DataFrame([
        {'seed': 'MEAN', 'roc_auc': auc_mean, 'n_top_features': np.nan},
        {'seed': 'STD', 'roc_auc': auc_std, 'n_top_features': np.nan},
        {'seed': 'MEAN_JACCARD', 'roc_auc': jaccard, 'n_top_features': np.nan},
    ])], ignore_index=True)
    # Persist the CV scheme so downstream (populate/KB) can tell an HONEST
    # lineage-CV AUC from a fallback holdout — auc_mean_seeds is only "lineage-CV"
    # when cv_method == lineage_group_kfold_* (M3-3).
    summary['cv_method'] = cv_method
    summary_path = EVAL_DIR / f"10_repeated_holdout_summary_{TARGET_ANTIBIOTIC}{OUT_SUFFIX}.csv"
    summary.to_csv(summary_path, index=False)

    n_seeds = len(seed_sets)
    counts = Counter()
    for s in seed_sets:
        counts.update(s)
    kmer_map = map_indices_to_kmers(set(counts))
    stab_rows = []
    for idx, c in counts.items():
        stab_rows.append({
            'feature_index': idx,
            'kmer': kmer_map.get(idx, 'UNKNOWN'),
            'selection_frequency': c / n_seeds,
            'mean_gain': float(np.mean(gain_accum.get(idx, [0.0]))),
            'stable': (c / n_seeds) >= 0.6,
        })
    # Fix the columns explicitly so an empty stab_rows (degenerate models with
    # no splits -> no selected features) still yields a valid, header-only CSV
    # instead of crashing on sort_values of a column-less frame.
    stab = pd.DataFrame(stab_rows, columns=['feature_index', 'kmer',
                                            'selection_frequency', 'mean_gain', 'stable'])
    if not stab.empty:
        stab = stab.sort_values(['selection_frequency', 'mean_gain'], ascending=False)
    stab_path = EXPLAIN_DIR / f"06_feature_stability_{TARGET_ANTIBIOTIC}{OUT_SUFFIX}.csv"
    stab.to_csv(stab_path, index=False)

    print("\n" + "=" * 80)
    print(f"CROSS-VALIDATION SUMMARY ({cv_method})")
    print("=" * 80)
    print(f"  ROC-AUC: {auc_mean:.4f} ± {auc_std:.4f}  (splits: {[f'{a:.3f}' for a in aucs]})")
    print(f"  Mean pairwise Jaccard (top-{TOP_N} sets): {jaccard:.4f}")
    print(f"  Stable k-mers (freq ≥ 0.6): {int(stab['stable'].sum()) if len(stab) else 0}")
    print(f"  Saved: {summary_path.name}, {stab_path.name}")
    print("=" * 80)


if __name__ == "__main__":
    import traceback
    try:
        main()
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
