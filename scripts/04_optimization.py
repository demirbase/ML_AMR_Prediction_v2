#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
XGBoost Hyperparameter Optimization Module for AMR Prediction

This script performs automated hyperparameter tuning using Optuna to optimize
XGBoost model performance for antimicrobial resistance (AMR) prediction.

Key Features:
    1. Smart Chunk Selection: Uses a subset of training data to reduce RAM usage
       while maintaining representative samples for hyperparameter optimization
    2. Bayesian Optimization: Employs Optuna's Tree-structured Parzen Estimator (TPE)
       to efficiently explore the hyperparameter space
    3. Early Stopping: Prevents overfitting during trial evaluations
    4. Cross-platform Compatibility: Uses pathlib for universal path handling

Smart Chunk Selection Strategy:
    Instead of loading all training data (which may exceed RAM), we select N chunks
    (optuna_chunk_count from config.yaml) using a stratified linspace approach:
    - Chunks are sorted from highest to lowest resistance ratio
    - N chunks are selected at evenly spaced intervals to represent the entire spectrum
    
    This strategy ensures:
    - Comprehensive representation of the resistance spectrum (hard, medium, easy)
    - Prevents RAM exhaustion by limiting optimization dataset size
    - Faster iteration during hyperparameter search

Hyperparameter Search Space:
    - learning_rate: [0.01, 0.3] - Controls gradient descent step size
    - max_depth: [3, 12] - Maximum tree depth to prevent overfitting
    - subsample: [0.6, 1.0] - Fraction of samples used per tree
    - colsample_bytree: log-scale window around 1/√p (the square-root heuristic
      from METHODOLOGY.md), derived dynamically from the feature count
    - min_child_weight: [1, 10] - Minimum sum of instance weight in child node
    - gamma: [0, 5] - Minimum loss reduction for further partition
    - n_estimators: determined by early stopping (not tuned directly)

    Class imbalance is handled in the final training stage (05) via Dynamic
    Instance Weighting, not via scale_pos_weight here (avoids double-counting).
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import pandas as pd
import numpy as np
import xgboost as xgb
import optuna
import yaml
import joblib
from pathlib import Path
from scipy.sparse import load_npz, vstack
from sklearn.model_selection import train_test_split
import os
import sys
import datetime
import shutil

# Shared label-slicing helper (single source of truth) — see scripts/utils.py
from utils import get_y_chunk
# MLOps run provenance (SCALE_MLOPS_PLAN.md §7.1) — additive, best-effort.
from lib import run_metadata as rm
from lib.config import resolve_path, get_target


# ============================================================================
# LOAD CONFIGURATION FROM YAML
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_DIR = PROJECT_ROOT / "config"
CONFIG_FILE = CONFIG_DIR / "config.yaml"

# Load configuration
if not CONFIG_FILE.exists():
    raise FileNotFoundError(
        f"Configuration file not found: {CONFIG_FILE}\n"
        f"Please ensure config.yaml exists in the config/ directory."
    )

with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

# Extract configuration values
TARGET_ANTIBIOTIC = get_target(config=config)[1]
CHUNK_SIZE = config['preprocessing']['chunk_size']

# Load fractional splitting parameters (fallback to defaults if undefined)
TEST_FRACTION = float(config['training'].get('test_fraction', 0.20))
OPTUNA_FRACTION = float(config['training'].get('optuna_fraction', 0.25))

RANDOM_SEED = config['training']['random_seed']
N_TRIALS = config['training']['n_trials']
# Optuna early stopping: stop once this many COMPLETED trials pass without
# improving the best value (TPE usually converges well before n_trials). null/0
# disables it (run all n_trials). Env AMR_OPTUNA_PATIENCE overrides config.
OPTUNA_PATIENCE = os.environ.get('AMR_OPTUNA_PATIENCE',
                                 config['training'].get('optuna_patience'))
OPTUNA_PATIENCE = int(OPTUNA_PATIENCE) if OPTUNA_PATIENCE not in (None, '') else 0
VAL_SPLIT = config['training']['validation_fraction']
EARLY_STOPPING_ROUNDS = config['xgboost_params']['early_stopping_rounds']

# Parallel HPO. A single trial on the small Optuna subset (a few hundred genomes)
# cannot saturate a many-core HPC allocation → very low SLURM CPU efficiency
# (TRUBA flags/kills such jobs). XGBoost releases the GIL during training, so we
# instead run several trials CONCURRENTLY (Optuna n_jobs threads), each XGBoost
# using only THREADS_PER_TRIAL threads — keeping all allocated cores busy.
TOTAL_CORES = int(os.environ.get('SLURM_CPUS_PER_TASK') or os.cpu_count() or 1)
THREADS_PER_TRIAL = max(1, int(config['training'].get('optuna_threads_per_trial', 2)))
OPTUNA_N_JOBS = max(1, min(N_TRIALS, TOTAL_CORES // THREADS_PER_TRIAL))

# XGBoost base parameters (fetched from config)
BASE_PARAMS = {
    'objective': config['xgboost_params'].get('objective', 'binary:logistic'),
    # Read the evaluation metric from config (was hardcoded 'aucpr', which
    # silently overrode config.yaml and mislabelled the reported score as
    # "ROC AUC"). With config eval_metric: "auc", study.best_value is a genuine
    # ROC AUC, consistent with the metadata field `best_auc_score`.
    'eval_metric': config['xgboost_params'].get('eval_metric', 'auc'),
    'tree_method': config['xgboost_params'].get('tree_method', 'hist'),
    # Threads PER trial — the remaining cores are used by running OPTUNA_N_JOBS
    # trials concurrently (see study.optimize(n_jobs=...) below).
    'nthread': THREADS_PER_TRIAL,
    'device': config['xgboost_params'].get('device', 'cpu'),
    'random_state': RANDOM_SEED,
    'verbosity': config['xgboost_params'].get('verbosity', 0),
    'max_bin': 2,  # CRITICAL RAM FIX: 0/1 binary data does not need 256 bins!
    # Pin base_score=0.5 (matches training stage 05): without it XGBoost derives
    # the intercept from the label mean, which is 0 or 1 for a pure-class chunk
    # and raises "base_score must be in (0,1)". Surfaces on small chunks.
    'base_score': 0.5,
}

# Organism-aware paths (SCALE_MLOPS_PLAN §4.2)
ORGANISM = get_target(config=config)[0]
MATRIX_DIR = resolve_path('matrix_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)

# Output directories (organism + antibiotic-specific)
MODELS_DIR = resolve_path('models_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
LOGS_DIR = resolve_path('logs_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)

# Create output directories
MODELS_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
# get_y_chunk() is imported from utils.py (shared with 05 and 06).


def count_features(matrix_dir):
    """
    Count the number of k-mer features (lines in features.txt).

    Used to derive the √p column-subsampling range so the Optuna search space
    is anchored to the actual feature dimensionality (see objective()).

    Returns:
        int: number of features, or 0 if features.txt is unavailable.
    """
    features_file = matrix_dir / "features.txt"
    if not features_file.exists():
        return 0
    with open(features_file, 'r', encoding='utf-8') as f:
        return sum(1 for _ in f)


def compute_colsample_range(n_features):
    """
    Derive a colsample_bytree search range centred on the √p heuristic.

    METHODOLOGY.md motivates colsample_bytree ≈ 1/√p (the random-forest square-
    root heuristic). Rather than hardcoding a fixed [0.05, 0.30] window — which
    is ~100x larger than 1/√p for p≈5M and invalidated the stated methodology —
    we build a log-scale window bracketing 1/√p so Optuna explores values
    consistent with the theory while retaining flexibility.

    Returns:
        tuple(float, float): (lower, upper) bounds for trial.suggest_float(log=True).
    """
    if n_features and n_features > 1:
        sqrt_p_ratio = 1.0 / np.sqrt(n_features)          # e.g. ≈ 4.5e-4 at p≈5M
        col_lower = max(sqrt_p_ratio * 0.5, 1e-5)
        col_upper = min(max(sqrt_p_ratio * 20.0, 1e-2), 1.0)
        return float(col_lower), float(col_upper)
    # Fallback when feature count is unknown: a conservative low-fraction range.
    return 1e-4, 1e-1


def analyze_and_stratify_all_chunks(y_all, all_files, chunk_size, test_fraction, optuna_fraction):
    """
    Analyze all chunks and perform stratified train/test split using resistance ratio.
    
    Dynamic Stratified Splitting Strategy (Scale-Invariant Fractional Method):
        This function implements a robust, data-driven splitting approach that ensures
        representative sampling across the resistance spectrum. Unlike static chunk counts,
        using fractions prevents crashes or statistical anomalies when chunk_size changes.
    
    Stratification Algorithm:
        Step 1: Calculate resistance ratio for each chunk (resistant / total)
        Step 2: Sort chunks from highest to lowest resistance (hard → easy)
        Step 3: Calculate dynamic chunk counts: N_test = max(1, floor(total * test_fraction))
        Step 4: Select N_test chunks using np.linspace for representative sampling
                - Ensures test set includes hard, medium, and easy cases
                - Maintains population distribution fidelity
        Step 5: Remaining chunks form training pool
        Step 6: Select Optuna subset (N_optuna = max(1, floor(train_total * optuna_fraction))) 
                from training pool for evaluation
    
    Why Fractional Stratification Matters:
        - Prevents bias: Random splits may over-represent easy/hard cases
        - Generalization: Model tested on full difficulty spectrum
        - Clinical validity: Real-world deployment faces all resistance levels
        - Scale-Invariant: Works seamlessly whether there are 10 chunks or 10,000 chunks
    
    Args:
        y_all: Complete label array for all samples
        all_files: List of ALL available chunk file paths
        chunk_size: Number of samples per chunk
        test_fraction: Fraction of all chunks to reserve for testing (e.g., 0.20)
        optuna_fraction: Fraction of training chunks to allocate for Optuna optimization (e.g., 0.25)
    
    Returns:
        tuple: (train_files, train_filenames, test_files, test_filenames, optuna_files, optuna_filenames)
            - train_files: Chunks for model training (excluding test)
            - test_files: Representative test set (stratified)
            - optuna_files: Subset for hyperparameter optimization (N representative chunks)
    """
    print("\n" + "=" * 80)
    print("DYNAMIC STRATIFIED DATA SPLITTING")
    print("=" * 80)
    print("Strategy: Representative sampling across resistance spectrum")
    print("=" * 80)
    
    # Analyze each chunk's resistance characteristics
    chunk_stats = []
    
    print("\n[Phase 1/4] Analyzing resistance distribution...")
    for f in all_files:
        chunk_num = int(f.stem.split('_')[-1])
        y_chunk = get_y_chunk(y_all, chunk_num, chunk_size, len(y_all))
        
        pos_count = sum(y_chunk)
        total_count = len(y_chunk)
        ratio = pos_count / total_count if total_count > 0 else 0
        
        chunk_stats.append({
            'file': f,
            'filename': f.name,
            'id': chunk_num,
            'pos_count': pos_count,
            'total': total_count,
            'ratio': ratio,
            'imbalance_score': abs(ratio - 0.5)  # Distance from perfect 50/50 balance
        })
    
    # Sort by imbalance_score (descending: Most Imbalanced (HARD) -> Most Balanced (EASY))
    df_stats = pd.DataFrame(chunk_stats).sort_values(by='imbalance_score', ascending=False).reset_index(drop=True)
    
    print("\n[Phase 2/4] Resistance spectrum analysis (sorted by Imbalance Score):")
    print("-" * 90)
    print(f"{'Rank':<6} | {'Chunk ID':<10} | {'Resistant':<12} | {'Total':<10} | {'Ratio':<10} | {'Imbalance Score':<15}")
    print("-" * 90)
    for idx, row in df_stats.head(5).iterrows():
        print(f"{idx+1:<6} | {row['id']:<10} | {row['pos_count']:<12} | {row['total']:<10} | {row['ratio']:<10.4f} | {row['imbalance_score']:<15.4f}")
    print("...")
    for idx, row in df_stats.tail(2).iterrows():
        print(f"{idx+1:<6} | {row['id']:<10} | {row['pos_count']:<12} | {row['total']:<10} | {row['ratio']:<10.4f} | {row['imbalance_score']:<15.4f}")
    print("-" * 90)
    
    # Validate and calculate fractional test chunks
    total_chunks = len(df_stats)
    test_chunks = max(1, int(np.floor(total_chunks * test_fraction)))
    
    if test_chunks >= total_chunks:
        print(f"\n⚠ WARNING: Calculated test chunks ({test_chunks}) >= total chunks ({total_chunks})")
        test_chunks = max(1, total_chunks - 2)  # Leave at least 2 for training
        print(f"  Auto-adjusted to: {test_chunks} test chunks")
    
    # Stratified test set selection using linspace (representative sampling)
    print(f"\n[Phase 3/4] Selecting {test_chunks} representative test chunks (Target fraction: {test_fraction:.1%})...")
    test_indices = np.linspace(0, total_chunks - 1, test_chunks, dtype=int)
    test_files = df_stats.iloc[test_indices]['file'].tolist()
    test_filenames = df_stats.iloc[test_indices]['filename'].tolist()
    
    print("Test set composition (stratified across difficulty):")
    print("-" * 80)
    for idx, test_idx in enumerate(test_indices):
        row = df_stats.iloc[test_idx]
        difficulty = "HARD (Imbalanced)" if row['imbalance_score'] > 0.3 else "MEDIUM" if row['imbalance_score'] > 0.15 else "EASY (Balanced)"
        print(f"  {idx+1}. {row['filename']:<30} | Ratio: {row['ratio']:.4f} | Imbalance: {row['imbalance_score']:.4f} | {difficulty}")
    print("-" * 80)
    
    # Training pool: all chunks NOT in test set
    train_mask = ~df_stats['file'].isin(test_files)
    train_files = df_stats[train_mask]['file'].tolist()
    train_filenames = df_stats[train_mask]['filename'].tolist()
    
    print(f"\n[Phase 4/4] Training pool: {len(train_files)} chunks")
    
    # Select Optuna subset from training pool (stratified by resistance)
    print("\nSelecting Optuna optimization subset using fractional logic...")
    train_df_sorted = df_stats[train_mask].reset_index(drop=True)
    
    final_optuna_count = max(1, int(np.floor(len(train_df_sorted) * optuna_fraction)))
    
    if final_optuna_count > len(train_df_sorted):
        print(f"  ⚠ WARNING: Calculated Optuna chunks ({final_optuna_count}) exceeds available training chunks ({len(train_df_sorted)}).")
        print(f"  Auto-adjusting to maximum available: {len(train_df_sorted)}")
        final_optuna_count = len(train_df_sorted)
        
    if final_optuna_count >= 1:
        # Stratified selection using np.linspace
        optuna_indices = np.linspace(0, len(train_df_sorted) - 1, final_optuna_count, dtype=int)
        
        optuna_files = train_df_sorted.iloc[optuna_indices]['file'].tolist()
        optuna_filenames = train_df_sorted.iloc[optuna_indices]['filename'].tolist()
        
        print(f"Loaded {final_optuna_count} chunks for Optuna based on config (stratified by resistance):")
        print("-" * 80)
        for idx_order, idx_val in enumerate(optuna_indices):
            row = train_df_sorted.iloc[idx_val]
            difficulty = "HARD (Imbalanced)" if row['imbalance_score'] > 0.3 else "MEDIUM" if row['imbalance_score'] > 0.15 else "EASY (Balanced)"
            print(f"  {idx_order+1}. {row['filename']:<30} | Ratio: {row['ratio']:.4f} | Imbalance: {row['imbalance_score']:.4f} | {difficulty}")
        print("-" * 80)
    else:
        # Fallback if no training chunks
        optuna_files = []
        optuna_filenames = []
        print("  ⚠ WARNING: No training chunks available for Optuna.")
    
    # Statistics
    total_train_samples = sum([len(get_y_chunk(y_all, int(f.stem.split('_')[-1]), chunk_size, len(y_all))) 
                               for f in train_files])
    total_test_samples = sum([len(get_y_chunk(y_all, int(f.stem.split('_')[-1]), chunk_size, len(y_all))) 
                              for f in test_files])
    total_optuna_samples = sum([len(get_y_chunk(y_all, int(f.stem.split('_')[-1]), chunk_size, len(y_all))) 
                                for f in optuna_files])
    
    print(f"\n{'='*80}")
    print("STRATIFIED SPLIT SUMMARY")
    print(f"{'='*80}")
    print(f"Total chunks:        {total_chunks}")
    print(f"Training chunks:     {len(train_files)} ({len(train_files)/total_chunks:.1%})")
    print(f"Test chunks:         {len(test_files)} ({len(test_files)/total_chunks:.1%})")
    print(f"Optuna chunks:       {len(optuna_files)} ({len(optuna_files)/len(train_files):.1%} of training)")
    print("")
    print(f"Training samples:    {total_train_samples:,}")
    print(f"Test samples:        {total_test_samples:,}")
    print(f"Optuna samples:      {total_optuna_samples:,}")
    print(f"{'='*80}")
    
    # Extract filenames for config storage
    optuna_filenames = [f.name for f in optuna_files]
    
    return train_files, train_filenames, test_files, test_filenames, optuna_files, optuna_filenames


def load_data_for_optuna(selected_files, y_all, chunk_size):
    """
    Load selected chunks and prepare data for Optuna optimization.
    
    Splits the selected data into training and validation sets for
    hyperparameter evaluation. Uses stratified splitting to maintain
    class balance in both sets.
    
    Args:
        selected_files: List of chunk file paths to load
        y_all: Complete label array
        chunk_size: Number of samples per chunk
    
    Returns:
        tuple: (dtrain, dval)
            - dtrain: XGBoost DMatrix for training
            - dval: XGBoost DMatrix for validation
    """
    print("\nLoading selected chunks into RAM...")
    
    X_list = []
    y_list = []
    
    for f in selected_files:
        try:
            X_chunk = load_npz(f)
            chunk_num = int(f.stem.split('_')[-1])
            y_chunk = get_y_chunk(y_all, chunk_num, chunk_size, len(y_all))
            
            X_list.append(X_chunk)
            y_list.append(y_chunk)
            
            print(f"  ✓ Loaded {f.name}: {X_chunk.shape[0]} samples")
            
        except Exception as e:
            print(f"ERROR: Failed to load {f.name}: {e}")
            sys.exit(1)
    
    # Stack all chunks
    X_opt = vstack(X_list)
    y_opt = np.concatenate(y_list)
    
    print(f"\nTotal optimization data: {X_opt.shape}")
    print(f"  Resistant: {np.sum(y_opt)} | Susceptible: {len(y_opt) - np.sum(y_opt)}")
    
    # Split into train/validation with stratification
    print(f"\nSplitting into train ({int((1-VAL_SPLIT)*100)}%) and validation ({int(VAL_SPLIT*100)}%)...")
    
    # Check if only one class exists (prevents stratification crash)
    unique_classes = np.unique(y_opt)
    use_stratify = y_opt if len(unique_classes) > 1 else None
    
    if len(unique_classes) == 1:
        print("  ⚠ WARNING: Only one class found in selected chunks. Stratification disabled.")
    
    X_train, X_val, y_train, y_val = train_test_split(
        X_opt, y_opt, 
        test_size=VAL_SPLIT, 
        random_state=RANDOM_SEED, 
        stratify=use_stratify
    )
    
    print(f"  Train: {X_train.shape} | Val: {X_val.shape}")

    # Build a compact, quantised QuantileDMatrix (max_bin=2 for binary data)
    # instead of a plain DMatrix. CRITICAL for parallel HPO: the quantised
    # gradient-index matrix is built ONCE and shared by all concurrent trials,
    # whereas a plain DMatrix makes every trial rebuild its own ~20 GB gmat →
    # running 20 trials at once OOM-killed the job (~550 GB). It also shrinks the
    # subset matrix itself from ~150 GB to a few GB. The validation matrix must
    # reference the training quantiles via ref=dtrain.
    max_bin = int(BASE_PARAMS.get('max_bin', 2))
    print(f"\nBuilding QuantileDMatrix (max_bin={max_bin})...")
    dtrain = xgb.QuantileDMatrix(X_train, label=y_train, max_bin=max_bin)
    dval = xgb.QuantileDMatrix(X_val, label=y_val, ref=dtrain, max_bin=max_bin)

    print("  ✓ Data preparation complete")

    return dtrain, dval


def save_antibiotic_specific_config(study, base_params, target_antibiotic, config_dir, 
                                    train_filenames, test_filenames, optuna_filenames):
    """
    Save optimization results to antibiotic-specific configuration file.
    
    Creates a dedicated YAML configuration file for each antibiotic containing:
    - Optimization metadata (date, best score, trials)
    - Fixed base XGBoost parameters
    - Optimized hyperparameters from Optuna study
    - Data split information (FULL train set, test set, optuna subset)
    
    CRITICAL DISTINCTION:
        - train_files: FULL training set (all chunks minus test) - for Script 05 final training
        - test_files: Stratified test set - for Script 06 evaluation
        - optuna_files: Small subset of train_files - ONLY for hyperparameter tuning (Script 04)
    
    This decoupling strategy enables:
    - Parallel optimization for multiple antibiotics without conflicts
    - Independent hyperparameter versioning per antibiotic
    - Reproducible data splits (scripts 05 and 06 use exact same files)
    - Clear separation between optimization subset and full training set
    - Clear audit trail of optimization results
    
    File Output: config/config_{antibiotic}.yaml
    
    Args:
        study: Completed Optuna study object
        base_params: Fixed XGBoost parameters (objective, eval_metric, etc.)
        target_antibiotic: Name of target antibiotic (e.g., 'ciprofloxacin')
        config_dir: Path to configuration directory
        train_filenames: List of FULL training chunk filenames (strings) - for Script 05
        test_filenames: List of test chunk filenames (strings) - for Script 06
        optuna_filenames: List of optuna subset filenames (strings) - for record only
    
    Returns:
        Path: Path to created configuration file
    
    Raises:
        IOError: If file write operation fails
    """
    # Organism-scoped experiment config (config/experiments/{organism}/...), so
    # two organisms that share an antibiotic name (e.g. ecoli/gentamicin vs
    # kpneumoniae/gentamicin) never overwrite each other's split + best params.
    antibiotic_config_path = resolve_path('experiment_config', organism=ORGANISM,
                                          antibiotic=target_antibiotic, config=config)
    antibiotic_config_path.parent.mkdir(parents=True, exist_ok=True)

    # Merge standard best params with the dynamically found n_estimators
    final_best_params = study.best_params.copy()
    if "n_estimators" in study.best_trial.user_attrs:
        final_best_params["n_estimators"] = study.best_trial.user_attrs["n_estimators"]

    # Construct configuration data structure
    model_config = {
        'antibiotic_metadata': {
            'target_antibiotic': target_antibiotic,
            'optimization_date': pd.Timestamp.now().isoformat(),
            'best_auc_score': float(study.best_value),
            'n_trials_completed': len(study.trials)
        },
        'data_split': {
            'description': 'Stratified split: FULL train (for 05_model_training.py), test (for 06_evaluation.py), optuna subset (record only)',
            'split_method': 'linspace_stratified_by_resistance_ratio',
            'train_files': train_filenames,
            'test_files': test_filenames,
            'optuna_files': optuna_filenames,
            'train_count': len(train_filenames),
            'test_count': len(test_filenames),
            'optuna_count': len(optuna_filenames)
        },
        'xgboost_params': base_params,
        'best_params': final_best_params
    }
    
    # Write to antibiotic-specific YAML file
    with open(antibiotic_config_path, 'w', encoding='utf-8') as f:
        # Add header comment
        f.write("# " + "=" * 78 + "\n")
        f.write(f"# AUTO-GENERATED CONFIGURATION: {target_antibiotic.upper()}\n")
        f.write("# " + "=" * 78 + "\n")
        f.write("# Generated by: 04_optimization.py\n")
        f.write(f"# Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"# Best ROC AUC: {study.best_value:.4f}\n")
        f.write(f"# Trials: {len(study.trials)}\n")
        f.write("#\n")
        f.write("# WARNING: This file is automatically generated. Manual edits will be\n")
        f.write("#          overwritten when re-running hyperparameter optimization.\n")
        f.write("# " + "=" * 78 + "\n\n")
        
        yaml.dump(model_config, f, default_flow_style=False, sort_keys=False)
    
    return antibiotic_config_path


def objective(trial, dtrain, dval, base_params, colsample_range):
    """
    Optuna objective function for hyperparameter optimization.

    This function is called by Optuna for each trial. It suggests a set of
    hyperparameters, trains an XGBoost model, and returns the validation
    performance metric (ROC AUC) for optimization.

    Hyperparameter Ranges Explained:
        - learning_rate [0.01, 0.3]: Lower values more stable but slower
        - max_depth [3, 12]: Deeper trees can model complex patterns but overfit
        - subsample [0.6, 1.0]: Using subset of data per tree prevents overfitting
        - colsample_bytree (log-scale window around 1/√p): aligns with the
          square-root heuristic in METHODOLOGY.md (see compute_colsample_range)
        - min_child_weight [1, 10]: Higher values more conservative (less overfitting)
        - gamma [0, 5]: Regularization parameter (minimum loss reduction)
        - n_estimators: determined by early stopping (saved as a user attr)

    NOTE on class imbalance: scale_pos_weight is intentionally NOT tuned here.
    The final model (05_model_training.py) handles imbalance via per-chunk
    Dynamic Instance Weighting (neg/pos ratio applied to the DMatrix). Tuning
    scale_pos_weight in HPO and then weighting again in training would
    double-count the imbalance correction.

    Args:
        trial: Optuna trial object
        dtrain: XGBoost DMatrix for training
        dval: XGBoost DMatrix for validation
        base_params: Fixed XGBoost parameters (objective, eval_metric, etc.)
        colsample_range: (lower, upper) bounds for colsample_bytree, derived
            from the feature count via compute_colsample_range().

    Returns:
        float: Best validation ROC AUC score achieved during training
    """
    # Create parameter dictionary starting with base parameters
    params = base_params.copy()

    col_lower, col_upper = colsample_range

    # Add trial-specific hyperparameters
    params.update({
        'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
        'max_depth': trial.suggest_int('max_depth', 3, 12),
        'subsample': trial.suggest_float('subsample', 0.6, 1.0),
        # √p-anchored, log-scale window (see compute_colsample_range / METHODOLOGY.md)
        'colsample_bytree': trial.suggest_float('colsample_bytree', col_lower, col_upper, log=True),
        'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
        'gamma': trial.suggest_float('gamma', 0, 5),
    })
    
    # Train model with early stopping (let XGBoost find the optimal trees)
    model = xgb.train(
        params,
        dtrain,
        num_boost_round=1000,  # Fixed high ceiling
        evals=[(dval, "validation")],
        early_stopping_rounds=EARLY_STOPPING_ROUNDS,
        verbose_eval=False
    )
    
    # Dynamically save the EXACT optimal tree count found by early stopping.
    # best_iteration is 0-indexed (the index of the best round), so the number
    # of trees to keep is best_iteration + 1. Storing best_iteration directly
    # undercounts by one and, when the very first tree is best, yields 0 — which
    # would make 05 train zero trees (model stays None). max(1, …) is a final guard.
    trial.set_user_attr("n_estimators", max(1, int(model.best_iteration) + 1))
    
    # Return best validation score (ROC AUC)
    return model.best_score


def generate_optuna_plots(study, target_antibiotic):
    """Generate and save Optuna optimization visualizations using matplotlib."""
    print("\nGenerating Optuna visualization plots...")
    
    # Derive output directory from centralised config (03_model_optimization)
    output_dir = resolve_path('dir_03_model_optimization', organism=ORGANISM,
                              antibiotic=target_antibiotic, config=config)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        from optuna.visualization.matplotlib import plot_optimization_history, plot_param_importances
        import warnings
        warnings.filterwarnings('ignore', category=UserWarning)  # Suppress optuna matplotlib experimental warnings
        
        # Set style
        sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
        
        # 1. Optimization History Plot
        fig1 = plt.figure(figsize=(10, 6))
        ax1 = plot_optimization_history(study)
        plt.title(f'Hyperparameter Optimization History ({target_antibiotic.upper()})', fontsize=14, pad=15)
        plt.tight_layout()
        hist_path = output_dir / f"01_optuna_history_{target_antibiotic}.png"
        plt.savefig(hist_path, dpi=300, bbox_inches='tight')
        
        # Save data to CSV
        history_df = pd.DataFrame([{
            'Trial_Number': t.number,
            'Objective_Value': t.value,
            'State': t.state.name
        } for t in study.trials if t.value is not None])
        hist_csv_path = output_dir / f"01_optuna_history_{target_antibiotic}.csv"
        history_df.to_csv(hist_csv_path, index=False)
        
        plt.close(fig1)
        print(f"  ✓ Saved history plot: {hist_path.name}")
        
        # 2. Hyperparameter Importance Plot
        fig2 = plt.figure(figsize=(10, 6))
        ax2 = plot_param_importances(study)
        plt.title(f'Hyperparameter Importance ({target_antibiotic.upper()})', fontsize=14, pad=15)
        plt.tight_layout()
        imp_path = output_dir / f"02_optuna_importance_{target_antibiotic}.png"
        plt.savefig(imp_path, dpi=300, bbox_inches='tight')
        
        # Save data to CSV
        try:
            importances = optuna.importance.get_param_importances(study)
            imp_df = pd.DataFrame({
                'Hyperparameter': list(importances.keys()),
                'Importance': list(importances.values())
            })
            imp_csv_path = output_dir / f"02_optuna_importance_{target_antibiotic}.csv"
            imp_df.to_csv(imp_csv_path, index=False)
        except Exception:
            pass
            
        plt.close(fig2)
        print(f"  ✓ Saved importance plot: {imp_path.name}")
        
    except ImportError:
        print("  ⚠ Could not generate plots. Ensure matplotlib is installed.")
    except Exception as e:
        print(f"  ⚠ Failed to generate Optuna plots: {e}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================
def main():
    """
    Execute hyperparameter optimization pipeline.
    
    Pipeline:
        1. Load label data and identify training chunks
        2. Perform smart chunk selection (memory-efficient subset)
        3. Load selected chunks into RAM
        4. Run Optuna optimization (N_TRIALS iterations)
        5. Save best hyperparameters to JSON file
        6. Export full study object for further analysis
    """
    print("=" * 80)
    print("XGBOOST HYPERPARAMETER OPTIMIZATION FOR AMR PREDICTION")
    print("=" * 80)
    print(f"Target Antibiotic: {TARGET_ANTIBIOTIC}")
    print(f"Optimization Trials: {N_TRIALS}")
    print(f"Validation Split: {int(VAL_SPLIT*100)}%")
    print("=" * 80)
    
    # ------------------------------------------------------------------------
    # STEP 1: Load Data Information
    # ------------------------------------------------------------------------
    print("\n[STEP 1/5] Loading data information...")
    
    try:
        # Load labels
        y_path = MATRIX_DIR / f"y_{TARGET_ANTIBIOTIC}.csv"
        if not y_path.exists():
            raise FileNotFoundError(f"Label file not found: {y_path}")
        
        y_all = pd.read_csv(y_path, encoding='utf-8')['label'].values
        print(f"  ✓ Loaded {len(y_all)} labels from {y_path.name}")
        
        # Find all matrix chunks
        all_files = sorted(
            list(MATRIX_DIR.glob(f"X_{TARGET_ANTIBIOTIC}_part_*.npz")),
            key=lambda x: int(x.stem.split('_')[-1])
        )
        
        if not all_files:
            raise FileNotFoundError(
                f"No matrix chunks found matching: X_{TARGET_ANTIBIOTIC}_part_*.npz"
            )
        
        print(f"  ✓ Found {len(all_files)} matrix chunks")
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        print("\nPlease ensure you have run:")
        print("  1. Data validation: 01_data_validation.py")
        print("  2. K-mer counting: 02_kmer_extraction.py")
        print("  3. Matrix creation: 03_matrix_construction.py")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Unexpected error: {e}")
        sys.exit(1)
    
    # ------------------------------------------------------------------------
    # STEP 2: Stratified Data Splitting
    # ------------------------------------------------------------------------
    # Split data using new fractional strategy
    print("\n[STEP 2/5] Dynamic Stratified Chunk Splitting...")
    train_files, train_filenames, test_files, test_filenames, optuna_files, optuna_filenames = analyze_and_stratify_all_chunks(
        y_all=y_all,
        all_files=all_files,
        chunk_size=CHUNK_SIZE,
        test_fraction=TEST_FRACTION,
        optuna_fraction=OPTUNA_FRACTION
    )
    
    # ------------------------------------------------------------------------
    # STEP 3: Prepare Optuna Data
    # ------------------------------------------------------------------------
    print("\n[STEP 3/5] Loading data for hyperparameter optimization...")

    dtrain_opt, dval_opt = load_data_for_optuna(optuna_files, y_all, CHUNK_SIZE)

    # Derive the √p-anchored colsample_bytree search range from the feature count
    n_features = count_features(MATRIX_DIR)
    colsample_range = compute_colsample_range(n_features)
    if n_features:
        print(f"  ✓ Features (p): {n_features:,} | 1/√p ≈ {1.0/np.sqrt(n_features):.2e}")
    else:
        print("  ⚠ features.txt not found; using fallback colsample_bytree range.")
    print(f"  ✓ colsample_bytree search range (log-scale): "
          f"[{colsample_range[0]:.2e}, {colsample_range[1]:.2e}]")
    
    # ------------------------------------------------------------------------
    # STEP 4: Run Optuna Optimization
    # ------------------------------------------------------------------------
    
    print("\n[STEP 4/5] Running Optuna hyperparameter optimization...")
    print("=" * 80)
    print(f"Starting {N_TRIALS} trials"
          + (f" (early stop after {OPTUNA_PATIENCE} no-improvement)" if OPTUNA_PATIENCE > 0
             else " (no early stopping)"))
    print(f"Parallelism: {OPTUNA_N_JOBS} concurrent trials × {THREADS_PER_TRIAL} threads "
          f"(allocated cores: {TOTAL_CORES})")
    print("=" * 80)
    
    # Create Optuna study with a SEEDED TPE sampler so the hyperparameter search
    # is reproducible (random_seed from config). NB: with parallel trials
    # (n_jobs>1) completion order still varies slightly, but the proposed
    # parameter sequence is fixed by the seed.
    study = optuna.create_study(
        direction='maximize',  # Maximize ROC AUC
        study_name=f"xgboost_{TARGET_ANTIBIOTIC}_optimization",
        sampler=optuna.samplers.TPESampler(seed=RANDOM_SEED),
    )
    
    # Early-stopping callback: stop dispatching new trials once OPTUNA_PATIENCE
    # completed trials pass without improving the best value (TPE typically
    # converges well before n_trials). With n_jobs>1, study.stop() lets in-flight
    # trials finish but launches no more. OPTUNA_PATIENCE<=0 disables it.
    _es = {'best': None, 'since': 0}

    def _early_stop_cb(study_, trial_):
        if OPTUNA_PATIENCE <= 0 or trial_.state != optuna.trial.TrialState.COMPLETE \
                or trial_.value is None:
            return
        best = study_.best_value
        if _es['best'] is None or best > _es['best'] + 1e-9:
            _es['best'], _es['since'] = best, 0
        else:
            _es['since'] += 1
            if _es['since'] >= OPTUNA_PATIENCE:
                print(f"\n  ⏹ HPO early stop: best {best:.4f} not improved for "
                      f"{OPTUNA_PATIENCE} trials ({len(study_.trials)} run); stopping.")
                study_.stop()

    # Run optimization with Graceful Shutdown (Fault Tolerance)
    try:
        study.optimize(
            lambda trial: objective(trial, dtrain_opt, dval_opt, BASE_PARAMS, colsample_range),
            n_trials=N_TRIALS,
            n_jobs=OPTUNA_N_JOBS,   # run trials concurrently to use all cores
            callbacks=[_early_stop_cb],
            show_progress_bar=False
        )
    except KeyboardInterrupt:
        print("\n\n  ⚠ Optimization manually interrupted by user (Ctrl+C).")
    except Exception as e:
        print(f"\n\n  ⚠ Optimization stopped unexpectedly: {e}")
    finally:
        # Graceful Shutdown: Save progress if at least one trial completed
        valid_trials = [t for t in study.trials if t.value is not None]
        
        if len(valid_trials) > 0:
            print("\n" + "=" * 80)
            print("SAVING PROGRESS (GRACEFUL SHUTDOWN)")
            print("=" * 80)
            print(f"  ✓ Found {len(valid_trials)} completed trials. Extracting best results so far...")
            
            # Generate Visualizations
            generate_optuna_plots(study, TARGET_ANTIBIOTIC)
            
            # Display best results found
            print(f"\nBest ROC AUC Score: {study.best_value:.4f}")
            print("\nBest Hyperparameters:")
            print("-" * 80)
            for param, value in study.best_params.items():
                print(f"  {param:<20} : {value}")
            print("-" * 80)
            
            # Save Results
            print("\nSaving configurations...")
            try:
                antibiotic_config_path = save_antibiotic_specific_config(
                    study=study,
                    base_params=BASE_PARAMS,
                    target_antibiotic=TARGET_ANTIBIOTIC,
                    config_dir=CONFIG_DIR,
                    train_filenames=train_filenames,
                    test_filenames=test_filenames,
                    optuna_filenames=optuna_filenames
                )
                print(f"  ✓ Antibiotic-specific config created: {antibiotic_config_path.name}")
            except Exception as e:
                print(f"ERROR: Failed to create antibiotic config: {e}")
            
            # Save complete study object
            study_path = MODELS_DIR / f"optuna_study_{TARGET_ANTIBIOTIC}.pkl"
            try:
                joblib.dump(study, study_path)
                print(f"  ✓ Full study object saved: {study_path}")
                
                # TIMESTAMP BACKUP
                ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                backup_study_path = MODELS_DIR / f"optuna_study_{TARGET_ANTIBIOTIC}_{ts}.pkl"
                shutil.copy2(study_path, backup_study_path)
                print(f"  ✓ Backup study object saved: {backup_study_path}")
            except Exception as e:
                print(f"WARNING: Failed to save study object: {e}")
            
            # ----------------------------------------------------------------
            # MLOps: write run_metadata.json (git hash, versions, seed, params,
            # data fingerprint). Best-effort — never break optimization.
            # ----------------------------------------------------------------
            try:
                organism = get_target(config=config)[0]
                run_id = rm.make_run_id(organism, TARGET_ANTIBIOTIC)
                run_dir = resolve_path('run_dir', organism=organism,
                                       antibiotic=TARGET_ANTIBIOTIC, run_id=run_id, config=config)
                # Provenance completeness (audit Issue 4): record the total genome
                # count + effective min_support so pipeline_runs isn't left NULL.
                n_genomes = None
                try:
                    _yp = (resolve_path('matrix_dir', organism=organism,
                                        antibiotic=TARGET_ANTIBIOTIC, config=config)
                           / f"y_{TARGET_ANTIBIOTIC}.csv")
                    if _yp.exists():
                        n_genomes = sum(1 for _ in open(_yp, encoding='utf-8')) - 1
                except Exception:
                    pass
                min_support = ((config.get('unitig', {}) or {}).get('min_support')
                               or (config.get('preprocessing', {}) or {}).get('min_support'))
                meta = rm.build_run_metadata(
                    organism=organism,
                    antibiotic=TARGET_ANTIBIOTIC,
                    run_id=run_id,
                    seed=RANDOM_SEED,
                    params=study.best_params,
                    data_files=[str(f) for f in optuna_files],
                    config=config,
                    extra={"stage": "04_optimization",
                           "best_score": float(study.best_value),
                           "n_trials": len(study.trials),
                           "n_genomes": n_genomes,
                           "min_support": min_support},
                )
                if rm.write_json(run_dir / "run_metadata.json", meta):
                    print(f"  ✓ Run metadata written: runs/.../{run_id}/run_metadata.json")
            except Exception as e:
                print(f"  ⚠ Could not write run metadata: {e}")

            # Display next steps
            print("\n" + "=" * 80)
            print("NEXT STEPS")
            print("=" * 80)
            print(f"1. Review optimized parameters in: config/experiments/{ORGANISM}/config_{TARGET_ANTIBIOTIC}.yaml")
            print("2. Run final training with optimized hyperparameters:")
            print("   python scripts/05_model_training.py")
            print("=" * 80)
        else:
            print("\n" + "=" * 80)
            print("  ⚠ No completed trials found. Cannot save configuration or plots.")
            print("=" * 80)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    main()
