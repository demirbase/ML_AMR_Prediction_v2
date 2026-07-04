#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Comprehensive Model Performance Analysis Module

This script performs in-depth evaluation of the trained AMR prediction model,
including threshold optimization, detailed performance metrics, error analysis,
and clinical performance interpretation.

Configuration Architecture:
    This script uses a two-stage configuration loading system:
    1. Global Config (config.yaml): Project settings, paths, and target antibiotic
    2. Specific Config (config_{antibiotic}.yaml): Exact test set from stratified
       split performed during optimization (04_optimization.py)

Evaluation Framework:
    1. Model loading from antibiotic-specific directory
    2. Test data preparation using exact test files from config
    3. Threshold optimization using multiple metrics
    4. Comprehensive performance assessment
    5. Error analysis (false positives/negatives)
    6. Clinical interpretation of results

Clinical Context:
    In antimicrobial resistance prediction, we must balance:
    - Sensitivity: Catching resistant bacteria (false negatives costly)
    - Specificity: Avoiding unnecessary antibiotics (false positives wasteful)
    
    Threshold optimization helps find the right balance for clinical deployment.
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import pandas as pd
import numpy as np
import xgboost as xgb
import yaml
from pathlib import Path
from scipy.sparse import load_npz, vstack
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    roc_auc_score,
    f1_score,
    matthews_corrcoef,
    precision_recall_curve,
    auc,
    cohen_kappa_score,
    balanced_accuracy_score,
    roc_curve,
    average_precision_score,
    precision_score,
    recall_score
)
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from sklearn.calibration import calibration_curve

# Shared label-slicing helper (single source of truth) — see scripts/utils.py
from utils import get_y_chunk
# MLOps run provenance (SCALE_MLOPS_PLAN.md §7) — additive, best-effort.
from lib import run_metadata as rm
from lib.config import resolve_path, get_target


# ============================================================================
# LOAD CONFIGURATION FROM YAML
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

# Load configuration
if not CONFIG_PATH.exists():
    raise FileNotFoundError(
        f"Configuration file not found: {CONFIG_PATH}\n"
        f"Please ensure config.yaml exists in the config/ directory."
    )

with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

# Extract configuration values
TARGET_ANTIBIOTIC = get_target(config=config)[1]
ORGANISM = get_target(config=config)[0]
CHUNK_SIZE = config['preprocessing']['chunk_size']
N_JOBS = config['xgboost_params'].get('n_jobs', -1)

# Organism-aware paths (SCALE_MLOPS_PLAN §4.2)
MATRIX_DIR = resolve_path('matrix_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
MODELS_DIR = resolve_path('models_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
# Output directory — derived from the centralised config (04_evaluation subdir)
OUTPUT_DIR = resolve_path('dir_04_evaluation', organism=ORGANISM,
                          antibiotic=TARGET_ANTIBIOTIC, config=config)

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)



# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
# get_y_chunk() is imported from utils.py (shared with 04 and 05).


def load_test_files_from_config():
    """
    Load test file list from antibiotic-specific configuration.
    
    This function ensures evaluation uses the EXACT test set defined during
    optimization (04_optimization.py). This maintains:
    - Reproducibility: Same test set across all evaluation runs
    - Stratification: Representative test set (hard, medium, easy cases)
    - Isolation: Test set never seen during training or hyperparameter tuning
    
    Why This Matters:
        Using a different test set than optimization would:
        1. Break reproducibility (results vary by random seed)
        2. Violate data leakage prevention (if test overlaps with train)
        3. Lose stratification benefits (may over/under-represent difficulty)
    
    Returns:
        list: Test chunk filenames (strings) from stratified split
    
    Raises:
        FileNotFoundError: If antibiotic-specific config doesn't exist
        KeyError: If data_split section is missing
    """
    # NOTE: do NOT pass config=config here — this function later rebinds a local
    # `config` (config = yaml.safe_load(...)), which makes `config` function-local
    # for the whole scope; referencing it here would raise UnboundLocalError.
    # resolve_path loads config.yaml itself when config is omitted.
    antibiotic_config_path = resolve_path('experiment_config', organism=ORGANISM,
                                          antibiotic=TARGET_ANTIBIOTIC)

    if not antibiotic_config_path.exists():
        raise FileNotFoundError(
            f"Antibiotic-specific configuration not found: {antibiotic_config_path}\n\n"
            f"This file is generated by hyperparameter optimization.\n"
            f"Please run optimization for '{TARGET_ANTIBIOTIC}' first:\n\n"
            f"  python scripts/04_optimization.py\n\n"
            f"This will create: {antibiotic_config_path}\n"
            f"containing the stratified test set for evaluation."
        )
    
    try:
        with open(antibiotic_config_path, 'r', encoding='utf-8') as f:
            config = yaml.safe_load(f)
        
        if 'data_split' not in config:
            raise KeyError(
                f"'data_split' section not found in {antibiotic_config_path.name}\n"
                f"The config file may be from an old version. Re-run optimization:\n"
                f"  python scripts/04_optimization.py"
            )
        
        test_files = config['data_split'].get('test_files', [])
        
        # Load the unbiased threshold calculated from training prevalence
        optimal_threshold = float(config.get('evaluation', {}).get('optimal_threshold', 0.5))
        
        if not test_files:
            raise ValueError(
                "No test files found in data_split section.\n"
                "Re-run optimization: python scripts/04_optimization.py"
            )
        
        print(f"✓ Loaded test set from: {antibiotic_config_path.name}")
        print(f"  Split method: {config['data_split'].get('split_method', 'N/A')}")
        print(f"  Test chunks: {len(test_files)}")
        print(f"  Unbiased Threshold: {optimal_threshold:.4f}")
        
        return test_files, optimal_threshold
        
    except yaml.YAMLError as e:
        print(f"ERROR: Malformed YAML file: {antibiotic_config_path}")
        print(f"  Details: {e}")
        sys.exit(1)
    except (KeyError, ValueError) as e:
        print(f"ERROR: {e}")
        sys.exit(1)


def load_test_data(y_all, all_chunk_files, test_filenames):
    """
    Load and aggregate test data from config-specified chunk files.
    
    Uses the exact test set defined during optimization (stratified split)
    to ensure evaluation is performed on the representative sample that
    includes hard, medium, and easy resistance cases.
    
    CRITICAL: This function uses the EXACT file list from config with NO
    randomization or reordering. The test set is fixed during optimization.
    
    Args:
        y_all: Complete label array
        all_chunk_files: List of ALL available matrix chunk file paths
        test_filenames: List of test chunk filenames from config (exact order preserved)
    
    Returns:
        tuple: (X_test, y_test, genome_ids_test)
            - X_test: Stacked sparse feature matrix
            - y_test: Corresponding labels
            - genome_ids_test: Genome identifiers (for error tracing)
    """
    print("\n" + "=" * 80)
    print("LOADING TEST DATA (FROM STRATIFIED SPLIT)")
    print("=" * 80)
    
    # Create filename -> path mapping
    file_map = {f.name: f for f in all_chunk_files}
    
    # Resolve test files from config (EXACT order preserved)
    test_files = []
    for fname in test_filenames:
        if fname in file_map:
            test_files.append(file_map[fname])
        else:
            print(f"WARNING: Test file '{fname}' not found in matrix directory")
    
    if not test_files:
        raise FileNotFoundError(
            "No test files could be resolved from config.\n"
            "Check that matrix chunks exist and match the config file list."
        )
    
    print(f"Test chunks: {[f.name for f in test_files]}")
    print(f"Total test chunks: {len(test_files)}")
    
    X_test_list = []
    y_test_list = []
    genome_ids_test = []
    
    # Load genome IDs
    try:
        genomes_file = MATRIX_DIR / f"genomes_{TARGET_ANTIBIOTIC}.csv"
        if not genomes_file.exists():
            raise FileNotFoundError(f"Genome ID file not found: {genomes_file}")
        
        # dtype=str (audit Issue 5): keep Genome IDs exact so 16_model_preds joins
        # against the phenotype/tools by the same string key (no float mangling).
        genomes_df = pd.read_csv(genomes_file, encoding='utf-8', dtype={'Genome ID': str})
        all_genome_ids = genomes_df['Genome ID'].values
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to load genome IDs: {e}")
        sys.exit(1)
    
    # Load each test chunk
    for chunk_file in test_files:
        try:
            X_chunk = load_npz(chunk_file)
            chunk_num = int(chunk_file.stem.split('_')[-1])
            
            # Get corresponding labels and genome IDs
            y_chunk = get_y_chunk(y_all, chunk_num, CHUNK_SIZE, len(y_all))
            genome_ids_chunk = get_y_chunk(all_genome_ids, chunk_num, CHUNK_SIZE, len(all_genome_ids))
            
            X_test_list.append(X_chunk)
            y_test_list.append(y_chunk)
            genome_ids_test.extend(genome_ids_chunk)
            
            print(f"  ✓ Loaded {chunk_file.name}: {X_chunk.shape[0]} samples")
            
        except Exception as e:
            print(f"ERROR: Failed to load {chunk_file.name}: {e}")
            sys.exit(1)
    
    # Stack all chunks
    X_test = vstack(X_test_list)
    y_test = np.concatenate(y_test_list)
    
    print(f"\n✓ Total test data: {X_test.shape}")
    print(f"  Resistant: {np.sum(y_test)} | Susceptible: {len(y_test) - np.sum(y_test)}")
    print("=" * 80)
    
    return X_test, y_test, np.array(genome_ids_test)



# ============================================================================
# Threshold logic is now handled in 05_model_training.py to prevent Data Leakage.
# We strictly use the Training Set prevalence threshold for unseen test data.
# ============================================================================



# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================
def plot_confusion_matrix_enhanced(y_true, y_pred, output_dir, antibiotic):
    """
    Generate publication-quality confusion matrix heatmap.
    
    Args:
        y_true (ndarray): Ground truth labels
        y_pred (ndarray): Predicted labels
        output_dir (Path): Directory to save the plot
        antibiotic (str): Antibiotic name for title
    """
    cm = confusion_matrix(y_true, y_pred)
    
    plt.figure(figsize=(8, 6))
    sns.heatmap(
        cm,
        annot=True,
        fmt='d',
        cmap='Blues',
        cbar_kws={'label': 'Count'},
        square=True,
        linewidths=1,
        linecolor='gray'
    )
    
    plt.title(f'Confusion Matrix: {antibiotic.title()} Resistance Prediction', 
              fontsize=14, fontweight='bold', pad=20)
    plt.suptitle("(Unbiased threshold from training; no test-set tuning)", fontsize=9, color='gray', y=0.01)
    plt.ylabel('True Label', fontsize=12)
    plt.xlabel('Predicted Label', fontsize=12)
    plt.xticks([0.5, 1.5], ['Susceptible', 'Resistant'], rotation=0)
    plt.yticks([0.5, 1.5], ['Susceptible', 'Resistant'], rotation=0)
    
    plt.tight_layout()
    output_path = output_dir / f'01_confusion_matrix_{antibiotic}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Save data to CSV
    tn, fp, fn, tp = cm.ravel()
    cm_df = pd.DataFrame({
        'Actual': ['Susceptible', 'Susceptible', 'Resistant', 'Resistant'],
        'Predicted': ['Susceptible', 'Resistant', 'Susceptible', 'Resistant'],
        'Count': [tn, fp, fn, tp]
    })
    cm_df.to_csv(output_dir / f'01_confusion_matrix_{antibiotic}.csv', index=False)
    
    plt.close()
    print(f"  ✓ Confusion matrix saved: {output_path.name}")


def plot_roc_curve_analysis(y_true, y_prob, output_dir, antibiotic):
    """
    Generate ROC curve with AUC score.
    
    Args:
        y_true (ndarray): Ground truth labels
        y_prob (ndarray): Predicted probabilities
        output_dir (Path): Directory to save the plot
        antibiotic (str): Antibiotic name for title
    """
    fpr, tpr, _ = roc_curve(y_true, y_prob)
    roc_auc = roc_auc_score(y_true, y_prob)
    
    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='#2E86AB', linewidth=2.5, label=f'ROC Curve (AUC = {roc_auc:.4f})')
    plt.plot([0, 1], [0, 1], color='gray', linestyle='--', linewidth=1.5, label='Random Classifier')
    
    plt.fill_between(fpr, tpr, alpha=0.2, color='#2E86AB')
    
    plt.xlabel('False Positive Rate', fontsize=12)
    plt.ylabel('True Positive Rate', fontsize=12)
    plt.title(f'ROC Curve: {antibiotic.title()} Resistance Prediction', 
              fontsize=14, fontweight='bold', pad=20)
    plt.legend(loc='lower right', fontsize=11)
    plt.grid(alpha=0.3, linestyle=':', linewidth=0.8)
    
    plt.tight_layout()
    output_path = output_dir / f'02_roc_curve_{antibiotic}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Save data to CSV
    fpr, tpr, roc_thresh = roc_curve(y_true, y_prob)
    roc_df = pd.DataFrame({
        'False_Positive_Rate': fpr,
        'True_Positive_Rate': tpr,
        'Threshold': roc_thresh
    })
    roc_df.to_csv(output_dir / f'02_roc_curve_{antibiotic}.csv', index=False)
    
    plt.close()
    print(f"  ✓ ROC curve saved: {output_path.name}")


def plot_precision_recall_curve_analysis(y_true, y_prob, output_dir, antibiotic):
    """
    Generate Precision-Recall curve with Average Precision score.
    
    Args:
        y_true (ndarray): Ground truth labels
        y_prob (ndarray): Predicted probabilities
        output_dir (Path): Directory to save the plot
        antibiotic (str): Antibiotic name for title
    """
    precision, recall, _ = precision_recall_curve(y_true, y_prob)
    avg_precision = average_precision_score(y_true, y_prob)
    
    # Calculate baseline (prevalence)
    baseline = np.sum(y_true) / len(y_true)
    
    plt.figure(figsize=(8, 6))
    plt.plot(recall, precision, color='#A23B72', linewidth=2.5, 
             label=f'PR Curve (AP = {avg_precision:.4f})')
    plt.axhline(y=baseline, color='gray', linestyle='--', linewidth=1.5, 
                label=f'Baseline (Prevalence = {baseline:.4f})')
    
    plt.fill_between(recall, precision, alpha=0.2, color='#A23B72')
    
    plt.xlabel('Recall', fontsize=12)
    plt.ylabel('Precision', fontsize=12)
    plt.title(f'Precision-Recall Curve: {antibiotic.title()} Resistance', 
              fontsize=14, fontweight='bold', pad=20)
    plt.legend(loc='best', fontsize=11)
    plt.grid(alpha=0.3, linestyle=':', linewidth=0.8)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    
    plt.tight_layout()
    output_path = output_dir / f'03_pr_curve_{antibiotic}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Save data to CSV
    precision_pts, recall_pts, pr_thresh = precision_recall_curve(y_true, y_prob)
    thresh_array = list(pr_thresh) + [1.0]
    pr_df = pd.DataFrame({
        'Recall': recall_pts,
        'Precision': precision_pts,
        'Threshold': thresh_array
    })
    pr_df.to_csv(output_dir / f'03_pr_curve_{antibiotic}.csv', index=False)
    
    plt.close()
    print(f"  ✓ Precision-Recall curve saved: {output_path.name}")


def plot_threshold_analysis(results_df, output_dir, antibiotic):
    """
    Generate multi-line plot showing how metrics vary with threshold.
    
    Args:
        results_df (DataFrame): Threshold analysis results from find_best_threshold
        output_dir (Path): Directory to save the plot
        antibiotic (str): Antibiotic name for title
    """
    plt.figure(figsize=(10, 6))
    
    plt.plot(results_df['threshold'], results_df['f1'], 
             marker='o', linewidth=2, markersize=4, label='F1 Score', color='#06A77D')
    plt.plot(results_df['threshold'], results_df['mcc'], 
             marker='s', linewidth=2, markersize=4, label='MCC', color='#D62246')
    plt.plot(results_df['threshold'], results_df['recall'], 
             marker='^', linewidth=2, markersize=4, label='Recall', color='#4E8098')
    plt.plot(results_df['threshold'], results_df['precision'], 
             marker='d', linewidth=2, markersize=4, label='Precision', color='#F18F01')
    
    # Mark optimal threshold
    best_idx = results_df['mcc'].idxmax()
    best_thresh = results_df.loc[best_idx, 'threshold']
    best_mcc = results_df.loc[best_idx, 'mcc']
    
    plt.axvline(x=best_thresh, color='red', linestyle='--', linewidth=1.5, 
                label=f'Optimal Threshold = {best_thresh:.2f}')
    
    plt.xlabel('Classification Threshold', fontsize=12)
    plt.ylabel('Metric Value', fontsize=12)
    plt.title(f'Threshold Sensitivity Analysis: {antibiotic.title()} Resistance', 
              fontsize=14, fontweight='bold', pad=20)
    plt.legend(loc='best', fontsize=10, ncol=2)
    plt.grid(alpha=0.3, linestyle=':', linewidth=0.8)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.0])
    
    plt.tight_layout()
    output_path = output_dir / f'06_threshold_analysis_{antibiotic}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  ✓ Threshold analysis saved: {output_path.name}")


def plot_probability_distribution(y_true, y_prob, output_dir, antibiotic, threshold=0.5):
    """Generate KDE plots for predicted probabilities grouped by true class."""
    plt.figure(figsize=(8, 6))

    sns.kdeplot(y_prob[y_true == 0], fill=True, color='#2E86AB', label='True: Susceptible (0)', alpha=0.5)
    sns.kdeplot(y_prob[y_true == 1], fill=True, color='#D62246', label='True: Resistant (1)', alpha=0.5)

    plt.axvline(x=threshold, color='black', linestyle='--', linewidth=1.5,
                label=f"Operating Threshold = {threshold:.4f} (training-derived)")

    plt.xlabel('Predicted Probability of Resistance', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.title(f'Probability Distribution: {antibiotic.title()} Resistance',
              fontsize=14, fontweight='bold', pad=20)
    plt.legend(loc='upper center', fontsize=11)
    plt.xlim([-0.05, 1.05])
    plt.grid(alpha=0.3, linestyle=':', linewidth=0.8)

    plt.tight_layout()
    output_path = output_dir / f'04_probability_distribution_{antibiotic}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Save data to CSV
    prob_df = pd.DataFrame({
        'True_Label': y_true,
        'Predicted_Probability': y_prob
    })
    prob_df.to_csv(output_dir / f'04_probability_distribution_{antibiotic}.csv', index=False)
    
    plt.close()
    print(f"  ✓ Probability distribution plot saved: {output_path.name}")


def plot_calibration_curve_analysis(y_true, y_prob, output_dir, antibiotic):
    """Generate a reliability diagram (calibration curve)."""
    prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=10, strategy='uniform')

    plt.figure(figsize=(8, 6))
    plt.plot(prob_pred, prob_true, marker='o', linewidth=2, color='#A23B72', label='XGBoost Model')
    plt.plot([0, 1], [0, 1], linestyle='--', color='gray', label='Perfectly Calibrated')

    plt.xlabel('Mean Predicted Probability', fontsize=12)
    plt.ylabel('Fraction of Positives (True Probability)', fontsize=12)
    plt.title(f'Calibration Curve (Reliability): {antibiotic.title()}',
              fontsize=14, fontweight='bold', pad=20)
    plt.legend(loc='best', fontsize=11)
    plt.grid(alpha=0.3, linestyle=':', linewidth=0.8)

    plt.tight_layout()
    output_path = output_dir / f'05_calibration_curve_{antibiotic}.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    # Save data to CSV
    prob_true, prob_pred = calibration_curve(y_true, y_prob, n_bins=10, strategy='uniform')
    cal_df = pd.DataFrame({
        'Mean_Predicted_Probability': prob_pred,
        'Fraction_Positives': prob_true
    })
    cal_df.to_csv(output_dir / f'05_calibration_curve_{antibiotic}.csv', index=False)
    
    plt.close()
    print(f"  ✓ Calibration curve saved: {output_path.name}")


def bootstrap_metric_ci(y_true, y_prob, metric_fn, n_bootstraps=1000,
                        random_state=42, alpha=0.05):
    """
    Estimate a percentile bootstrap confidence interval for a probability-based
    metric (e.g. ROC-AUC, average precision).

    A single train/test split yields a point estimate with no variance
    information — reviewers rightly object to "AUC = 0.99" reported without an
    interval (B15/P-05). Resampling the test set with replacement gives an
    empirical sampling distribution and a 95% CI, which is the minimal,
    tractable robustness check for an out-of-core pipeline (full k-fold CV
    remains future work).

    Args:
        y_true:        Ground-truth labels (array).
        y_prob:        Predicted probabilities (array).
        metric_fn:     Callable (y_true, y_prob) -> float.
        n_bootstraps:  Number of resamples.
        random_state:  Seed for reproducibility.
        alpha:         Significance level (0.05 -> 95% CI).

    Returns:
        tuple(point, lower, upper): point estimate and CI bounds (np.nan if the
        metric is undefined, e.g. a resample with a single class).
    """
    y_true = np.asarray(y_true)
    y_prob = np.asarray(y_prob)
    rng = np.random.RandomState(random_state)
    n = len(y_true)

    point = metric_fn(y_true, y_prob)
    scores = []
    for _ in range(n_bootstraps):
        idx = rng.randint(0, n, n)
        # Skip degenerate resamples that contain only one class
        if len(np.unique(y_true[idx])) < 2:
            continue
        scores.append(metric_fn(y_true[idx], y_prob[idx]))

    if not scores:
        return point, np.nan, np.nan
    lower = float(np.percentile(scores, 100 * (alpha / 2)))
    upper = float(np.percentile(scores, 100 * (1 - alpha / 2)))
    return float(point), lower, upper


def save_comprehensive_metrics(y_true, y_pred, y_prob, best_thresh, output_dir, antibiotic):
    """
    Export comprehensive evaluation metrics to CSV.
    
    Args:
        y_true (ndarray): Ground truth labels
        y_pred (ndarray): Predicted labels (at optimal threshold)
        y_prob (ndarray): Predicted probabilities
        best_thresh (float): Optimal classification threshold
        output_dir (Path): Directory to save the CSV
        antibiotic (str): Antibiotic name
    """
    # Confusion matrix components
    cm = confusion_matrix(y_true, y_pred)
    tn, fp, fn, tp = cm.ravel()
    
    # Compute all metrics
    metrics = {
        'Antibiotic': antibiotic.title(),
        'Optimal_Threshold': best_thresh,
        'Accuracy': accuracy_score(y_true, y_pred),
        'Balanced_Accuracy': balanced_accuracy_score(y_true, y_pred),
        'Sensitivity_Recall_TPR': recall_score(y_true, y_pred),
        'Specificity_TNR': tn / (tn + fp) if (tn + fp) > 0 else 0,
        'Precision_PPV': precision_score(y_true, y_pred),
        'NPV': tn / (tn + fn) if (tn + fn) > 0 else 0,
        'F1_Score': f1_score(y_true, y_pred),
        'MCC': matthews_corrcoef(y_true, y_pred),
        'Cohen_Kappa': cohen_kappa_score(y_true, y_pred),
        'ROC_AUC': roc_auc_score(y_true, y_prob),
        'PR_AUC': average_precision_score(y_true, y_prob),
        'True_Positives': int(tp),
        'True_Negatives': int(tn),
        'False_Positives': int(fp),
        'False_Negatives': int(fn),
        'Total_Samples': len(y_true),
        'Resistant_Count': int(np.sum(y_true)),
        'Susceptible_Count': int(len(y_true) - np.sum(y_true))
    }
    
    metrics_df = pd.DataFrame([metrics])
    output_path = output_dir / f'06_comprehensive_metrics_{antibiotic}.csv'
    metrics_df.to_csv(output_path, index=False, encoding='utf-8')
    
    print(f"\n  ✓ Comprehensive metrics saved: {output_path.name}")
    print("\nPublication-Ready Metrics Summary:")
    print("=" * 80)
    for key, value in metrics.items():
        if isinstance(value, float):
            print(f"{key:30s}: {value:.4f}")
        else:
            print(f"{key:30s}: {value}")
    print("=" * 80)


# ============================================================================
# MAIN ANALYSIS FUNCTION
# ============================================================================
def main():
    """
    Execute comprehensive model performance analysis.
    
    Pipeline:
        1. Load test data and trained model
        2. Generate predictions
        3. Optimize classification threshold
        4. Compute comprehensive performance metrics
        5. Perform error analysis
        6. Export results
    """
    # Set publication-quality matplotlib/seaborn styling
    sns.set_style('whitegrid')
    sns.set_context('paper', font_scale=1.2)
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
    
    print("=" * 80)
    print(f"COMPREHENSIVE MODEL ANALYSIS: {config['project']['model'].upper()}")
    print("=" * 80)
    
    # ------------------------------------------------------------------------
    # STEP 1: Validate and Load Files
    # ------------------------------------------------------------------------
    print("\n[STEP 1/6] Validating data files...")
    print("=" * 80)
    
    try:
        # Check label file
        y_path = MATRIX_DIR / f"y_{TARGET_ANTIBIOTIC}.csv"
        if not y_path.exists():
            raise FileNotFoundError(f"Label file not found: {y_path}")
        
        # Check model file - try possible model name variations
        model_paths_to_try = [
            MODELS_DIR / f"xgboost_{TARGET_ANTIBIOTIC}_final_v2.json",
            MODELS_DIR / f"xgboost_{TARGET_ANTIBIOTIC}.json",
            MODELS_DIR / f"xgboost_{TARGET_ANTIBIOTIC}_incremental.json"
        ]
        
        model_path = None
        for path in model_paths_to_try:
            if path.exists():
                model_path = path
                break
        
        if model_path is None:
            raise FileNotFoundError(
                "Model file not found. Tried:\n" +
                "\n".join([f"  - {p}" for p in model_paths_to_try])
            )
        
        print(f"  ✓ Label file: {y_path.name}")
        print(f"  ✓ Model file: {model_path.name}")
        
        # Load labels
        y_all = pd.read_csv(y_path, encoding='utf-8')['label'].values
        
        # Find matrix chunks
        chunk_files = sorted(
            list(MATRIX_DIR.glob(f"X_{TARGET_ANTIBIOTIC}_part_*.npz")),
            key=lambda x: int(x.stem.split('_')[-1])
        )
        
        if not chunk_files:
            raise FileNotFoundError(
                f"No matrix chunks found matching: X_{TARGET_ANTIBIOTIC}_part_*.npz"
            )
        
        print(f"  ✓ Matrix chunks: {len(chunk_files)} files")
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        print("\nPlease ensure you have run:")
        print("  1. Matrix creation: 03_matrix_construction.py")
        print("  2. Model training: 05_model_training.py")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Unexpected error: {e}")
        sys.exit(1)
    
    # ------------------------------------------------------------------------
    # STEP 2: Load Data and Model
    # ------------------------------------------------------------------------
    print("\n[STEP 2/6] Loading data and model...")

    # Load test filenames and threshold from config
    test_filenames, best_thresh = load_test_files_from_config()

    # Load test data
    X_test, y_test, ids_test = load_test_data(y_all, chunk_files, test_filenames)

    # Load model natively as Booster
    try:
        print(f"\nLoading model: {model_path.name}")
        model = xgb.Booster()
        model.load_model(str(model_path))
        print("  ✓ Model loaded successfully")
        print(f"  ✓ Number of features: {model.num_features()}")
        print(f"  ✓ Number of trees: {model.num_boosted_rounds()}")
        
    except Exception as e:
        print(f"ERROR: Failed to load model: {e}")
        sys.exit(1)
    
    # ------------------------------------------------------------------------
    # STEP 4: Generate Predictions
    # ------------------------------------------------------------------------
    print("\n[STEP 3/6] Generating predictions...")

    try:
        dtest = xgb.DMatrix(X_test, nthread=N_JOBS)
        y_prob = model.predict(dtest)
        print(f"  ✓ Predictions generated for {len(y_prob)} samples")
        print(f"  ✓ Probability range: [{y_prob.min():.4f}, {y_prob.max():.4f}]")
        
    except Exception as e:
        print(f"ERROR: Prediction failed: {e}")
        sys.exit(1)
    
    # ------------------------------------------------------------------------
    # STEP 3b: Operating threshold (NO DATA LEAKAGE)
    # ------------------------------------------------------------------------
    # CRITICAL FIX (was P-01 / B02): the previous version computed Youden's J on
    # the TEST set and OVERWROTE the config threshold. Optimising the decision
    # threshold on the same data you report metrics on is textbook data leakage
    # and inflates Sensitivity/Specificity — a guaranteed reviewer rejection.
    #
    # We now apply the threshold that was fixed BEFORE seeing the test set
    # (training-derived, written by 05_model_training.py and loaded from config).
    # Youden's J on the test set is still computed and printed for INFORMATION
    # ONLY — it is never applied and never written back to config.
    fpr_roc, tpr_roc, roc_thresholds = roc_curve(y_test, y_prob)
    youden_j = tpr_roc - fpr_roc
    optimal_idx = np.argmax(youden_j)
    youden_threshold_test = float(roc_thresholds[optimal_idx])

    print(f"\n[INFO] Operating threshold (from training, applied): {best_thresh:.4f}")
    print(f"[INFO] Test-set Youden's J threshold (reported only, NOT applied): {youden_threshold_test:.4f}")
    print("[INFO] Config is NOT modified by this script (data-leakage prevention).")

    # ------------------------------------------------------------------------
    # STEP 5: Applying the unbiased threshold and Performance Evaluation
    # ------------------------------------------------------------------------
    print("\n[STEP 4/6] Evaluating performance with the unbiased (training-derived) threshold...")
    print(f"  Applying Threshold: {best_thresh:.4f}")
    
    y_pred_opt = (y_prob >= best_thresh).astype(int)

    # M13: persist held-out test predictions for external-validation concordance
    # (16_external_concordance.py). These are the config test split — never seen
    # in training — so the model-vs-tool head-to-head stays leakage-free.
    try:
        ev_dir = OUTPUT_DIR.parent.parent / "external_validation"
        ev_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"Genome ID": ids_test, "y_true": y_test.astype(int),
                      "model_pred": y_pred_opt.astype(int), "model_prob": y_prob}).to_csv(
            ev_dir / f"16_model_preds_{TARGET_ANTIBIOTIC}.csv", index=False, encoding="utf-8")
        print(f"  ✓ Held-out test predictions saved for M13: 16_model_preds_{TARGET_ANTIBIOTIC}.csv")
    except Exception as e:
        print(f"  ⚠ Could not save M13 test predictions: {e}")

    # Compute comprehensive metrics
    print("\n" + "=" * 80)
    print("FINAL PERFORMANCE METRICS")
    print("=" * 80)
    
    acc = accuracy_score(y_test, y_pred_opt)
    balanced_acc = balanced_accuracy_score(y_test, y_pred_opt)
    roc_auc = roc_auc_score(y_test, y_prob)

    # Use average_precision_score for PR-AUC (consistent with
    # save_comprehensive_metrics at the bottom of this file). The trapezoidal
    # auc(recall, precision) over precision_recall_curve points can disagree
    # with AP due to interpolation, producing two different "PR-AUC" numbers in
    # the same report.
    pr_auc = average_precision_score(y_test, y_prob)

    kappa = cohen_kappa_score(y_test, y_pred_opt)
    mcc = matthews_corrcoef(y_test, y_pred_opt)
    
    # Bootstrap 95% confidence intervals (variance estimate for the single split)
    roc_auc_pt, roc_lo, roc_hi = bootstrap_metric_ci(y_test, y_prob, roc_auc_score)
    pr_auc_pt, pr_lo, pr_hi = bootstrap_metric_ci(y_test, y_prob, average_precision_score)

    print(f"Accuracy                   : {acc:.4f}")
    print(f"Balanced Accuracy          : {balanced_acc:.4f}")
    print(f"ROC AUC Score              : {roc_auc:.4f}  (95% CI: {roc_lo:.4f}–{roc_hi:.4f})")
    print(f"Precision-Recall AUC       : {pr_auc:.4f}  (95% CI: {pr_lo:.4f}–{pr_hi:.4f})")
    print(f"Cohen's Kappa              : {kappa:.4f}")
    print(f"Matthews Correlation Coef  : {mcc:.4f}")
    print(f"Operating Threshold        : {best_thresh:.4f}")
    print("=" * 80)

    # Persist the bootstrap CIs alongside the other evaluation artefacts
    pd.DataFrame([{
        'Antibiotic': TARGET_ANTIBIOTIC,
        'ROC_AUC': roc_auc_pt, 'ROC_AUC_CI_low': roc_lo, 'ROC_AUC_CI_high': roc_hi,
        'PR_AUC': pr_auc_pt, 'PR_AUC_CI_low': pr_lo, 'PR_AUC_CI_high': pr_hi,
        'n_bootstraps': 1000
    }]).to_csv(OUTPUT_DIR / f'08_bootstrap_ci_{TARGET_ANTIBIOTIC}.csv', index=False)

    # ------------------------------------------------------------------------
    # MLOps: write metrics.json (run provenance). Best-effort — never break eval.
    # ------------------------------------------------------------------------
    try:
        organism = get_target(config=config)[0]
        run_id = rm.make_run_id(organism, TARGET_ANTIBIOTIC)
        metrics_payload = {
            "run_id": run_id,
            "organism": organism,
            "antibiotic": TARGET_ANTIBIOTIC,
            "git_commit": rm.git_commit_hash(),
            "operating_threshold": float(best_thresh),
            "metrics": {
                "accuracy": float(acc),
                "balanced_accuracy": float(balanced_acc),
                "roc_auc": float(roc_auc), "roc_auc_ci": [roc_lo, roc_hi],
                "pr_auc": float(pr_auc), "pr_auc_ci": [pr_lo, pr_hi],
                "mcc": float(mcc), "cohen_kappa": float(kappa),
            },
        }
        run_dir = resolve_path('run_dir', organism=organism,
                               antibiotic=TARGET_ANTIBIOTIC, run_id=run_id, config=config)
        rm.write_json(run_dir / "metrics.json", metrics_payload)
        # also drop a copy next to the evaluation outputs for convenience
        rm.write_json(OUTPUT_DIR / f"09_metrics_{TARGET_ANTIBIOTIC}.json", metrics_payload)
    except Exception as e:
        print(f"  ⚠ Could not write metrics.json: {e}")
    
    print("\nDetailed Classification Report:")
    print("-" * 80)
    print(classification_report(
        y_test,
        y_pred_opt,
        target_names=['Susceptible (0)', 'Resistant (1)'],
        digits=4
    ))
    
    print("Confusion Matrix:")
    print("-" * 80)
    cm = confusion_matrix(y_test, y_pred_opt)
    print("                     Predicted")
    print("                   Susceptible  Resistant")
    print(f"Actual Susceptible     {cm[0,0]:<8}     {cm[0,1]:<8}")
    print(f"       Resistant       {cm[1,0]:<8}     {cm[1,1]:<8}")
    print("-" * 80)
    
    # Clinical interpretation
    tn, fp, fn, tp = cm.ravel()
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0
    
    print("\nClinical Performance Metrics:")
    print("-" * 80)
    print(f"Sensitivity (True Positive Rate):  {sensitivity:.4f}")
    print("  → Proportion of resistant bacteria correctly identified")
    print(f"Specificity (True Negative Rate):  {specificity:.4f}")
    print("  → Proportion of susceptible bacteria correctly identified")
    print(f"Positive Predictive Value (PPV):   {ppv:.4f}")
    print("  → When model predicts resistant, probability it's correct")
    print(f"Negative Predictive Value (NPV):   {npv:.4f}")
    print("  → When model predicts susceptible, probability it's correct")
    print("-" * 80)
    
    # ------------------------------------------------------------------------
    # STEP 6: Generate Publication-Quality Visualizations
    # ------------------------------------------------------------------------
    print("\n[STEP 5/6] Generating publication-quality outputs...")
    print("=" * 80)

    # Save comprehensive metrics CSV
    save_comprehensive_metrics(
        y_test, 
        y_pred_opt, 
        y_prob, 
        best_thresh, 
        OUTPUT_DIR, 
        TARGET_ANTIBIOTIC
    )
    
    # Generate all visualizations
    print("\nGenerating visualizations...")
    plot_confusion_matrix_enhanced(y_test, y_pred_opt, OUTPUT_DIR, TARGET_ANTIBIOTIC)
    plot_roc_curve_analysis(y_test, y_prob, OUTPUT_DIR, TARGET_ANTIBIOTIC)
    plot_precision_recall_curve_analysis(y_test, y_prob, OUTPUT_DIR, TARGET_ANTIBIOTIC)
    plot_probability_distribution(y_test, y_prob, OUTPUT_DIR, TARGET_ANTIBIOTIC, threshold=best_thresh)
    plot_calibration_curve_analysis(y_test, y_prob, OUTPUT_DIR, TARGET_ANTIBIOTIC)
    
    print("\n" + "=" * 80)
    print("ALL PUBLICATION-READY OUTPUTS GENERATED")
    print("=" * 80)
    
    # ------------------------------------------------------------------------
    # STEP 6: Error Analysis
    # ------------------------------------------------------------------------
    print("\n[STEP 6/6] Performing error analysis...")
    print("=" * 80)
    
    # False Negatives (most critical in medical context)
    fn_mask = (y_test == 1) & (y_pred_opt == 0)
    fn_ids = ids_test[fn_mask]
    fn_probs = y_prob[fn_mask]
    
    print(f"False Negatives (Missed Resistant): {len(fn_ids)}")
    
    if len(fn_ids) > 0:
        error_df = pd.DataFrame({
            'Genome_ID': fn_ids,
            'Predicted_Probability': fn_probs,
            'True_Label': 'Resistant',
            'Predicted_Label': 'Susceptible',
            'Error_Type': 'False_Negative'
        })
        
        print("\nTop 10 False Negatives (by probability):")
        print(error_df.sort_values('Predicted_Probability').head(10).to_string(index=False))
        
        error_file = OUTPUT_DIR / f"07_error_analysis_false_negatives_{TARGET_ANTIBIOTIC}.csv"
        error_df.to_csv(error_file, index=False, encoding='utf-8')
        print(f"\n  ✓ Full error analysis saved: {error_file}")
    
    # False Positives
    fp_mask = (y_test == 0) & (y_pred_opt == 1)
    fp_count = np.sum(fp_mask)
    print(f"\nFalse Positives (Overestimated): {fp_count}")
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print("\nResults saved to:")
    print(f"  {OUTPUT_DIR}")
    print("\nNext steps:")
    print("  1. Review false negatives - these are high-risk predictions")
    print("  2. Consider adjusting threshold if clinical context demands higher sensitivity")
    print("  3. Run feature extraction: 07_explainability.py")
    print("=" * 80)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    main()
