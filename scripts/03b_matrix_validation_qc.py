#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Feature Matrix Validation & Visualization Module (03_vis.py)

This script provides statistical validation and visualization of the final sparse binary 
feature matrix before it is fed into machine learning models (XGBoost).

Scientific Purpose:
- To visually confirm that the matrix is properly constructed (Sparsity check).
- To analyze the chunk-based memory distributions.
- To verify that the ML algorithm will receive a scientifically balanced dataset.

Visualizations generated:
1. Matrix Sparsity Distribution (How much of the DNA matrix is just "0"s)
2. Chunk Size & Feature Space (Memory/Processing map)
3. Cumulative Class Balance (Is the target variable ready for ML?)
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scipy.sparse as sp
from tqdm import tqdm
import warnings
import gc
from sklearn.decomposition import TruncatedSVD

warnings.filterwarnings('ignore')

# Set publication-ready seaborn style
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# ============================================================================
# CONFIGURATION & PATHS
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

try:
    with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    TARGET_ANTIBIOTIC = config['project']['target_antibiotic']

    MATRIX_DIR = PROJECT_ROOT / config['paths']['matrix_dir'].format(antibiotic=TARGET_ANTIBIOTIC)

    # Derive output directory from the centralised config key (02_matrix_qc)
    OUTPUT_DIR = PROJECT_ROOT / config['paths']['dir_02_matrix_qc'].format(
        antibiotic=TARGET_ANTIBIOTIC
    )

except Exception as e:
    print(f"ERROR loading config: {e}")
    sys.exit(1)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# DATA LOADING & STATISTICAL EXTRACTION
# ============================================================================

def analyze_matrix_structure():
    """
    Load metadata and iteratively analyze sparse matrix chunks to extract 
    statistical characteristics without loading everything into RAM.
    """
    print(f"Analyzing constructed feature matrices in: {MATRIX_DIR}")
    
    # 1. Load Labels & Validate Sizes
    y_file = MATRIX_DIR / f"y_{TARGET_ANTIBIOTIC}.csv"
    if not y_file.exists():
        print(f"ERROR: Label file not found: {y_file}")
        sys.exit(1)
        
    y_df = pd.read_csv(y_file)
    total_samples = len(y_df)
    res_count = y_df['label'].sum()
    sus_count = total_samples - res_count
    
    print(f"✓ Found {total_samples} Valid Genomes for ML.")
    
    # 2. Analyze Sparse Chunks (Memory efficient reading)
    chunk_files = sorted(list(MATRIX_DIR.glob(f"X_{TARGET_ANTIBIOTIC}_part_*.npz")), 
                         key=lambda x: int(x.stem.split('_part_')[1]))
                         
    if not chunk_files:
        print("ERROR: No matrix chunks (.npz) found.")
        sys.exit(1)
        
    print(f"✓ Found {len(chunk_files)} Matrix Chunks. Extracting statistics...")
    
    chunk_stats = []
    
    for chunk_file in tqdm(chunk_files, desc="Parsing Sparsity"):
        # Load sparse matrix
        X_chunk = sp.load_npz(chunk_file)
        
        # Calculate statistics
        rows, cols = X_chunk.shape
        non_zero_elements = X_chunk.nnz
        total_elements = rows * cols
        
        # Sparsity = (1.0 - (non-zero elements) / (total elements)) * 100
        sparsity_pct = (1.0 - (non_zero_elements / total_elements)) * 100
        
        # Memory tracking
        file_size_mb = chunk_file.stat().st_size / (1024 * 1024)
        
        chunk_stats.append({
            'Chunk': int(chunk_file.stem.split('_part_')[1]),
            'Genomes': rows,
            'Features': cols,
            'Sparsity_Pct': sparsity_pct,
            'File_Size_MB': file_size_mb,
            'Non_Zeros': non_zero_elements
        })
        
    return y_df, pd.DataFrame(chunk_stats)


# ============================================================================
# VISUALIZATION PLATFORM
# ============================================================================

def plot_class_balance(y_df):
    """
    Scientific Validation: Visualizes the final class balance going into XGBoost.
    Helps justify if SMOTE or class weighting is needed later.
    """
    print("Generating Final ML Class Balance...")
    output_path = OUTPUT_DIR / f"01_class_balance_{TARGET_ANTIBIOTIC}.png"
    if output_path.exists():
        print(f" -> Skipping: {output_path.name} already exists.")
        return
    
    total = len(y_df)
    res = int(y_df['label'].sum())
    sus = total - res
    
    plt.figure(figsize=(8, 6))
    
    # Assign x to `hue` + legend=False: required by seaborn >=0.14 to use a
    # per-category palette without triggering the deprecated-palette warning.
    ax = sns.barplot(x=['Susceptible (0)', 'Resistant (1)'], y=[sus, res],
                     hue=['Susceptible (0)', 'Resistant (1)'],
                     palette=['#1f78b4', '#d95f02'], legend=False)


    # Annotate with exact numbers and percentages
    for i, v in enumerate([sus, res]):
        pct = (v / total) * 100
        ax.text(i, v + (total * 0.02), f"n={v}\n({pct:.1f}%)", 
                ha='center', va='bottom', fontweight='bold', fontsize=12)
                
    # Calculate Shannon Entropy
    p_sus = sus / total
    p_res = res / total
    entropy = 0
    if p_sus > 0: entropy -= p_sus * np.log2(p_sus)
    if p_res > 0: entropy -= p_res * np.log2(p_res)
    
    # Add Statistical power context with formulas
    imbalance_ratio = max(sus, res) / min(sus, res) if min(sus, res) > 0 else float('inf')
    
    info_text = (
        f"Total N = {total}\n"
        f"Imbalance Ratio = {imbalance_ratio:.2f}:1\n"
        f"Shannon Entropy (H) = {entropy:.3f} bits\n"
        r"Formula: $H = - \sum P(x) \log_2 P(x)$"
    )
    if entropy < 0.5:
        info_text += "\n\nCRITICAL WARNING:\nLow Information Content (H < 0.5)"
        bbox_color = "#ffcccc"
        edge_color = "red"
    else:
        bbox_color = "#f8f9fa"
        edge_color = "gray"
        
    plt.annotate(info_text, xy=(0.95, 0.95), xycoords='axes fraction',
                 fontsize=11, ha='right', va='top', 
                 bbox=dict(boxstyle="round,pad=0.5", fc=bbox_color, ec=edge_color, alpha=0.9))
                 
    plt.title(f'Final ML Input: {TARGET_ANTIBIOTIC.upper()} Phenotypes', fontsize=16, pad=20)
    plt.ylabel('Number of Genomes', fontsize=12)
    plt.ylim(0, max(res, sus) * 1.2)
    
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f" -> Saved: {output_path.name}")


def plot_matrix_sparsity(chunk_df):
    """
    Scientific Validation: A perfectly engineered K-mer matrix should be extremely sparse (>90% zeros).
    If it's dense (e.g. 10% sparsity), the k-mer length or filtering logic is flawed.
    """
    print("Generating Matrix Sparsity Distribution...")
    output_path = OUTPUT_DIR / f"02_sparsity_check_{TARGET_ANTIBIOTIC}.png"
    if output_path.exists():
        print(f" -> Skipping: {output_path.name} already exists.")
        return
    
    plt.figure(figsize=(10, 6))
    
    # Line plot showing sparsity across all chunks
    ax1 = sns.lineplot(data=chunk_df, x='Chunk', y='Sparsity_Pct', marker='o', 
                       color='#9b59b6', linewidth=2, markersize=8)
                       
    mean_sparsity = chunk_df['Sparsity_Pct'].mean()
    plt.axhline(y=mean_sparsity, color='red', linestyle='--', alpha=0.7, 
                label=f'Avg Sparsity: {mean_sparsity:.2f}%')
                
    # Format Y-axis to clearly show it's near 100%
    y_min = max(0, min(chunk_df['Sparsity_Pct']) - 2)
    y_max = min(100, max(chunk_df['Sparsity_Pct']) + 2)
    plt.ylim(y_min, 100)
    
    plt.title('Feature Matrix Sparsity Control (QC Check)', fontsize=16, pad=20)
    plt.xlabel('Matrix Data Chunk Index', fontsize=12)
    plt.ylabel('Percentage of Zeros (%)', fontsize=12)
    
    # Calculate Theoretical Dense vs Actual Sparse Size
    first_chunk = chunk_df.iloc[0]
    theoretical_dense_gb = (first_chunk['Genomes'] * first_chunk['Features'] * 8) / (1024**3)
    actual_sparse_mb = first_chunk['File_Size_MB']
    
    # Add theoretical explanation with formula
    explanation = (
        "Biyolojik Geçerlilik (Biological Validity):\n"
        "Yüksek Sparsity (>%90) k-mer verisinin doğasında vardır.\n"
        r"$Sparsity = \left(1 - \frac{NonZeros}{Rows \times Cols}\right) \times 100$" "\n\n"
        f"Memory Proof (Chunk 0):\n"
        f"Theoretical Dense RAM: {theoretical_dense_gb:.2f} GB\n"
        f"Actual Sparse File: {actual_sparse_mb:.2f} MB"
    )
    plt.annotate(explanation, xy=(0.02, 0.05), xycoords='axes fraction',
                 fontsize=10, ha='left', va='bottom', 
                 bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="gray", alpha=0.8))
    
    plt.legend(loc='lower right')
    sns.despine()
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f" -> Saved: {output_path.name}")


def plot_chunk_memory_footprint(chunk_df):
    """
    Performance Validation: Shows how the data was split and physical memory sizes.
    Important for reproducibility and hardware requirements logic.
    """
    print("Generating Chunk Memory Profile...")
    output_path = OUTPUT_DIR / f"03_memory_profile_{TARGET_ANTIBIOTIC}.png"
    if output_path.exists():
        print(f" -> Skipping: {output_path.name} already exists.")
        return
    
    fig, ax1 = plt.subplots(figsize=(12, 6))
    
    color1 = '#3498db'
    ax1.set_xlabel('Matrix Data Chunk Index', fontsize=12)
    ax1.set_ylabel('Included Genomes (N)', color=color1, fontsize=12)
    # Bar plot for genomes per chunk
    sns.barplot(data=chunk_df, x='Chunk', y='Genomes', color=color1, ax=ax1, alpha=0.7)
    ax1.tick_params(axis='y', labelcolor=color1)
    
    # Hide some x labels if too many chunks
    for ind, label in enumerate(ax1.get_xticklabels()):
        if len(chunk_df) > 20 and ind % 5 != 0:
            label.set_visible(False)
            
    # Second y-axis for File Size
    ax2 = ax1.twinx()  
    color2 = '#e74c3c'
    ax2.set_ylabel('Disk Storage Size (MB)', color=color2, fontsize=12)
    sns.lineplot(data=chunk_df, x='Chunk', y='File_Size_MB', ax=ax2, 
                 color=color2, marker='s', linewidth=2, markersize=8)
    ax2.tick_params(axis='y', labelcolor=color2)
    
    total_gb = chunk_df['File_Size_MB'].sum() / 1024
    total_features = chunk_df['Features'].iloc[0] # Features are same for all chunks
    
    plt.title(f'Data Pipeline Storage & Memory Architecture ({TARGET_ANTIBIOTIC})', fontsize=16, pad=20)
    
    info_text = f"Total Features (K-mers): {total_features:,.0f}\nTotal Storage: {total_gb:.2f} GB"
    plt.annotate(info_text, xy=(0.5, 0.95), xycoords='axes fraction',
                 fontsize=11, ha='center', va='top', 
                 bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="gray", alpha=0.9))
                 
    fig.tight_layout()  
    plt.savefig(output_path)
    plt.close()
    print(f" -> Saved: {output_path.name}")


def plot_feature_prevalence(MATRIX_DIR, TARGET_ANTIBIOTIC, OUTPUT_DIR):
    """
    Scientific Validation: Analyzes the distribution of k-mer frequencies (prevalence)
    across genomes in the first matrix chunk. Reveals the balance between core and accessory genes.
    """
    print("\nGenerating Global Feature Prevalence Distribution...")
    output_path = OUTPUT_DIR / f"04_feature_prevalence_{TARGET_ANTIBIOTIC}.png"
    if output_path.exists():
        print(f" -> Skipping: {output_path.name} already exists.")
        return
    
    chunk_files = sorted(list(MATRIX_DIR.glob(f"X_{TARGET_ANTIBIOTIC}_part_*.npz")), 
                         key=lambda x: int(x.stem.split('_part_')[1]))
    if not chunk_files:
        return
        
    try:
        global_counts = None
        
        for f in tqdm(chunk_files, desc="Calculating Global Prevalence"):
            X_chunk = sp.load_npz(f)
            if global_counts is None:
                global_counts = np.zeros(X_chunk.shape[1], dtype=np.int32)
                
            # Sum down the columns (axis=0) to get frequency of each k-mer
            global_counts += np.array(X_chunk.sum(axis=0))[0].astype(np.int32)
            del X_chunk
            gc.collect()
            
        # We only want to plot non-zero counts (features actually present)
        active_features = global_counts[global_counts > 0]
        
        if len(active_features) == 0:
            print("  ⚠ No active features found across all chunks.")
            return

        plt.figure(figsize=(10, 6))
        
        # Plot distribution using log scale on x-axis to handle the massive range cleanly
        sns.histplot(active_features, bins=50, color='#e67e22', log_scale=(True, False))
        
        plt.title('Global K-mer Feature Prevalence (Across All Chunks)', fontsize=16, pad=20)
        plt.xlabel('Number of Genomes Sharing the K-mer (Log Scale)', fontsize=12)
        plt.ylabel('Count of Features', fontsize=12)
        
        # Add Statistical info
        median_prev = np.median(active_features)
        mean_prev = np.mean(active_features)
        
        explanation = (
            f"Mean Sub-Prevalence: {mean_prev:.1f} genomes\n"
            f"Median Sub-Prevalence: {median_prev:.1f} genomes\n\n"
            "İstatistiksel Yorum (Statistical Insight):\n"
            "Sol taraf: Sadece birkac suşta görülen nadir varyantlar (Accessory).\n"
            "Sağ kuyruk: Tüm türlerde ortak olan çekirdek diziler (Core Genome)."
        )
        plt.annotate(explanation, xy=(0.95, 0.95), xycoords='axes fraction',
                     fontsize=10, ha='right', va='top', 
                     bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="gray", alpha=0.9))
        
        sns.despine()
        plt.tight_layout()
        
        plt.savefig(output_path, dpi=300)
        plt.close()
        
        print(f" -> Saved: {output_path.name}")
        
    except Exception as e:
        print(f"  ⚠ Failed to generate Global Feature Prevalence: {e}")
        sys.exit(1)


def plot_svd_separability(MATRIX_DIR, TARGET_ANTIBIOTIC, y_data):
    print("\nGenerating Global SVD Separability Proof (Truncated SVD, all k-mers)...")
    output_path_2d = OUTPUT_DIR / f"05_svd_2d_separability_{TARGET_ANTIBIOTIC}.png"
    output_path_3d = OUTPUT_DIR / f"05_svd_3d_separability_{TARGET_ANTIBIOTIC}.png"

    if output_path_2d.exists() and output_path_3d.exists():
        print(" -> Skipping SVD: 2D and 3D plots already exist.")
        return

    chunk_files = sorted(list(MATRIX_DIR.glob(f"X_{TARGET_ANTIBIOTIC}_part_*.npz")),
                         key=lambda x: int(x.stem.split('_part_')[1]))
    if not chunk_files: return

    try:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers 3d projection)

        # ------------------------------------------------------------------
        # SCALABILITY FIX (was B08/P-12):
        # The previous implementation materialised the full N x N Gram matrix
        # (XX^T) as a dense float32 array, which is O(N^2) RAM — ~10 GB at
        # N=50,000 — and could OOM the machine. We instead run sklearn's
        # TruncatedSVD directly on the stacked SPARSE matrix. TruncatedSVD uses
        # randomized SVD and never densifies the feature space, so peak memory
        # scales with the (sparse) data, not N^2.
        # ------------------------------------------------------------------
        print("  [Pass 1/2] Stacking sparse chunks (kept sparse, no densification)...")
        X_stacked = sp.vstack([sp.load_npz(f) for f in tqdm(chunk_files, desc="Loading chunks")])
        X_stacked = X_stacked.astype(np.float32).tocsr()

        n_components = 3
        if min(X_stacked.shape) <= n_components:
            print(f"  ⚠ Matrix too small for {n_components} components; skipping SVD.")
            return

        print(f"  [Pass 2/2] Computing TruncatedSVD ({n_components} components)...")
        svd = TruncatedSVD(n_components=n_components, random_state=42)
        X_proj = svd.fit_transform(X_stacked)
        var_ratio = svd.explained_variance_ratio_

        del X_stacked
        gc.collect()

        N = X_proj.shape[0]
        y_global = y_data['label'].iloc[:N].values

        print("  Generating Plots...")
        # 2D Plot
        plt.figure(figsize=(10, 8))
        scatter = sns.scatterplot(
            x=X_proj[:, 0], y=X_proj[:, 1], hue=y_global, 
            palette=['#1f78b4', '#d95f02'], alpha=0.7, s=80, edgecolor='k'
        )
        plt.title('Exact Global SVD Separability (100% Data, All K-mers)', fontsize=16, pad=20)
        plt.xlabel(f'Principal Component 1 ({var_ratio[0]*100:.2f}%)')
        plt.ylabel(f'Principal Component 2 ({var_ratio[1]*100:.2f}%)')
        handles, labels = scatter.get_legend_handles_labels()
        plt.legend(handles, ['Susceptible (0)', 'Resistant (1)'], title="Phenotype")
        sns.despine()
        plt.tight_layout()
        plt.savefig(output_path_2d, dpi=300)
        plt.close()
        print(f" -> Saved 2D plot: {output_path_2d.name}")
        
        # 3D Plot
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        mask_0, mask_1 = (y_global == 0), (y_global == 1)
        ax.scatter(X_proj[mask_0, 0], X_proj[mask_0, 1], X_proj[mask_0, 2], 
                   c='#1f78b4', alpha=0.7, s=50, edgecolor='k', label='Susceptible (0)')
        ax.scatter(X_proj[mask_1, 0], X_proj[mask_1, 1], X_proj[mask_1, 2], 
                   c='#d95f02', alpha=0.7, s=50, edgecolor='k', label='Resistant (1)')
        ax.set_title('Exact Global SVD 3D Separability (100% Data)', fontsize=16)
        ax.set_xlabel(f'PC 1 ({var_ratio[0]*100:.2f}%)')
        ax.set_ylabel(f'PC 2 ({var_ratio[1]*100:.2f}%)')
        ax.set_zlabel(f'PC 3 ({var_ratio[2]*100:.2f}%)')
        ax.legend(title="Phenotype", loc='best')
        plt.tight_layout()
        plt.savefig(output_path_3d, dpi=300)
        plt.close()
        print(f" -> Saved 3D plot: {output_path_3d.name}")
        
    except Exception as e:
        print(f"  ⚠ Failed to generate Truncated SVD: {e}")
    finally:
        # Locals (X_stacked, X_proj, svd) are released when the function
        # returns; an explicit collection here reclaims the large arrays
        # promptly on memory-constrained machines.
        gc.collect()

# ============================================================================
# ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("FEATURE MATRIX (03) STATISTICAL VALIDATION & VISUALIZATION")
    print("=" * 60)
    
    y_data, chunk_data = analyze_matrix_structure()
    
    print("\n[Running Visualizations]")
    plot_class_balance(y_data)
    plot_matrix_sparsity(chunk_data)
    plot_chunk_memory_footprint(chunk_data)
    plot_feature_prevalence(MATRIX_DIR, TARGET_ANTIBIOTIC, OUTPUT_DIR)
    plot_svd_separability(MATRIX_DIR, TARGET_ANTIBIOTIC, y_data)

    print("\n" + "=" * 60)
    print(f"All visualizations saved to: {OUTPUT_DIR}")
    print("=" * 60)

