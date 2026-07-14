#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Exploratory Data Analysis (EDA) & Visualization Module for AMR Data

This script generates publication-quality visualizations for the raw genomic 
metadata (amr_phenotypes.csv) validated in 01_data_validation.py.
It provides visual confirmation of data integrity, class distributions and
missing data patterns before proceeding to computationally expensive k-mer
extraction and model training.

Visualizations generated:
1. Antibiotic Resistance Distribution (Horizontal Bar Plot)
2. Missing Data Heatmap (Sparsity of phenotypes)
3. Antibiotic Class Representation (Pie Chart)
4. Target Antibiotic Deep Dive (Class Imbalance visual)
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
import warnings
import gc

# Suppress minor seaborn/matplotlib warnings for clean output
warnings.filterwarnings('ignore')

# Set publication-ready seaborn style
sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.family'] = 'sans-serif'

# ============================================================================
# CONFIGURATION & PATHS
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

# Shared antibiotic classification — single source of truth in
# config/registry/antibiotics.yaml, accessed via the registry (SCALE_MLOPS_PLAN §3).
from lib.registry import load_antibiotic_classes
from lib.config import resolve_path
ANTIBIOTIC_CLASSES = load_antibiotic_classes()

try:
    with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    TARGET_ANTIBIOTIC = config.get('project', {}).get('target_antibiotic', 'unknown')
    ORGANISM = config.get('project', {}).get('organism', 'ecoli')
    MATRIX_FILE = resolve_path('metadata_file', organism=ORGANISM, config=config)
    OUTPUT_DIR = resolve_path('dir_01_data_exploration', organism=ORGANISM,
                              antibiotic=TARGET_ANTIBIOTIC, config=config)
except Exception as e:
    print(f"ERROR loading config: {e}")
    sys.exit(1)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# MEMORY & DATA LOADERS
# ============================================================================

def free_memory():
    """Forces garbage collection and closes any dangling matplot objects to avoid OOM on 8GB machines."""
    gc.collect()

def load_and_preprocess_data():
    """Load metadata and filter out unused categorical columns."""
    if not MATRIX_FILE.exists():
        print(f"ERROR: Metadata file not found at {MATRIX_FILE}")
        sys.exit(1)

    print(f"Loading metadata from: {MATRIX_FILE.name}")
    df = pd.read_csv(MATRIX_FILE)
    df['Genome ID'] = df['Genome ID'].astype(str)

    # Isolate antibiotic columns (drop Genome ID)
    antibiotics_df = df.drop(columns=['Genome ID'])

    # Robustly convert strictly to float representing resistance (0.0/1.0), setting errors to NaN
    for col in antibiotics_df.columns:
        antibiotics_df[col] = pd.to_numeric(antibiotics_df[col], errors='coerce')

    # Filter to antibiotics that actually have valid binary data
    valid_cols = [col for col in antibiotics_df.columns if antibiotics_df[col].notna().sum() > 0]
    df_clean = antibiotics_df[valid_cols]

    print(f"Loaded {len(df)} genomes and {len(valid_cols)} valid antibiotics.")
    return df, df_clean

# ============================================================================
# MAIN VISUALIZATION FUNCTIONS
# ============================================================================

def plot_resistance_distribution(df_clean):
    """Plot the number of Resistant (1) vs Susceptible (0) samples per antibiotic."""
    print("Generating Resistance Distribution Plot...")

    stats = []
    for col in df_clean.columns:
        counts = df_clean[col].value_counts()
        res = int(counts.get(1.0, 0))
        sus = int(counts.get(0.0, 0))
        total = res + sus
        if total > 0:
            stats.append({'Antibiotic': col, 'Resistant': res, 'Susceptible': sus, 'Total': total})

    stats_df = pd.DataFrame(stats).sort_values('Total', ascending=False).head(25)

    melted_df = pd.melt(stats_df, id_vars=['Antibiotic'], value_vars=['Resistant', 'Susceptible'],
                        var_name='Phenotype', value_name='Count')

    plt.figure(figsize=(14, 10))
    # Using 'hue' as mandated by latest seaborn stack syntax specifications
    ax = sns.histplot(data=melted_df, y='Antibiotic', hue='Phenotype', weights='Count',
                      multiple='stack', palette={'Resistant': '#d95f02', 'Susceptible': '#1f78b4'},
                      shrink=0.8)

    plt.title('Top 25 Antibiotics by Sample Size: Resistance vs Susceptibility', fontsize=16, pad=20)
    plt.xlabel('Number of Genomes', fontsize=12)
    plt.ylabel('Antibiotic', fontsize=12)
    sns.despine(left=True, bottom=True)

    # Highlight target antibiotic if it's in the top 25
    for i, label in enumerate(ax.get_yticklabels()):
        if label.get_text() == TARGET_ANTIBIOTIC:
            label.set_fontweight("bold")
            label.set_color("black")

    plt.tight_layout()
    output_path = OUTPUT_DIR / "01_resistance_distribution.png"
    plt.savefig(output_path)
    
    # Save underlying data to CSV
    csv_path = OUTPUT_DIR / "01_resistance_distribution.csv"
    stats_df.to_csv(csv_path, index=False)
    
    plt.close()

    del stats, stats_df, melted_df, ax
    free_memory()
    print(f" -> Saved: {output_path.name}")

def plot_missing_data_heatmap(df_clean):
    """Generate a heatmap showing the sparsity/missingness of AMR phenotypes."""
    print("Generating Missing Data Heatmap...")

    plot_df = df_clean.sample(min(5000, len(df_clean)), random_state=42)

    cols_sorted = plot_df.notna().sum().sort_values(ascending=False).index
    plot_df = plot_df[cols_sorted]
    plot_df = plot_df.iloc[:, :30]

    plt.figure(figsize=(16, 8))
    cmap = sns.color_palette("viridis", as_cmap=True)

    sns.heatmap(plot_df.notna(), cmap=cmap, cbar=False, yticklabels=False)

    plt.title('AMR Phenotype Data Density (Top 30 Antibiotics)', fontsize=16, pad=20)
    plt.xlabel('Antibiotics', fontsize=12)
    plt.ylabel(f'Genomes (n={len(plot_df)})', fontsize=12)
    plt.xticks(rotation=45, ha='right')

    plt.tight_layout()
    output_path = OUTPUT_DIR / "02_missing_data_heatmap.png"
    plt.savefig(output_path)
    
    # Save underlying data to CSV
    csv_path = OUTPUT_DIR / "02_missing_data_heatmap.csv"
    plot_df.notna().astype(int).to_csv(csv_path, index=True, index_label='Genome_Index')
    
    plt.close()

    del plot_df, cmap, cols_sorted
    free_memory()
    print(f" -> Saved: {output_path.name}")

def plot_antibiotic_classes(df_clean):
    """Visualize the distribution of data across major antibiotic classes."""
    print("Generating Antibiotic Class Representation Plot...")

    class_counts = {}
    for class_name, antibiotics in ANTIBIOTIC_CLASSES.items():
        total = 0
        for ab in antibiotics:
            if ab in df_clean.columns:
                total += df_clean[ab].notna().sum()
        if total > 0:
            class_counts[class_name] = total

    if not class_counts:
        print(" -> Warning: No class mappings matched the dataset. Skipping pie chart.")
        return

    labels = list(class_counts.keys())
    sizes = list(class_counts.values())

    colors = sns.color_palette('Set3', len(labels))

    plt.figure(figsize=(10, 10))
    explode = [0.05 if i == np.argmax(sizes) else 0 for i in range(len(sizes))]

    plt.pie(sizes, explode=explode, labels=None, colors=colors, autopct='%1.1f%%',
            shadow=False, startangle=140, pctdistance=0.85,
            textprops={'fontsize': 11, 'weight': 'bold'})

    centre_circle = plt.Circle((0,0), 0.70, fc='white')
    fig = plt.gcf()
    fig.gca().add_artist(centre_circle)

    plt.title('Data Point Distribution by Antibiotic Class', fontsize=16, pad=20)
    plt.legend(labels, loc="center left", bbox_to_anchor=(1, 0, 0.5, 1))

    plt.tight_layout()
    output_path = OUTPUT_DIR / "03_class_distribution_donut.png"
    plt.savefig(output_path)
    
    # Save underlying data to CSV
    csv_path = OUTPUT_DIR / "03_class_distribution_donut.csv"
    pd.DataFrame({'Antibiotic_Class': labels, 'Genome_Count': sizes}).to_csv(csv_path, index=False)
    
    plt.close()

    del class_counts, labels, sizes, colors, explode, fig
    free_memory()
    print(f" -> Saved: {output_path.name}")

def plot_target_antibiotic_deepdive(df_clean, target):
    """Create a focused visualization of the specific target antibiotic's balance."""
    print(f"Generating Deep-Dive for Target: {target}...")

    if target not in df_clean.columns:
        print(f" -> Warning: Target {target} not found in metadata. Skipping deep-dive.")
        return

    target_data = df_clean[target].dropna()
    res_count = (target_data == 1.0).sum()
    sus_count = (target_data == 0.0).sum()
    total = len(target_data)

    if total == 0:
        return

    plt.figure(figsize=(8, 6))

    ax = sns.barplot(
        x=['Susceptible (0)', 'Resistant (1)'],
        y=[sus_count, res_count],
        palette=['#1f78b4', '#d95f02'],
        hue=['Susceptible (0)', 'Resistant (1)'],
        legend=False
    )

    for i, v in enumerate([sus_count, res_count]):
        ax.text(i, v + (total*0.02), f"{v}\n({(v/total)*100:.1f}%)",
                ha='center', va='bottom', fontweight='bold', fontsize=12)

    plt.title(f'Class Balance: {target.upper()}', fontsize=16, pad=20)
    plt.ylabel('Number of Genomes', fontsize=12)

    min_ratio = 40.0
    if total >= 2000: min_ratio = 2.0
    elif total >= 1000: min_ratio = 5.0
    elif total >= 500: min_ratio = 10.0
    elif total >= 200: min_ratio = 20.0

    minority_pct = (min(res_count, sus_count) / total) * 100
    status = "VALID" if minority_pct >= min_ratio and min(res_count, sus_count) >= 40 else "INVALID"

    info_text = (
        f"Total Samples: {total}\n"
        f"Minority Class: {minority_pct:.1f}%\n"
        f"Required Threshold: ≥ {min_ratio:.1f}%\n"
        f"Status: {status}"
    )

    plt.annotate(info_text, xy=(0.95, 0.95), xycoords='axes fraction',
                 fontsize=11, ha='right', va='top', bbox=dict(boxstyle="round,pad=0.5", fc="white", ec="gray", alpha=0.9))

    plt.ylim(0, max(res_count, sus_count) * 1.2)
    sns.despine()
    plt.tight_layout()
    output_path = OUTPUT_DIR / f"04_target_deepdive_{target.replace('/', '_')}.png"
    plt.savefig(output_path)
    
    # Save underlying data to CSV
    csv_path = OUTPUT_DIR / f"04_target_deepdive_{target.replace('/', '_')}.csv"
    pd.DataFrame({'Phenotype': ['Susceptible (0)', 'Resistant (1)'], 'Count': [sus_count, res_count]}).to_csv(csv_path, index=False)
    
    plt.close()

    del target_data, ax
    free_memory()
    print(f" -> Saved: {output_path.name}")


# ============================================================================
# ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    print("=" * 60)
    print("AMR DATA EXPLORATION AND VISUALIZATION")
    print("=" * 60)

    full_df, clean_df = load_and_preprocess_data()

    print("\n[Running Visualizations]")
    plot_resistance_distribution(clean_df)
    plot_missing_data_heatmap(clean_df)
    plot_antibiotic_classes(clean_df)
    plot_target_antibiotic_deepdive(clean_df, TARGET_ANTIBIOTIC)

    print("\n" + "=" * 60)
    print(f"All visualizations saved to: {OUTPUT_DIR}")
    print("=" * 60)
