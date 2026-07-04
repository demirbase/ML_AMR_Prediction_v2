#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Feature Importance Extraction Module

This script extracts and analyzes the most important k-mer features from a
trained XGBoost model for AMR prediction. Feature importance is measured using
Gain metric, which represents the average improvement in accuracy brought by
a feature across all splits where it was used.

Scientific Value:
    Identifying important k-mers helps:
    1. Understand genetic basis of resistance (biologically interpretable)
    2. Validate model (do important features correspond to known resistance genes?)
    3. Design targeted diagnostic assays
    4. Generate hypotheses for experimental validation

Gain vs Other Metrics:
    - Gain: Total improvement in accuracy (preferred for interpretation)
    - Weight: Number of times feature is used (may overweight redundant features)
    - Cover: Number of samples affected (biased toward common features)
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import xgboost as xgb
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
import sys


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

from lib.config import get_target  # early: env>config target before module globals

# Extract configuration values
TARGET_ANTIBIOTIC = get_target(config=config)[1]
ORGANISM = get_target(config=config)[0]
TOP_N = config['analysis']['top_n_features']
STABILITY_THRESHOLD = config.get('analysis', {}).get('stability_threshold', 0.6)

# Model filename to analyze
MODEL_FILE = f"xgboost_{TARGET_ANTIBIOTIC}_final_v2.json"

# Organism-aware paths (SCALE_MLOPS_PLAN §4.2)
from lib.config import resolve_path
MATRIX_DIR = resolve_path('matrix_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
MODELS_DIR = resolve_path('models_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
OUTPUT_DIR = resolve_path('dir_05_explainability', organism=ORGANISM,
                          antibiotic=TARGET_ANTIBIOTIC, config=config)

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ============================================================================
# STABILITY INTEGRATION (07b → 07)
# ============================================================================
def load_stability_map():
    """
    Load 07b's per-k-mer stability table if it exists.

    Returns {feature_index: (selection_frequency, mean_gain, kmer)}. Empty dict
    if 07b has not been run — in that case 07 falls back to gain-only output
    (fully backward compatible). For the stability-aware thesis run the order is
    05 → 07b → 07 → 08 → 09 so this file is present.
    """
    stab_path = OUTPUT_DIR / f"06_feature_stability_{TARGET_ANTIBIOTIC}.csv"
    if not stab_path.exists() or stab_path.stat().st_size == 0:
        return {}
    try:
        df = pd.read_csv(stab_path)
    except Exception:
        return {}
    out = {}
    for _, r in df.iterrows():
        try:
            out[int(r['feature_index'])] = (
                float(r['selection_frequency']),
                float(r.get('mean_gain', 0.0)),
                str(r['kmer']),
            )
        except (KeyError, TypeError, ValueError):
            continue
    return out


# ============================================================================
# FEATURE EXTRACTION FUNCTIONS
# ============================================================================
def extract_top_features():
    """
    Extract and export top k-mer features from trained XGBoost model.
    
    This function:
    1. Loads the trained XGBoost model
    2. Extracts feature importance scores (Gain metric)
    3. Maps feature indices to actual k-mer sequences
    4. Exports results in CSV and FASTA formats
    
    The Gain importance metric represents the average improvement in accuracy
    that each feature provides when used for splitting. Features with high
    gain are most crucial for resistance prediction.
    
    Output Files:
        - CSV: Ranked list with scores (for analysis)
        - FASTA: Sequences only (for BLAST/alignment tools)
    
    Returns:
        None: Saves results to disk and prints summary
    
    Raises:
        SystemExit: If model or feature dictionary files are missing
    """
    print("=" * 80)
    print(f"FEATURE IMPORTANCE ANALYSIS: {TARGET_ANTIBIOTIC.upper()}")
    print("=" * 80)
    
    # ------------------------------------------------------------------------
    # SECTION 1: Load Trained Model
    # ------------------------------------------------------------------------
    print("\n[STEP 1/4] Loading trained model...")
    
    try:
        model_path = MODELS_DIR / MODEL_FILE
        
        if not model_path.exists():
            raise FileNotFoundError(
                f"Model file not found: {model_path}\n"
                f"Please train the model first (05_model_training.py)"
            )
        
        # Load model natively as Booster
        model = xgb.Booster()
        model.load_model(str(model_path))
        
        print(f"  ✓ Model loaded: {MODEL_FILE}")
        print(f"  ✓ Total features in model: {model.num_features()}")
        print(f"  ✓ Total trees: {model.num_boosted_rounds()}")
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to load model: {e}")
        sys.exit(1)
    
    
    # ------------------------------------------------------------------------
    # SECTION 2: Extract Feature Importance Scores
    # ------------------------------------------------------------------------
    print(f"\n[STEP 2/4] Extracting top {TOP_N} features...")
    
    try:
        # Get feature importance using 'gain' metric
        # Gain = average improvement in accuracy when feature is used for splitting
        # This is the most interpretable metric for understanding feature contribution
        importance_dict = model.get_score(importance_type='gain')
        
        print(f"  ✓ Extracted importance scores for {len(importance_dict)} features")
        
        # Sort features by importance (descending) and take top N
        sorted_importance = sorted(
            importance_dict.items(),
            key=lambda x: x[1],
            reverse=True
        )[:TOP_N]
        
        print(f"  ✓ Selected top {len(sorted_importance)} features")

        # A model with no splits (e.g. a single-leaf stump on trivial data, or a
        # fully-regularised model) yields an EMPTY importance dict. Handle this
        # gracefully instead of crashing on sorted_importance[0]: write empty
        # CSV + FASTA so downstream steps (08/09) still find their inputs.
        if not sorted_importance:
            print("  ⚠ WARNING: the model exposes no feature importances (no tree splits).")
            print("    Writing empty top-feature outputs and stopping cleanly.")
            empty_cols = ['Rank', 'Feature_ID', 'Feature_Index', 'Gain_Score',
                          'Kmer_Sequence', 'Kmer_Length',
                          'in_gain_topN', 'selection_frequency', 'stable']
            csv_path = OUTPUT_DIR / f"01_top_{TOP_N}_features_{TARGET_ANTIBIOTIC}.csv"
            pd.DataFrame(columns=empty_cols).to_csv(csv_path, index=False, encoding='utf-8')
            fasta_path = OUTPUT_DIR / f"02_top_{TOP_N}_features_{TARGET_ANTIBIOTIC}.fasta"
            fasta_path.write_text("", encoding='utf-8')
            print(f"  ✓ Wrote empty outputs: {csv_path.name}, {fasta_path.name}")
            return

        print(f"  ✓ Importance range: {sorted_importance[0][1]:.2f} (max) to {sorted_importance[-1][1]:.2f} (min)")

    except Exception as e:
        print(f"ERROR: Failed to extract feature importance: {e}")
        sys.exit(1)
    
    
    # ------------------------------------------------------------------------
    # SECTION 3: Map Feature Indices to K-mer Sequences
    # ------------------------------------------------------------------------
    print("\n[STEP 3/4] Mapping features to k-mer sequences...")
    
    try:
        # XGBoost features are named as 'f0', 'f1', 'f2', etc.
        # These correspond to line numbers in features.txt (0-indexed)
        # Extract the indices we need
        needed_indices = set()
        for feat_name, score in sorted_importance:
            # Parse feature index from name (e.g., 'f123' -> 123).
            # Use [1:] to strip ONLY the leading 'f' prefix. replace('f','')
            # would also delete any 'f' inside the token — harmless for the
            # all-digit indices XGBoost emits, but fragile; [1:] is exact.
            idx = int(feat_name[1:])
            needed_indices.add(idx)
        
        print(f"  Feature indices to map: {len(needed_indices)}")
        
        # Load k-mer dictionary
        features_file = MATRIX_DIR / "features.txt"
        
        if not features_file.exists():
            raise FileNotFoundError(
                f"Feature dictionary not found: {features_file}\n"
                f"Please run matrix creation script first (03_matrix_construction.py)"
            )
        
        # Map indices to k-mer sequences
        # Read file line by line (memory efficient for large dictionaries)
        features_map = {}
        
        with open(features_file, 'r', encoding='utf-8') as f:
            for line_idx, line in enumerate(f):
                if line_idx in needed_indices:
                    # Format: "KMER_SEQUENCE COUNT"
                    kmer_sequence = line.split()[0]
                    features_map[line_idx] = kmer_sequence
                    
                    # Early exit if all features found (optimization)
                    if len(features_map) == len(needed_indices):
                        break
        
        print(f"  ✓ Successfully mapped {len(features_map)} k-mer sequences")
        
        # Verify all features were found
        if len(features_map) < len(needed_indices):
            missing = needed_indices - set(features_map.keys())
            print(f"  WARNING: {len(missing)} features not found in dictionary")
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to map features: {e}")
        sys.exit(1)
    
    
    # ------------------------------------------------------------------------
    # SECTION 4: Export Results
    # ------------------------------------------------------------------------
    print("\n[STEP 4/4] Exporting results...")
    
    try:
        # 07b stability (if present) — lets us flag/extend the candidate set with
        # k-mers that are reproducible across seeds (ROADMAP H1/H2). Both the
        # single-model gain top-N and the stable set are carried forward so the
        # biological validation (08/09) covers BOTH.
        stability_map = load_stability_map()
        gain_indices = {int(name[1:]) for name, _ in sorted_importance}

        # Prepare results data
        results_data = []
        fasta_lines = []

        for rank, (feat_name, score) in enumerate(sorted_importance, 1):
            # Extract feature index (strip only the leading 'f' prefix)
            idx = int(feat_name[1:])

            # Get k-mer sequence (use 'UNKNOWN' if not found)
            kmer_seq = features_map.get(idx, "UNKNOWN")
            sel_freq = stability_map.get(idx, (np.nan, None, None))[0]

            # Add to results table
            results_data.append({
                'Rank': rank,
                'Feature_ID': feat_name,
                'Feature_Index': idx,
                'Gain_Score': score,
                'Kmer_Sequence': kmer_seq,
                'Kmer_Length': len(kmer_seq) if kmer_seq != "UNKNOWN" else 0,
                'in_gain_topN': True,
                'selection_frequency': sel_freq,
                'stable': bool(sel_freq >= STABILITY_THRESHOLD) if pd.notna(sel_freq) else False,
            })

            # Format for FASTA output
            # FASTA header includes rank and importance score for reference
            fasta_header = f">Rank_{rank}|Score_{score:.4f}|Feature_{feat_name}"
            fasta_lines.append(f"{fasta_header}\n{kmer_seq}")

        # ---- Append STABLE k-mers (07b) not already in the gain top-N --------
        # These are reproducible across seeds but did not make the single model's
        # top-N; the thesis's H2 BLAST validation must include them.
        rank = len(results_data)
        n_stable_added = 0
        for idx, (sel_freq, mean_gain, kmer_seq) in sorted(
                stability_map.items(), key=lambda kv: kv[1][0], reverse=True):
            if idx in gain_indices or sel_freq < STABILITY_THRESHOLD:
                continue
            if not kmer_seq or kmer_seq == "UNKNOWN":
                continue
            rank += 1
            n_stable_added += 1
            feat_name = f"f{idx}"
            score = float(mean_gain or 0.0)
            results_data.append({
                'Rank': rank,
                'Feature_ID': feat_name,
                'Feature_Index': idx,
                'Gain_Score': score,
                'Kmer_Sequence': kmer_seq,
                'Kmer_Length': len(kmer_seq),
                'in_gain_topN': False,
                'selection_frequency': sel_freq,
                'stable': True,
            })
            fasta_lines.append(f">Rank_{rank}|Score_{score:.4f}|Feature_{feat_name}\n{kmer_seq}")

        if stability_map:
            print(f"  ✓ Stability merged: {n_stable_added} stable k-mer(s) "
                  f"added beyond the gain top-{TOP_N} (threshold ≥ {STABILITY_THRESHOLD})")

        # Create DataFrame
        results_df = pd.DataFrame(results_data)
        
        # Export to CSV (sequentially prefixed for consistent folder sort)
        csv_path = OUTPUT_DIR / f"01_top_{TOP_N}_features_{TARGET_ANTIBIOTIC}.csv"
        results_df.to_csv(csv_path, index=False, encoding='utf-8')
        print(f"  ✓ CSV saved: {csv_path}")
        
        # Export to FASTA (sequentially prefixed for consistent folder sort)
        fasta_path = OUTPUT_DIR / f"02_top_{TOP_N}_features_{TARGET_ANTIBIOTIC}.fasta"
        with open(fasta_path, 'w', encoding='utf-8') as f:
            f.write("\n".join(fasta_lines))
        print(f"  ✓ FASTA saved: {fasta_path}")
        
    except Exception as e:
        print(f"ERROR: Failed to export results: {e}")
        sys.exit(1)
    
    # ------------------------------------------------------------------------
    # Display Summary
    # ------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("FEATURE EXTRACTION COMPLETE")
    print("=" * 80)
    print("\nTop 10 Most Important Features:")
    print("-" * 80)
    print(results_df[['Rank', 'Feature_ID', 'Gain_Score', 'Kmer_Sequence']].head(10).to_string(index=False))
    print("-" * 80)
    
    print("\nOutput Files:")
    print(f"  1. {csv_path}")
    print("     → CSV format for analysis and visualization")
    print(f"  2. {fasta_path}")
    print("     → FASTA format for BLAST analysis")
    
    print("\nNext Steps for Biological Validation:")
    print("  1. BLAST Analysis:")
    print("     → Upload FASTA file to NCBI BLAST")
    print("     → Search against bacterial genomes database")
    print("     → Identify if k-mers match known resistance genes")
    print("  2. Literature Review:")
    print(f"     → Check if identified genes are documented for {TARGET_ANTIBIOTIC} resistance")
    print(f"     → Look for known {TARGET_ANTIBIOTIC} resistance mechanisms / target genes")
    print("  3. Experimental Validation:")
    print("     → Design primers targeting these k-mers")
    print("     → Validate presence in resistant vs susceptible isolates")
    print("=" * 80)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    extract_top_features()
