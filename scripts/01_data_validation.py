#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data Quality Assessment Module for AMR Prediction Project

This script performs comprehensive validation of genomic data and antimicrobial
resistance (AMR) metadata before the machine learning pipeline. It ensures data
integrity by cross-referencing genome files with metadata records and provides
statistical insights into antibiotic resistance distribution.

Academic Purpose:
    - Validates data availability and consistency
    - Identifies class imbalance issues in AMR phenotypes
    - Provides actionable recommendations for model training

Multi-Antibiotic Architecture:
    This script validates data for ALL antibiotics. Outputs are shared across
    all targets since it only performs validation without creating antibiotic-specific files.
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import sys
import traceback
import pandas as pd
import yaml
from pathlib import Path

# Shared antibiotic classification — single source of truth in
# config/registry/antibiotics.yaml, accessed via the registry (SCALE_MLOPS_PLAN §3).
from lib.registry import load_antibiotic_classes
from lib.config import resolve_path
ANTIBIOTIC_CLASSES = load_antibiotic_classes()

# ============================================================================
# CONFIGURATION: CROSS-PLATFORM COMPATIBLE PATHS
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

try:
    with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    CURRENT_TARGET = config.get('project', {}).get('target_antibiotic', 'unknown')
except Exception:
    CURRENT_TARGET = "unknown"
    config = None

if config:
    ORGANISM = config.get('project', {}).get('organism', 'ecoli')
    # Organism-aware path resolution (SCALE_MLOPS_PLAN §4.2). Reads the
    # {organism} templates from config 'paths_organism:' (falls back to 'paths:').
    MATRIX_FILE = resolve_path('metadata_file', organism=ORGANISM, config=config)
    GENOMES_DIR = resolve_path('raw_genomes_dir', organism=ORGANISM, config=config)
    REPORT_PATH = resolve_path('dir_global_exploration', organism=ORGANISM, config=config) / "validation_report.txt"
else:
    # Config failed to load — fall back to a best-effort E. coli path so the
    # script can still be imported; a real run requires config.yaml.
    BASE_DIR = PROJECT_ROOT / "data"
    MATRIX_FILE = BASE_DIR / "external" / "ecoli" / "metadata" / "amr_phenotypes.csv"
    GENOMES_DIR = BASE_DIR / "raw" / "ecoli" / "genomes"
    REPORT_PATH = BASE_DIR / "external" / "ecoli" / "metadata" / "validation_report.txt"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================
_log_file = None

def log(message: str = "") -> None:
    """
    Output message to both console and log file.

    Args:
        message (str): Message to output. Default is an empty string for blank lines.
    """
    print(message)
    if _log_file is not None and not _log_file.closed:
        try:
            _log_file.write(message + "\n")
            _log_file.flush()
        except IOError as e:
            print(f"WARNING: Failed to write to log file: {e}")

def validate_dataset_scientific(resistant_count: int, susceptible_count: int) -> tuple:
    """
    Validate dataset suitability for XGBoost training using scientific criteria.

    This function implements rigorous validation based on machine learning best practices
    for imbalanced genomic datasets. It ensures statistical reliability by enforcing:
    1. Class existence (both resistant and susceptible must be present)
    2. Absolute minority threshold (minimum 40 samples for hold-out validation)
    3. Dynamic imbalance tolerance (scales with dataset size)

    Args:
        resistant_count (int): Number of resistant (1) samples.
        susceptible_count (int): Number of susceptible (0) samples.

    Returns:
        tuple: (is_valid: bool, status_message: str)
            - is_valid: True if dataset meets all scientific criteria.
            - status_message: "VALID" or specific rejection reason.
    """
    if resistant_count == 0 or susceptible_count == 0:
        return (False, "MISSING CLASS")

    minority_count = min(resistant_count, susceptible_count)
    if minority_count < 40:
        return (False, "INSUFFICIENT MINORITY (<40)")

    total_count = resistant_count + susceptible_count
    minority_ratio = (minority_count / total_count) * 100

    if total_count >= 2000:
        min_ratio = 2.0
    elif total_count >= 1000:
        min_ratio = 5.0
    elif total_count >= 500:
        min_ratio = 10.0
    elif total_count >= 200:
        min_ratio = 20.0
    else:
        min_ratio = 40.0

    if minority_ratio < min_ratio:
        return (False, f"IMBALANCED (<{min_ratio:.1f}%)")

    return (True, "VALID")

# ============================================================================
# MAIN DATA VALIDATION FUNCTION
# ============================================================================
def check_data() -> None:
    """
    Perform comprehensive data quality assessment for AMR prediction project.

    This function orchestrates the complete data validation workflow:
    1. Loads antimicrobial resistance metadata from CSV
    2. Scans physical genome files (.fna format) in the filesystem
    3. Identifies the intersection of metadata and available genome files
    4. Computes and reports resistance statistics for each antibiotic
    5. Provides data balance recommendations for model training
    6. Saves complete report to file for permanent record

    Raises:
        SystemExit: If critical data files are missing or no usable data is found.
    """
    global _log_file

    try:
        REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
        _log_file = open(REPORT_PATH, 'w', encoding='utf-8')
    except Exception as e:
        print(f"WARNING: Could not open log file {REPORT_PATH}: {e}")
        print("Continuing with console output only...")
        _log_file = None

    try:
        log("=" * 80)
        log("DATA QUALITY ASSESSMENT - AMR PREDICTION PROJECT")
        log("=" * 80)
        log(f"Current config target: {CURRENT_TARGET}")
        log(f"Metadata file: {MATRIX_FILE}")
        log(f"Genome directory: {GENOMES_DIR}")
        log(f"Report saved to: {REPORT_PATH}")
        log("=" * 80)

        # 1. Load and Validate Metadata
        try:
            if not MATRIX_FILE.exists():
                raise FileNotFoundError(f"Metadata file not found: {MATRIX_FILE}")

            df = pd.read_csv(MATRIX_FILE, encoding='utf-8')
            df['Genome ID'] = df['Genome ID'].astype(str)
            log(f"✓ Metadata loaded successfully: {len(df)} records found")
        except Exception as e:
            log(f"ERROR: Failed loading metadata: {e}")
            sys.exit(1)

        # 2. Scan Physical Genome Files
        try:
            if not GENOMES_DIR.exists():
                raise FileNotFoundError(f"Genome directory not found: {GENOMES_DIR}")

            log("\nScanning genome files...")
            files = list(GENOMES_DIR.glob("*.fna"))
            file_ids = {f.stem for f in files}
            log(f"✓ Physical genome files found: {len(file_ids)}")
        except Exception as e:
            log(f"ERROR: Failed scanning genome files: {e}")
            sys.exit(1)

        # 3. Find Intersection
        common_ids = set(df['Genome ID']).intersection(file_ids)
        log("\n" + "=" * 80)
        log(f"USABLE DATA (intersection): {len(common_ids)} genomes")
        log("=" * 80)

        if not common_ids:
            log("CRITICAL ERROR: No matching genome IDs found!")
            sys.exit(1)

        df_final = df[df['Genome ID'].isin(common_ids)].copy()

        # 4. Compute Resistance Statistics
        antibiotics = [col for col in df_final.columns if col != 'Genome ID']
        stats_dict = {}
        for antibiotic in antibiotics:
            try:
                # Robust extraction: convert to numeric, coercing unparsable strings to NaN
                col_data = pd.to_numeric(df_final[antibiotic], errors='coerce')
                counts = col_data.value_counts()

                # After pd.to_numeric the index is float64, so 1.0 and 1 reference
                # the SAME bucket. Summing counts.get(1.0) + counts.get(1) would
                # double-count every resistant/susceptible isolate. Use the float
                # keys only (matches the correct logic in 01b/03).
                resistant = int(counts.get(1.0, 0))
                susceptible = int(counts.get(0.0, 0))
                total = resistant + susceptible

                if total > 0:
                    resistance_ratio = (resistant / total) * 100
                    is_valid, status = validate_dataset_scientific(resistant, susceptible)

                    stats_dict[antibiotic] = {
                        'total': total,
                        'resistant': resistant,
                        'susceptible': susceptible,
                        'ratio': resistance_ratio,
                        'status': status,
                        'is_valid': is_valid
                    }
            except Exception as e:
                log(f"Warning: Error processing antibiotic '{antibiotic}': {e}")

        # Display Categorized Reporting
        log("\n" + "=" * 110)
        log(f"{'ANTIBIOTIC':<35} | {'TOTAL':<8} | {'RES (1)':<10} | {'SUS (0)':<10} | {'RATIO (%)':<11} | {'STATUS':<18}")
        log("=" * 110)

        for class_name, antibiotic_list in ANTIBIOTIC_CLASSES.items():
            class_antibiotics = [(ab, stats_dict[ab]) for ab in antibiotic_list if ab in stats_dict]
            if class_antibiotics:
                log(f"\n[ {class_name.upper()} ]")
                log("-" * 110)
                class_antibiotics.sort(key=lambda x: x[1]['total'], reverse=True)
                for antibiotic, data in class_antibiotics:
                    log(f"{antibiotic:<35} | {data['total']:<8} | {data['resistant']:<10} | "
                        f"{data['susceptible']:<10} | {data['ratio']:>6.2f}%     | {data['status']:<18}")
        log("\n" + "=" * 110)

        # 5. Scientific Target Recommendation
        log("\n" + "=" * 80)
        log("ML TARGET RECOMMENDATION SYSTEM (SCIENTIFIC CRITERIA)")
        log("=" * 80)

        candidates = [
            {'name': ab, **data}
            for ab, data in stats_dict.items() if data.get('is_valid', False)
        ]

        if not candidates:
            log("CRITICAL: No scientifically valid ML targets found.")
            all_sorted = sorted(stats_dict.items(), key=lambda x: x[1]['total'], reverse=True)
            if all_sorted:
                best_fallback_name, best_fallback_data = all_sorted[0]
                log(f"  Best Fallback (highest count): {best_fallback_name} (Total: {best_fallback_data['total']}, Status: {best_fallback_data['status']})")
        else:
            candidates.sort(key=lambda x: -x['total'])
            best_target = candidates[0]
            minority_count = min(best_target['resistant'], best_target['susceptible'])

            log(f"Analysis identified {len(candidates)} scientifically valid candidates for Machine Learning.")
            log("\n" + "=" * 80)
            log("GLOBAL TOP RECOMMENDATION")
            log("=" * 80)
            log(f"TARGET: {best_target['name'].upper()}")
            log("   Reasoning: (Highest sample count among scientifically valid datasets)")
            log(f"   - Total Samples: {best_target['total']} (Maximum Statistical Power)")
            log(f"   - Resistant: {best_target['resistant']} | Susceptible: {best_target['susceptible']}")
            log(f"   - Minority Class: {minority_count} samples ({(minority_count/best_target['total'])*100:.1f}%)")
            log("   - Validation: PASSED")

            log("\n" + "=" * 80)
            log("TOP 2 VIABLE CANDIDATES PER ANTIBIOTIC CLASS")
            log("=" * 80)
            for class_name, antibiotic_list in ANTIBIOTIC_CLASSES.items():
                class_candidates = [c for c in candidates if c['name'] in antibiotic_list]
                if class_candidates:
                    class_candidates.sort(key=lambda x: -x['total'])
                    log(f"\n[ {class_name.upper()} ]")
                    for idx, candidate in enumerate(class_candidates[:2], 1):
                        minority = min(candidate['resistant'], candidate['susceptible'])
                        log(f"  {idx}. {candidate['name'].upper():<30} | Total: {candidate['total']:<5} | "
                            f"Minority: {minority:<4} ({(minority/candidate['total'])*100:>4.1f}%)")

            log("\n" + "=" * 80)
            log("ACTION REQUIRED")
            log("=" * 80)
            log(f"Current target in config.yaml: {CURRENT_TARGET.upper()}")
            if best_target['name'].lower() != CURRENT_TARGET.lower():
                log(f"RECOMMENDATION: Consider changing to {best_target['name'].upper()} for better performance")
            else:
                log("Current target matches recommendation.")

        log("\n" + "=" * 80)
        log("DATA QUALITY ASSESSMENT COMPLETE")
        log("=" * 80)

    except Exception as e:
        log(f"UNEXPECTED ERROR: {e}")
        log(traceback.format_exc())
    finally:
        if _log_file is not None and not _log_file.closed:
            _log_file.close()
            print(f"\n[Report saved to: {REPORT_PATH}]")

if __name__ == "__main__":
    check_data()
