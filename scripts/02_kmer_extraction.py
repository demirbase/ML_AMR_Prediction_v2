#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
K-mer Counting Module for Genomic Feature Extraction

This script performs k-mer counting on bacterial genome assemblies using KMC
(https://github.com/refresh-bio/KMC) tool. K-mers are short DNA subsequences used as genomic features
for machine learning-based antimicrobial resistance (AMR) prediction.

K-mer Analysis Background:
    K-mers capture local genomic patterns that may be associated with resistance
    genes or regulatory elements. By representing each genome as a k-mer frequency
    vector, we can apply machine learning algorithms to identify resistance patterns.
    
    K-mer length (k=21, set via config.yaml `preprocessing.k_length`) is chosen
    to balance (see METHODOLOGY.md for the full justification):
    - Specificity: Longer k-mers are more unique
    - Computational efficiency: Shorter k-mers are faster to process
    - Biological relevance: 21-mers capture gene-level features while keeping
      the feature space tractable

Multi-Antibiotic Architecture:
    K-mer databases are stored in antibiotic-specific directories to prevent
    collisions when analyzing multiple antibiotics simultaneously:
    - data/{antibiotic}/kmc_outputs/ for each target antibiotic
    - Changing 'target_antibiotic' in config.yaml automatically redirects outputs

Technical Notes:
    - Uses KMC in multi-FASTA mode for genome assemblies
    - Processes each genome independently to enable parallel computing
    - Outputs binary databases (.kmc_pre, .kmc_suf) for efficient storage
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import subprocess
import yaml
from pathlib import Path
from tqdm import tqdm
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

# Organism-aware path resolution (SCALE_MLOPS_PLAN §4.2)
from lib.config import get_target, resolve_path, resolve_tool
# AMR_ORGANISM env overrides config (parallel per-organism QC without a config
# edit — same mechanism as 03u). Organism-level step: antibiotic not needed.
ORGANISM = get_target(config=config)[0]

# Extract configuration values
K_LENGTH = config['preprocessing']['k_length']
MEMORY_GB = config['preprocessing']['kmc_mem']
THREADS = config['preprocessing']['threads']
# Minimum k-mer count threshold WITHIN a single genome. Fixed at 1 because
# inputs are assembled genomes (not raw reads): every k-mer that exists should
# be retained at this stage. Cross-genome rarity filtering (min_support) is
# applied later, globally, in 03_matrix_construction.py. Falls back to config
# if a `preprocessing.min_count` key is ever added.
MIN_COUNT = config['preprocessing'].get('min_count', 1)

# ============================================================================
# CROSS-PLATFORM COMPATIBLE PATHS (ANTIBIOTIC-SPECIFIC)
# ============================================================================
RAW_GENOMES_DIR = resolve_path('raw_genomes_dir', organism=ORGANISM, config=config)
AMR_MATRIX_PATH = resolve_path('metadata_file', organism=ORGANISM, config=config)

# Organism-scoped KMC output directory
KMC_OUTPUTS_DIR = resolve_path('kmc_outputs_dir', organism=ORGANISM, config=config)
TEMP_DIR = KMC_OUTPUTS_DIR / "tmp"

# KMC binary location — cross-platform: prefer the project-bundled binary if it
# exists and is executable for this OS, otherwise fall back to `kmc` on PATH
# (the normal case on Linux/HPC where KMC is installed via conda). Override with
# the AMR_KMC_BIN environment variable.
KMC_BIN = resolve_tool('kmc_bin', 'kmc', config=config)

# Create necessary directories if they don't exist
KMC_OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# MAIN K-MER COUNTING FUNCTION
# ============================================================================
def count_kmers():
    """
    Execute k-mer counting for all genome assemblies in the input directory.
    
    This function orchestrates the k-mer extraction pipeline:
    1. Validates KMC binary availability
    2. Discovers all genome assembly files (.fna format)
    3. Processes each genome independently using KMC
    4. Tracks progress and reports success/failure statistics
    
    KMC Workflow:
        For each genome:
        - Create a temporary file list (required by KMC list mode)
        - Execute KMC with optimized parameters
        - Clean up temporary files
        - Store results as binary database (.kmc_pre/.kmc_suf)
    
    Returns:
        None: Outputs progress to console and saves k-mer databases to disk
    
    Raises:
        SystemExit: If KMC binary is not found or no genome files exist
    """
    print("=" * 80)
    print("K-MER COUNTING OPERATION - GENOMIC FEATURE EXTRACTION")
    print("=" * 80)
    print(f"K-mer length: {K_LENGTH}")
    print(f"Minimum count threshold: {MIN_COUNT}")
    print(f"Memory allocation: {MEMORY_GB} GB")
    print(f"CPU threads: {THREADS}")
    print("=" * 80)
    
    # ------------------------------------------------------------------------
    # SECTION 1: Validate KMC Binary
    # ------------------------------------------------------------------------
    try:
        if not KMC_BIN:
            raise FileNotFoundError(
                "KMC executable not found.\n"
                "Install KMC (conda: `conda install -c bioconda kmc`) so `kmc` is on PATH,\n"
                "or set the AMR_KMC_BIN environment variable to its full path,\n"
                "or place the binary at bin/bin/kmc.\n"
                "Docs: https://github.com/refresh-bio/KMC"
            )
        print(f"✓ KMC binary located: {KMC_BIN}")
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # ------------------------------------------------------------------------
    # SECTION 2: Discover Genome Files
    # ------------------------------------------------------------------------
    try:
        # Find all FASTA nucleotide assembly files (.fna extension)
        # These represent bacterial genome assemblies from NCBI or similar databases
        genome_files = list(RAW_GENOMES_DIR.glob("*.fna"))
        total_files = len(genome_files)
        
        if total_files == 0:
            raise FileNotFoundError(
                f"No genome files (.fna) found in: {RAW_GENOMES_DIR}\n"
                f"Please ensure genome assemblies are placed in this directory."
            )
        
        print(f"✓ Genome files detected: {total_files}")
        print(f"  Input directory: {RAW_GENOMES_DIR}")
        print(f"  Output directory: {KMC_OUTPUTS_DIR}")
        print("=" * 80)
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    
    # ------------------------------------------------------------------------
    # SECTION 3: Process Each Genome
    # ------------------------------------------------------------------------
    # Track processing statistics
    success_count = 0
    failure_count = 0
    skipped_count = 0
    
    # Use tqdm for visual progress tracking
    for genome_file in tqdm(genome_files, desc="K-mer counting progress", unit="genome"):
        
        # Extract genome identifier from filename (e.g., "562.1234" from "562.1234.fna")
        genome_id = genome_file.stem
        
        output_db = KMC_OUTPUTS_DIR / genome_id
        
        # Skip if k-mer database already exists (resume capability).
        # A complete KMC database is BOTH .kmc_pre and .kmc_suf. Requiring both
        # prevents reusing a corrupt/partial database left behind by an
        # interrupted run (which would silently produce truncated features).
        kmc_pre = KMC_OUTPUTS_DIR / f"{genome_id}.kmc_pre"
        kmc_suf = KMC_OUTPUTS_DIR / f"{genome_id}.kmc_suf"
        if kmc_pre.exists() and kmc_suf.exists():
            skipped_count += 1
            continue
        
        # Create temporary input file list for KMC
        # KMC requires a text file containing input file paths (list mode)
        temp_list_file = TEMP_DIR / f"{genome_id}_input_list.txt"
        
        try:
            # Write genome file path to temporary list
            with open(temp_list_file, 'w', encoding='utf-8') as f:
                f.write(str(genome_file) + "\n")
            
            # Construct KMC command with optimized parameters
            kmc_command = [
                str(KMC_BIN),                     # KMC executable path
                f"-k{K_LENGTH}",                  # K-mer length
                f"-m{MEMORY_GB}",                 # Memory limit (GB)
                f"-t{THREADS}",                   # Number of threads
                f"-ci{MIN_COUNT}",                # Minimum k-mer count threshold
                "-fm",                            # Multi-FASTA mode (for assemblies)
                f"@{str(temp_list_file)}",        # Input file list (@ prefix for list mode)
                str(output_db),                   # Output database prefix
                str(TEMP_DIR)                     # Temporary directory for KMC workspace
            ]
            
            # Execute KMC as a subprocess
            # capture_output=True captures stdout/stderr for error reporting
            result = subprocess.run(
                kmc_command,
                capture_output=True,
                text=True,
                check=False  # Don't raise exception, handle errors manually
            )
            
            # Clean up temporary list file
            if temp_list_file.exists():
                temp_list_file.unlink()
            
            # Check execution status
            if result.returncode == 0:
                success_count += 1
            else:
                failure_count += 1
                # Extract first line of error message for concise reporting
                error_msg = result.stderr.splitlines()[0] if result.stderr else 'Unknown error'
                tqdm.write(f"✗ ERROR ({genome_id}): {error_msg}")
        
        except Exception as e:
            failure_count += 1
            tqdm.write(f"✗ SYSTEM ERROR ({genome_id}): {str(e)}")
            # Ensure cleanup even if error occurs
            if temp_list_file.exists():
                temp_list_file.unlink()
    
    # ------------------------------------------------------------------------
    # SECTION 4: Report Processing Statistics
    # ------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("K-MER COUNTING OPERATION COMPLETE")
    print("=" * 80)
    print(f"Successfully processed: {success_count} genomes")
    print(f"Skipped (already done): {skipped_count} genomes")
    print(f"Failed: {failure_count} genomes")
    print(f"Output location: {KMC_OUTPUTS_DIR}")
    print("=" * 80)
    
    # Suggest next steps
    if failure_count > 0:
        print("\nWARNING: Some genomes failed to process.")
        print("Check error messages above for details.")
        print("Common issues: corrupted files, insufficient memory, invalid FASTA format")
    
    if success_count > 0:
        print("\n✓ K-mer databases created successfully.")
        print(f"  Output directory: {KMC_OUTPUTS_DIR}")
        print("  Next step: Run matrix creation script (03_matrix_construction.py)")


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    count_kmers()
