#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Feature Matrix Construction Module for AMR Prediction

This script transforms k-mer count data into a sparse binary feature matrix
suitable for machine learning. The matrix represents presence/absence of k-mers
across all genomes, where rows are samples (genomes) and columns are features (k-mers).

Scientific Rationale:
    Binary encoding (presence/absence) is used instead of frequency counts because:
    1. Genome assemblies have uniform coverage (unlike raw reads)
    2. Binary features reduce noise and computational complexity
    3. Presence of resistance-associated k-mers is more important than count
    
    Feature selection via minimum support threshold removes rare k-mers that:
    - May represent sequencing errors or strain-specific variations
    - Don't contribute to generalizable resistance patterns
    - Increase dimensionality without improving predictive power

Multi-Antibiotic Architecture:
    Outputs are organized by antibiotic name to prevent collisions:
    - data/{antibiotic}/matrix/ for feature matrices
    - Filenames include antibiotic-specific prefix (e.g., y_{prefix}.csv)
    - Changing 'target_antibiotic' in config.yaml automatically redirects outputs

Memory Management Strategy:
    To handle large datasets (thousands of genomes, millions of k-mers), we:
    1. Process genomes in chunks to avoid RAM overflow
    2. Use sparse matrix format (CSR) to store binary data efficiently
    3. Perform garbage collection after each chunk
    4. Save intermediate results to disk for fault tolerance
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from scipy.sparse import csr_matrix, save_npz
from tqdm import tqdm
import sys
import gc  # Garbage collector for explicit memory management
import array

# Shared safe subprocess wrapper (shlex-based, NO shell=True) — see scripts/utils.py
from utils import run_command


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

# Extract target via get_target: CLI-arg > AMR_ORGANISM/AMR_ANTIBIOTIC env >
# config (parallel per-(organism,antibiotic) runs without a config edit, like 03u).
ORGANISM, TARGET_ANTIBIOTIC = get_target(config=config)
K_LENGTH          = config['preprocessing']['k_length']   # Must match 02_kmer_extraction.py
# Data-adaptive minimum support (see config comments). MIN_SUPPORT is an optional
# hard override; when null the effective value is derived from the genome count
# at runtime via resolve_min_support().
MIN_SUPPORT       = config['preprocessing'].get('min_support', None)
MIN_SUPPORT_FLOOR = int(config['preprocessing'].get('min_support_floor', 5))
MIN_PREVALENCE    = float(config['preprocessing'].get('min_prevalence', 0.0))
CHUNK_SIZE        = config['preprocessing']['chunk_size']
KMC_MEMORY_GB     = config['preprocessing']['kmc_mem']
THREADS           = config['preprocessing']['threads']

# ============================================================================
# CROSS-PLATFORM COMPATIBLE PATHS (ANTIBIOTIC-SPECIFIC)
# ============================================================================
# Organism-scoped paths (resolved via lib.config)
RAW_GENOMES_DIR = resolve_path('raw_genomes_dir', organism=ORGANISM, config=config)
METADATA_FILE = resolve_path('metadata_file', organism=ORGANISM, config=config)

# K-mer specific directories
KMC_OUTPUTS_DIR = resolve_path('kmc_outputs_dir', organism=ORGANISM, config=config)
MATRIX_OUTPUT_DIR = resolve_path('matrix_dir', organism=ORGANISM, antibiotic=TARGET_ANTIBIOTIC, config=config)
TEMP_DIR = KMC_OUTPUTS_DIR / "tmp"

# KMC binaries — PATH-aware (conda/module on Linux/HPC; bundled macOS binary is
# only a Darwin fallback). Override with AMR_KMC_BIN / AMR_KMC_TOOLS_BIN.
KMC_BIN = resolve_tool('kmc_bin', 'kmc', config=config)
KMC_TOOLS_BIN = resolve_tool('kmc_tools_bin', 'kmc_tools', config=config)
if not KMC_BIN or not KMC_TOOLS_BIN:
    sys.exit("ERROR: KMC not found. Install KMC (conda install -c bioconda kmc) so "
             "`kmc`/`kmc_tools` are on PATH, or set AMR_KMC_BIN / AMR_KMC_TOOLS_BIN.")

# Create output directories
MATRIX_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
TEMP_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================
# run_command() is imported from utils.py. It tokenises with shlex and runs
# with shell=False, eliminating the previous shell-injection vector (a genome
# ID / path containing shell metacharacters can no longer execute commands).


# ============================================================================
# MAIN FEATURE MATRIX CONSTRUCTION
# ============================================================================
def create_feature_matrix():
    """
    Construct sparse binary feature matrix from k-mer data.
    
    This function orchestrates the complete feature engineering pipeline:
    
    STEP 1: Data Validation
        - Load AMR metadata and filter for target antibiotic
        - Verify k-mer databases exist for all genomes
        
    STEP 2: Global Feature Dictionary Construction
        - Create a unified k-mer vocabulary across ALL genomes
        - Apply minimum support threshold to filter rare k-mers
        - This ensures consistent feature space for all samples
        
    STEP 3: Chunked Matrix Construction
        - Process genomes in batches to manage memory
        - For each genome, extract k-mers and map to feature indices
        - Store as sparse binary matrix (CSR format)
        - Save chunk to disk and free memory before next chunk
        
    STEP 4: Save Metadata
        - Store genome IDs and resistance labels separately
        - These will be used to reconstruct full dataset during training
    
    Returns:
        None: Saves matrix chunks and metadata to disk
    
    Raises:
        SystemExit: If critical files are missing or processing fails
    """
    print("=" * 80)
    print("FEATURE MATRIX CONSTRUCTION - CHUNKED PROCESSING")
    print("=" * 80)
    print(f"Target antibiotic: {TARGET_ANTIBIOTIC}")
    if MIN_SUPPORT is not None:
        print(f"Minimum support: {MIN_SUPPORT} genomes (config override)")
    else:
        print(f"Minimum support: adaptive = max({MIN_SUPPORT_FLOOR}, "
              f"ceil({MIN_PREVALENCE} x n_genomes)) — resolved after QC")
    print(f"Chunk size: {CHUNK_SIZE} genomes")
    print("=" * 80)

    # ------------------------------------------------------------------------
    # SECTION 1: Data Validation and Preparation
    # ------------------------------------------------------------------------
    print("\n[STEP 1/4] Loading and validating metadata...")
    
    try:
        if not METADATA_FILE.exists():
            raise FileNotFoundError(f"Metadata file not found: {METADATA_FILE}")
        
        # Load resistance metadata
        metadata_df = pd.read_csv(METADATA_FILE, encoding='utf-8')
        metadata_df['Genome ID'] = metadata_df['Genome ID'].astype(str)
        
        # Filter for genomes with known resistance status for target antibiotic
        # NaN values are removed to ensure clean training data
        metadata_filtered = metadata_df.dropna(subset=[TARGET_ANTIBIOTIC]).copy()
        
        print(f"  ✓ Total genomes in metadata: {len(metadata_df)}")
        print(f"  ✓ Genomes with {TARGET_ANTIBIOTIC} labels: {len(metadata_filtered)}")
        
        # --- LOAD QC OUTLIERS ---
        outlier_file = resolve_path('dir_global_exploration', organism=ORGANISM, config=config) / "global_qc_outliers.csv"
        outlier_ids = set()
        if outlier_file.exists():
            outliers_df = pd.read_csv(outlier_file)
            if 'Genome' in outliers_df.columns:
                outlier_ids = set(outliers_df['Genome'].astype(str))
            else:
                outlier_ids = set(outliers_df.iloc[:, 0].astype(str))
            print(f"  ✓ Loaded {len(outlier_ids)} QC outlier genomes to exclude.")
        else:
            print(f"  ⚠ Warning: Outlier file not found at {outlier_file}")
            
        # Extract genome IDs and resistance labels
        all_genome_ids = metadata_filtered['Genome ID'].values
        all_labels = metadata_filtered[TARGET_ANTIBIOTIC].astype(int).values
        
    except FileNotFoundError as e:
        print(f"ERROR: {e}")
        sys.exit(1)
    except KeyError:
        print(f"ERROR: Column '{TARGET_ANTIBIOTIC}' not found in metadata.")
        print(f"Available columns: {metadata_df.columns.tolist()}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Unexpected error while loading metadata: {e}")
        sys.exit(1)
    
    # Verify k-mer databases AND raw genome FASTA files exist for each genome.
    # Both are required: .kmc_pre for matrix population (Step 3) and .fna for
    # the global vocabulary construction pass (Step 2).
    valid_genomes = []
    valid_labels = []
    missing_fna_count = 0
    missing_kmc_count = 0
    skipped_outliers_count = 0
    
    print("  Verifying k-mer databases and raw genome files...")
    for i, genome_id in enumerate(all_genome_ids):
        # Explicitly REMOVE any genome ID that is present in the outlier set
        if genome_id in outlier_ids:
            skipped_outliers_count += 1
            continue
            
        kmc_db_file = KMC_OUTPUTS_DIR / f"{genome_id}.kmc_pre"
        fna_file    = RAW_GENOMES_DIR  / f"{genome_id}.fna"
        
        if not kmc_db_file.exists():
            missing_kmc_count += 1
            continue  # Cannot populate matrix without k-mer database
        if not fna_file.exists():
            missing_fna_count += 1
            continue  # Cannot include in global vocabulary pass without FASTA
        
        valid_genomes.append(genome_id)
        valid_labels.append(all_labels[i])
    
    if skipped_outliers_count > 0:
        print(f"  ✓ Skipped {skipped_outliers_count} genomes: Explicitly removed due to QC outlier blacklist.")
    if missing_kmc_count > 0:
        print(f"  ⚠ Skipped {missing_kmc_count} genomes: .kmc_pre database missing (run 02_kmer_extraction.py).")
    if missing_fna_count > 0:
        print(f"  ⚠ Skipped {missing_fna_count} genomes: raw .fna assembly missing in {RAW_GENOMES_DIR}.")
    
    if len(valid_genomes) == 0:
        print("ERROR: No genomes passed validation!")
        print(f"  K-mer output directory: {KMC_OUTPUTS_DIR}")
        print(f"  Raw genome directory:   {RAW_GENOMES_DIR}")
        print("Please run k-mer counting script first (02_kmer_extraction.py)")
        sys.exit(1)
    
    print(f"  ✓ Valid genomes with k-mer data: {len(valid_genomes)}")
    
    # Report class distribution
    resistant_count = sum(valid_labels)
    susceptible_count = len(valid_labels) - resistant_count
    print(f"  ✓ Class distribution: {resistant_count} resistant, {susceptible_count} susceptible")
    print(f"    (Ratio: {resistant_count/len(valid_labels)*100:.1f}% resistant)")

    # ------------------------------------------------------------------------
    # SECTION 2: Global Feature Dictionary Construction
    # ------------------------------------------------------------------------
    print("\n[STEP 2/4] Constructing global k-mer vocabulary...")
    
    features_file = MATRIX_OUTPUT_DIR / "features.txt"
    
    if features_file.exists():
        print("  ✓ Feature dictionary already exists, loading from disk...")
    else:
        print("  Creating unified k-mer vocabulary across all genomes...")
        print(f"  (This may take several minutes for {len(valid_genomes)} genomes)")
        
        # Create a text file listing all genome FASTA files.
        # All genomes in valid_genomes are guaranteed to have their .fna files
        # present (verified in Section 1), so no existence re-check is needed.
        # This allows KMC to process all genomes in a single, efficient pass.
        global_genome_list = TEMP_DIR / "global_genome_list.txt"
        global_db = TEMP_DIR / "global_features_db"
        
        try:
            with open(global_genome_list, 'w', encoding='utf-8') as f:
                for genome_id in valid_genomes:
                    genome_file = RAW_GENOMES_DIR / f"{genome_id}.fna"
                    f.write(str(genome_file) + "\n")
            
            # Filter out zero-variance features (core genome k-mers present in 100% of samples)
            # If a k-mer is in all genomes, its frequency is len(valid_genomes). We set max allowed to len - 1.
            max_support = len(valid_genomes) - 1

            # Data-adaptive minimum support: scale with the genome count so the
            # same config works across antibiotics/organisms of very different
            # sizes (see config comments). Explicit MIN_SUPPORT overrides.
            if MIN_SUPPORT is not None:
                eff_min_support = int(MIN_SUPPORT)
                _ms_src = "config override"
            else:
                from math import ceil
                eff_min_support = max(MIN_SUPPORT_FLOOR,
                                      ceil(MIN_PREVALENCE * len(valid_genomes)))
                _ms_src = (f"adaptive: max(floor={MIN_SUPPORT_FLOOR}, "
                           f"ceil({MIN_PREVALENCE}*{len(valid_genomes)}))")
            eff_min_support = min(eff_min_support, max(1, max_support))  # never exceed max

            # Run KMC on ALL genomes to identify k-mers meeting minimum support.
            # K-mer length uses K_LENGTH from config (consistent with 02_kmer_extraction.py).
            #
            # IMPORTANT — semantics caveat:
            #   -ci / -cx filter on the TOTAL OCCURRENCE COUNT of a k-mer across the
            #   concatenated multi-FASTA input, NOT on the strict number of distinct
            #   genomes it appears in (document frequency). For assembled bacterial
            #   genomes most informative k-mers occur ~once per genome, so this
            #   occurrence count is a close PROXY for genome prevalence. It can,
            #   however, over-count k-mers that repeat within a single genome
            #   (e.g. multi-copy elements, rRNA operons). This is an accepted
            #   approximation here; a strict document-frequency filter would require
            #   per-genome canonicalisation and is noted as future work.
            # -ci{MIN_SUPPORT} filters rare k-mers (occurrence count < MIN_SUPPORT).
            kmc_cmd = (
                f"{KMC_BIN} -k{K_LENGTH} -m{KMC_MEMORY_GB} -t{THREADS} "
                f"-ci{eff_min_support} -cx{max_support} -fm @{global_genome_list} "
                f"{global_db} {TEMP_DIR}"
            )

            print("  Running KMC to extract global k-mer vocabulary...")
            print(f"  Filtering parameters: min_support={eff_min_support} ({_ms_src}), "
                  f"max_support={max_support} (zero-variance core)")
            run_command(kmc_cmd)
            
            # Export k-mer database to text format (k-mer sequence + count)
            print("  Exporting k-mer dictionary to text format...")
            dump_cmd = f"{KMC_TOOLS_BIN} transform {global_db} dump {features_file}"
            run_command(dump_cmd)
            
        finally:
            # Guarantee cleanup of temporary files regardless of success or failure.
            # global_genome_list is no longer needed after KMC completes.
            # global_db.* (binary KMC database) is superseded by the text dump.
            if global_genome_list.exists():
                global_genome_list.unlink()
            for suffix in [".kmc_pre", ".kmc_suf"]:
                db_part = global_db.parent / (global_db.name + suffix)
                if db_part.exists():
                    db_part.unlink()
        
        print(f"  ✓ Feature dictionary saved to: {features_file}")
    
    # Load feature dictionary into memory
    # This creates a mapping from k-mer sequence to column index in the matrix
    print("  Loading feature dictionary into memory...")
    kmer_to_index = {}
    
    try:
        with open(features_file, 'r', encoding='utf-8') as f:
            for index, line in enumerate(f):
                # Format: "KMER_SEQUENCE COUNT"
                kmer_sequence = line.split()[0]
                kmer_to_index[kmer_sequence] = index
        
        num_features = len(kmer_to_index)
        print(f"  ✓ Feature dictionary loaded: {num_features:,} k-mers")
        
        if num_features == 0:
            print("ERROR: Feature dictionary is empty!")
            sys.exit(1)
            
    except Exception as e:
        print(f"ERROR: Failed to load feature dictionary: {e}")
        sys.exit(1)

    # ------------------------------------------------------------------------
    # SECTION 3: Chunked Matrix Construction
    # ------------------------------------------------------------------------
    print("\n[STEP 3/4] Constructing feature matrix (chunked processing)...")
    
    # Save labels and genome IDs (these are small, load into memory entirely)
    labels_df = pd.DataFrame(valid_labels, columns=['label'])
    labels_file = MATRIX_OUTPUT_DIR / f"y_{TARGET_ANTIBIOTIC}.csv"
    labels_df.to_csv(labels_file, index=False, encoding='utf-8')
    print(f"  ✓ Labels saved: {labels_file}")
    
    genomes_df = pd.DataFrame(valid_genomes, columns=['Genome ID'])
    genomes_file = MATRIX_OUTPUT_DIR / f"genomes_{TARGET_ANTIBIOTIC}.csv"
    genomes_df.to_csv(genomes_file, index=False, encoding='utf-8')
    print(f"  ✓ Genome IDs saved: {genomes_file}")
    
    # Calculate number of chunks needed
    total_genomes = len(valid_genomes)
    num_chunks = (total_genomes + CHUNK_SIZE - 1) // CHUNK_SIZE
    
    print(f"  Processing {total_genomes} genomes in {num_chunks} chunks...")
    print(f"  Chunk size: {CHUNK_SIZE} genomes")
    print("=" * 80)
    
    # Process each chunk
    for chunk_id in range(num_chunks):
        start_idx = chunk_id * CHUNK_SIZE
        end_idx = min((chunk_id + 1) * CHUNK_SIZE, total_genomes)
        
        chunk_genomes = valid_genomes[start_idx:end_idx]
        chunk_output_file = MATRIX_OUTPUT_DIR / f"X_{TARGET_ANTIBIOTIC}_part_{chunk_id}.npz"
        
        # Resume capability: skip if chunk already processed
        if chunk_output_file.exists():
            print(f"  [Chunk {chunk_id+1}/{num_chunks}] Already processed, skipping...")
            continue
        
        print(f"  [Chunk {chunk_id+1}/{num_chunks}] Processing {len(chunk_genomes)} genomes...")
        
        # Build each genome's column indices in parallel — kmc_tools dump per
        # genome is the bottleneck; a thread pool (workers = preprocessing.threads)
        # saturates the core allocation. ThreadPoolExecutor.map preserves genome
        # order, so the CSR indptr stays correct.
        from concurrent.futures import ThreadPoolExecutor

        def _dump_cols(genome_id):
            kmc_db = KMC_OUTPUTS_DIR / genome_id
            temp_dump_file = TEMP_DIR / f"{genome_id}_dump.txt"
            cols = []
            try:
                run_command(f"{KMC_TOOLS_BIN} transform {kmc_db} dump {temp_dump_file}")
                with open(temp_dump_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        ci = kmer_to_index.get(line.split()[0])
                        if ci is not None:
                            cols.append(ci)
            except Exception as e:
                print(f"    WARNING: Failed to process {genome_id}: {e}")
            finally:
                if temp_dump_file.exists():
                    temp_dump_file.unlink()
            return cols

        col_indices = array.array('i')
        indptr = array.array('q', [0])
        with ThreadPoolExecutor(max_workers=max(1, int(THREADS))) as _ex:
            for cols in _ex.map(_dump_cols, chunk_genomes):
                col_indices.extend(cols)
                indptr.append(len(col_indices))
        
        # Construct sparse matrix for this chunk
        # shape = (number of genomes in chunk, total number of features)
        # dtype=np.int8 saves memory (only need 0 or 1)
        chunk_matrix = csr_matrix(
            (np.ones(len(col_indices), dtype=np.int8),
             np.frombuffer(col_indices, dtype=np.int32),
             np.frombuffer(indptr, dtype=np.int64)),
            shape=(len(chunk_genomes), num_features)
        )

        # Safety: guarantee strict binary (presence/absence) encoding. Each
        # k-mer is dumped once per genome by kmc_tools, so duplicates should not
        # occur — but if any (genome, k-mer) pair appeared twice, csr_matrix
        # SUMS the values, yielding a 2 and silently breaking the binary
        # assumption. Clamp any accumulated value back to 1.
        if chunk_matrix.nnz > 0 and chunk_matrix.data.max() > 1:
            np.clip(chunk_matrix.data, 0, 1, out=chunk_matrix.data)

        # Save chunk to disk in compressed sparse format
        save_npz(chunk_output_file, chunk_matrix)
        print(f"    ✓ Saved: {chunk_output_file}")
        print(f"      Matrix shape: {chunk_matrix.shape}")
        print(f"      Sparsity: {(1 - chunk_matrix.nnz / (chunk_matrix.shape[0] * chunk_matrix.shape[1])) * 100:.2f}%")
        
        # Explicit memory cleanup
        del chunk_matrix, col_indices, indptr
        gc.collect()

    # ------------------------------------------------------------------------
    # SECTION 4: Completion Report
    # ------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("FEATURE MATRIX CONSTRUCTION COMPLETE")
    print("=" * 80)
    print(f"Output directory: {MATRIX_OUTPUT_DIR}")
    print(f"Total chunks created: {num_chunks}")
    print(f"Features per genome: {num_features:,} k-mers")
    print(f"Total samples: {total_genomes} genomes")
    print("\nGenerated files:")
    print("  - Feature dictionary: features.txt")
    print(f"  - Labels: y_{TARGET_ANTIBIOTIC}.csv")
    print(f"  - Genome IDs: genomes_{TARGET_ANTIBIOTIC}.csv")
    print(f"  - Matrix chunks: X_{TARGET_ANTIBIOTIC}_part_*.npz (×{num_chunks})")
    print("\n✓ Ready for hyperparameter optimization and model training")
    print("  Next: python scripts/04_optimization.py")
    print("=" * 80)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    # Global guard: any unhandled exception (disk full, OOM, etc.) is reported
    # with a full traceback and forces a non-zero exit, so a crash cannot leave
    # partial outputs that a later run silently treats as complete.
    import traceback
    try:
        create_feature_matrix()
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
