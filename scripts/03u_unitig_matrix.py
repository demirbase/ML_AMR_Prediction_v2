#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unitig Feature Matrix Construction (ROADMAP §0.1 M12 — replaces raw k-mers).

This is the unitig-based alternative to 03_matrix_construction.py. Instead of a
genome×k-mer presence/absence matrix it builds a genome×UNITIG presence/absence
matrix, where unitigs are the maximal non-branching paths of the population
compacted de Bruijn graph (unitig-caller, Bifrost backend). Unitigs are longer,
fewer, BLAST-mappable and GWAS-standard, which dissolves the raw-k-mer
speed/memory/min_support pressure while keeping the *downstream XGBoost unchanged*.

Drop-in output contract (identical to 03 so 03b/04/05/06/07/07b read it as-is):
    <out_subdir>/
        features.txt                 one line per unitig: "<unitig_seq>\\t1"
                                     (line index == matrix column index)
        y_{antibiotic}.csv           column 'label' (genome order == rows)
        genomes_{antibiotic}.csv     column 'Genome ID' (same order)
        X_{antibiotic}_part_{c}.npz  CSR int8 binary chunks of chunk_size genomes

The ONLY semantic difference vs 03: features.txt rows are variable-length unitig
sequences, not fixed 21-mers. Downstream steps that hard-assume k=21 (08 BLAST
task, 09 coverage = aln_len/k, 11 SNP codon mapping) are handled in later
ROADMAP §0 steps, NOT here — this script only produces the matrix.

Output goes to a SEPARATE sibling dir (default 'matrix_unitig') so the working
raw-k-mer 'matrix' (the baseline) is never overwritten.

KMC (02/02b) stays for QC/spectra only; this step needs the raw .fna assemblies.
"""

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.sparse import csc_matrix, save_npz

from utils import run_command

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

from lib.config import resolve_path, resolve_tool  # noqa: E402


def _load_config():
    if not CONFIG_PATH.exists():
        raise FileNotFoundError(f"Configuration file not found: {CONFIG_PATH}")
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def select_genomes(config, organism, antibiotic):
    """Genomes with a label for `antibiotic` AND a present .fna, minus QC outliers.

    Returns (valid_genomes, valid_labels) in metadata order. Unlike 03 this does
    NOT require a KMC database — unitig-caller consumes the .fna assemblies.
    """
    metadata_file = resolve_path("metadata_file", organism=organism, config=config)
    raw_genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)

    if not metadata_file.exists():
        raise FileNotFoundError(f"Metadata file not found: {metadata_file}")

    meta = pd.read_csv(metadata_file, encoding="utf-8")
    meta["Genome ID"] = meta["Genome ID"].astype(str)
    if antibiotic not in meta.columns:
        raise KeyError(
            f"Column '{antibiotic}' not found in metadata. "
            f"Available: {meta.columns.tolist()}"
        )
    meta = meta.dropna(subset=[antibiotic]).copy()

    # QC outlier blacklist (same source as 03)
    outlier_ids = set()
    outlier_file = (
        resolve_path("dir_global_exploration", organism=organism, config=config)
        / "global_qc_outliers.csv"
    )
    if outlier_file.exists():
        odf = pd.read_csv(outlier_file)
        col = "Genome" if "Genome" in odf.columns else odf.columns[0]
        outlier_ids = set(odf[col].astype(str))
        print(f"  ✓ Loaded {len(outlier_ids)} QC outlier genomes to exclude.")
    else:
        print(f"  ⚠ Outlier file not found at {outlier_file} (none excluded).")

    valid_genomes, valid_labels = [], []
    missing_fna = skipped_outliers = 0
    for gid, label in zip(meta["Genome ID"].values, meta[antibiotic].astype(int).values):
        if gid in outlier_ids:
            skipped_outliers += 1
            continue
        if not (raw_genomes_dir / f"{gid}.fna").exists():
            missing_fna += 1
            continue
        valid_genomes.append(gid)
        valid_labels.append(int(label))

    if skipped_outliers:
        print(f"  ✓ Skipped {skipped_outliers} QC-outlier genomes.")
    if missing_fna:
        print(f"  ⚠ Skipped {missing_fna} genomes: .fna missing in {raw_genomes_dir}.")
    if not valid_genomes:
        raise SystemExit(f"ERROR: no genomes passed validation (looked in {raw_genomes_dir}).")

    pos = sum(valid_labels)
    print(f"  ✓ Valid genomes: {len(valid_genomes)} "
          f"({pos} resistant / {len(valid_labels) - pos} susceptible, "
          f"{pos / len(valid_labels) * 100:.1f}% R)")
    return valid_genomes, valid_labels


def run_unitig_caller(valid_genomes, raw_genomes_dir, out_dir, threads, config):
    """Call unitigs across all valid genomes -> <out_dir>/unitigs.rtab.

    Refs file = one absolute .fna path per line (unitig-caller v1.3.x format;
    sample names are derived from the file basenames == genome_id). Returns the
    rtab path. Skips the call if the rtab already exists (resume-safe).
    """
    rtab = out_dir / "unitigs.rtab"
    if rtab.exists():
        print(f"  ✓ Unitig rtab already exists, reusing: {rtab}")
        return rtab

    unitig_caller = resolve_tool(
        "unitig_caller_bin", "unitig-caller", config=config,
        env_var="AMR_UNITIG_CALLER_BIN",
    )
    if not unitig_caller:
        sys.exit(
            "ERROR: unitig-caller not found. Install it (conda install -c bioconda "
            "unitig-caller) so it is on PATH, or set AMR_UNITIG_CALLER_BIN."
        )

    refs_file = out_dir / "unitig_refs.txt"
    with open(refs_file, "w", encoding="utf-8") as f:
        for gid in valid_genomes:
            f.write(str((raw_genomes_dir / f"{gid}.fna").resolve()) + "\n")

    out_prefix = out_dir / "unitigs"
    print(f"  Running unitig-caller on {len(valid_genomes)} genomes "
          f"(threads={threads})...")
    run_command(
        f"{unitig_caller} --call --refs {refs_file} --rtab "
        f"--out {out_prefix} --threads {int(threads)}"
    )
    if not rtab.exists():
        sys.exit(f"ERROR: unitig-caller did not produce expected rtab: {rtab}")
    print(f"  ✓ Unitig rtab written: {rtab}")
    return rtab


def rtab_to_chunks(rtab, valid_genomes, valid_labels, out_dir, antibiotic,
                   chunk_size, min_support):
    """Transpose unitig×genome rtab -> genome×unitig CSR chunks (03 contract).

    Filtering (absolute, NOT proportional — ROADMAP §0.7 risk 4): keep a unitig
    only if min_support <= (#genomes carrying it) <= n_genomes-1. The upper bound
    drops zero-variance core unitigs (present in every genome), mirroring 03's
    -cx max_support. Singletons / ultra-rare unitigs below min_support are dropped.
    """
    genome_to_row = {gid: i for i, gid in enumerate(valid_genomes)}
    n_genomes = len(valid_genomes)
    max_support = n_genomes - 1

    with open(rtab, "r", encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sample_ids = header[1:]  # first field is the 'Unitig_sequence' label

        # Map each rtab sample COLUMN -> our output ROW index. Do NOT assume the
        # rtab column order matches valid_genomes; derive it from the header.
        missing = [s for s in sample_ids if s not in genome_to_row]
        if missing:
            sys.exit(
                f"ERROR: rtab has {len(missing)} sample(s) not in the genome set "
                f"(e.g. {missing[:3]}). Refs/metadata mismatch."
            )
        if len(sample_ids) != n_genomes:
            print(f"  ⚠ rtab has {len(sample_ids)} samples but {n_genomes} were "
                  f"requested; building rows only for the {len(sample_ids)} present.")
        row_of_rtabcol = np.array([genome_to_row[s] for s in sample_ids], dtype=np.int32)

        features_file = out_dir / "features.txt"
        # Global CSC accumulators (column = unitig, row = genome).
        indices = []          # row indices, concatenated column by column
        indptr = [0]
        n_unitigs_kept = 0
        n_seen = 0

        with open(features_file, "w", encoding="utf-8") as feat:
            for line in fh:
                n_seen += 1
                tab = line.find("\t")
                if tab < 0:
                    continue
                seq = line[:tab]
                # Presence values are single-char 0/1, tab-separated, one per sample.
                vals = np.fromstring(line[tab + 1:], dtype=np.int8, sep="\t")
                if vals.size != len(sample_ids):
                    sys.exit(
                        f"ERROR: unitig '{seq[:20]}...' has {vals.size} values but "
                        f"{len(sample_ids)} samples expected (malformed rtab)."
                    )
                support = int(vals.sum())
                if support < min_support or support > max_support:
                    continue  # rare/singleton or zero-variance core -> drop
                rows = row_of_rtabcol[np.nonzero(vals)[0]]
                indices.append(rows)
                indptr.append(indptr[-1] + rows.size)
                feat.write(f"{seq}\t1\n")
                n_unitigs_kept += 1

    if n_unitigs_kept == 0:
        sys.exit(
            f"ERROR: no unitigs passed the support filter "
            f"(min_support={min_support}, max_support={max_support}, "
            f"{n_seen} candidates). Lower --min-support?"
        )

    print(f"  ✓ Kept {n_unitigs_kept:,} / {n_seen:,} unitigs "
          f"(min_support={min_support}, max_support={max_support}).")

    # Assemble the global sparse matrix once (column-major from the stream), then
    # slice row-blocks into the 03-style chunk files.
    indices_arr = (np.concatenate(indices) if indices
                   else np.empty(0, dtype=np.int32))
    data_arr = np.ones(indices_arr.size, dtype=np.int8)
    indptr_arr = np.asarray(indptr, dtype=np.int64)
    full = csc_matrix(
        (data_arr, indices_arr, indptr_arr),
        shape=(n_genomes, n_unitigs_kept),
    )
    del indices, indices_arr, data_arr, indptr, indptr_arr
    gc.collect()

    # Labels + genome IDs (row order)
    pd.DataFrame(valid_labels, columns=["label"]).to_csv(
        out_dir / f"y_{antibiotic}.csv", index=False, encoding="utf-8")
    pd.DataFrame(valid_genomes, columns=["Genome ID"]).to_csv(
        out_dir / f"genomes_{antibiotic}.csv", index=False, encoding="utf-8")

    num_chunks = (n_genomes + chunk_size - 1) // chunk_size
    print(f"  Writing {num_chunks} chunk(s) of up to {chunk_size} genomes...")
    for c in range(num_chunks):
        start, end = c * chunk_size, min((c + 1) * chunk_size, n_genomes)
        chunk = full[start:end].tocsr()
        if chunk.nnz and chunk.data.max() > 1:  # safety: enforce strict binary
            np.clip(chunk.data, 0, 1, out=chunk.data)
        out_npz = out_dir / f"X_{antibiotic}_part_{c}.npz"
        save_npz(out_npz, chunk)
        sparsity = (1 - chunk.nnz / (chunk.shape[0] * chunk.shape[1])) * 100 \
            if chunk.shape[1] else 0.0
        print(f"    ✓ {out_npz.name}  shape={chunk.shape}  sparsity={sparsity:.2f}%")
        del chunk
        gc.collect()

    return n_unitigs_kept, num_chunks


def main():
    config = _load_config()
    unitig_cfg = config.get("unitig", {}) or {}
    default_org = config.get("project", {}).get("organism", "ecoli")
    default_ab = config["project"]["target_antibiotic"]
    # unitig.threads overrides preprocessing.threads when set (null -> fall back).
    default_threads = unitig_cfg.get("threads") or config["preprocessing"].get("threads", 8)
    # chunk_size MUST match the value 04/05/06 use (they slice y_{ab}.csv by it via
    # get_y_chunk), so it is sourced from the same preprocessing.chunk_size key.
    default_chunk = config["preprocessing"].get("chunk_size", 200)
    default_out_subdir = unitig_cfg.get("out_subdir", "matrix_unitig")
    default_min_support = int(unitig_cfg.get("min_support", 1))

    ap = argparse.ArgumentParser(description="Build the genome×unitig binary matrix.")
    ap.add_argument("--organism", default=default_org)
    ap.add_argument("--antibiotic", default=default_ab)
    ap.add_argument("--threads", type=int, default=default_threads)
    ap.add_argument("--chunk-size", type=int, default=default_chunk)
    ap.add_argument("--out-subdir", default=default_out_subdir,
                    help="Sibling of the raw-k-mer 'matrix' dir (kept separate so "
                         "the baseline matrix is never overwritten). "
                         "Default from config unitig.out_subdir.")
    ap.add_argument("--min-support", type=int, default=default_min_support,
                    help="Drop unitigs carried by fewer than this many genomes "
                         "(absolute count; ROADMAP §0.7 recommends >=10 for the "
                         "full run, config unitig.min_support). Zero-variance core "
                         "(present in all genomes) is always dropped.")
    ap.add_argument("--rtab", default=None,
                    help="Use an existing unitigs rtab instead of running "
                         "unitig-caller (debug/resume).")
    args = ap.parse_args()

    organism, antibiotic = args.organism, args.antibiotic
    print("=" * 80)
    print("UNITIG FEATURE MATRIX CONSTRUCTION")
    print("=" * 80)
    print(f"Organism: {organism} | Antibiotic: {antibiotic}")
    print(f"min_support: {args.min_support} (absolute) | chunk_size: {args.chunk_size}")
    print("=" * 80)

    raw_genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    matrix_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic,
                              config=config)
    out_dir = matrix_dir.parent / args.out_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output dir: {out_dir}")

    print("\n[1/3] Selecting genomes...")
    valid_genomes, valid_labels = select_genomes(config, organism, antibiotic)

    print("\n[2/3] Calling unitigs...")
    if args.rtab:
        rtab = Path(args.rtab)
        if not rtab.exists():
            sys.exit(f"ERROR: --rtab path not found: {rtab}")
        print(f"  ✓ Using provided rtab: {rtab}")
    else:
        rtab = run_unitig_caller(valid_genomes, raw_genomes_dir, out_dir,
                                 args.threads, config)

    print("\n[3/3] Building genome×unitig matrix...")
    n_unitigs, num_chunks = rtab_to_chunks(
        rtab, valid_genomes, valid_labels, out_dir, antibiotic,
        args.chunk_size, args.min_support,
    )

    print("\n" + "=" * 80)
    print("UNITIG MATRIX CONSTRUCTION COMPLETE")
    print("=" * 80)
    print(f"Output directory: {out_dir}")
    print(f"Unitigs (features): {n_unitigs:,} | Genomes: {len(valid_genomes)} | "
          f"Chunks: {num_chunks}")
    print("Downstream: point 03b/04/05 at this dir (matrix_unitig) to train on unitigs.")
    print("=" * 80)


if __name__ == "__main__":
    import traceback
    try:
        main()
    except SystemExit:
        raise
    except Exception as e:
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
