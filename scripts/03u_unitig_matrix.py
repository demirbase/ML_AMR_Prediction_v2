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

Two modes (unitigs are sequence features, independent of antibiotic):
  • --build-db : ORGANISM-LEVEL — run unitig-caller ONCE over ALL the organism's
    assemblies and write a reusable store (processed/{organism}/unitig_all/).
  • default    : PER-ANTIBIOTIC — if the store exists, just SUBSET it (select the
    antibiotic's genome rows + re-filter min_support; NO unitig-caller re-run);
    otherwise fall back to calling unitig-caller on the antibiotic subset.
This makes the 2nd..Nth antibiotic/organism near-instant instead of a fresh
multi-hour unitig-caller run each time.
"""

import argparse
import gc
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.sparse import csc_matrix, load_npz, save_npz, vstack

from utils import run_command

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = PROJECT_ROOT / "config" / "config.yaml"

from lib.config import resolve_path, resolve_tool, get_target  # noqa: E402


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

    # Lineage-cluster intersection. PopPUNK's --qc-db (02c) drops distance/length
    # outliers from the sketch DB, so those genomes carry NO lineage label. They are
    # low quality by PopPUNK's own QC and must not train the model: if they stay,
    # 07b finds genomes with no cluster and silently falls back from lineage-CV to a
    # 5-seed holdout (lineage-INFLATED AUC). Dropping them here keeps the matrix,
    # pyseer (14) and the CV (07b) on ONE consistent population — the QC-passed,
    # clustered, labelled genomes. Skipped when no cluster file exists yet (e.g. the
    # synthetic integration test, or a first run before 02c) so prior behaviour holds
    # and the run still completes (07b then warns and falls back, as before).
    clustered_ids = None
    try:
        lineage_dir = resolve_path("lineage_dir", organism=organism, config=config)
    except KeyError:
        lineage_dir = resolve_path("data_dir", config=config) / "processed" / organism / "lineage"
    cluster_file = lineage_dir / "poppunk_clusters.csv"
    if cluster_file.exists():
        cdf = pd.read_csv(cluster_file, encoding="utf-8")
        id_col = "Genome ID" if "Genome ID" in cdf.columns else cdf.columns[0]
        clustered_ids = set(cdf[id_col].astype(str))
        print(f"  ✓ Loaded {len(clustered_ids)} clustered genomes (lineage intersection).")
    else:
        print(f"  ⚠ No lineage clusters at {cluster_file}; skipping intersection "
              f"(07b will fall back to holdout CV — run 02c first for lineage-CV).")

    valid_genomes, valid_labels = [], []
    missing_fna = skipped_outliers = skipped_unclustered = 0
    for gid, label in zip(meta["Genome ID"].values, meta[antibiotic].astype(int).values):
        if gid in outlier_ids:
            skipped_outliers += 1
            continue
        if clustered_ids is not None and gid not in clustered_ids:
            skipped_unclustered += 1
            continue
        if not (raw_genomes_dir / f"{gid}.fna").exists():
            missing_fna += 1
            continue
        valid_genomes.append(gid)
        valid_labels.append(int(label))

    if skipped_outliers:
        print(f"  ✓ Skipped {skipped_outliers} QC-outlier genomes.")
    if skipped_unclustered:
        print(f"  ✓ Skipped {skipped_unclustered} genomes absent from the lineage "
              f"clusters (PopPUNK --qc-db failures; excluded from the model).")
    if missing_fna:
        print(f"  ⚠ Skipped {missing_fna} genomes: .fna missing in {raw_genomes_dir}.")
    if not valid_genomes:
        raise SystemExit(f"ERROR: no genomes passed validation (looked in {raw_genomes_dir}).")

    pos = sum(valid_labels)
    print(f"  ✓ Valid genomes: {len(valid_genomes)} "
          f"({pos} resistant / {len(valid_labels) - pos} susceptible, "
          f"{pos / len(valid_labels) * 100:.1f}% R)")

    # Deterministic (seed-42) shuffle to break any phenotype-BLOCKED ordering in
    # the metadata (e.g. clonal MRSA / A. baumannii, where amr_phenotypes.csv can
    # list all-R then all-S). Without it, the downstream chunk split (03/04) can
    # produce a single-class chunk -> single-class CV fold -> XGBoost NaN. The
    # fixed seed keeps the genome/chunk assignment fully reproducible.
    perm = np.random.RandomState(42).permutation(len(valid_genomes))
    valid_genomes = [valid_genomes[i] for i in perm]
    valid_labels = [valid_labels[i] for i in perm]
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
                # Presence values are single-char 0/1, tab-separated, one per
                # sample. Parse via frombuffer (fast, C-level) instead of the
                # deprecated np.fromstring(sep=...): strip tabs/newline -> a run of
                # '0'/'1' chars -> ASCII bytes -> subtract ord('0') to get 0/1.
                vals = (np.frombuffer(line[tab + 1:].replace("\t", "").strip().encode("ascii"),
                                      dtype=np.int8) - ord("0"))
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


def store_dir_for(organism, config):
    """Organism-level unitig store (one unitig-caller run, reused by every
    antibiotic). Derived from data_dir so no extra config key is required."""
    return resolve_path("data_dir", config=config) / "processed" / organism / "unitig_all"


def subset_store_to_antibiotic(store_dir, out_dir, antibiotic, valid_genomes,
                               valid_labels, chunk_size, min_support):
    """Build matrix_unitig for one antibiotic by SUBSETTING the organism-level
    store (no unitig-caller re-run). Selects the antibiotic's genome rows from the
    store, re-applies the absolute min_support + zero-variance-core filter over the
    SUBSET (so support is genome-count within this antibiotic), re-indexes unitigs,
    and writes the same chunk contract as rtab_to_chunks.
    """
    store_genomes = pd.read_csv(store_dir / "genomes_all.csv",
                                encoding="utf-8")["Genome ID"].astype(str).tolist()
    g2row = {g: i for i, g in enumerate(store_genomes)}
    missing = [g for g in valid_genomes if g not in g2row]
    if missing:
        sys.exit(f"ERROR: {len(missing)} antibiotic genome(s) not in the unitig store "
                 f"(e.g. {missing[:3]}). Rebuild the store with --build-db.")
    sel_rows = [g2row[g] for g in valid_genomes]   # store row indices, in valid_genomes order

    chunks = sorted(store_dir.glob("X_all_part_*.npz"),
                    key=lambda p: int(p.stem.split("_")[-1]))
    if not chunks:
        sys.exit(f"ERROR: no store chunks (X_all_part_*.npz) in {store_dir}")
    print(f"  Loading organism store ({len(store_genomes)} genomes, {len(chunks)} chunks)...")
    X_all = vstack([load_npz(f) for f in chunks]).tocsr()
    X_sub = X_all[sel_rows]                          # (n_sel × n_unitigs_store), valid_genomes order
    del X_all
    gc.collect()

    n_sel = len(valid_genomes)
    max_support = n_sel - 1
    support = np.asarray(X_sub.sum(axis=0)).ravel()
    keep = np.where((support >= min_support) & (support <= max_support))[0]
    if keep.size == 0:
        sys.exit(f"ERROR: no unitigs pass support in the {antibiotic} subset "
                 f"(min_support={min_support}, n={n_sel}). Lower --min-support?")
    X_kept = X_sub[:, keep].tocsr()
    if X_kept.nnz and X_kept.data.max() > 1:
        np.clip(X_kept.data, 0, 1, out=X_kept.data)
    del X_sub
    gc.collect()
    print(f"  ✓ Kept {keep.size:,} / {len(support):,} store unitigs for {antibiotic} "
          f"(min_support={min_support}, max_support={max_support}).")

    # features.txt: keep the store's unitig sequences at the surviving column indices.
    all_feats = [ln.split("\t")[0] for ln in
                 (store_dir / "features.txt").read_text(encoding="utf-8").splitlines()]
    kept_feats = [all_feats[ci] for ci in keep]
    with open(out_dir / "features.txt", "w", encoding="utf-8") as f:
        for seq in kept_feats:
            f.write(f"{seq}\t1\n")

    # unitigs.rtab — pyseer (14) consumes the unitig×genome presence table in
    # unitig-caller's Rtab format (header 'Unitig_sequence<TAB>samples', then one
    # 0/1 row per unitig). The fallback path gets it straight from unitig-caller, but
    # the store-subset path built ONLY the npz chunks, so pyseer failed with 'unitig
    # Rtab not found'. Emit it here from the subset so both 03u paths leave the same
    # artefacts. (Big text file — the per-antibiotic chain rm's it after populate.)
    Xt = X_kept.T.tocsr()                       # (n_unitigs × n_genomes), 0/1
    n_g = Xt.shape[1]
    with open(out_dir / "unitigs.rtab", "w", encoding="utf-8") as rf:
        rf.write("Unitig_sequence\t" + "\t".join(valid_genomes) + "\n")
        indptr, indices = Xt.indptr, Xt.indices
        for j in range(Xt.shape[0]):
            row = np.zeros(n_g, dtype=np.int8)
            row[indices[indptr[j]:indptr[j + 1]]] = 1
            rf.write(kept_feats[j])
            rf.write("\t")
            row.tofile(rf, sep="\t")
            rf.write("\n")
    print(f"  ✓ Wrote unitigs.rtab ({Xt.shape[0]:,} unitigs × {n_g} genomes) for pyseer (14).")

    pd.DataFrame(valid_labels, columns=["label"]).to_csv(
        out_dir / f"y_{antibiotic}.csv", index=False, encoding="utf-8")
    pd.DataFrame(valid_genomes, columns=["Genome ID"]).to_csv(
        out_dir / f"genomes_{antibiotic}.csv", index=False, encoding="utf-8")

    num_chunks = (n_sel + chunk_size - 1) // chunk_size
    print(f"  Writing {num_chunks} chunk(s)...")
    for c in range(num_chunks):
        start, end = c * chunk_size, min((c + 1) * chunk_size, n_sel)
        chunk = X_kept[start:end]
        out_npz = out_dir / f"X_{antibiotic}_part_{c}.npz"
        save_npz(out_npz, chunk)
        sparsity = (1 - chunk.nnz / (chunk.shape[0] * chunk.shape[1])) * 100 \
            if chunk.shape[1] else 0.0
        print(f"    ✓ {out_npz.name}  shape={chunk.shape}  sparsity={sparsity:.2f}%")
        del chunk
        gc.collect()
    return int(keep.size), num_chunks


def main():
    config = _load_config()
    unitig_cfg = config.get("unitig", {}) or {}
    default_org = get_target(config=config)[0]
    default_ab = get_target(config=config)[1]
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
    ap.add_argument("--build-db", action="store_true",
                    help="ORGANISM-LEVEL mode: run unitig-caller ONCE over ALL of the "
                         "organism's assemblies and write a reusable store "
                         "(processed/{organism}/unitig_all/). Then per-antibiotic runs "
                         "just subset this store — no unitig-caller re-run.")
    ap.add_argument("--db-min-support", type=int,
                    default=int(unitig_cfg.get("db_min_support", 2)),
                    help="Absolute support floor for the organism store (--build-db). "
                         "Keep <= every antibiotic's --min-support so nothing a "
                         "per-antibiotic filter would keep is pre-dropped (default 2: "
                         "drops only singletons, always safe).")
    args = ap.parse_args()

    organism, antibiotic = args.organism, args.antibiotic
    raw_genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    store_dir = store_dir_for(organism, config)

    # ----- ORGANISM-LEVEL store build (once, reused by every antibiotic) ---------
    if args.build_db:
        print("=" * 80)
        print(f"UNITIG STORE BUILD (organism-level) — {organism}")
        print("=" * 80)
        store_dir.mkdir(parents=True, exist_ok=True)
        all_genomes = sorted(p.stem for p in raw_genomes_dir.glob("*.fna"))
        if not all_genomes:
            sys.exit(f"ERROR: no .fna assemblies in {raw_genomes_dir}")
        print(f"  Genomes: {len(all_genomes)} | db_min_support: {args.db_min_support}")
        if args.rtab:
            rtab = Path(args.rtab)
            if not rtab.exists():
                sys.exit(f"ERROR: --rtab path not found: {rtab}")
            print(f"  ✓ Using provided rtab: {rtab}")
        else:
            rtab = run_unitig_caller(all_genomes, raw_genomes_dir, store_dir,
                                     args.threads, config)
        # Reuse rtab_to_chunks with name 'all' + dummy labels: writes features.txt,
        # genomes_all.csv, X_all_part_*.npz (y_all.csv is a harmless dummy).
        n_unitigs, num_chunks = rtab_to_chunks(
            rtab, all_genomes, [0] * len(all_genomes), store_dir, "all",
            args.chunk_size, args.db_min_support)
        print("\n" + "=" * 80)
        print("UNITIG STORE BUILD COMPLETE")
        print(f"  Store: {store_dir} | unitigs: {n_unitigs:,} | "
              f"genomes: {len(all_genomes)} | chunks: {num_chunks}")
        print("  Per-antibiotic: run 03u normally — it will SUBSET this store.")
        print("=" * 80)
        return

    # ----- PER-ANTIBIOTIC matrix --------------------------------------------------
    matrix_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic,
                              config=config)
    out_dir = matrix_dir.parent / args.out_subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    print("=" * 80)
    print("UNITIG FEATURE MATRIX CONSTRUCTION")
    print("=" * 80)
    print(f"Organism: {organism} | Antibiotic: {antibiotic}")
    print(f"min_support: {args.min_support} (absolute) | chunk_size: {args.chunk_size}")
    print(f"Output dir: {out_dir}")
    print("=" * 80)

    print("\n[1/3] Selecting genomes...")
    valid_genomes, valid_labels = select_genomes(config, organism, antibiotic)

    store_ready = (store_dir / "genomes_all.csv").exists() and not args.rtab
    if store_ready:
        # Fast path: subset the organism store (no unitig-caller re-run).
        print(f"\n[2/3] Subsetting organism store: {store_dir}")
        print("[3/3] Building genome×unitig matrix from store...")
        n_unitigs, num_chunks = subset_store_to_antibiotic(
            store_dir, out_dir, antibiotic, valid_genomes, valid_labels,
            args.chunk_size, args.min_support)
    else:
        # Fallback: call unitig-caller on just this antibiotic's genomes.
        print("\n[2/3] Calling unitigs (no organism store; per-antibiotic run)...")
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
            args.chunk_size, args.min_support)

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
