#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 00 — Build the binary phenotype matrix from the cleaned AMR table.

Reads the cleaned long table produced by 00a_download_bvbrc.py, pivots it into
the wide binary matrix the rest of the pipeline expects, and keeps only the
genomes whose assembly (.fna) is actually present on disk.

    input :  data/external/{org}/metadata/amr_cleaned_long.csv   (genome_id, antibiotic, label)
    output:  data/external/{org}/metadata/amr_phenotypes.csv      (Genome ID + antibiotic 0/1 columns)

This step is RE-RUNNABLE: after retrying failed downloads (00a --retry-failed),
run it again and the matrix is refreshed to include the newly arrived genomes.
The labels (y_*.csv) are still materialised later, in 03_matrix_construction.py.
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.bvbrc import pivot_binary               # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402
from lib.registry import normalize_antibiotic    # noqa: E402

log = logging.getLogger("prepare")


def setup_logging(log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    fh = logging.FileHandler(log_path, encoding="utf-8"); fh.setFormatter(fmt)
    sh = logging.StreamHandler(); sh.setFormatter(fmt)
    log.addHandler(fh); log.addHandler(sh)


def main():
    p = argparse.ArgumentParser(description="Build amr_phenotypes.csv from the cleaned AMR table.")
    p.add_argument("--organism", default=None, help="registry slug (default: config project.organism)")
    args = p.parse_args()

    config = load_config()
    organism = args.organism or config.get("project", {}).get("organism", "ecoli")

    metadata_file = resolve_path("metadata_file", organism=organism, config=config)  # amr_phenotypes.csv
    meta_dir = metadata_file.parent
    genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    logs_dir = resolve_path("logs_dir", organism=organism, antibiotic="_global", config=config).parent
    cleaned_csv = meta_dir / "amr_cleaned_long.csv"

    setup_logging(logs_dir / "00_prepare.log")
    log.info("=" * 70)
    log.info(f"Prepare phenotype matrix — organism={organism}")
    log.info("=" * 70)

    if not cleaned_csv.exists():
        log.error(f"Cleaned table not found: {cleaned_csv}\n"
                  f"  Run 00a_download_bvbrc.py first.")
        sys.exit(1)

    cleaned = pd.read_csv(cleaned_csv, dtype={"genome_id": str})
    # Re-normalise antibiotic names to canonical registry spelling BEFORE pivoting.
    # 00a's clean_amr_table already normalises, but an amr_cleaned_long.csv written
    # before an alias was added keeps the raw name — e.g. 'ampicillin/sulbactam'
    # instead of 'ampicillin_sulbactam' — which then becomes a phenotype column the
    # registry panel (and 03u) can't match, failing that model. Applying the current
    # alias here makes the matrix self-healing against a stale cleaned table.
    cleaned["antibiotic"] = cleaned["antibiotic"].map(normalize_antibiotic)
    wide = pivot_binary(cleaned)
    n_candidates = len(wide)

    # keep only genomes whose assembly is present on disk
    present = {f.stem for f in genomes_dir.glob("*.fna")} if genomes_dir.exists() else set()
    log.info(f"Candidate genomes (cleaned): {n_candidates} | assemblies present: {len(present)}")
    if present:
        wide = wide[wide["Genome ID"].astype(str).isin(present)].copy()
    else:
        log.warning("No .fna assemblies found — writing the full cleaned matrix anyway. "
                    "Run 00a to download genomes, then re-run this step.")

    meta_dir.mkdir(parents=True, exist_ok=True)
    wide.to_csv(metadata_file, index=False, encoding="utf-8")

    n_written = len(wide)
    ab_cols = [c for c in wide.columns if c != "Genome ID"]
    log.info(f"Wrote {metadata_file}")
    log.info(f"  genomes: {n_written} (dropped {n_candidates - n_written} without an assembly)")
    log.info(f"  antibiotics: {len(ab_cols)} -> {ab_cols}")
    # per-antibiotic tested counts (non-NaN)
    if n_written:
        counts = {c: int(wide[c].notna().sum()) for c in ab_cols}
        log.info(f"  tested counts per antibiotic: {counts}")
    log.info("Done. Next: 01_data_validation.py / 02_kmer_extraction.py")


if __name__ == "__main__":
    main()
