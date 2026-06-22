#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Step 02c — PopPUNK lineage clustering (ROADMAP §0.1 M2).

Clusters every assembly of an organism into a **lineage** (PopPUNK population
cluster) ONCE, antibiotic-independent, so the lineage-aware cross-validation
(lib/lineage.py + 06/07b) can keep a whole lineage on one side of the train/test
split — removing the ~20-30 % AUC inflation from lineage leakage (Yu 2024).

Workflow (validated on a 100-genome smoke, PopPUNK 2.7.8):
    poppunk --create-db --r-files refs.txt --output <db>      # sketch + distances
    poppunk --fit-model <model> --ref-db <db> --output <fit>  # cluster
    -> <fit>/<fit>_clusters.csv  (columns: Taxon, Cluster)

PopPUNK rewrites '.'→'_' in sample names (562.100036 -> 562_100036), so its raw
Taxon column will NOT match the pipeline's genome ids. This script reverses that
and writes the CANONICAL table the rest of the pipeline reads:
    data/processed/{organism}/lineage/poppunk_clusters.csv  (Genome ID, Cluster)

Organism-level: clusters ALL assemblies in the organism's genomes dir (not the
per-antibiotic subset) so one clustering is reused by every antibiotic. KMC/unitig
feature steps are unaffected — this only produces the CV split labels.
"""

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

from utils import run_command

PROJECT_ROOT = Path(__file__).resolve().parent.parent

from lib.config import load_config, resolve_path, resolve_tool  # noqa: E402
from lib.lineage import lineage_summary  # noqa: E402


def write_refs(genome_ids, genomes_dir: Path, refs_path: Path):
    """PopPUNK --r-files: one 'name<TAB>/abs/path.fna' per line."""
    with open(refs_path, "w", encoding="utf-8") as f:
        for gid in genome_ids:
            f.write(f"{gid}\t{(genomes_dir / f'{gid}.fna').resolve()}\n")


def normalize_clusters(clusters_csv, genome_ids, *,
                       taxon_col="Taxon", cluster_col="Cluster") -> pd.DataFrame:
    """Map PopPUNK's mangled Taxon back to the real genome id.

    PopPUNK sanitises sample names by replacing '.' with '_'. Since we control the
    input ids, we build a reverse map (mangled -> original) and translate. Errors
    loudly if any Taxon fails to map (so a silent mismatch can't corrupt the CV).
    Returns a DataFrame with columns ['Genome ID', 'Cluster'] (Cluster as str).
    """
    df = pd.read_csv(clusters_csv, encoding="utf-8")
    for c in (taxon_col, cluster_col):
        if c not in df.columns:
            raise KeyError(f"Column '{c}' not in {clusters_csv} (have {df.columns.tolist()}).")

    rev = {}
    for gid in genome_ids:
        rev.setdefault(str(gid).replace(".", "_"), str(gid))
        rev.setdefault(str(gid), str(gid))   # in case a name was not mangled

    rows, unmatched = [], []
    for taxon, cluster in zip(df[taxon_col].astype(str), df[cluster_col].astype(str)):
        gid = rev.get(taxon)
        if gid is None:
            unmatched.append(taxon)
            continue
        rows.append({"Genome ID": gid, "Cluster": cluster})
    if unmatched:
        raise ValueError(
            f"{len(unmatched)} PopPUNK taxa did not map back to a genome id "
            f"(e.g. {unmatched[:3]}). Name-mangling beyond '.'→'_'?"
        )
    return pd.DataFrame(rows)


def run_poppunk(poppunk, refs_path: Path, work_dir: Path, model: str, threads: int):
    """create-db -> fit-model; returns the path to PopPUNK's raw clusters CSV."""
    db_dir = work_dir / "db"
    fit_dir = work_dir / "fit"
    # PopPUNK refuses to overwrite silently; clear stale dirs for a clean re-run.
    for d in (db_dir, fit_dir):
        shutil.rmtree(d, ignore_errors=True)

    run_command(f"{poppunk} --create-db --r-files {refs_path} "
                f"--output {db_dir} --threads {int(threads)}")
    run_command(f"{poppunk} --fit-model {model} --ref-db {db_dir} "
                f"--output {fit_dir} --threads {int(threads)}")

    clusters = fit_dir / f"{fit_dir.name}_clusters.csv"   # = fit/fit_clusters.csv
    if not clusters.exists():
        found = list(fit_dir.glob("*_clusters.csv"))
        if not found:
            sys.exit(f"ERROR: PopPUNK produced no clusters CSV in {fit_dir}")
        # Avoid the *_unword_clusters.csv variant (different schema).
        clusters = next((p for p in found if not p.name.endswith("_unword_clusters.csv")),
                        found[0])
    return clusters


def main():
    config = load_config()
    lin_cfg = config.get("lineage", {}) or {}
    default_org = config.get("project", {}).get("organism", "ecoli")
    default_threads = lin_cfg.get("threads") or config["preprocessing"].get("threads", 8)

    ap = argparse.ArgumentParser(description="PopPUNK lineage clustering (organism-level).")
    ap.add_argument("--organism", default=default_org)
    ap.add_argument("--threads", type=int, default=default_threads)
    ap.add_argument("--model", default=lin_cfg.get("model", "bgmm"),
                    help="PopPUNK --fit-model (bgmm | dbscan). Default from config lineage.model.")
    ap.add_argument("--clusters-csv", default=None,
                    help="Skip PopPUNK and normalize an existing raw clusters CSV "
                         "(debug/resume; e.g. a smoke's pp_fit_clusters.csv).")
    args = ap.parse_args()

    organism = args.organism
    print("=" * 80)
    print(f"LINEAGE CLUSTERING (PopPUNK) — {organism}")
    print("=" * 80)

    genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    # lineage_dir from config when present; else derive from data_dir so this runs
    # on configs that predate the 'lineage_dir' key (e.g. a manually-tuned HPC
    # config.yaml) without needing a config edit — same canonical location.
    try:
        lineage_dir = resolve_path("lineage_dir", organism=organism, config=config)
    except KeyError:
        lineage_dir = resolve_path("data_dir", config=config) / "processed" / organism / "lineage"
    lineage_dir.mkdir(parents=True, exist_ok=True)

    genome_ids = sorted(p.stem for p in genomes_dir.glob("*.fna"))
    if not genome_ids:
        sys.exit(f"ERROR: no .fna assemblies in {genomes_dir}")
    print(f"  Genomes to cluster: {len(genome_ids)} (organism-level, all antibiotics)")

    if args.clusters_csv:
        raw_clusters = Path(args.clusters_csv)
        if not raw_clusters.exists():
            sys.exit(f"ERROR: --clusters-csv not found: {raw_clusters}")
        print(f"  Using existing raw clusters: {raw_clusters}")
    else:
        poppunk = resolve_tool("poppunk_bin", "poppunk", config=config,
                               env_var="AMR_POPPUNK_BIN")
        if not poppunk:
            sys.exit("ERROR: poppunk not found. Install it (conda install -c bioconda "
                     "poppunk) so it is on PATH, or set AMR_POPPUNK_BIN.")
        work_dir = lineage_dir / "_poppunk_work"
        work_dir.mkdir(parents=True, exist_ok=True)
        refs_path = work_dir / "refs.txt"
        write_refs(genome_ids, genomes_dir, refs_path)
        print(f"  Running PopPUNK (model={args.model}, threads={args.threads})...")
        raw_clusters = run_poppunk(poppunk, refs_path, work_dir, args.model, args.threads)

    print(f"  Normalizing clusters (un-mangling '.'→'_') from: {raw_clusters.name}")
    out_df = normalize_clusters(raw_clusters, genome_ids)

    out_path = lineage_dir / "poppunk_clusters.csv"
    out_df.to_csv(out_path, index=False, encoding="utf-8")

    groups = out_df["Cluster"].to_numpy()
    summ = lineage_summary(groups)
    print("\n" + "=" * 80)
    print("LINEAGE CLUSTERING COMPLETE")
    print("=" * 80)
    print(f"  Output: {out_path}")
    print(f"  Genomes clustered: {summ['n_genomes']} | lineages: {summ['n_clusters']} | "
          f"singletons: {summ['n_singletons']}")
    print(f"  Largest lineage: {summ['largest_cluster']} "
          f"({summ['largest_cluster_size']} genomes, {summ['largest_cluster_frac']*100:.1f}%)")
    if summ["n_clusters"] < int(lin_cfg.get("n_splits", 5)):
        print(f"  ⚠ Only {summ['n_clusters']} lineage(s) — fewer than n_splits "
              f"({lin_cfg.get('n_splits', 5)}). Try --model dbscan or PopPUNK refine "
              f"for finer resolution before GroupKFold.")
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
