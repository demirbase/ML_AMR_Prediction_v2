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

from lib import registry  # noqa: E402
from lib.config import load_config, resolve_path, resolve_tool  # noqa: E402
from lib.lineage import lineage_summary  # noqa: E402


def lineage_params(organism, config):
    """Global ``lineage:`` defaults from config.yaml, overlaid with this
    organism's ``lineage:`` override in the registry (organisms.yaml).

    ESKAPEE population structures differ enough that one setting cannot serve all
    seven (docs/literature/E3.md): S. aureus is highly clonal with very low core
    diversity and needs a bigger sketch to separate strains at all, while
    A. baumannii mainly needs a tighter length filter. Everything not overridden
    stays on the shared default, so the panel remains comparable by construction
    and each deviation is explicit and justified in the registry.
    """
    params = dict(config.get("lineage", {}) or {})
    try:
        override = registry.get_organism(organism).get("lineage") or {}
    except KeyError:
        override = {}   # organism not in the registry: fall back to globals
    params.update(override)
    return params


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


def _sketch_args(p):
    """PopPUNK --create-db sketching flags from the merged lineage params.

    Passed EXPLICITLY rather than relying on PopPUNK's defaults (13/29/4/10000):
    the sketch settings decide the clustering, the clustering is the CV grouping,
    so they are scientific parameters and belong in the recorded config — not in
    whatever the installed PopPUNK happens to default to.
    """
    return (f"--min-k {int(p['min_k'])} --max-k {int(p['max_k'])} "
            f"--k-step {int(p['k_step'])} --sketch-size {int(p['sketch_size'])}")


def _qc_args(p):
    """PopPUNK --qc-db flags from the merged lineage params (E3 §3.4)."""
    args = (f"--max-a-dist {float(p['max_a_dist'])} "
            f"--length-sigma {int(p['length_sigma'])} "
            f"--prop-n {float(p['prop_n'])}")
    lr = p.get("length_range")
    if lr:
        lo, hi = lr
        args += f" --length-range {int(lo)} {int(hi)}"   # takes TWO values
    return args


def run_poppunk(poppunk, refs_path: Path, work_dir: Path, model: str, threads: int,
                *, params: dict, refine: bool = True, reuse_db: bool = False):
    """create-db [-> qc-db] -> fit-model [-> refine]; returns the raw clusters CSV.

    Uses the canonical single-dir PopPUNK pattern (--output == --ref-db, with
    --overwrite so a batch job never blocks on a prompt). ``refine`` chains a
    boundary-refinement step after the initial ``model`` fit. ``reuse_db`` skips
    the expensive re-sketching when a prior run's database is still present.

    MODEL CHOICE — we deviate from the literature, deliberately. E3 §3.2 (and the
    iterative-PopPUNK papers) recommend BGMM K=2 + refine and warn that HDBSCAN
    (PopPUNK's ``dbscan``) can fail to converge on highly recombinant species.
    On OUR data the opposite happened: dbscan gives good strain-level resolution
    while bgmm collapsed into a single ~94 % mega-cluster whose degenerate
    2-Gaussian fit then made refine die on a NaN boundary. Note the literature's
    recommendation assumes iterative-PopPUNK's multi-boundary fit, which is not
    what we ran. So: ``model: dbscan``, ``refine: false``, justified empirically
    rather than by convention — and the Methods should say exactly this.

    Because E3's warning targets high-recombination organisms, dbscan may yet
    fail on an ESKAPEE member we have not clustered (Enterobacter is the prime
    suspect). It fails LOUDLY: main() checks n_clusters >= n_splits and tells you
    to try ``--model bgmm``. Verify the resolution for every new organism — the
    clusters ARE the CV folds, so a bad clustering silently degrades the split
    rather than raising.
    """
    db = work_dir / "db"
    has_sketch = db.exists() and any(db.glob("*.h5"))
    if reuse_db and has_sketch:
        print("  ✓ Reusing existing PopPUNK sketch database (skipping --create-db).")
    else:
        shutil.rmtree(db, ignore_errors=True)
        print(f"  Sketching (k {params['min_k']}-{params['max_k']} step "
              f"{params['k_step']}, sketch {params['sketch_size']})...")
        run_command(f"{poppunk} --create-db --r-files {refs_path} "
                    f"--output {db} --threads {int(threads)} {_sketch_args(params)}")
        if params.get("qc", True):
            print("  Quality control on the sketch database (--qc-db)...")
            run_command(f"{poppunk} --qc-db --ref-db {db} --output {db} "
                        f"--overwrite --threads {int(threads)} {_qc_args(params)}")

    run_command(f"{poppunk} --fit-model {model} --ref-db {db} --output {db} "
                f"--overwrite --threads {int(threads)}")
    if refine:
        print("  Refining model boundary for strain-level resolution (bgmm/dbscan "
              "alone under-cluster)...")
        run_command(f"{poppunk} --fit-model refine --ref-db {db} --output {db} "
                    f"--overwrite --threads {int(threads)}")

    clusters = db / f"{db.name}_clusters.csv"   # = db/db_clusters.csv
    if not clusters.exists():
        found = [p for p in db.glob("*_clusters.csv")
                 if not p.name.endswith("_unword_clusters.csv")]
        if not found:
            sys.exit(f"ERROR: PopPUNK produced no clusters CSV in {db}")
        clusters = found[0]
    return clusters


def main():
    config = load_config()
    lin_cfg = config.get("lineage", {}) or {}
    default_org = config.get("project", {}).get("organism", "ecoli")
    default_threads = lin_cfg.get("threads") or config["preprocessing"].get("threads", 8)

    ap = argparse.ArgumentParser(description="PopPUNK lineage clustering (organism-level).")
    ap.add_argument("--organism", default=default_org)
    ap.add_argument("--threads", type=int, default=default_threads)
    ap.add_argument("--model", default=lin_cfg.get("model", "dbscan"),
                    help="PopPUNK initial --fit-model (dbscan | bgmm). Default from "
                         "config lineage.model. Try bgmm if dbscan under-clusters a "
                         "new organism (n_clusters < n_splits).")
    ap.add_argument("--refine", dest="refine", action="store_true", default=None,
                    help="Chain a refine step after the initial fit (default: config "
                         "lineage.refine, which is OFF — refine failed on a degenerate "
                         "bgmm fit here; see run_poppunk).")
    ap.add_argument("--no-refine", dest="refine", action="store_false",
                    help="Use only the initial model (the validated setting).")
    ap.add_argument("--reuse-db", action="store_true",
                    help="Reuse an existing PopPUNK sketch database from a prior run "
                         "(skip the expensive --create-db re-sketching).")
    ap.add_argument("--clusters-csv", default=None,
                    help="Skip PopPUNK and normalize an existing raw clusters CSV "
                         "(debug/resume; e.g. a smoke's pp_fit_clusters.csv).")
    # Sketch/QC overrides. These exist so a CONTROLLED comparison never requires
    # hand-editing config.yaml: to test whether a rebuilt container changed the
    # clustering you must re-run the OLD settings against the OLD labels, varying
    # exactly one thing. Everything unset falls through to config + registry.
    ap.add_argument("--min-k", type=int, default=None)
    ap.add_argument("--max-k", type=int, default=None)
    ap.add_argument("--k-step", type=int, default=None)
    ap.add_argument("--sketch-size", type=int, default=None)
    ap.add_argument("--no-qc", dest="qc", action="store_false", default=None,
                    help="Skip the PopPUNK --qc-db pass (e.g. to reproduce a "
                         "pre-QC clustering for comparison).")
    ap.add_argument("--out-name", default="poppunk_clusters.csv",
                    help="Output filename under the lineage dir. Use a different "
                         "name for a comparison run so it does not clobber the "
                         "canonical labels the pipeline reads.")
    args = ap.parse_args()
    refine = lin_cfg.get("refine", True) if args.refine is None else args.refine

    organism = args.organism
    params = lineage_params(organism, config)   # globals + this organism's override
    # CLI wins over config + registry (same precedence as get_target).
    for key in ("min_k", "max_k", "k_step", "sketch_size", "qc"):
        val = getattr(args, key, None)
        if val is not None:
            params[key] = val
    # model/refine are resolved separately (args.model, and `refine` computed
    # above), so fold the ACTUAL values back into params — otherwise the printed
    # provenance line shows config's model/refine while a CLI override ran, i.e.
    # the log lies about what clustered. Must run after `refine` is computed.
    params["model"] = args.model
    params["refine"] = refine
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
        print(f"  Running PopPUNK (model={args.model}, refine={refine}, "
              f"reuse_db={args.reuse_db}, threads={args.threads})...")
        print(f"  Lineage params: {params}")   # provenance: what actually clustered
        raw_clusters = run_poppunk(poppunk, refs_path, work_dir, args.model, args.threads,
                                   params=params, refine=refine, reuse_db=args.reuse_db)

    print(f"  Normalizing clusters (un-mangling '.'→'_') from: {raw_clusters.name}")
    out_df = normalize_clusters(raw_clusters, genome_ids)

    out_path = lineage_dir / args.out_name
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
