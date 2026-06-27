#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Population-structure-corrected unitig association (pyseer LMM) — Step 14

Why
---
The model + CPSS tell us which unitigs are *predictive* and *reproducibly
selected*. A reviewer will still ask: are they associated with resistance once
the clonal population structure is accounted for, or do they just track a
lineage? pyseer's linear mixed model (FaST-LMM) answers this — it fits the unitig
presence/absence against the R/S phenotype with a random effect (a genome-genome
kinship/similarity matrix) that absorbs population structure (Lees 2018; Jaillard
2018). Unitigs passing the **Bonferroni** threshold (0.05 / #patterns) are an
independent, lineage-corrected cross-check of the CPSS selection.

This script is an orchestrator (pyseer/similarity_pyseer/count_patterns are CLI
tools in amr-tools.sif) + post-processing:
  1. write the phenotype TSV (samples + R/S) from the matrix labels;
  2. compute the kinship/similarity matrix from the genome-wide unitig presence
     (``similarity_pyseer``) unless one is supplied;
  3. run ``pyseer --lmm`` over the unitig Rtab (``--pres unitigs.rtab``);
  4. Bonferroni threshold via ``count_patterns.py``; flag significant unitigs and
     mark which ones are CPSS-stable (step 13).

Inputs default to the unitig matrix dir (``unitigs.rtab``, ``y_{ab}.csv``,
``genomes_{ab}.csv``). Heavy steps are skipped when their output already exists,
so the (long) LMM can be resumed.

Output (results/{org}/{ab}/05_explainability/)
    14_pyseer_assoc_{ab}.txt          — raw pyseer association table
    14_pyseer_significant_{ab}.csv    — Bonferroni-significant unitigs (+ is_cpss_stable)
    14_pyseer_summary_{ab}.json       — threshold, #tested, #significant, #CPSS-confirmed
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config, resolve_path  # noqa: E402


def _tool(name):
    p = shutil.which(name)
    if not p:
        print(f"ERROR: '{name}' not on PATH (run inside amr-tools.sif).")
        sys.exit(1)
    return p


def write_phenotype(genomes_csv, y_csv, out_tsv, samples_txt):
    """pyseer phenotype TSV ('samples<TAB>resistant') + a plain sample-name list
    (one genome id per line) for similarity_pyseer's positional argument."""
    g = pd.read_csv(genomes_csv, encoding="utf-8")
    gid = g[g.columns[0]].astype(str)          # 'Genome ID' (first col)
    y = pd.read_csv(y_csv, encoding="utf-8")["label"].astype(int)
    pd.DataFrame({"samples": gid, "resistant": y.values}).to_csv(
        out_tsv, sep="\t", index=False)
    Path(samples_txt).write_text("\n".join(gid) + "\n", encoding="utf-8")
    return int((y == 1).sum()), int((y == 0).sum())


def bonferroni_threshold(patterns_file):
    """0.05 / (number of unique presence/absence patterns) — pyseer's standard
    multiple-testing correction (count_patterns.py isn't in this container, so we
    count unique pattern hashes from --output-patterns directly)."""
    pats = set()
    with open(patterns_file, encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if s:
                pats.add(s)
    n = len(pats)
    return (0.05 / n if n else float("nan")), n


def parse_and_flag(assoc_file, threshold, cpss_kmers):
    """Read pyseer output, keep Bonferroni-significant variants, flag CPSS-stable."""
    df = pd.read_csv(assoc_file, sep="\t")
    pcol = "lrt-pvalue" if "lrt-pvalue" in df.columns else "filter-pvalue"
    df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
    sig = df[df[pcol] <= threshold].copy()
    sig["is_cpss_stable"] = sig["variant"].astype(str).isin(cpss_kmers).astype(int)
    sig = sig.sort_values(pcol)
    return df, sig, pcol


def main():
    ap = argparse.ArgumentParser(description="pyseer LMM unitig association (M14).")
    ap.add_argument("--organism", default=None)
    ap.add_argument("--antibiotic", default=None)
    ap.add_argument("--pres", default=None, help="unitig Rtab (default: matrix_dir/unitigs.rtab)")
    ap.add_argument("--similarity", default=None, help="precomputed kinship (else computed)")
    ap.add_argument("--cpus", type=int, default=8)
    args = ap.parse_args()

    config = load_config()
    organism = args.organism or config["project"].get("organism", "ecoli")
    antibiotic = args.antibiotic or config["project"]["target_antibiotic"]
    matrix_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic, config=config)
    out_dir = resolve_path("dir_05_explainability", organism=organism,
                           antibiotic=antibiotic, config=config)
    out_dir.mkdir(parents=True, exist_ok=True)

    pres = Path(args.pres) if args.pres else (matrix_dir / "unitigs.rtab")
    if not pres.exists():
        print(f"ERROR: unitig Rtab not found: {pres}"); sys.exit(1)

    pyseer = _tool("pyseer")
    sim_bin = _tool("similarity_pyseer")

    print("=" * 74)
    print(f"PYSEER LMM  —  {organism} / {antibiotic}")
    print(f"  pres: {pres.name}")
    print("=" * 74)

    # 1) phenotype -----------------------------------------------------------
    pheno = out_dir / f"14_phenotype_{antibiotic}.tsv"
    samples_txt = out_dir / f"14_samples_{antibiotic}.txt"
    nR, nS = write_phenotype(matrix_dir / f"genomes_{antibiotic}.csv",
                             matrix_dir / f"y_{antibiotic}.csv", pheno, samples_txt)
    print(f"  [1/4] phenotype: {nR} R / {nS} S  -> {pheno.name}")

    # 2) kinship / similarity (genome-wide population structure) --------------
    sim = Path(args.similarity) if args.similarity else (out_dir / f"14_similarity_{antibiotic}.tsv")
    if not sim.exists():
        print("  [2/4] computing similarity (similarity_pyseer, streamed)...", flush=True)
        with open(sim, "w") as fh:
            subprocess.run([sim_bin, str(samples_txt), "--pres", str(pres)],
                           stdout=fh, check=True)
    else:
        print(f"  [2/4] using existing similarity: {sim.name}")

    # 3) LMM -----------------------------------------------------------------
    assoc = out_dir / f"14_pyseer_assoc_{antibiotic}.txt"
    patterns = out_dir / f"14_patterns_{antibiotic}.txt"
    if not assoc.exists() or assoc.stat().st_size == 0:
        print("  [3/4] pyseer --lmm (streamed over unitigs)...", flush=True)
        with open(assoc, "w") as fh:
            subprocess.run([pyseer, "--lmm", "--phenotypes", str(pheno),
                            "--pres", str(pres), "--similarity", str(sim),
                            "--output-patterns", str(patterns),
                            "--cpu", str(args.cpus)], stdout=fh, check=True)
    else:
        print(f"  [3/4] using existing association: {assoc.name}")

    # 4) Bonferroni + CPSS cross-check ---------------------------------------
    threshold, n_pat = bonferroni_threshold(patterns)
    cpss_csv = out_dir / f"13_stability_selection_{antibiotic}.csv"
    cpss_kmers = set()
    if cpss_csv.exists():
        c = pd.read_csv(cpss_csv, encoding="utf-8")
        cpss_kmers = set(c[c["stable"] == 1]["kmer"].astype(str))

    df, sig, pcol = parse_and_flag(assoc, threshold, cpss_kmers)
    sig.to_csv(out_dir / f"14_pyseer_significant_{antibiotic}.csv", index=False)

    n_cpss_sig = int(sig["is_cpss_stable"].sum()) if not sig.empty else 0
    summary = {
        "antibiotic": antibiotic, "organism": organism,
        "n_patterns": n_pat, "bonferroni_threshold": threshold,
        "n_variants_tested": int(len(df)),
        "n_significant": int(len(sig)),
        "n_cpss_stable_significant": n_cpss_sig,
        "n_cpss_stable_total": len(cpss_kmers),
        "pvalue_column": pcol,
    }
    (out_dir / f"14_pyseer_summary_{antibiotic}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 74)
    print(f"  threshold {threshold:.2e} ({n_pat} patterns) | tested {len(df)} | "
          f"significant {len(sig)} | CPSS-stable & significant {n_cpss_sig}/{len(cpss_kmers)}")
    print(f"  ✓ 14_pyseer_significant_{antibiotic}.csv  ✓ 14_pyseer_summary_{antibiotic}.json")
    print("=" * 74)


if __name__ == "__main__":
    main()
