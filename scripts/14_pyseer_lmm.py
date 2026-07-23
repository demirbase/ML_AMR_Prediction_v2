#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Population-structure-corrected unitig association (pyseer LMM) — Step 14

Why
---
The model + CPSS tell us which unitigs are *predictive* and *reproducibly
selected*. A reviewer will still ask: are they associated with resistance once
the clonal population structure is accounted for, or do they just track a
lineage? pyseer's linear mixed model (FaST-LMM) fits unitig presence/absence
against the R/S phenotype with a genome-genome kinship random effect that absorbs
population structure (Lees 2018; Jaillard 2018). Unitigs passing the **Bonferroni**
threshold (0.05 / #patterns) are an independent, lineage-corrected cross-check of
the CPSS selection.

Two-container split (each tool in the container that has it)
-----------------------------------------------------------
`pyseer`/`similarity_pyseer` live in ``amr-tools.sif`` (which has NO PyYAML), the
config-driven Python prep/post needs ``amr.sif`` (yaml+pandas). So this script
runs only the Python halves and the SLURM job chains the pyseer CLIs between them:

    amr.sif:        14_pyseer_lmm.py --mode prep   # phenotype + samples + paths.sh
    amr-tools.sif:  similarity_pyseer $SAMPLES --pres $PRES > $SIM
    amr-tools.sif:  pyseer --lmm --phenotypes $PHENO --pres $PRES --similarity $SIM \
                           --output-patterns $PATTERNS --cpu N > $ASSOC
    amr.sif:        14_pyseer_lmm.py --mode post    # Bonferroni + CPSS cross-check

``prep`` writes ``14_pyseer_paths_{ab}.sh`` (PRES/PHENO/SAMPLES/SIM/ASSOC/PATTERNS)
so the SLURM script sources the config-resolved paths instead of hardcoding them.

Output (results/{org}/{ab}/05_explainability/)
    14_pyseer_assoc_{ab}.txt          — raw pyseer association table (from the LMM step)
    14_pyseer_significant_{ab}.csv    — Bonferroni-significant unitigs (+ is_cpss_stable)
    14_pyseer_summary_{ab}.json       — threshold, #tested, #significant, #CPSS-confirmed
"""

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config, resolve_path, get_target  # noqa: E402


def _pyseer_version():
    """The pyseer version that produced this run's associations, or None.

    `--mode post` runs in amr.sif, where pyseer is ABSENT — it ships in
    amr-tools.sif and is invoked separately (run_pyseer_env.slurm). So the direct
    probe below returns None here; the step that ACTUALLY ran pyseer exports
    AMR_PYSEER_VERSION, which we prefer. (The old code claimed "only this step can
    see the binary it called" and then probed a binary that isn't there, leaving
    pyseer_version null in every summary.) Best-effort: never fails the run.
    """
    env_ver = os.environ.get("AMR_PYSEER_VERSION")
    if env_ver and env_ver.strip():
        return env_ver.strip()
    import subprocess
    try:
        r = subprocess.run(["pyseer", "--version"], capture_output=True,
                           text=True, check=False, timeout=30)
        out = (r.stdout or r.stderr or "").strip().splitlines()
        return out[0] if out else None
    except Exception:
        return None


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
    multiple-testing correction (count_patterns.py isn't in the container, so we
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
    return df, sig.sort_values(pcol), pcol


def main():
    ap = argparse.ArgumentParser(description="pyseer LMM unitig association (M14).")
    ap.add_argument("--mode", choices=["prep", "post"], required=True)
    ap.add_argument("--organism", default=None)
    ap.add_argument("--antibiotic", default=None)
    ap.add_argument("--pres", default=None, help="unitig Rtab (default matrix_dir/unitigs.rtab)")
    ap.add_argument("--similarity", default=None, help="kinship matrix path")
    args = ap.parse_args()

    config = load_config()
    organism, antibiotic = get_target(args, config=config)
    matrix_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic, config=config)
    out_dir = resolve_path("dir_05_explainability", organism=organism,
                           antibiotic=antibiotic, config=config)
    out_dir.mkdir(parents=True, exist_ok=True)

    pres = Path(args.pres) if args.pres else (matrix_dir / "unitigs.rtab")
    pheno = out_dir / f"14_phenotype_{antibiotic}.tsv"
    samples_txt = out_dir / f"14_samples_{antibiotic}.txt"
    sim = Path(args.similarity) if args.similarity else (out_dir / f"14_similarity_{antibiotic}.tsv")
    assoc = out_dir / f"14_pyseer_assoc_{antibiotic}.txt"
    patterns = out_dir / f"14_patterns_{antibiotic}.txt"

    if args.mode == "prep":
        if not pres.exists():
            print(f"ERROR: unitig Rtab not found: {pres}"); sys.exit(1)
        nR, nS = write_phenotype(matrix_dir / f"genomes_{antibiotic}.csv",
                                 matrix_dir / f"y_{antibiotic}.csv", pheno, samples_txt)
        paths_sh = out_dir / f"14_pyseer_paths_{antibiotic}.sh"
        paths_sh.write_text("\n".join([
            f'PRES="{pres}"', f'PHENO="{pheno}"', f'SAMPLES="{samples_txt}"',
            f'SIM="{sim}"', f'ASSOC="{assoc}"', f'PATTERNS="{patterns}"']) + "\n",
            encoding="utf-8")
        print(f"  prep: phenotype {nR} R / {nS} S -> {pheno.name}")
        print(f"  paths for SLURM -> {paths_sh.name}")
        return

    # mode == post
    if not assoc.exists() or assoc.stat().st_size == 0:
        print(f"ERROR: association table missing/empty: {assoc}\n"
              f"  Run the pyseer --lmm step (amr-tools.sif) first."); sys.exit(1)
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
        "n_variants_tested": int(len(df)), "n_significant": int(len(sig)),
        "n_cpss_stable_significant": n_cpss_sig, "n_cpss_stable_total": len(cpss_kmers),
        "pvalue_column": pcol,
        # pyseer lives in amr-tools.sif, populate runs in amr.sif — so the step
        # that actually invokes the tool is the only one that can honestly report
        # its version. populate reads it back from here into pipeline_runs.
        "pyseer_version": _pyseer_version(),
    }
    (out_dir / f"14_pyseer_summary_{antibiotic}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print("=" * 74)
    print(f"  threshold {threshold:.2e} ({n_pat} patterns) | tested {len(df)} | "
          f"significant {len(sig)} | CPSS-stable & significant {n_cpss_sig}/{len(cpss_kmers)}")
    print(f"  ✓ 14_pyseer_significant_{antibiotic}.csv  ✓ 14_pyseer_summary_{antibiotic}.json")
    print("=" * 74)


if __name__ == "__main__":
    main()
