#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pipeline orchestrator — run the numbered steps in order for one target.

Thin, dependency-free wrapper that invokes each numbered script as a subprocess
(with the SAME Python interpreter running this file, so the conda environment is
preserved) and logs progress + per-step timing to the console and a log file.
It does not re-implement any step; it just sequences them and fails fast.

Examples:
    python scripts/run_pipeline.py --organism ecoli --antibiotic ampicillin
    python scripts/run_pipeline.py --from 02 --to 07            # a sub-range
    python scripts/run_pipeline.py --only 09 10 11             # specific steps
    python scripts/run_pipeline.py --list                      # show the step plan

Notes:
  - Steps 00a/00 (data acquisition) and 08/11 (BLAST/CARD) are OPTIONAL and only
    run when explicitly selected, since they need network / external databases.
  - The default plan is the analysis core 01 → 10.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.logging_utils import get_logger  # noqa: E402

# Ordered (step-id, script) plan. Step ids are strings to match the file prefixes.
# NOTE: the canonical (unitig, publication-grade) pipeline is HPC/SLURM-driven
# and multi-container — see docs/ROADMAP.md §0 and the TRUBA guide. Several steps
# need a specific container, an env override, --mode prep/post, or internet, so
# they are NOT runnable by this plain single-container orchestrator (they are
# listed in HPC_SLURM_STEPS below for visibility). This orchestrator is a local
# convenience for the single-container Python core; on TRUBA each step is a
# targeted SLURM job.
ALL_STEPS: list[tuple[str, str]] = [
    ("00a", "00a_download_bvbrc.py"),
    ("00",  "00_prepare_metadata.py"),
    ("01",  "01_data_validation.py"),
    ("02",  "02_kmer_extraction.py"),        # KMC (QC + k-mer baseline)
    ("02b", "02b_global_qc_analysis.py"),
    ("02c", "02c_lineage_poppunk.py"),       # PopPUNK lineage (amr-pp.sif)
    ("02d", "02d_genome_qc.py"),             # CheckM2+QUAST QC (M15; --mode, multi-container)
    ("03",  "03_matrix_construction.py"),    # raw k-mer matrix (baseline)
    ("03u", "03u_unitig_matrix.py"),         # unitig matrix (CANONICAL; AMR_FEATURE_REPR=unitig)
    ("03b", "03b_matrix_validation_qc.py"),
    ("04",  "04_optimization.py"),
    ("05",  "05_model_training.py"),
    ("06",  "06_evaluation.py"),
    ("07b", "07b_feature_stability.py"),
    ("07",  "07_explainability.py"),
    ("08",  "08_blast_annotation.py"),       # CARD local + NCBI remote (internet)
    ("09",  "09_biological_summary.py"),
    ("10",  "10_kmer_background_frequency.py"),
    ("11",  "11_variant_snp_check.py"),
    ("12",  "12_permutation_test.py"),       # MDA permutation (M9)
    ("12b", "12b_label_permutation_test.py"),# label-permutation null (M9)
    ("13",  "13_stability_selection.py"),    # CPSS + SHAP (M4)
    ("13b", "13b_stable_annotation.py"),
    ("14",  "14_pyseer_lmm.py"),             # pyseer LMM (M14; --mode, amr-tools.sif)
    ("15",  "15_cross_antibiotic.py"),       # cross-antibiotic overlap / H3 (S1)
    ("16",  "16_external_concordance.py"),   # AMRFinderPlus/ResFinder concordance (M13; --mode)
]
# Steps that need SLURM + a specific container / env / --mode / internet and so
# cannot be launched by this plain orchestrator (run them as SLURM jobs):
#   02c (amr-pp.sif) · 02d (--mode + amr-checkm2/amr-tools) · 03u (AMR_FEATURE_REPR=unitig)
#   08-NCBI (internet) · 14 (--mode + amr-tools) · 16 (--mode + amr-tools) · populate_database.py
HPC_SLURM_STEPS = {"02c", "02d", "14", "16"}

# Default plan: the local single-container analysis core (raw-k-mer baseline).
# The canonical unitig run uses 03u (+ AMR_FEATURE_REPR=unitig) and runs on HPC.
DEFAULT_PLAN = ["01", "02", "02b", "03", "04", "05", "06", "07b", "07", "09", "10"]


def _index(step_id: str) -> int:
    for i, (sid, _) in enumerate(ALL_STEPS):
        if sid == step_id:
            return i
    raise SystemExit(f"Unknown step id: {step_id} (valid: {[s for s, _ in ALL_STEPS]})")


def select_steps(args) -> list[tuple[str, str]]:
    if args.only:
        wanted = set(args.only)
        return [(s, f) for s, f in ALL_STEPS if s in wanted]
    if args.from_ or args.to:
        lo = _index(args.from_) if args.from_ else 0
        hi = _index(args.to) if args.to else len(ALL_STEPS) - 1
        return ALL_STEPS[lo:hi + 1]
    return [(s, f) for s, f in ALL_STEPS if s in DEFAULT_PLAN]


def main() -> None:
    p = argparse.ArgumentParser(description="Run the AMR pipeline steps in order.")
    p.add_argument("--organism", default=None, help="registry slug (default: config)")
    p.add_argument("--antibiotic", default=None, help="antibiotic id (default: config)")
    p.add_argument("--from", dest="from_", default=None, help="first step id (e.g. 02)")
    p.add_argument("--to", default=None, help="last step id (e.g. 09)")
    p.add_argument("--only", nargs="+", default=None, help="run only these step ids")
    p.add_argument("--list", action="store_true", help="print the step plan and exit")
    p.add_argument("--continue-on-error", action="store_true",
                   help="keep going if a step fails (default: stop)")
    args = p.parse_args()

    steps = select_steps(args)
    if args.list:
        print("Planned steps:")
        for sid, script in steps:
            print(f"  [{sid:>3}] {script}")
        return

    # Per-invocation overrides consumed by lib.config.get_target in each script.
    env = os.environ.copy()
    if args.organism:
        env["AMR_ORGANISM"] = args.organism
    if args.antibiotic:
        env["AMR_ANTIBIOTIC"] = args.antibiotic

    org = args.organism or "default"
    log = get_logger("pipeline", logfile=PROJECT_ROOT / "logs" / f"run_pipeline_{org}.log")
    log.info("Pipeline plan: %s", " -> ".join(s for s, _ in steps))
    if args.organism or args.antibiotic:
        log.info("Target override via env: AMR_ORGANISM=%s AMR_ANTIBIOTIC=%s",
                 args.organism, args.antibiotic)

    # The target is taken from config.yaml; --organism/--antibiotic are exported
    # as env vars (honoured by scripts that resolve via lib.config.get_target).
    # We do NOT forward them as CLI flags, because not every script defines those
    # arguments and argparse would reject unknown options.
    failures = []
    for sid, script in steps:
        if sid in HPC_SLURM_STEPS and not (args.only and sid in args.only):
            log.warning("SKIP %s (%s): SLURM / multi-container / --mode step — run it "
                        "as a SLURM job (see the TRUBA guide), not via this plain "
                        "orchestrator. (Force by naming it explicitly in --only.)",
                        sid, script)
            continue
        # Just-in-time resolution for --antibiotic auto (after metadata is prepared)
        if sid >= "01" and env.get("AMR_ANTIBIOTIC") == "auto":
            try:
                from lib.config import get_target
                import pandas as pd
                import yaml
                
                org, _ = get_target()
                reg_path = PROJECT_ROOT / "config/registry/organisms.yaml"
                with open(reg_path) as f:
                    registry = yaml.safe_load(f)
                
                org_conf = registry["organisms"].get(org, {})
                candidates = org_conf.get("antibiotics", [])
                meta_file = PROJECT_ROOT / org_conf.get("metadata_file", f"data/external/{org}/metadata/amr_phenotypes.csv")
                
                if not meta_file.exists():
                    best_ab = candidates[0] if candidates else "ampicillin"
                    log.warning("Metadata %s not found. Auto-selecting '%s' fallback.", meta_file.name, best_ab)
                else:
                    df = pd.read_csv(meta_file)
                    best_ab = candidates[0] if candidates else "ampicillin"
                    best_score = -1
                    for ab in candidates:
                        if ab in df.columns:
                            # amr_phenotypes.csv is binary 0/1 (00_prepare_metadata),
                            # NOT the strings "Resistant"/"Susceptible" (audit Issue 6:
                            # the old value_counts lookup always scored 0 -> auto always
                            # picked candidates[0]).
                            col = df[df[ab].notna()][ab]
                            r = int((col == 1).sum())
                            s = int((col == 0).sum())
                            score = min(r, s)  # Maximize the minority class size
                            if score > best_score:
                                best_score = score
                                best_ab = ab
                    log.info("Auto-selected ideal antibiotic '%s' based on class balance (minority class size: %d).", best_ab, best_score)
                env["AMR_ANTIBIOTIC"] = best_ab
            except Exception as e:
                log.warning("Failed to auto-select antibiotic, using fallback 'ampicillin': %s", e)
                env["AMR_ANTIBIOTIC"] = "ampicillin"

        script_path = PROJECT_ROOT / "scripts" / script
        log.info("=== STEP %s : %s ===", sid, script)
        t0 = time.time()
        rc = subprocess.run([sys.executable, str(script_path)], env=env).returncode
        dt = time.time() - t0
        if rc == 0:
            log.info("STEP %s done in %.1fs", sid, dt)
        else:
            log.error("STEP %s FAILED (exit %s) after %.1fs", sid, rc, dt)
            failures.append(sid)
            if not args.continue_on_error:
                sys.exit(f"Pipeline aborted at step {sid}.")

    if failures:
        sys.exit(f"Pipeline finished with failures: {failures}")
    log.info("Pipeline complete: %d steps OK.", len(steps))


if __name__ == "__main__":
    main()
