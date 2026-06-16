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
ALL_STEPS: list[tuple[str, str]] = [
    ("00a", "00a_download_bvbrc.py"),
    ("00",  "00_prepare_metadata.py"),
    ("01",  "01_data_validation.py"),
    ("02",  "02_kmer_extraction.py"),
    ("02b", "02b_global_qc_analysis.py"),
    ("03",  "03_matrix_construction.py"),
    ("04",  "04_optimization.py"),
    ("05",  "05_model_training.py"),
    ("06",  "06_evaluation.py"),
    ("07b", "07b_feature_stability.py"),
    ("07",  "07_explainability.py"),
    ("08",  "08_blast_annotation.py"),
    ("09",  "09_biological_summary.py"),
    ("10",  "10_kmer_background_frequency.py"),
    ("11",  "11_variant_snp_check.py"),
]
# Default plan: the analysis core (skip network/DB-heavy data + BLAST steps).
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
