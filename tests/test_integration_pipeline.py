#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
End-to-end integration test on a TINY synthetic dataset.

This is the answer to "to test the pipeline I rerun everything for days": it
drives the REAL scripts (01 → 07) as subprocesses on ~16 synthetic genomes
under a throwaway organism `testorg`, so the whole chain finishes in minutes.
The 9 numbered scripts are NOT modified — the test only:
  * stages synthetic inputs at the organism's resolved paths,
  * temporarily swaps in a small test config (restored in a finally block),
  * runs each step and asserts its key outputs appear.

Opt-in (marked integration + slow); a plain `pytest` never runs it and never
touches your real config. Run it with:

    pytest -m integration -s

Requirements (skipped automatically if missing):
  * xgboost + optuna (steps 04–07)
  * the bundled KMC binary at bin/bin/kmc (steps 02–03)

Steps 08/09 (BLAST/Nextflow/NCBI network) are intentionally out of scope here.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = PROJECT_ROOT / "scripts"
CONFIG = PROJECT_ROOT / "config" / "config.yaml"

pytestmark = [pytest.mark.integration, pytest.mark.slow]


def _require_tools():
    """Skip the whole module unless the heavy deps / KMC binary are present."""
    for mod in ("xgboost", "optuna"):
        try:
            __import__(mod)
        except ImportError:
            pytest.skip(f"integration test needs '{mod}' (not installed)")
    kmc = PROJECT_ROOT / "bin" / "bin" / "kmc"
    if not kmc.exists():
        pytest.skip(f"integration test needs the KMC binary at {kmc}")


def _run_step(script, env_note=""):
    """Run `python scripts/<script>` from the repo root; assert success."""
    proc = subprocess.run(
        [sys.executable, str(SCRIPTS / script)],
        cwd=str(PROJECT_ROOT), capture_output=True, text=True,
    )
    if proc.returncode != 0:
        print(f"\n----- {script} STDOUT -----\n{proc.stdout[-3000:]}")
        print(f"\n----- {script} STDERR -----\n{proc.stderr[-3000:]}")
        pytest.fail(f"{script} failed (rc={proc.returncode}) {env_note}")
    return proc


def _cleanup(organism, antibiotic=None):
    """Remove every artifact the test created for `organism`/`antibiotic`."""
    targets = [
        PROJECT_ROOT / "data" / "raw" / organism,
        PROJECT_ROOT / "data" / "external" / organism,
        PROJECT_ROOT / "data" / "interim" / organism,
        PROJECT_ROOT / "data" / "processed" / organism,
        PROJECT_ROOT / "models" / organism,
        PROJECT_ROOT / "results" / organism,
        PROJECT_ROOT / "logs" / organism,
        PROJECT_ROOT / "runs" / organism,
        PROJECT_ROOT / "config" / "experiments" / organism,
        PROJECT_ROOT / "config" / f"config_{organism}.yaml",
    ]
    if antibiotic:
        # Step 04 writes config/config_{antibiotic}.yaml (flat, not organism-scoped).
        targets.append(PROJECT_ROOT / "config" / f"config_{antibiotic}.yaml")
    for t in targets:
        if t.is_dir():
            shutil.rmtree(t, ignore_errors=True)
        elif t.exists():
            t.unlink()


def test_pipeline_end_to_end(synthetic_dataset):
    _require_tools()

    organism = synthetic_dataset["organism"]
    antibiotic = synthetic_dataset["antibiotic"]
    k_length = synthetic_dataset["k_length"]

    # Resolve the organism-scoped input locations the scripts will read.
    sys.path.insert(0, str(SCRIPTS))
    from lib.config import resolve_path  # noqa: E402

    base_cfg = yaml.safe_load(CONFIG.read_text())
    genomes_dst = resolve_path("raw_genomes_dir", organism=organism, config=base_cfg)
    meta_dst = resolve_path("metadata_file", organism=organism, config=base_cfg)

    config_backup = CONFIG.read_text()
    try:
        # --- stage synthetic inputs -------------------------------------------
        _cleanup(organism, antibiotic)
        genomes_dst.mkdir(parents=True, exist_ok=True)
        for fna in synthetic_dataset["genomes_dir"].glob("*.fna"):
            shutil.copy2(fna, genomes_dst / fna.name)
        meta_dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(synthetic_dataset["metadata_file"], meta_dst)

        # --- write a small test config (restored in finally) ------------------
        cfg = yaml.safe_load(CONFIG.read_text())
        cfg["project"]["organism"] = organism
        cfg["project"]["target_antibiotic"] = antibiotic
        cfg["preprocessing"]["k_length"] = k_length
        cfg["preprocessing"]["min_support"] = 2
        cfg["preprocessing"]["chunk_size"] = 4
        cfg["preprocessing"]["threads"] = 2
        cfg["preprocessing"]["kmc_mem"] = 4
        cfg["training"]["n_trials"] = 2
        cfg["training"]["test_fraction"] = 0.25
        cfg["training"]["validation_fraction"] = 0.5
        cfg["training"]["optuna_fraction"] = 0.5
        CONFIG.write_text(yaml.safe_dump(cfg, sort_keys=False))

        matrix_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic, config=cfg)
        models_dir = resolve_path("models_dir", organism=organism, antibiotic=antibiotic, config=cfg)
        explain_dir = resolve_path("dir_05_explainability", organism=organism, antibiotic=antibiotic, config=cfg)

        # --- 02: k-mer counting ----------------------------------------------
        _run_step("02_kmer_extraction.py")
        kmc_out = resolve_path("kmc_outputs_dir", organism=organism, config=cfg)
        assert any(kmc_out.glob("*.kmc_pre")), "02 produced no KMC databases"

        # --- 03: matrix construction -----------------------------------------
        _run_step("03_matrix_construction.py")
        assert (matrix_dir / "features.txt").exists(), "03 produced no features.txt"
        assert (matrix_dir / f"y_{antibiotic}.csv").exists()
        assert any(matrix_dir.glob(f"X_{antibiotic}_part_*.npz")), "03 produced no matrix chunks"

        # --- 04: optimization (writes the organism-scoped experiment config) --
        _run_step("04_optimization.py")
        ab_cfg = resolve_path("experiment_config", organism=organism,
                              antibiotic=antibiotic, config=cfg)
        assert ab_cfg.exists(), "04 did not write the organism-scoped antibiotic config"

        # --- 05: training (writes model + manifest.json) ---------------------
        _run_step("05_model_training.py")
        assert any(models_dir.glob(f"xgboost_{antibiotic}_final_v2.json")), "05 produced no model"
        assert (models_dir / "manifest.json").exists(), "05 produced no manifest.json"

        # --- 06: evaluation (writes metrics + does NOT mutate the threshold) --
        thresh_before = yaml.safe_load(ab_cfg.read_text()).get("evaluation", {}).get("optimal_threshold")
        _run_step("06_evaluation.py")
        eval_dir = resolve_path("dir_04_evaluation", organism=organism, antibiotic=antibiotic, config=cfg)
        assert any(eval_dir.glob(f"06_comprehensive_metrics_{antibiotic}.csv")), "06 produced no metrics CSV"
        thresh_after = yaml.safe_load(ab_cfg.read_text()).get("evaluation", {}).get("optimal_threshold")
        # Data-leakage guard: 06 must NOT overwrite the training-derived threshold.
        assert thresh_before == thresh_after, "06 overwrote the config threshold (data leakage regression!)"

        # --- 07: explainability (top features CSV + FASTA) -------------------
        _run_step("07_explainability.py")
        top_n = cfg["analysis"]["top_n_features"]
        assert (explain_dir / f"01_top_{top_n}_features_{antibiotic}.csv").exists(), "07 produced no top-features CSV"

        # --- 07b: 5-seed repeated holdout (AUC mean/std, stability, Jaccard) --
        _run_step("07b_feature_stability.py")
        assert (eval_dir / f"10_repeated_holdout_summary_{antibiotic}.csv").exists(), "07b produced no holdout summary"
        assert (explain_dir / f"06_feature_stability_{antibiotic}.csv").exists(), "07b produced no stability CSV"

    finally:
        CONFIG.write_text(config_backup)
        _cleanup(organism, antibiotic)
