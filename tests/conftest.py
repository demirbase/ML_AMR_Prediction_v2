#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared pytest fixtures & helpers for the AMR pipeline test suite.

Goal: let you validate the pipeline in SECONDS/MINUTES instead of re-running the
full multi-day job. Nothing here modifies the 9 numbered scripts — the suite is
purely additive (SCALE_MLOPS_PLAN.md §7.5).
"""

import importlib.util
import random
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

# Make `lib`, `utils`, `constants`, and the numbered scripts importable exactly
# as they are when launched via `python scripts/0X_*.py`.
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))


@pytest.fixture(scope="session")
def repo_root():
    return PROJECT_ROOT


@pytest.fixture
def load_script():
    """
    Return a loader that imports a numbered script module by filename.

    The numbered scripts can't be imported with a normal `import` (names start
    with a digit), so we load them via importlib from their path. If the module
    needs a heavy dependency that isn't installed in this environment
    (e.g. xgboost, optuna), the loader raises pytest.skip with a clear reason —
    so the same test file runs fully in your real env and skips gracefully here.

    Usage:
        def test_x(load_script):
            mod = load_script("01_data_validation.py")
            assert mod.validate_dataset_scientific(50, 50)[0] is True
    """
    def _load(filename):
        path = SCRIPTS_DIR / filename
        if not path.exists():
            pytest.fail(f"Script not found: {path}")
        mod_name = "amrtest_" + path.stem.replace(".", "_")
        spec = importlib.util.spec_from_file_location(mod_name, path)
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
        except ImportError as e:
            pytest.skip(f"{filename} needs an unavailable dependency: {e}")
        except SystemExit as e:
            pytest.skip(f"{filename} called sys.exit at import (config/data missing): {e}")
        return module

    return _load


# ---------------------------------------------------------------------------
# Synthetic mini-dataset
# ---------------------------------------------------------------------------
def _random_dna(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


@pytest.fixture
def synthetic_dataset(tmp_path):
    """
    Build a tiny, deterministic genomic dataset for end-to-end testing.

    16 genomes (~4 kb each), one binary phenotype ('gentamicin'). A short
    "resistance motif" is planted into the 8 resistant genomes so a model can
    learn real signal — letting the full 02→07 chain run in minutes instead of
    days.

    Returns a dict with: genomes_dir, metadata_file, organism, antibiotic,
    ids, labels, k_length, motif.
    """
    rng = random.Random(1234)
    organism = "testorg"
    # Deliberately NOT a real antibiotic name: step 04 writes
    # config/config_{antibiotic}.yaml, so using e.g. "gentamicin" would clobber
    # the user's real config_gentamicin.yaml. "testdrug" cannot collide.
    antibiotic = "testdrug"
    k_length = 13            # small k for tiny genomes (config-driven downstream)
    motif = "ACGTACGTACGTACGT"  # planted resistance signal (>= k_length)

    genomes_dir = tmp_path / "genomes"
    genomes_dir.mkdir(parents=True, exist_ok=True)

    n = 16
    ids, labels = [], []
    for i in range(n):
        gid = f"test.{i+1}"
        # Interleave classes (R,S,R,S,…) so that every contiguous chunk the
        # pipeline builds is class-balanced. This keeps the test set from being
        # single-class (roc_auc needs both classes) and exercises the realistic
        # mixed-chunk path. The pure-chunk edge case is covered separately by the
        # base_score fix in 05_model_training.py.
        resistant = 1 if i % 2 == 0 else 0
        seq = _random_dna(rng, 4000)
        if resistant:
            # Insert the motif a few times so it clears any min_support filter.
            for pos in (500, 1500, 2500):
                seq = seq[:pos] + motif + seq[pos + len(motif):]
        fna = genomes_dir / f"{gid}.fna"
        with open(fna, "w") as f:
            f.write(f">{gid} synthetic test genome\n")
            for j in range(0, len(seq), 70):
                f.write(seq[j:j + 70] + "\n")
        ids.append(gid)
        labels.append(resistant)

    metadata_file = tmp_path / "amr_phenotypes.csv"
    with open(metadata_file, "w") as f:
        f.write(f"Genome ID,{antibiotic}\n")
        for gid, lab in zip(ids, labels):
            f.write(f"{gid},{lab}\n")

    return {
        "genomes_dir": genomes_dir,
        "metadata_file": metadata_file,
        "organism": organism,
        "antibiotic": antibiotic,
        "ids": ids,
        "labels": labels,
        "k_length": k_length,
        "motif": motif,
        "tmp": tmp_path,
    }
