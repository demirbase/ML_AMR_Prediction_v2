#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Smoke tests — load every numbered script and assert it imports cleanly.

This is the cheapest, highest-value safety net: it catches syntax errors, bad
imports, broken config wiring, and unresolved paths in SECONDS — instead of
discovering them hours into a multi-day pipeline run. Scripts that need a heavy
dependency missing from this environment (xgboost/optuna) are skipped, so the
same suite runs fully in your real environment.
"""

import pytest

PIPELINE_SCRIPTS = [
    "00a_download_bvbrc.py",
    "00_prepare_metadata.py",
    "01_data_validation.py",
    "01b_data_validation.py",
    "02_kmer_extraction.py",
    "02p_kmer_parallel.py",
    "02b_global_qc_analysis.py",
    "02c_lineage_poppunk.py",
    "03_matrix_construction.py",
    "03u_unitig_matrix.py",
    "03b_matrix_validation_qc.py",
    "04_optimization.py",
    "05_model_training.py",
    "06_evaluation.py",
    "07_explainability.py",
    "07b_feature_stability.py",
    "08_blast_annotation.py",
    "09_biological_summary.py",
    "10_kmer_background_frequency.py",
    "11_variant_snp_check.py",
    "migrate_to_organism_layout.py",
]


@pytest.mark.smoke
@pytest.mark.parametrize("script", PIPELINE_SCRIPTS)
def test_script_imports(load_script, script):
    """Each script loads without syntax/import/config errors (or skips on missing dep)."""
    mod = load_script(script)
    assert mod is not None


@pytest.mark.smoke
def test_lib_package_imports():
    """The shared lib package and its public API import cleanly."""
    from lib import registry, config, chunking, io_utils, run_metadata  # noqa: F401
    assert callable(config.resolve_path)
    assert callable(registry.load_antibiotic_classes)
    assert callable(chunking.get_y_chunk)
    assert callable(io_utils.run_command)
    assert callable(run_metadata.make_run_id)


@pytest.mark.smoke
def test_compat_shims_import():
    """Backward-compat shims still expose the legacy names."""
    import constants
    import utils
    assert isinstance(constants.ANTIBIOTIC_CLASSES, dict)
    assert callable(utils.get_y_chunk)
    assert callable(utils.run_command)
