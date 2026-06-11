#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unit tests for the bug-prone pure functions inside the numbered scripts.

These pin down the exact logic that the three audit reports flagged (R/S
counting, data-leakage threshold, BLAST confidence tiers, √p colsample range,
bootstrap CI), so a regression is caught in seconds rather than after a
multi-day run. Functions living in modules that need xgboost/optuna are loaded
via the `load_script` fixture, which skips them when that dep is unavailable.
"""

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# 01_data_validation.validate_dataset_scientific
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_validate_dataset_missing_class(load_script):
    m = load_script("01_data_validation.py")
    assert m.validate_dataset_scientific(0, 50) == (False, "MISSING CLASS")
    assert m.validate_dataset_scientific(50, 0) == (False, "MISSING CLASS")


@pytest.mark.unit
def test_validate_dataset_insufficient_minority(load_script):
    m = load_script("01_data_validation.py")
    ok, msg = m.validate_dataset_scientific(10, 1000)  # minority 10 < 40
    assert ok is False and "INSUFFICIENT" in msg


@pytest.mark.unit
def test_validate_dataset_valid_and_imbalanced(load_script):
    m = load_script("01_data_validation.py")
    # 500/500: total>=1000 -> min_ratio 5%, minority 50% -> VALID
    assert m.validate_dataset_scientific(500, 500) == (True, "VALID")
    # 45/4955: total>=2000 -> min_ratio 2%, minority ~0.9% -> IMBALANCED
    ok, msg = m.validate_dataset_scientific(45, 4955)
    assert ok is False and "IMBALANCED" in msg


# ---------------------------------------------------------------------------
# 09_biological_summary — confidence tiers & sseqid parsing
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_confidence_tiers(load_script):
    m = load_script("09_biological_summary.py")
    assert m.classify_confidence(99.0, 1e-5) == "confirmed"
    assert m.classify_confidence(92.0, 0.5) == "candidate"
    assert m.classify_confidence(80.0, 0.5) == "weak"     # identity too low
    assert m.classify_confidence(99.0, 5.0) == "weak"     # 21-mer E=5 is noise


@pytest.mark.unit
def test_card_gene_and_accession_parsing(load_script):
    m = load_script("09_biological_summary.py")
    assert m.extract_card_gene("gb|NG_068181.1|+|100-925|ARO:3006096|OXA-909") == "OXA-909"
    assert m._extract_accession("ref|NZ_CP012345.1|") == "NZ_CP012345.1"
    assert m._extract_accession("NZ_CP012345.1") == "NZ_CP012345.1"


@pytest.mark.unit
def test_clean_ncbi_stitle(load_script):
    m = load_script("09_biological_summary.py")
    cleaned = m.clean_ncbi_stitle("Escherichia coli strain K12, complete genome")
    assert "complete genome" not in cleaned.lower()


# ---------------------------------------------------------------------------
# 04_optimization — √p colsample range & feature counting
# (module imports xgboost -> load_script skips when unavailable)
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_colsample_range_brackets_sqrt_p(load_script):
    m = load_script("04_optimization.py")
    lower, upper = m.compute_colsample_range(1_000_000)   # 1/sqrt(p) = 1e-3
    assert lower < 1e-3 < upper
    assert 0 < lower < upper <= 1.0


@pytest.mark.unit
def test_colsample_range_fallback(load_script):
    m = load_script("04_optimization.py")
    assert m.compute_colsample_range(0) == (1e-4, 1e-1)


@pytest.mark.unit
def test_count_features(load_script, tmp_path):
    m = load_script("04_optimization.py")
    feats = tmp_path / "features.txt"
    feats.write_text("AAA 3\nCCC 5\nGGG 1\n")
    assert m.count_features(tmp_path) == 3
    assert m.count_features(tmp_path / "nope") == 0


# ---------------------------------------------------------------------------
# 06_evaluation — bootstrap CI (module imports xgboost -> may skip)
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_bootstrap_ci_perfect_separation(load_script):
    from sklearn.metrics import roc_auc_score
    m = load_script("06_evaluation.py")
    y = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    prob = np.array([0.1, 0.2, 0.15, 0.05, 0.9, 0.85, 0.95, 0.8])
    point, lo, hi = m.bootstrap_metric_ci(y, prob, roc_auc_score, n_bootstraps=200)
    assert point == 1.0
    assert 0.0 <= lo <= hi <= 1.0
    assert hi == pytest.approx(1.0, abs=1e-9)
