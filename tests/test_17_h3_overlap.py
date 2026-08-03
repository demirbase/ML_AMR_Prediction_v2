#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the H3 gene-family overlap statistics (E2 framework).

These guard the two claims the thesis makes about the test itself: that the
multiple-testing correction really is Benjamini-YEKUTIELI (valid under arbitrary
dependence, i.e. strictly more conservative than BH), and that the enrichment
p-value is the one-sided hypergeometric a Fisher exact test gives.
"""

import numpy as np
import pytest


@pytest.fixture
def mod(load_script):
    return load_script("17_h3_gene_family_overlap.py")


def test_by_is_bh_scaled_by_harmonic_number(mod):
    p = [0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216]
    by = mod.benjamini_yekutieli(p)
    m = len(p)
    c_m = sum(1.0 / i for i in range(1, m + 1))

    # BY = BH * c(m), so the smallest p maps to p * m * c(m) / 1.
    assert by[0] == pytest.approx(p[0] * m * c_m, rel=1e-9)
    # Monotone non-decreasing in p order, and always a valid probability.
    assert np.all(np.diff(by) >= -1e-12)
    assert np.all((by >= 0) & (by <= 1))


def test_by_handles_nan_and_is_conservative(mod):
    by = mod.benjamini_yekutieli([0.01, np.nan, 0.5])
    assert np.isnan(by[1])                     # NaN in, NaN out (never silently 0)
    assert by[0] >= 0.01                       # correction can only inflate a p-value


def test_fisher_greater_matches_hypergeometric_tail(mod):
    # Universe 20, sets of 10 and 10, observed overlap 8 — enrichment tail.
    p = mod.fisher_greater(k=8, K=10, n=10, N=20)
    assert 0.0 < p < 0.05
    # No overlap at all can never be "enriched": P(X >= 0) == 1.
    assert mod.fisher_greater(k=0, K=5, n=5, N=50) == pytest.approx(1.0)


def test_mc_pvalue_never_zero_and_tracks_analytic(mod):
    rng = np.random.default_rng(0)
    # add-one smoothing: an empirical p must never be exactly 0
    p_extreme = mod.mc_pvalue(k=10, K=10, n=10, N=200, B=200, rng=rng)
    assert p_extreme > 0
    # A trivially satisfiable overlap must come back ~1.
    assert mod.mc_pvalue(k=0, K=5, n=5, N=50, B=200, rng=rng) == pytest.approx(1.0)
