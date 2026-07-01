#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for lib/concordance.py (M13 external-validation metrics)."""

import pytest

from lib import concordance as C


def test_confusion_and_basic_metrics():
    # 3 R + 3 S; predictor gets one FN and one FP.
    yt = [1, 1, 1, 0, 0, 0]
    yp = [1, 1, 0, 0, 0, 1]
    cm = C.confusion(yt, yp)
    assert cm == {"TP": 2, "FP": 1, "TN": 2, "FN": 1, "n": 6}
    assert C.sensitivity(cm) == pytest.approx(2 / 3)
    assert C.specificity(cm) == pytest.approx(2 / 3)
    assert C.balanced_accuracy(cm) == pytest.approx(2 / 3)
    fe = C.fda_errors(cm)
    assert fe["major_error_rate"] == pytest.approx(1 / 3)
    assert fe["very_major_error_rate"] == pytest.approx(1 / 3)
    assert C.cohen_kappa(yt, yp) == pytest.approx(1 / 3)


def test_perfect_agreement():
    yt = [1, 0, 1, 0, 1]
    cm = C.confusion(yt, yt)
    assert C.sensitivity(cm) == 1.0 and C.specificity(cm) == 1.0
    assert C.cohen_kappa(yt, yt) == 1.0
    fe = C.fda_errors(cm)
    assert fe["major_error_rate"] == 0.0 and fe["very_major_error_rate"] == 0.0
    assert fe["me_within_fda_band"] and fe["vme_within_fda_band"]


def test_none_values_dropped():
    # positions with a missing call on either side are ignored
    yt = [1, 1, None, 0, 0]
    yp = [1, None, 1, 0, 1]
    cm = C.confusion(yt, yp)   # only indices 0,3,4 survive -> TP1, TN1, FP1
    assert cm["n"] == 3
    assert cm == {"TP": 1, "FP": 1, "TN": 1, "FN": 0, "n": 3}


def test_fda_bands_pass_and_fail():
    # 200 R (2 missed -> VME 1%), 200 S (10 false-R -> ME 5%)
    yt = [1] * 200 + [0] * 200
    yp = [1] * 198 + [0] * 2 + [0] * 190 + [1] * 10
    fe = C.fda_errors(C.confusion(yt, yp))
    assert fe["very_major_error_rate"] == pytest.approx(0.01)
    assert fe["vme_within_fda_band"] is True          # 1% <= 7.5%
    assert fe["major_error_rate"] == pytest.approx(0.05)
    assert fe["me_within_fda_band"] is False           # 5% > 3%


def test_mcnemar_exact_small():
    # A calls 3 resistant that B calls susceptible; no reverse discordance.
    ya = [1, 1, 1, 0, 0]
    yb = [0, 0, 0, 0, 0]
    r = C.mcnemar(ya, yb)
    assert r["b"] == 3 and r["c"] == 0
    assert r["statistic"] is None                      # exact branch (b+c<25)
    assert r["p_value"] == pytest.approx(0.25)         # 2 * P(X<=0 | n=3)


def test_mcnemar_no_discordance():
    r = C.mcnemar([1, 0, 1], [1, 0, 1])
    assert r["b"] == 0 and r["c"] == 0 and r["p_value"] == 1.0


def test_score_pair_bundle():
    yt = [1, 1, 0, 0]
    yp = [1, 0, 0, 0]
    s = C.score_pair(yt, yp)
    assert set(s) >= {"sensitivity", "specificity", "balanced_accuracy",
                      "cohen_kappa", "major_error_rate", "very_major_error_rate", "n"}
    assert s["n"] == 4 and s["sensitivity"] == pytest.approx(0.5)
