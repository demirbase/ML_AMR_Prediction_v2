#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Binary concordance / clinical-error metrics for external validation (M13).

Pure, dependency-light functions (scipy only for p-values, with a stdlib
fallback) shared by ``16_external_concordance.py`` to compare, on the same set of
genomes, any two of: our k-mer model's calls, AMRFinderPlus/ResFinder genotypic
calls, and the EUCAST/CLSI phenotype (ground truth).

Convention: label 1 = **resistant** (the positive/clinically-important class),
0 = susceptible. ``y_true`` is the reference (usually the phenotype); ``y_pred``
is the predictor being graded (a tool or the model).

FDA/CLSI genotype-vs-phenotype error bands (ROADMAP §0.1): Major Error (ME) =
predict R when truly S (false resistant), band ME ≤ 3%; Very Major Error (VME) =
predict S when truly R (missed resistance — the dangerous one), band ≤ 1.5–7.5%.
"""

from __future__ import annotations

VME_BAND = 0.075   # FDA upper bound for Very Major Error rate (1.5–7.5%)
ME_BAND = 0.03     # FDA upper bound for Major Error rate


def _pairs(y_true, y_pred):
    """Zip to aligned int pairs, dropping any position where either side is None
    (a genome the tool/phenotype has no call for)."""
    out = []
    for t, p in zip(y_true, y_pred):
        if t is None or p is None:
            continue
        out.append((int(t), int(p)))
    return out


def confusion(y_true, y_pred):
    """2×2 counts with 1=resistant as positive. Returns dict TP/FP/TN/FN/n."""
    tp = fp = tn = fn = 0
    for t, p in _pairs(y_true, y_pred):
        if t == 1 and p == 1:
            tp += 1
        elif t == 0 and p == 1:
            fp += 1
        elif t == 0 and p == 0:
            tn += 1
        else:
            fn += 1
    return {"TP": tp, "FP": fp, "TN": tn, "FN": fn, "n": tp + fp + tn + fn}


def _safe_div(a, b):
    return (a / b) if b else None


def sensitivity(cm):
    """TP / (TP+FN) — recall of resistance."""
    return _safe_div(cm["TP"], cm["TP"] + cm["FN"])


def specificity(cm):
    """TN / (TN+FP)."""
    return _safe_div(cm["TN"], cm["TN"] + cm["FP"])


def balanced_accuracy(cm):
    se, sp = sensitivity(cm), specificity(cm)
    return (se + sp) / 2 if (se is not None and sp is not None) else None


def fda_errors(cm):
    """Major-error and very-major-error rates + FDA-band pass flags.

    ME = FP / (all truly susceptible); VME = FN / (all truly resistant)."""
    me = _safe_div(cm["FP"], cm["TN"] + cm["FP"])
    vme = _safe_div(cm["FN"], cm["TP"] + cm["FN"])
    return {
        "major_error_rate": me,
        "very_major_error_rate": vme,
        "me_within_fda_band": (me is not None and me <= ME_BAND),
        "vme_within_fda_band": (vme is not None and vme <= VME_BAND),
    }


def cohen_kappa(y_true, y_pred):
    """Cohen's κ for two binary raters (chance-corrected agreement)."""
    cm = confusion(y_true, y_pred)
    n = cm["n"]
    if not n:
        return None
    po = (cm["TP"] + cm["TN"]) / n
    p_pos = (cm["TP"] + cm["FN"]) / n * (cm["TP"] + cm["FP"]) / n
    p_neg = (cm["TN"] + cm["FP"]) / n * (cm["TN"] + cm["FN"]) / n
    pe = p_pos + p_neg
    return (po - pe) / (1 - pe) if pe != 1 else 1.0


def mcnemar(y_a, y_b):
    """McNemar's test of marginal homogeneity between two predictors on the same
    samples. Discordant cells b=(A=R,B=S), c=(A=S,B=R). Uses an exact two-sided
    binomial p when b+c is small (<25), else the continuity-corrected χ² approx.
    Returns {b, c, statistic, p_value} (statistic is None for the exact branch)."""
    b = c = 0
    for a, d in _pairs(y_a, y_b):
        if a == 1 and d == 0:
            b += 1
        elif a == 0 and d == 1:
            c += 1
    nd = b + c
    if nd == 0:
        return {"b": 0, "c": 0, "statistic": None, "p_value": 1.0}
    if nd < 25:
        p = _binom_two_sided(min(b, c), nd)
        return {"b": b, "c": c, "statistic": None, "p_value": p}
    stat = (abs(b - c) - 1) ** 2 / nd
    try:
        from scipy.stats import chi2
        p = float(chi2.sf(stat, 1))
    except Exception:
        p = _chi2_sf_1df(stat)
    return {"b": b, "c": c, "statistic": stat, "p_value": p}


def _binom_two_sided(k, n, prob=0.5):
    """Exact two-sided binomial p for k successes in n at p=0.5 (McNemar exact)."""
    try:
        from scipy.stats import binomtest
        return float(binomtest(k, n, prob, alternative="two-sided").pvalue)
    except Exception:
        from math import comb
        # two-sided at p=0.5 is symmetric: 2 * P(X <= k), capped at 1.
        tail = sum(comb(n, i) for i in range(0, k + 1)) / (2 ** n)
        return min(1.0, 2 * tail)


def _chi2_sf_1df(x):
    """Survival function of χ²₁ via erfc (stdlib fallback, no scipy)."""
    import math
    if x <= 0:
        return 1.0
    return math.erfc(math.sqrt(x / 2.0))


def score_pair(y_true, y_pred):
    """Full metric bundle grading ``y_pred`` against reference ``y_true``."""
    cm = confusion(y_true, y_pred)
    out = {
        "n": cm["n"], "confusion": cm,
        "sensitivity": sensitivity(cm), "specificity": specificity(cm),
        "balanced_accuracy": balanced_accuracy(cm),
        "cohen_kappa": cohen_kappa(y_true, y_pred),
    }
    out.update(fda_errors(cm))
    return out
