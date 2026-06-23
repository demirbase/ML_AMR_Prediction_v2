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
    t, k = m.DEFAULT_TIERS, 21
    # identity + coverage are primary; E-value is a loose secondary gate
    assert m.classify_confidence(99.0, 1e-5, 21, k, t) == "confirmed"  # full-length near-exact
    assert m.classify_confidence(92.0, 0.5, 18, k, t) == "candidate"   # high id, partial coverage
    assert m.classify_confidence(99.0, 1e-5, 18, k, t) == "candidate"  # 99% but cov<0.95 → demoted
    assert m.classify_confidence(85.0, 5.0, 14, k, t) == "weak"        # kept + flagged
    assert m.classify_confidence(70.0, 5.0, 14, k, t) == "none"        # identity below weak floor
    assert m.classify_confidence(99.0, 1e-5, 10, k, t) == "none"       # coverage below weak floor


# ---------------------------------------------------------------------------
# 09_biological_summary.composite_score / build_kb_candidates (M7 / H4)
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_composite_score(load_script):
    m = load_script("09_biological_summary.py")
    # stability × log10(1/E) × (identity/100): 0.8 × 5 × 1.0 = 4.0
    assert abs(m.composite_score(0.8, 100.0, 1e-5) - 4.0) < 1e-9
    # weak hit (E>1) → log10(1/E) negative → clamped to 0
    assert m.composite_score(0.8, 100.0, 1.5) == 0.0
    # missing stability → NaN
    assert m.composite_score(float("nan"), 100.0, 1e-5) != m.composite_score(float("nan"), 100.0, 1e-5)


@pytest.mark.unit
def test_kb_recovery_and_novel(load_script):
    import pandas as pd
    m = load_script("09_biological_summary.py")
    feats = pd.DataFrame([
        {"Rank": 1, "Gain_Score": 9.0, "Feature_ID": "f1", "Kmer_Sequence": "AAA",
         "in_gain_topN": True, "selection_frequency": 0.8, "stable": True},
        {"Rank": 2, "Gain_Score": 5.0, "Feature_ID": "f2", "Kmer_Sequence": "CCC",
         "in_gain_topN": True, "selection_frequency": 0.6, "stable": True},
        {"Rank": 3, "Gain_Score": 4.0, "Feature_ID": "f3", "Kmer_Sequence": "GGG",
         "in_gain_topN": False, "selection_frequency": 0.8, "stable": True},
    ])
    card = pd.DataFrame([
        {"qseqid": "Rank_1|Score_9.0000|Feature_f1", "pident": 100.0, "evalue": 1e-5,
         "Gene_Match": "TEM-1", "Confidence": "confirmed"},
        {"qseqid": "Rank_2|Score_5.0000|Feature_f2", "pident": 92.0, "evalue": 1.5,
         "Gene_Match": "gyrA", "Confidence": "candidate"},
    ])  # f3 has no CARD hit → novel
    kb, met = m.build_kb_candidates(feats, card, 0.6)
    assert met["n_stable"] == 3
    assert abs(met["known_mechanism_recovery_rate"] - 1 / 3) < 1e-9   # 1 confirmed / 3 stable
    assert met["H2_pass"] is False                                    # < 0.40
    assert abs(met["novel_candidate_fraction"] - 1 / 3) < 1e-9        # 1 novel / 3 stable
    assert met["tier_counts_all"]["confirmed"] == 1


# ---------------------------------------------------------------------------
# 10_kmer_background_frequency.compute_kmer_stats (discriminativeness)
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_kmer_background_frequency(load_script):
    m = load_script("10_kmer_background_frequency.py")
    # Strong resistance marker: present in most R, few S -> discriminative
    s = m.compute_kmer_stats(present_r=90, n_r=100, present_s=10, n_s=100)
    assert abs(s["prevalence_resistant"] - 0.9) < 1e-9
    assert abs(s["prevalence_susceptible"] - 0.1) < 1e-9
    assert s["enriched_in"] == "resistant"
    assert s["fisher_p"] < 0.05
    assert s["discriminative"] is True
    # Ubiquitous: present everywhere, no R/S difference -> NOT discriminative
    u = m.compute_kmer_stats(present_r=95, n_r=100, present_s=95, n_s=100)
    assert u["discriminative"] is False
    assert u["enriched_in"] == "equal"


@pytest.mark.unit
def test_benjamini_hochberg(load_script):
    import numpy as np
    m = load_script("10_kmer_background_frequency.py")
    # Known BH example: p=[0.01,0.02,0.03,0.04,0.05], m=5
    # q_i = min over k>=i of (p_k * 5 / k): largest stays 0.05, others adjust up.
    q = m.benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05])
    assert np.all(np.isfinite(q))
    assert np.all(q >= np.array([0.01, 0.02, 0.03, 0.04, 0.05]) - 1e-12)  # q >= p
    assert np.all(np.diff(q) >= -1e-12)                                   # monotone non-decreasing
    assert q[-1] == pytest.approx(0.05, abs=1e-9)
    assert q[0] == pytest.approx(0.05, abs=1e-9)                          # 0.01*5/1
    # NaN passes through; finite values still corrected
    q2 = m.benjamini_hochberg([0.01, np.nan, 0.04])
    assert np.isnan(q2[1]) and np.all(np.isfinite(q2[[0, 2]]))
    # all q-values clipped to <= 1
    assert np.all(m.benjamini_hochberg([0.9, 0.95, 0.99]) <= 1.0)


# ---------------------------------------------------------------------------
# 11_variant_snp_check — SNP parsing / codon mapping / allele classification
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_variant_snp_helpers(load_script):
    m = load_script("11_variant_snp_check.py")
    assert m.parse_snp_token("S83L") == ("S", 83, "L")
    assert m.parse_snp_token("not-a-snp") is None
    assert m.protein_pos_to_codon_nt(83) == (247, 249)   # (83-1)*3+1
    assert m.translate_codon("CTG") == "L"
    assert m.aro_from_sseqid("gb|AF469609.1|+|566-3353|ARO:3003297|gyrA") == "3003297"

    # Plus strand, ungapped: subject CDS 245..253, codon 247-249.
    # k-mer carries CTG (Leu) at the codon -> resistant for S83L.
    qc = m.query_codon_from_alignment(245, 253, "plus",
                                      "AG" + "CTG" + "AAGT", "AG" + "AGC" + "AAGT",
                                      (247, 248, 249))
    assert qc == "CTG"
    assert m.classify_allele(qc, "S", "L") == "resistant_allele"
    # Same codon but wildtype AGC (Ser)
    qc_wt = m.query_codon_from_alignment(245, 253, "plus",
                                         "AG" + "AGC" + "AAGT", "AG" + "AGC" + "AAGT",
                                         (247, 248, 249))
    assert m.classify_allele(qc_wt, "S", "L") == "wildtype"
    # Minus strand: CDS-sense base = complement of aligned query char
    qcm = m.query_codon_from_alignment(253, 245, "minus",
                                       "AAAA" + "C" + "A" + "G" + "CT",
                                       "TTAA" + "C" + "A" + "G" + "CT",
                                       (247, 248, 249))
    assert qcm == "CTG"
    # Codon not fully covered -> None -> ambiguous
    assert m.query_codon_from_alignment(248, 253, "plus", "CTGAAG", "AGCAAG", (247, 248, 249)) is None
    assert m.classify_allele(None, "S", "L") == "ambiguous"


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
