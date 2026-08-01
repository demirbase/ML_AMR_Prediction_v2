#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for 07b build_cv_splits — lineage-aware CV with 5-seed fallback (ROADMAP §0.1 M2).

Verifies the scheme selection: lineage StratifiedGroupKFold when PopPUNK labels
exist (no lineage spans train/test), else the legacy 5-seed repeated holdout.
07b imports xgboost; the test skips cleanly where that's unavailable.
"""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def mod(load_script):
    return load_script("07b_feature_stability.py")


def test_fallback_when_no_lineage_file(mod, tmp_path):
    y = np.array([0, 1] * 50)
    splits, method, labels = mod.build_cv_splits(
        y, 100, tmp_path / "missing_genomes.csv", tmp_path / "missing_lineage.csv", 5)
    assert method == "repeated_holdout_5seed"
    assert len(splits) == 5 and labels == mod.SEEDS
    for tr, te in splits:
        assert not (tr & te).any()                  # disjoint
        assert tr.sum() + te.sum() == 100


def test_lineage_group_kfold_when_labels_present(mod, tmp_path):
    n = 100
    genomes = [f"g{i}" for i in range(n)]
    pd.DataFrame({"Genome ID": genomes}).to_csv(tmp_path / "genomes_x.csv", index=False)
    pd.DataFrame({"Genome ID": genomes,
                  "Cluster": np.repeat(np.arange(20), 5)}).to_csv(tmp_path / "lin.csv", index=False)
    y = np.repeat(np.random.default_rng(0).integers(0, 2, 20), 5)

    splits, method, labels = mod.build_cv_splits(
        y, n, tmp_path / "genomes_x.csv", tmp_path / "lin.csv", 5)
    assert method == "lineage_group_kfold_5fold"
    assert len(splits) == 5 and labels == [0, 1, 2, 3, 4]

    from lib.lineage import load_lineage, no_group_leakage
    groups = load_lineage(tmp_path / "genomes_x.csv", tmp_path / "lin.csv")
    test_union = np.zeros(n, dtype=bool)
    for tr, te in splits:
        assert no_group_leakage(tr, te, groups)     # the key guarantee
        test_union |= te
    assert test_union.all()                         # every genome tested once


def test_lineage_falls_back_when_too_few_clusters(mod, tmp_path):
    n = 100
    genomes = [f"g{i}" for i in range(n)]
    pd.DataFrame({"Genome ID": genomes}).to_csv(tmp_path / "g.csv", index=False)
    # only 3 lineages -> cannot do 5 folds -> fall back to 5-seed holdout
    pd.DataFrame({"Genome ID": genomes,
                  "Cluster": np.repeat(np.arange(3), [40, 30, 30])}).to_csv(tmp_path / "l.csv", index=False)
    y = np.array([0, 1] * 50)
    splits, method, labels = mod.build_cv_splits(y, n, tmp_path / "g.csv", tmp_path / "l.csv", 5)
    assert method == "repeated_holdout_5seed"


def test_random_cv_mode_is_lineage_blind_and_suffixes_outputs(load_script, tmp_path, monkeypatch):
    """AMR_CV_MODE=random must ignore lineage labels AND write to suffixed paths.

    The suffix is the safety property: this mode exists to produce the inflated
    comparison baseline, so if it ever wrote to the canonical filenames it would
    overwrite the lineage-CV numbers the KB is populated from.
    """
    monkeypatch.setenv("AMR_CV_MODE", "random")
    mod = load_script("07b_feature_stability.py")
    assert mod.RANDOM_CV is True
    assert mod.OUT_SUFFIX == "_randomcv"

    n = 100
    genomes = [f"g{i}" for i in range(n)]
    pd.DataFrame({"Genome ID": genomes}).to_csv(tmp_path / "genomes_x.csv", index=False)
    # One lineage per half: a lineage-aware split would never mix them; random does.
    pd.DataFrame({"Genome ID": genomes,
                  "Cluster": ["A"] * (n // 2) + ["B"] * (n // 2)}).to_csv(
        tmp_path / "poppunk_clusters.csv", index=False)

    y = np.array([0, 1] * (n // 2))
    splits, method, labels = mod.build_cv_splits(
        y, n, tmp_path / "genomes_x.csv", tmp_path / "poppunk_clusters.csv", 5)

    assert method == "random_stratified_kfold_5fold"
    assert len(splits) == 5 and labels == list(range(5))
    seen_test = np.zeros(n, dtype=bool)
    for tr, te in splits:
        assert not (tr & te).any()                 # disjoint
        assert tr.sum() + te.sum() == n            # a partition
        seen_test |= te
    assert seen_test.all()                         # every sample tested exactly once
