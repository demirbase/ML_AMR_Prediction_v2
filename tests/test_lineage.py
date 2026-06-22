#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for lib.lineage (lineage-aware grouped K-fold, ROADMAP §0.1 M2).

The defining property: a lineage (PopPUNK cluster) never spans the train and
test side of a fold. These tests verify that on synthetic groups/labels, plus
the label alignment, rare-cluster pooling, and summary helpers — no PopPUNK or
container needed.
"""

import numpy as np
import pandas as pd
import pytest

from lib import lineage


def test_group_kfold_no_leakage_and_full_coverage():
    rng = np.random.default_rng(0)
    # 20 lineages, 5 genomes each -> 100 genomes; label correlated with lineage.
    groups = np.repeat(np.arange(20), 5).astype(object)
    y = (np.repeat(rng.integers(0, 2, size=20), 5)).astype(int)

    masks = lineage.group_kfold_masks(y, groups, n_splits=5, stratified=True, seed=42)
    assert len(masks) == 5

    test_union = np.zeros(len(y), dtype=bool)
    for tr, te in masks:
        # train/test partition the rows within a fold
        assert not (tr & te).any()
        assert (tr | te).all()
        # the key guarantee: no lineage on both sides
        assert lineage.no_group_leakage(tr, te, groups)
        test_union |= te
    # every genome is in the test side of exactly one fold
    assert test_union.all()


def test_group_kfold_non_stratified_also_leak_free():
    groups = np.repeat(np.arange(12), 4).astype(object)
    y = np.tile([0, 1], 24).astype(int)
    masks = lineage.group_kfold_masks(y, groups, n_splits=4, stratified=False)
    assert len(masks) == 4
    for tr, te in masks:
        assert lineage.no_group_leakage(tr, te, groups)


def test_group_kfold_rejects_too_many_splits():
    groups = np.array([0, 0, 1, 1], dtype=object)   # only 2 lineages
    y = np.array([0, 1, 0, 1])
    with pytest.raises(ValueError):
        lineage.group_kfold_masks(y, groups, n_splits=5)


def test_load_lineage_aligns_to_genome_order(tmp_path):
    genomes = tmp_path / "genomes_testdrug.csv"
    pd.DataFrame({"Genome ID": ["562.3", "562.1", "562.2"]}).to_csv(genomes, index=False)
    clusters = tmp_path / "poppunk_clusters.csv"
    # PopPUNK-style columns (Taxon, Cluster); deliberately different row order.
    pd.DataFrame({"Taxon": ["562.1", "562.2", "562.3"],
                  "Cluster": [7, 7, 42]}).to_csv(clusters, index=False)

    groups = lineage.load_lineage(genomes, clusters)
    # aligned to genomes_csv order: 562.3 -> 42, 562.1 -> 7, 562.2 -> 7
    assert groups.tolist() == ["42", "7", "7"]


def test_load_lineage_missing_genome_raises(tmp_path):
    genomes = tmp_path / "g.csv"
    pd.DataFrame({"Genome ID": ["a", "b", "c"]}).to_csv(genomes, index=False)
    clusters = tmp_path / "c.csv"
    pd.DataFrame({"Taxon": ["a", "b"], "Cluster": [1, 1]}).to_csv(clusters, index=False)
    with pytest.raises(ValueError):
        lineage.load_lineage(genomes, clusters)
    # ...but allow_missing pools the gap
    groups = lineage.load_lineage(genomes, clusters, allow_missing=True)
    assert groups.tolist() == ["1", "1", "UNCLUSTERED"]


def test_collapse_rare_clusters():
    groups = np.array(["A", "A", "A", "B", "C", "C"], dtype=object)  # A=3, B=1, C=2
    out = lineage.collapse_rare_clusters(groups, min_size=2)
    # B (size 1) pooled into OTHER; A and C kept
    assert out.tolist() == ["A", "A", "A", "OTHER", "C", "C"]
    # no-op when min_size <= 1
    assert lineage.collapse_rare_clusters(groups, min_size=1).tolist() == groups.tolist()


def test_lineage_summary():
    groups = np.array(["A", "A", "A", "B", "C"], dtype=object)
    s = lineage.lineage_summary(groups)
    assert s["n_genomes"] == 5
    assert s["n_clusters"] == 3
    assert s["n_singletons"] == 2          # B and C
    assert s["largest_cluster"] == "A" and s["largest_cluster_size"] == 3
    assert abs(s["largest_cluster_frac"] - 0.6) < 1e-9
