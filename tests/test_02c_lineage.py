#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for 02c_lineage_poppunk.normalize_clusters (PopPUNK name un-mangling).

PopPUNK rewrites '.'→'_' in sample names, so its raw Taxon column won't match the
pipeline's genome ids. These verify the reverse mapping against a synthetic
PopPUNK clusters CSV — no PopPUNK / container needed.
"""

import pandas as pd
import pytest


@pytest.fixture
def mod(load_script):
    return load_script("02c_lineage_poppunk.py")


def test_normalize_unmangles_dots(mod, tmp_path):
    # PopPUNK output: dots already turned into underscores in Taxon.
    raw = tmp_path / "pp_fit_clusters.csv"
    pd.DataFrame({"Taxon": ["562_100036", "562_100039", "562_100004"],
                  "Cluster": [1, 1, 7]}).to_csv(raw, index=False)
    genome_ids = ["562.100004", "562.100036", "562.100039"]  # real ids (with dots)

    out = mod.normalize_clusters(raw, genome_ids)
    assert list(out.columns) == ["Genome ID", "Cluster"]
    mapping = dict(zip(out["Genome ID"], out["Cluster"]))
    assert mapping == {"562.100036": "1", "562.100039": "1", "562.100004": "7"}


def test_normalize_raises_on_unmatched(mod, tmp_path):
    raw = tmp_path / "c.csv"
    pd.DataFrame({"Taxon": ["562_1", "ZZZ_9"], "Cluster": [1, 2]}).to_csv(raw, index=False)
    with pytest.raises(ValueError):
        mod.normalize_clusters(raw, ["562.1"])   # ZZZ_9 has no matching id


def test_normalize_missing_column_raises(mod, tmp_path):
    raw = tmp_path / "c.csv"
    pd.DataFrame({"Sample": ["562_1"], "Cluster": [1]}).to_csv(raw, index=False)
    with pytest.raises(KeyError):
        mod.normalize_clusters(raw, ["562.1"])
