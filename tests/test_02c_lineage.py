#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for 02c_lineage_poppunk.normalize_clusters (PopPUNK name un-mangling).

PopPUNK rewrites '.'→'_' in sample names, so its raw Taxon column won't match the
pipeline's genome ids. These verify the reverse mapping against a synthetic
PopPUNK clusters CSV — no PopPUNK / container needed.
"""

from pathlib import Path

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent


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


# ---- lineage params: global E3 defaults + per-organism registry override ----
# The sketch/k settings decide the clustering, and the clustering IS the CV
# grouping — so these are scientific parameters and must be explicit, not left
# to whatever the installed PopPUNK defaults to.

def test_params_apply_e3_range_to_every_organism(mod):
    """One k-range across the panel: mixing PopPUNK defaults for some species and
    E3 settings for others would undermine the cross-organism comparison."""
    from lib.config import load_config
    cfg = load_config()
    for org in ("ecoli", "kpneumoniae", "staphylococcus_aureus", "acinetobacter_baumannii"):
        p = mod.lineage_params(org, cfg)
        assert (p["min_k"], p["max_k"], p["k_step"]) == (15, 35, 2)


def test_staph_overrides_sketch_size(mod):
    """S. aureus is the panel's most clonal species — E3 §4 requires a 10^5 sketch
    to separate strains the 10^4 default would collapse."""
    from lib.config import load_config
    cfg = load_config()
    assert mod.lineage_params("staphylococcus_aureus", cfg)["sketch_size"] == 100_000
    assert mod.lineage_params("ecoli", cfg)["sketch_size"] == 10_000   # global default


def test_unknown_organism_falls_back_to_globals(mod):
    from lib.config import load_config
    cfg = load_config()
    assert mod.lineage_params("not_an_organism", cfg)["sketch_size"] == 10_000


def test_sketch_args_are_passed_explicitly(mod):
    args = mod._sketch_args({"min_k": 15, "max_k": 35, "k_step": 2, "sketch_size": 100000})
    assert args == "--min-k 15 --max-k 35 --k-step 2 --sketch-size 100000"


def test_qc_args_and_length_range_pair(mod):
    base = {"max_a_dist": 0.5, "length_sigma": 5, "prop_n": 0.1, "length_range": None}
    assert mod._qc_args(base) == "--max-a-dist 0.5 --length-sigma 5 --prop-n 0.1"
    # PopPUNK's --length-range takes TWO values (lower upper)
    assert mod._qc_args({**base, "length_range": [3_500_000, 4_200_000]}).endswith(
        "--length-range 3500000 4200000")


def test_only_acinetobacter_overrides_refine(mod):
    """Refinement is NOT a universal improvement — measured 2026-07-16 it HARMED
    K. pneumoniae (22.3% -> 58.6% largest lineage, merging its high-risk clones)
    while it rescued A. baumannii (76.3% -> 52.1%). The global default stays off;
    only A. baumannii overrides. Pins the rule so a future 'let's just enable
    refine everywhere' cannot land silently."""
    from lib.config import load_config
    cfg = load_config()
    assert mod.lineage_params("acinetobacter_baumannii", cfg)["refine"] is True
    for org in ("ecoli", "kpneumoniae", "staphylococcus_aureus",
                "pseudomonas_aeruginosa", "enterococcus_faecium"):
        assert mod.lineage_params(org, cfg)["refine"] is False, org


def test_registry_lineage_override_is_not_shadowed_by_argparse(mod):
    """--model/--refine default to None so lineage_params (config + REGISTRY) wins.
    They used to default to config's value, which silently shadowed the registry —
    A. baumannii's refine override would have been ignored."""
    import argparse
    src = (PROJECT_ROOT / "scripts" / "02c_lineage_poppunk.py").read_text()
    assert 'ap.add_argument("--model", default=None' in src
    assert 'default=lin_cfg.get("model"' not in src
