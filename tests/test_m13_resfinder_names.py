#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ResFinder antibiotic-name matching (M13, step 16).

ResFinder writes "amoxicillin+clavulanic acid"; this project keys the same agent
as "amoxicillin_clavulanic_acid". The original literal comparison matched
neither, so every combination agent silently produced no ResFinder call — 7 of
the 45 models scored against an empty ResFinder result. See parse_resfinder.
"""

import importlib.util
from pathlib import Path

import pytest

SPEC = importlib.util.spec_from_file_location(
    "m13", Path(__file__).resolve().parents[1] / "scripts" / "16_external_concordance.py")
m13 = importlib.util.module_from_spec(SPEC)
try:
    SPEC.loader.exec_module(m13)
except SystemExit:                                   # module guards on config
    pass

HEADER = ("# ResFinder phenotype results.\n"
          "# Antimicrobial\tClass\tWGS-predicted phenotype\tMatch\n")


def _table(tmp_path, rows):
    p = tmp_path / "pheno_table_escherichia_coli.txt"
    p.write_text(HEADER + "".join(f"{a}\tcls\t{v}\t1\n" for a, v in rows), encoding="utf-8")
    return p


def test_separator_style_does_not_block_a_match(tmp_path):
    calls = m13.parse_resfinder(_table(tmp_path, [
        ("amoxicillin+clavulanic acid", "Resistant"),
        ("piperacillin+tazobactam", "No resistance"),
    ]))
    assert calls["amoxicillin_clavulanic_acid"] == 1
    assert calls["piperacillin_tazobactam"] == 0


def test_component_reported_combination_is_assembled(tmp_path):
    """ResFinder publishes no row for TMP-SMX, only its two components."""
    r = m13.parse_resfinder(_table(tmp_path, [
        ("trimethoprim", "Resistant"), ("sulfamethoxazole", "No resistance")]))
    assert r["trimethoprim_sulfamethoxazole"] == 1        # either component suffices
    s = m13.parse_resfinder(_table(tmp_path, [
        ("trimethoprim", "No resistance"), ("sulfamethoxazole", "No resistance")]))
    assert s["trimethoprim_sulfamethoxazole"] == 0


def test_inhibitor_combination_never_matches_the_bare_drug(tmp_path):
    """ampicillin_sulbactam must not be read off ResFinder's plain 'ampicillin':
    the inhibitor changes the phenotype, so a missing call is the honest result."""
    calls = m13.parse_resfinder(_table(tmp_path, [("ampicillin", "Resistant")]))
    assert calls["ampicillin"] == 1
    assert "ampicillin_sulbactam" not in calls
    assert "oxacillin" not in calls


@pytest.mark.parametrize("name", ["Ciprofloxacin", "CIPROFLOXACIN", " ciprofloxacin "])
def test_case_and_whitespace_are_ignored(tmp_path, name):
    assert m13.parse_resfinder(_table(tmp_path, [(name, "Resistant")]))["ciprofloxacin"] == 1


def test_phenotype_columns_match_across_separator_styles(tmp_path):
    """BV-BRC writes "trimethoprim/sulfamethoxazole"; the project keys it with
    underscores. A literal lookup returned None, so every combination agent was
    scored against no phenotype at all and reported n=0."""
    md = tmp_path / "amr_phenotypes.csv"
    md.write_text("Genome ID,ampicillin,trimethoprim/sulfamethoxazole,"
                  "amoxicillin/clavulanic acid\n"
                  "562.1,1,0,1\n562.2,0,1,\n", encoding="utf-8")
    pheno = m13.load_phenotype(md, ["ampicillin", "trimethoprim_sulfamethoxazole",
                                    "amoxicillin_clavulanic_acid", "oxacillin"])
    assert pheno["562.1"]["trimethoprim_sulfamethoxazole"] == 0
    assert pheno["562.2"]["trimethoprim_sulfamethoxazole"] == 1
    assert pheno["562.1"]["amoxicillin_clavulanic_acid"] == 1
    assert pheno["562.2"]["amoxicillin_clavulanic_acid"] is None   # blank stays unknown
    assert pheno["562.1"]["oxacillin"] is None                     # no such column
