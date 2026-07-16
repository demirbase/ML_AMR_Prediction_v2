#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for BV-BRC AMR cleaning + name normalisation (step 00)."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import registry                      # noqa: E402
from lib.bvbrc import clean_amr_table, pivot_binary  # noqa: E402


# ---------------------------------------------------------------------------
# antibiotic name normalisation (registry)
# ---------------------------------------------------------------------------
@pytest.mark.unit
def test_normalize_aliases_and_typos():
    n = registry.normalize_antibiotic
    # ampicillin/sulbactam canonicalises to the UNDERSCORE form (slash-safe for
    # file paths), like the other combo drugs. This test used to assert the slash
    # form — it was encoding the path-breaking bug fixed 2026-07-16.
    assert n("amipicillin_sulbactam") == "ampicillin_sulbactam"
    assert n("ampicillin/sulbactam") == "ampicillin_sulbactam"     # BV-BRC spelling -> canonical
    assert n("rifampicin") == "rifampin"
    assert n("cefalotin") == "cephalothin"
    assert n("tigecyklin") == "tigecycline"
    assert n("amoxicillin_clavulanat") == "amoxicillin_clavulanic_acid"


@pytest.mark.unit
def test_normalize_cotrimoxazole_merged():
    n = registry.normalize_antibiotic
    # co-trimoxazole and its variants all fold to one canonical (path-safe underscore)
    assert n("co-trimoxazole") == "trimethoprim_sulfamethoxazole"
    assert n("sulfamethoxazole/trimethoprim") == "trimethoprim_sulfamethoxazole"


@pytest.mark.unit
def test_normalize_case_trim_and_passthrough():
    n = registry.normalize_antibiotic
    assert n("  GENTAMICIN ") == "gentamicin"
    assert n("some_unlisted_drug") == "some_unlisted_drug"   # unknown kept, not dropped
    assert n(None) is None


# ---------------------------------------------------------------------------
# cleaning
# ---------------------------------------------------------------------------
def _df(rows):
    return pd.DataFrame(rows, columns=["Genome ID", "Antibiotic", "Resistant Phenotype",
                                       "Testing standard", "Testing standard year"])


@pytest.mark.unit
def test_clean_standard_and_phenotype_filters():
    df = _df([
        ("562.1", "ampicillin", "Resistant", "EUCAST", 2020),
        ("562.1", "tetracycline", "Resistant", "NARMS", 2021),       # standard dropped
        ("562.1", "ciprofloxacin", "Intermediate", "EUCAST", 2021),  # phenotype dropped
        ("562.1", "gentamicin", "Susceptible", "EUCAST, CLSI", 2019),  # combined std kept
    ])
    cleaned, rep = clean_amr_table(df)
    got = {(r.genome_id, r.antibiotic): r.label for r in cleaned.itertuples()}
    assert got == {("562.1", "ampicillin"): 1, ("562.1", "gentamicin"): 0}
    assert rep["rows_after_standard"] == 3
    assert rep["rows_after_phenotype"] == 2


@pytest.mark.unit
def test_clean_conflict_majority():
    df = _df([
        ("562.3", "ampicillin", "Resistant", "EUCAST", 2019),
        ("562.3", "ampicillin", "Resistant", "CLSI", 2020),
        ("562.3", "ampicillin", "Susceptible", "EUCAST", 2018),
    ])
    cleaned, rep = clean_amr_table(df)
    assert cleaned.iloc[0]["label"] == 1          # 2R vs 1S -> R
    assert rep["pairs_conflicted"] == 1


@pytest.mark.unit
def test_clean_conflict_tie_newest_year():
    df = _df([
        ("562.4", "gentamicin", "Resistant", "EUCAST", 2015),
        ("562.4", "gentamicin", "Susceptible", "CLSI", 2022),
    ])
    cleaned, _ = clean_amr_table(df)
    assert cleaned.iloc[0]["label"] == 0          # tie -> 2022 Susceptible


@pytest.mark.unit
def test_clean_conflict_tie_no_year_dropped():
    df = _df([
        ("562.5", "cefotaxime", "Resistant", "EUCAST", np.nan),
        ("562.5", "cefotaxime", "Susceptible", "CLSI", np.nan),
    ])
    cleaned, rep = clean_amr_table(df)
    assert cleaned.empty
    assert rep["pairs_unresolved_dropped"] == 1


@pytest.mark.unit
def test_clean_conflict_tie_partial_nan_year():
    # tie with a partially-missing year: the row that HAS a year must win
    # (regression for the np.argmax->np.nanargmax fix; argmax would pick the NaN row).
    df = _df([
        ("562.6", "gentamicin", "Resistant", "EUCAST", 2015),
        ("562.6", "gentamicin", "Susceptible", "CLSI", np.nan),
    ])
    cleaned, _ = clean_amr_table(df)
    assert cleaned.iloc[0]["label"] == 1          # 2015 Resistant wins, not the NaN-year S
    # and the reverse orientation
    df2 = _df([
        ("562.7", "gentamicin", "Resistant", "EUCAST", np.nan),
        ("562.7", "gentamicin", "Susceptible", "CLSI", 2018),
    ])
    cleaned2, _ = clean_amr_table(df2)
    assert cleaned2.iloc[0]["label"] == 0          # 2018 Susceptible wins


@pytest.mark.unit
def test_clean_unknown_antibiotic_reported_and_optionally_dropped():
    df = _df([
        ("562.8", "ampicillin", "Resistant", "EUCAST", 2020),
        ("562.8", "fluoroquinolones", "Resistant", "EUCAST", 2020),  # class label, not a drug
    ])
    cleaned, rep = clean_amr_table(df)                 # default: keep + report
    assert "fluoroquinolones" in rep["unknown_antibiotics"]
    assert set(cleaned["antibiotic"]) == {"ampicillin", "fluoroquinolones"}
    cleaned_s, rep_s = clean_amr_table(df, strict_antibiotics=True)   # strict: drop it
    assert set(cleaned_s["antibiotic"]) == {"ampicillin"}
    assert rep_s["n_unknown_antibiotic_names"] == 1


@pytest.mark.unit
def test_clean_intermediate_policy():
    df = _df([
        ("562.9", "ciprofloxacin", "Intermediate", "EUCAST", 2021),
    ])
    # default drop -> nothing survives
    cleaned_drop, rep_d = clean_amr_table(df)
    assert cleaned_drop.empty
    assert rep_d["phenotype_dropped"].get("intermediate") == 1
    # policy 'resistant' -> Intermediate folds into R (label 1)
    cleaned_r, _ = clean_amr_table(df, intermediate_policy="resistant")
    assert cleaned_r.iloc[0]["label"] == 1


@pytest.mark.unit
def test_clean_cli_prefixed_headers_and_evidence_filter():
    # Simulates BV-BRC CLI output: table-prefixed headers + an evidence column.
    df = pd.DataFrame({
        "genome.genome_id": ["562.1", "562.1", "562.2"],
        "genome_drug.antibiotic": ["ampicillin", "gentamicin", "rifampicin"],
        "genome_drug.resistant_phenotype": ["Resistant", "Susceptible", "Resistant"],
        "genome_drug.testing_standard": ["EUCAST", "CLSI", "EUCAST"],
        "genome_drug.testing_standard_year": [2020, 2019, 2021],
        "genome_drug.evidence": ["Laboratory Method", "Laboratory Method", "Computational Method"],
    })
    cleaned, rep = clean_amr_table(df)
    assert rep["rows_after_evidence"] == 2           # computational row dropped
    assert set(cleaned["antibiotic"]) == {"ampicillin", "gentamicin"}


@pytest.mark.unit
def test_pivot_binary_shape_and_nan():
    cleaned = pd.DataFrame({
        "genome_id": ["562.1", "562.1", "562.2"],
        "antibiotic": ["ampicillin", "gentamicin", "ampicillin"],
        "label": [1, 0, 1],
    })
    wide = pivot_binary(cleaned)
    assert list(wide.columns) == ["Genome ID", "ampicillin", "gentamicin"]
    row2 = wide[wide["Genome ID"] == "562.2"].iloc[0]
    assert row2["ampicillin"] == 1 and pd.isna(row2["gentamicin"])


# ---- evidence polarity: drop computational, don't demand "Laboratory Method" --
# BV-BRC leaves `evidence` EMPTY on many real CLSI/EUCAST measurements. The old
# filter required "laborator" and silently discarded them (26 608 rows for
# K. pneumoniae alone). The rescue must not, however, wave through rows with no
# proof of real AST — the testing_standard filter still has to earn its keep.

def _frame(rows):
    import pandas as pd
    return pd.DataFrame(rows, columns=["genome_id", "antibiotic", "resistant_phenotype",
                                       "evidence", "testing_standard"])


def test_empty_evidence_with_clsi_is_kept():
    from lib.bvbrc import clean_amr_table
    df = _frame([["1.1", "meropenem", "Resistant", None, "CLSI"],
                 ["1.2", "meropenem", "Susceptible", "Laboratory Method", "EUCAST"]])
    out, _ = clean_amr_table(df)
    assert set(out["genome_id"]) == {"1.1", "1.2"}   # the empty-evidence CLSI row survives


def test_computational_rows_are_dropped_even_with_a_standard():
    from lib.bvbrc import clean_amr_table
    df = _frame([["2.1", "meropenem", "Resistant", "Computational Method", "CLSI"],
                 ["2.2", "meropenem", "Resistant", "Laboratory Method", "CLSI"]])
    out, _ = clean_amr_table(df)
    assert set(out["genome_id"]) == {"2.2"}


def test_empty_evidence_without_a_standard_is_still_dropped():
    """The rescue must not become a hole: no evidence AND no standard = no proof."""
    from lib.bvbrc import clean_amr_table
    df = _frame([["3.1", "meropenem", "Resistant", None, None],
                 ["3.2", "meropenem", "Resistant", None, "CLSI"]])
    out, _ = clean_amr_table(df)
    assert set(out["genome_id"]) == {"3.2"}


def test_computational_hidden_in_typing_method_is_dropped():
    """"Computational" lives in TWO columns. Rows with an empty evidence but
    laboratory_typing_method="Computational Prediction" (9814 of them for
    K. pneumoniae) slip past an evidence-only filter. They currently carry no
    testing_standard so step 1 catches them anyway — this pins the explicit
    defence so that stays true if step 1 is ever loosened."""
    import pandas as pd
    from lib.bvbrc import clean_amr_table
    df = pd.DataFrame(
        [["4.1", "meropenem", "Resistant", None, "Computational Prediction", "CLSI"],
         ["4.2", "meropenem", "Resistant", None, "Broth dilution", "CLSI"]],
        columns=["genome_id", "antibiotic", "resistant_phenotype", "evidence",
                 "laboratory_typing_method", "testing_standard"])
    out, _ = clean_amr_table(df)
    assert set(out["genome_id"]) == {"4.2"}   # the predicted row must not become a label
