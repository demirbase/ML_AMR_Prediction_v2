#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BV-BRC AMR table cleaning + binary pivoting (used by steps 00a / 00).

Pure data logic — no network, import-light — so it is fully unit-testable on
synthetic frames. The download itself lives in 00a_download_bvbrc.py.

Cleaning rules (agreed):
  1. Keep only EUCAST / CLSI testing standards (case-insensitive; combined forms
     like "EUCAST, CLSI" / "EUCAST and CLSI" kept; NARMS / SFM / BSAC and blank
     dropped).
  2. Keep only Resistant / Susceptible phenotypes -> label 1 / 0
     (Intermediate, Non-susceptible, Susceptible-dose dependent, undefined dropped).
  3. Normalise antibiotic names to canonical registry spelling.
  4. Resolve duplicate / conflicting (genome, antibiotic) cells:
       majority vote -> on a tie, the most recent (max testing_standard_year)
       -> still tied / no year: drop the cell (NaN) and count it.
"""

import numpy as np
import pandas as pd

from lib.registry import normalize_antibiotic as _default_normalize

# Accepted testing-standard substrings (case-insensitive)
_ALLOWED_STANDARD_SUBSTRINGS = ("eucast", "clsi")

# Phenotype text -> binary label
_PHENOTYPE_MAP = {"resistant": 1, "susceptible": 0}

# Map raw column headers (API snake_case OR web-export Title Case) -> canonical
_COLUMN_ALIASES = {
    "genome id": "genome_id",
    "genome_id": "genome_id",
    "genome name": "genome_name",
    "genome_name": "genome_name",
    "antibiotic": "antibiotic",
    "resistant phenotype": "resistant_phenotype",
    "resistant_phenotype": "resistant_phenotype",
    "testing standard": "testing_standard",
    "testing_standard": "testing_standard",
    "testing standard year": "testing_standard_year",
    "testing_standard_year": "testing_standard_year",
    "taxon id": "taxon_id",
    "taxon_id": "taxon_id",
    "evidence": "evidence",
}


def standardise_columns(df):
    """
    Return a copy with column names mapped to the canonical snake_case set.

    Handles both the HTTP API / web-export headers and the BV-BRC CLI headers,
    which prefix fields with their table name (e.g. ``genome_drug.antibiotic``,
    ``genome.genome_id``) — the prefix before the last '.' is stripped first.
    """
    rename = {}
    for col in df.columns:
        key = str(col).strip().lower()
        if "." in key:                       # CLI prefix, e.g. genome_drug.antibiotic
            key = key.rsplit(".", 1)[-1]
        if key in _COLUMN_ALIASES:
            rename[col] = _COLUMN_ALIASES[key]
    return df.rename(columns=rename)


def _resolve_group(labels, years):
    """
    Resolve one (genome, antibiotic) group to a single label.

    Returns (label_or_nan, conflicted_bool). conflicted is True when the group
    contained both classes (regardless of whether the tie-break succeeded).
    """
    labels = pd.Series(labels).dropna().astype(int)
    if labels.empty:
        return np.nan, False
    counts = labels.value_counts()
    if len(counts) == 1:
        return int(counts.index[0]), False
    # both classes present -> conflict
    n1, n0 = int(counts.get(1, 0)), int(counts.get(0, 0))
    if n1 != n0:
        return (1 if n1 > n0 else 0), True
    # tie -> most recent year
    yr = pd.to_numeric(pd.Series(years), errors="coerce")
    if yr.notna().any():
        return int(pd.Series(labels.values)[int(np.argmax(yr.values))]), True
    return np.nan, True  # unresolved -> drop the cell


def clean_amr_table(df, normalize_fn=None):
    """
    Clean a raw BV-BRC genome_amr frame into a one-row-per-(genome, antibiotic)
    long table with a binary `label`.

    Args:
        df:            raw frame (API or web-export columns).
        normalize_fn:  antibiotic name normaliser (default: registry).

    Returns:
        (cleaned_long_df, report_dict)
        cleaned_long_df columns: ['genome_id', 'antibiotic', 'label']
        report_dict: row/pair counts at each step (for the cleaning report).
    """
    normalize_fn = normalize_fn or _default_normalize
    report = {}
    df = standardise_columns(df)
    report["rows_raw"] = len(df)

    if "genome_id" not in df or "antibiotic" not in df or "resistant_phenotype" not in df:
        raise ValueError("Input must have genome_id, antibiotic, resistant_phenotype columns")

    df = df.drop_duplicates()
    report["rows_dedup"] = len(df)

    # 0) evidence filter — keep laboratory-measured phenotypes only (drop
    #    computational predictions). No-op if the column is absent.
    if "evidence" in df.columns:
        # fillna("") before astype(str): newer pandas can leave NaN as a float
        # after astype(str).str.lower(), which then breaks substring tests.
        ev = df["evidence"].fillna("").astype(str).str.lower()
        df = df[ev.str.contains("laborator", na=False)]
        report["rows_after_evidence"] = len(df)

    # 1) testing standard filter (EUCAST / CLSI only). Vectorised + NaN-safe:
    #    empty / missing testing_standard rows are simply dropped (don't match).
    if "testing_standard" in df.columns:
        std = df["testing_standard"].fillna("").astype(str).str.lower()
        keep = std.str.contains("|".join(_ALLOWED_STANDARD_SUBSTRINGS), na=False)
        df = df[keep]
    else:
        report["warning"] = "no testing_standard column — standard filter skipped"
    report["rows_after_standard"] = len(df)

    # 2) phenotype filter + label
    pheno = df["resistant_phenotype"].astype(str).str.strip().str.lower()
    df = df.assign(label=pheno.map(_PHENOTYPE_MAP))
    df = df[df["label"].notna()].copy()
    df["label"] = df["label"].astype(int)
    report["rows_after_phenotype"] = len(df)

    # 3) normalise antibiotic
    df["antibiotic"] = df["antibiotic"].apply(normalize_fn)
    df = df[df["antibiotic"].notna() & (df["antibiotic"].astype(str).str.len() > 0)]

    # 4) conflict resolution per (genome_id, antibiotic)
    if "testing_standard_year" not in df.columns:
        df["testing_standard_year"] = np.nan

    resolved, n_conflict, n_unresolved = [], 0, 0
    for (gid, ab), grp in df.groupby(["genome_id", "antibiotic"], sort=False):
        label, conflicted = _resolve_group(grp["label"].values, grp["testing_standard_year"].values)
        if conflicted:
            n_conflict += 1
        if pd.isna(label):
            n_unresolved += 1
            continue
        resolved.append({"genome_id": str(gid), "antibiotic": ab, "label": int(label)})

    cleaned = pd.DataFrame(resolved, columns=["genome_id", "antibiotic", "label"])
    report["pairs_resolved"] = len(cleaned)
    report["pairs_conflicted"] = n_conflict
    report["pairs_unresolved_dropped"] = n_unresolved
    report["n_genomes"] = cleaned["genome_id"].nunique()
    report["n_antibiotics"] = cleaned["antibiotic"].nunique()
    return cleaned, report


def pivot_binary(cleaned_long):
    """
    Pivot the cleaned long table into a wide binary phenotype matrix.

    Returns a frame with a 'Genome ID' column followed by one column per
    antibiotic (values in {0, 1}; NaN = untested for that genome/antibiotic).
    """
    if cleaned_long.empty:
        return pd.DataFrame(columns=["Genome ID"])
    wide = cleaned_long.pivot_table(
        index="genome_id", columns="antibiotic", values="label", aggfunc="first"
    )
    wide = wide.reindex(sorted(wide.columns), axis=1)
    wide = wide.reset_index().rename(columns={"genome_id": "Genome ID"})
    wide.columns.name = None
    return wide
