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

from lib.registry import antibiotic_to_class
from lib.registry import normalize_antibiotic as _default_normalize

# Accepted testing-standard substrings (case-insensitive)
_ALLOWED_STANDARD_SUBSTRINGS = ("eucast", "clsi")

# Phenotype text -> binary label. The Intermediate / Non-susceptible family is
# handled by the `intermediate_policy` argument to clean_amr_table (default
# 'drop' = binary R/S only).
_PHENOTYPE_MAP = {"resistant": 1, "susceptible": 0}
_INTERMEDIATE_TERMS = ("intermediate", "non-susceptible", "nonsusceptible",
                       "susceptible-dose dependent", "susceptible dose dependent", "sdd")

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
    lab = pd.Series(labels)
    yr = pd.to_numeric(pd.Series(years).reset_index(drop=True), errors="coerce")
    lab = lab.reset_index(drop=True)
    mask = lab.notna()
    lab = lab[mask].astype(int)
    yr = yr[mask]                       # keep labels and years positionally aligned
    if lab.empty:
        return np.nan, False
    counts = lab.value_counts()
    if len(counts) == 1:
        return int(counts.index[0]), False
    # both classes present -> conflict
    n1, n0 = int(counts.get(1, 0)), int(counts.get(0, 0))
    if n1 != n0:
        return (1 if n1 > n0 else 0), True
    # tie -> most recent year. nanargmax (NOT argmax): a partially-missing year
    # column would otherwise make np.argmax return the NaN row and pick the wrong
    # label. Ignore NaN years; drop the cell if no year is available at all.
    if yr.notna().any():
        idx = int(np.nanargmax(yr.to_numpy(dtype=float)))
        return int(lab.to_numpy()[idx]), True
    return np.nan, True  # unresolved -> drop the cell


def clean_amr_table(df, normalize_fn=None, intermediate_policy="drop",
                    strict_antibiotics=False):
    """
    Clean a raw BV-BRC genome_amr frame into a one-row-per-(genome, antibiotic)
    long table with a binary `label`.

    Args:
        df:                 raw frame (API or web-export columns).
        normalize_fn:       antibiotic name normaliser (default: registry).
        intermediate_policy: how to treat the Intermediate / Non-susceptible / SDD
            family — 'drop' (default, binary R/S only), 'resistant' (fold into R,
            CLSI-cautious), or 'susceptible'.
        strict_antibiotics: if True, drop rows whose (normalised) antibiotic is not
            a known registry drug (default False: keep + report, so genuinely new
            drugs are never silently lost).

    Returns:
        (cleaned_long_df, report_dict)
        cleaned_long_df columns: ['genome_id', 'antibiotic', 'label']
        report_dict: row/pair counts at each step + intermediate_policy,
        phenotype_dropped, unknown_antibiotics.
    """
    normalize_fn = normalize_fn or _default_normalize
    if intermediate_policy not in ("drop", "resistant", "susceptible"):
        raise ValueError(
            f"intermediate_policy must be drop|resistant|susceptible, got {intermediate_policy!r}")
    report = {"intermediate_policy": intermediate_policy}
    df = standardise_columns(df)
    report["rows_raw"] = len(df)

    if "genome_id" not in df or "antibiotic" not in df or "resistant_phenotype" not in df:
        raise ValueError("Input must have genome_id, antibiotic, resistant_phenotype columns")

    df = df.drop_duplicates()
    report["rows_dedup"] = len(df)

    # 0) evidence filter — DROP software predictions. Note the polarity: this
    #    excludes "Computational Method" rather than requiring "Laboratory
    #    Method", because BV-BRC leaves `evidence` EMPTY on many real CLSI/EUCAST
    #    measurements. Requiring "laborator" discarded those (measured 2026-07-15:
    #    26 608 such rows for K. pneumoniae alone), and dropping a CLSI-standard
    #    MIC because a neighbouring column was blank is not defensible.
    #
    #    Empty-evidence rows are NOT waved through: step 1 below keeps only
    #    EUCAST/CLSI testing standards, so a row with neither evidence nor a
    #    standard still goes. The two filters compose — this one removes what is
    #    known to be computational, that one demands positive proof of real AST.
    #    (Computational rows carry no testing_standard either, so step 1 would
    #    catch them anyway; this stays explicit rather than relying on that.)
    #    "Computational" hides in TWO columns, not one: besides evidence=
    #    "Computational Method", BV-BRC also has rows with an empty evidence but
    #    laboratory_typing_method="Computational Prediction" (9 814 for
    #    K. pneumoniae, 209 for A. baumannii — measured 2026-07-15). Those slip
    #    past an evidence-only filter. They happen to carry no testing_standard,
    #    so step 1 currently catches every one of them — but that is luck, not
    #    design: loosen step 1 some day and software predictions would quietly
    #    become training labels. Check both columns explicitly.
    _comp = None
    for col in ("evidence", "laboratory_typing_method"):
        if col in df.columns:
            # fillna("") before astype(str): newer pandas can leave NaN as a float
            # after astype(str).str.lower(), which then breaks substring tests.
            hit = df[col].fillna("").astype(str).str.lower().str.contains(
                "computational", na=False)
            _comp = hit if _comp is None else (_comp | hit)
    if _comp is not None:
        df = df[~_comp]
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

    # 2) phenotype filter + label. intermediate_policy folds the Intermediate /
    #    Non-susceptible / SDD family into R (or S), or drops it (default).
    pheno = df["resistant_phenotype"].astype(str).str.strip().str.lower()
    pmap = dict(_PHENOTYPE_MAP)
    if intermediate_policy in ("resistant", "susceptible"):
        tgt = 1 if intermediate_policy == "resistant" else 0
        for term in _INTERMEDIATE_TERMS:
            pmap[term] = tgt
    # transparency (Methods): what phenotype strings get dropped, and how many
    report["phenotype_dropped"] = pheno[~pheno.isin(pmap)].value_counts().to_dict()
    df = df.assign(label=pheno.map(pmap))
    df = df[df["label"].notna()].copy()
    df["label"] = df["label"].astype(int)
    report["rows_after_phenotype"] = len(df)

    # 3) normalise antibiotic names to the registry's canonical spelling.
    df["antibiotic"] = df["antibiotic"].apply(normalize_fn)
    df = df[df["antibiotic"].notna() & (df["antibiotic"].astype(str).str.len() > 0)]
    # data hygiene: surface (and optionally drop) names that are NOT a known
    # registry drug — usually phenotype labels ("fluoroquinolones", "extended
    # spectrum beta lactamase"), never real ML targets.
    unknown = sorted({a for a in df["antibiotic"].unique() if antibiotic_to_class(a) is None})
    report["unknown_antibiotics"] = unknown
    report["n_unknown_antibiotic_names"] = len(unknown)
    if strict_antibiotics and unknown:
        df = df[df["antibiotic"].apply(lambda a: antibiotic_to_class(a) is not None)]

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
