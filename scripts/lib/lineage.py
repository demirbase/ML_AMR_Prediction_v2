#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lineage-aware cross-validation helpers (ROADMAP §0.1 M2).

Random / chunk-level train-test splits let genomes from the SAME lineage (PopPUNK
cluster — near-genetic-twins, e.g. E. coli ST131) fall on both sides, so the model
memorises the lineage instead of the resistance mechanism and the reported AUC is
inflated ~20-30 % (Yu/Barquist 2024). This module replaces that with **grouped
K-fold**: an entire lineage goes to either train OR test, never both — so the
held-out AUC honestly answers "can the model predict resistance in a lineage it
has never seen?".

It is feature-agnostic and out-of-core friendly: it emits **sample-level boolean
masks** (one per fold) that plug straight into ``lib.xgb_data`` (``row_mask=``),
the same mechanism 07b already uses. Lineage labels come from PopPUNK (organism-
level, antibiotic-independent — computed once, reused for every antibiotic).

PopPUNK clusters are used AS-IS: every cluster, including singletons, stays its
own group. Pooling rare clusters into one bucket was considered and rejected —
unrelated singletons would become a single group and therefore all land in the
same fold, skewing it, whereas separate singletons cost nothing (each still sits
wholly on one side of a split, so ``no_group_leakage`` holds either way).

Public API:
    load_lineage(genomes_csv, clusters_csv)        -> np.ndarray (groups, genome order)
    group_kfold_masks(y, groups, n_splits)         -> [(train_mask, test_mask), ...]
    lineage_summary(groups)                         -> dict (n_clusters, sizes, ...)
    no_group_leakage(train_mask, test_mask, groups)-> bool
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def load_lineage(genomes_csv: str | Path, clusters_csv: str | Path, *,
                 genome_col: str = "Genome ID",
                 clusters_genome_col: str = "Genome ID",
                 cluster_col: str = "Cluster",
                 allow_missing: bool = False,
                 missing_label: str = "UNCLUSTERED") -> np.ndarray:
    """Return per-genome lineage labels aligned to ``genomes_csv`` row order.

    ``genomes_csv`` is the pipeline's ``genomes_{ab}.csv`` (column ``genome_col``),
    which fixes the row order of the matrix; ``clusters_csv`` is the CANONICAL
    lineage table written by 02c (``poppunk_clusters.csv``: ``Genome ID,Cluster``,
    already un-mangled — PopPUNK rewrites '.'→'_' in its raw Taxon column, which
    02c reverses, so here both files key on the same ``Genome ID``).

    Every genome in ``genomes_csv`` must have a cluster unless ``allow_missing`` —
    a missing genome means PopPUNK never saw it (re-run clustering). With
    ``allow_missing`` the gaps get ``missing_label`` (each treated as one group).
    """
    # dtype=str at read time (audit Issue 5): a Genome ID like "562.10" parsed as
    # float becomes 562.1 -> "562.1"; astype(str) after the fact cannot recover it.
    genomes = pd.read_csv(genomes_csv, encoding="utf-8",
                          dtype={genome_col: str})[genome_col].astype(str).tolist()
    cl = pd.read_csv(clusters_csv, encoding="utf-8", dtype={clusters_genome_col: str})
    for c in (clusters_genome_col, cluster_col):
        if c not in cl.columns:
            raise KeyError(
                f"Column '{c}' not in clusters file {clusters_csv} "
                f"(have: {cl.columns.tolist()})."
            )
    mapping = dict(zip(cl[clusters_genome_col].astype(str),
                       cl[cluster_col].astype(str)))

    missing = [g for g in genomes if g not in mapping]
    if missing and not allow_missing:
        raise ValueError(
            f"{len(missing)} genome(s) in {Path(genomes_csv).name} have no PopPUNK "
            f"cluster (e.g. {missing[:3]}). Re-run lineage clustering, or pass "
            f"allow_missing=True to pool them as '{missing_label}'."
        )
    groups = np.array([mapping.get(g, missing_label) for g in genomes], dtype=object)
    return groups


def group_kfold_masks(y: np.ndarray, groups: np.ndarray, n_splits: int = 5, *,
                      stratified: bool = True, seed: int = 42
                      ) -> list[tuple[np.ndarray, np.ndarray]]:
    """Grouped K-fold -> list of ``(train_mask, test_mask)`` boolean arrays.

    With ``stratified`` (default) uses ``StratifiedGroupKFold`` so folds keep
    groups intact AND balance the R/S label; otherwise plain ``GroupKFold``. The
    masks are sample-level (length == len(y)) and feed ``lib.xgb_data`` row_mask
    directly, so training/eval stay out-of-core (07b regime).
    """
    from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

    y = np.asarray(y)
    groups = np.asarray(groups)
    n = len(y)
    if len(groups) != n:
        raise ValueError(f"y ({n}) and groups ({len(groups)}) length mismatch.")
    n_groups = len(set(groups.tolist()))
    if n_splits > n_groups:
        raise ValueError(
            f"n_splits ({n_splits}) cannot exceed the number of lineages "
            f"({n_groups}); lower n_splits or pool rare clusters."
        )

    if stratified:
        splitter = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        gen = splitter.split(np.zeros(n), y, groups)
    else:
        splitter = GroupKFold(n_splits=n_splits)  # deterministic; no shuffle param
        gen = splitter.split(np.zeros(n), y, groups)

    masks = []
    for train_idx, test_idx in gen:
        tr = np.zeros(n, dtype=bool); tr[train_idx] = True
        te = np.zeros(n, dtype=bool); te[test_idx] = True
        masks.append((tr, te))
    return masks


def no_group_leakage(train_mask: np.ndarray, test_mask: np.ndarray,
                     groups: np.ndarray) -> bool:
    """True iff no lineage appears in both the train and test side of a fold."""
    groups = np.asarray(groups)
    return not (set(groups[train_mask].tolist()) & set(groups[test_mask].tolist()))


def lineage_summary(groups: np.ndarray) -> dict[str, Any]:
    """Descriptive stats for reporting / QC (n clusters, sizes, dominant clade)."""
    groups = np.asarray(groups)
    counts = Counter(groups.tolist())
    sizes = sorted(counts.values(), reverse=True)
    n = len(groups)
    largest_label, largest_size = (max(counts.items(), key=lambda kv: kv[1])
                                   if counts else (None, 0))
    return {
        "n_genomes": int(n),
        "n_clusters": int(len(counts)),
        "n_singletons": int(sum(1 for s in sizes if s == 1)),
        "largest_cluster": largest_label,
        "largest_cluster_size": int(largest_size),
        "largest_cluster_frac": float(largest_size / n) if n else 0.0,
        "median_cluster_size": float(np.median(sizes)) if sizes else 0.0,
    }
