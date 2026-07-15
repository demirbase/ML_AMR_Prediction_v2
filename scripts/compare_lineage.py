#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare two PopPUNK clusterings of the same genomes — the container-rebuild check.

Rebuilding the pipeline container re-solves transitive dependencies. PopPUNK's
own version can be pinned while `graph-tool` (its network-analysis backend) jumps
a major version underneath it — which is exactly what happened on 2026-07-15
(2.98 -> 3.0). PopPUNK clusters ARE the cross-validation groups, so if they shift,
every lineage-aware AUC shifts with them and results stop being comparable to the
published ones. "The build succeeded" says nothing about that; only re-running the
OLD settings and diffing the labels does.

Cluster IDs are arbitrary: PopPUNK may number the same partition differently. So
identity is judged on the PARTITION (which genomes share a cluster), via the
Adjusted Rand Index, not on the labels.

Usage:
    python scripts/compare_lineage.py <old_clusters.csv> <new_clusters.csv>

Exit code 0 == identical partition (ARI 1.0), 1 == they differ.
"""

import sys
from pathlib import Path

import pandas as pd


def load(path):
    df = pd.read_csv(path, dtype=str)
    missing = {"Genome ID", "Cluster"} - set(df.columns)
    if missing:
        sys.exit(f"ERROR: {path} lacks column(s): {sorted(missing)}")
    return df.set_index("Genome ID")["Cluster"]


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    old_path, new_path = Path(sys.argv[1]), Path(sys.argv[2])
    for p in (old_path, new_path):
        if not p.exists():
            sys.exit(f"ERROR: not found: {p}")

    old, new = load(old_path), load(new_path)

    print("=" * 72)
    print("POPPUNK CLUSTERING COMPARISON")
    print("=" * 72)
    print(f"  old: {old_path}  ({len(old)} genomes, {old.nunique()} clusters)")
    print(f"  new: {new_path}  ({len(new)} genomes, {new.nunique()} clusters)")

    only_old = set(old.index) - set(new.index)
    only_new = set(new.index) - set(old.index)
    if only_old or only_new:
        print(f"\n  ⚠ genome sets differ — old-only {len(only_old)}, new-only {len(only_new)}")
        print("    (a QC pass that drops genomes will do this; comparing the overlap)")
    shared = sorted(set(old.index) & set(new.index))
    if not shared:
        sys.exit("ERROR: no genomes in common — are these the same organism?")

    from sklearn.metrics import adjusted_rand_score
    ari = adjusted_rand_score(old.loc[shared].to_numpy(), new.loc[shared].to_numpy())

    print(f"\n  shared genomes    : {len(shared)}")
    print(f"  Adjusted Rand Index: {ari:.6f}")

    identical = (ari == 1.0) and not only_old and not only_new
    if identical:
        print("\n  ✓ IDENTICAL partition — the rebuilt container reproduces the clustering.")
        print("    Lineage-aware results stay comparable to the published runs.")
        return 0

    print("\n  ✗ CLUSTERING CHANGED.")
    print("    The CV groups are not the same, so lineage-aware AUCs from this")
    print("    container are NOT comparable to earlier ones. Before going further:")
    print("      - if only the container changed, a dependency altered PopPUNK's")
    print("        behaviour (check graph-tool / pp-sketchlib in environment.lock.yml)")
    print("        -> pin it to the previous version and rebuild, or accept the new")
    print("           clustering knowingly and re-run EVERYTHING on it;")
    print("      - if sketch/QC settings also changed, this comparison proves nothing:")
    print("        re-run with the OLD settings so only the container varies.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
