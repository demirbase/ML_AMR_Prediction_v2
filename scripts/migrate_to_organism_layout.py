#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
One-time, REVERSIBLE migration: single-organism layout -> {organism} layout
(SCALE_MLOPS_PLAN.md §8, Faz 2).

Moves the existing (E. coli) data/model/result/log directories into an
organism-scoped layout so the pipeline can scale to multiple organisms:

    data/raw/raw_genomes                         -> data/raw/{org}/genomes
    data/external/metadata/genome_amr_matrix.csv -> data/external/{org}/metadata/amr_phenotypes.csv
    data/interim/global_kmc_outputs              -> data/interim/{org}/kmc_outputs
    data/processed/{ab}                          -> data/processed/{org}/{ab}
    models/{ab}                                  -> models/{org}/{ab}
    results/{ab}                                 -> results/{org}/{ab}
    logs/{ab}                                    -> logs/{org}/{ab}

SAFETY:
    * Default mode is a DRY RUN — nothing is moved until you pass --apply.
    * --apply records every move in a manifest (runs/migration_manifest.json).
    * --revert reads that manifest and moves everything back.
    * A destination that already exists is skipped (never overwritten).
    * Only antibiotic directories listed in the registry for {organism} are
      moved, so global folders (e.g. results/global_exploration) are left alone.

Usage:
    python scripts/migrate_to_organism_layout.py                 # dry run (ecoli)
    python scripts/migrate_to_organism_layout.py --apply
    python scripts/migrate_to_organism_layout.py --organism ecoli --apply
    python scripts/migrate_to_organism_layout.py --revert
"""

import argparse
import json
import shutil
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import registry  # noqa: E402
from lib.config import load_config  # noqa: E402

MANIFEST_PATH = PROJECT_ROOT / "runs" / "migration_manifest.json"


def planned_moves(organism):
    """Build the list of (src, dst) Path pairs for the migration."""
    moves = []

    # 1. Global single-file / single-dir moves
    moves.append((
        PROJECT_ROOT / "data" / "raw" / "raw_genomes",
        PROJECT_ROOT / "data" / "raw" / organism / "genomes",
    ))
    moves.append((
        PROJECT_ROOT / "data" / "external" / "metadata" / "genome_amr_matrix.csv",
        PROJECT_ROOT / "data" / "external" / organism / "metadata" / "amr_phenotypes.csv",
    ))
    moves.append((
        PROJECT_ROOT / "data" / "interim" / "global_kmc_outputs",
        PROJECT_ROOT / "data" / "interim" / organism / "kmc_outputs",
    ))

    # 2. Per-antibiotic moves under processed/models/results/logs
    try:
        antibiotics = registry.get_organism(organism).get("antibiotics", [])
    except KeyError:
        antibiotics = []

    for base in ("data/processed", "models", "results", "logs"):
        base_dir = PROJECT_ROOT / base
        for ab in antibiotics:
            src = base_dir / ab
            dst = base_dir / organism / ab
            moves.append((src, dst))

    return moves


def do_apply(organism, moves):
    recorded = []
    for src, dst in moves:
        if not src.exists():
            print(f"  · skip (no source): {src.relative_to(PROJECT_ROOT)}")
            continue
        if dst.exists():
            print(f"  ⚠ skip (dest exists): {dst.relative_to(PROJECT_ROOT)}")
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(src), str(dst))
        print(f"  ✓ moved: {src.relative_to(PROJECT_ROOT)}  ->  {dst.relative_to(PROJECT_ROOT)}")
        recorded.append({"src": str(src), "dst": str(dst)})

    MANIFEST_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(MANIFEST_PATH, "w", encoding="utf-8") as f:
        json.dump({"organism": organism, "moves": recorded}, f, indent=2)
    print(f"\n✓ Migration applied. {len(recorded)} item(s) moved.")
    print(f"  Manifest (for --revert): {MANIFEST_PATH.relative_to(PROJECT_ROOT)}")
    print("\nNext: point scripts at the new layout via lib/config.resolve_path(),")
    print("      or update config.yaml `paths:` keys to the {organism} templates.")


def do_revert():
    if not MANIFEST_PATH.exists():
        print(f"ERROR: No migration manifest found at {MANIFEST_PATH}")
        sys.exit(1)
    with open(MANIFEST_PATH, "r", encoding="utf-8") as f:
        manifest = json.load(f)

    moves = manifest.get("moves", [])
    reverted = 0
    for item in reversed(moves):
        src, dst = Path(item["src"]), Path(item["dst"])
        if not dst.exists():
            print(f"  · skip (nothing at dest): {dst}")
            continue
        if src.exists():
            print(f"  ⚠ skip (original path occupied): {src}")
            continue
        src.parent.mkdir(parents=True, exist_ok=True)
        shutil.move(str(dst), str(src))
        print(f"  ✓ reverted: {dst}  ->  {src}")
        reverted += 1
    print(f"\n✓ Revert complete. {reverted} item(s) moved back.")


def main():
    parser = argparse.ArgumentParser(description="Migrate to {organism} data layout (reversible).")
    parser.add_argument("--organism", default=None,
                        help="Organism slug (default: config.yaml project.organism).")
    parser.add_argument("--apply", action="store_true", help="Execute the moves (default is dry run).")
    parser.add_argument("--revert", action="store_true", help="Undo a previous migration using the manifest.")
    args = parser.parse_args()

    if args.revert:
        do_revert()
        return

    organism = args.organism or load_config().get("project", {}).get("organism", "ecoli")
    print("=" * 78)
    print(f"ORGANISM LAYOUT MIGRATION — organism = '{organism}'  "
          f"({'APPLY' if args.apply else 'DRY RUN'})")
    print("=" * 78)

    moves = planned_moves(organism)

    if not args.apply:
        print("Planned moves (no changes made — pass --apply to execute):\n")
        for src, dst in moves:
            tag = "ok   " if src.exists() and not dst.exists() else "skip "
            print(f"  [{tag}] {src.relative_to(PROJECT_ROOT)}")
            print(f"          -> {dst.relative_to(PROJECT_ROOT)}")
        print("\nRun again with --apply to perform the migration.")
        return

    do_apply(organism, moves)


if __name__ == "__main__":
    main()
