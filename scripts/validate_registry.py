#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Registry consistency guard (M1 — "watertight registry").

Checks the registry files (config/registry/{organisms,antibiotics}.yaml) for the
invariants the pipeline relies on, and — with --db — that the KB agrees with the
registry (the registry is the single source of truth; the KB is its derivative).

    python scripts/validate_registry.py
    python scripts/validate_registry.py --db results/kb/amrk.db

Exit code 0 = all checks pass, 1 = one or more violations (suitable for CI).
"""
from __future__ import annotations

import argparse
import re
import sqlite3
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import registry  # noqa: E402

_SLUG_RE = re.compile(r"^[a-z][a-z0-9_]*$")
_AWARE_VALUES = {"Access", "Watch", "Reserve"}
_MECH_VALUES = {"acquired", "target_snp"}


def _check_registry(errors: list[str], warnings: list[str]) -> None:
    doc = registry._antibiotics_doc()
    classes = doc.get("classes", {})

    # 1) every member belongs to exactly one class
    seen: dict[str, str] = {}
    for cid, block in classes.items():
        for m in block.get("members", []):
            key = str(m).lower()
            if key in seen:
                errors.append(f"antibiotic '{m}' is in two classes: {seen[key]} and {cid}")
            seen[key] = cid

    # 2) every alias canonical is a member
    for canonical in (doc.get("aliases", {}) or {}):
        if str(canonical).lower() not in seen:
            errors.append(f"alias canonical '{canonical}' is not a member of any class")

    # 3) class_mechanism_type keys valid, values in {acquired,target_snp}
    for cid, val in (doc.get("class_mechanism_type", {}) or {}).items():
        if cid not in classes:
            errors.append(f"class_mechanism_type references unknown class '{cid}'")
        if val not in _MECH_VALUES:
            errors.append(f"class_mechanism_type['{cid}'] = '{val}' not in {sorted(_MECH_VALUES)}")

    # 4) who_aware antibiotics are members; values valid
    for ab, val in (doc.get("who_aware", {}) or {}).items():
        if str(ab).lower() not in seen:
            warnings.append(f"who_aware lists '{ab}' which is not a class member")
        if val not in _AWARE_VALUES:
            errors.append(f"who_aware['{ab}'] = '{val}' not in {sorted(_AWARE_VALUES)}")

    # 5) amrfinder_keywords antibiotics are members
    for ab in registry.load_amrfinder_keywords():
        if str(ab).lower() not in seen:
            warnings.append(f"amrfinder_keywords lists '{ab}' which is not a class member")

    # 6) organisms: lowercase slugs, valid status, priority_classes exist,
    #    antibiotics classified
    for slug, block in registry.load_organisms().items():
        if not _SLUG_RE.match(slug):
            errors.append(f"organism slug '{slug}' is not lowercase snake_case")
        # An unknown status silently makes the organism INACTIVE (is_active() just
        # tests membership), so a typo would quietly drop it from the panel with
        # nothing failing. Pin the vocabulary shut.
        status = block.get("status")
        if status is not None and status not in registry.VALID_STATUS:
            errors.append(f"organism '{slug}' status '{status}' not in "
                          f"{sorted(registry.VALID_STATUS)}")
        for cid in block.get("priority_classes", []) or []:
            if cid not in classes:
                errors.append(f"organism '{slug}' priority_class '{cid}' is not a known class")
        for ab in block.get("antibiotics", []) or []:
            if registry.antibiotic_to_class(ab) is None:
                warnings.append(f"organism '{slug}' target '{ab}' has no registry class")


def _check_kb(db: Path, errors: list[str], warnings: list[str]) -> None:
    if not db.exists():
        warnings.append(f"KB not found at {db} — skipped KB↔registry checks")
        return
    conn = sqlite3.connect(str(db))
    known_orgs = set(registry.load_organisms())

    for (org,) in conn.execute("SELECT DISTINCT organism FROM pipeline_runs"):
        if org not in known_orgs:
            errors.append(f"KB model organism '{org}' is not in the registry")

    for ab, dc, mech, aware in conn.execute(
            "SELECT antibiotic, drug_class, mechanism_type, who_aware FROM antibiotics"):
        exp_dc = registry.antibiotic_to_class(ab)
        if dc != exp_dc:
            errors.append(f"KB drug_class['{ab}'] = {dc!r} but registry says {exp_dc!r}")
        exp_mech = registry.antibiotic_mechanism_type(ab)
        if mech != exp_mech:
            errors.append(f"KB mechanism_type['{ab}'] = {mech!r} but registry says {exp_mech!r}")
        exp_aware = registry.antibiotic_who_aware(ab)
        if aware != exp_aware:
            errors.append(f"KB who_aware['{ab}'] = {aware!r} but registry says {exp_aware!r}")
    conn.close()


def main() -> int:
    ap = argparse.ArgumentParser(description="Validate the AMR registry (and optionally the KB).")
    ap.add_argument("--db", default=None, help="KB SQLite path for registry↔KB consistency checks")
    args = ap.parse_args()

    registry.clear_cache()
    errors: list[str] = []
    warnings: list[str] = []
    _check_registry(errors, warnings)
    if args.db:
        _check_kb(Path(args.db), errors, warnings)

    for w in warnings:
        print(f"  WARN  {w}")
    for e in errors:
        print(f"  FAIL  {e}")
    print(f"\nregistry validation: {len(errors)} error(s), {len(warnings)} warning(s)")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
