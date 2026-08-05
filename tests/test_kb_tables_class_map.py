#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Guard: kb_tables' drug-class -> ARO keyword map must cover the antibiotic registry.

`mechanisms.csv` marks a CARD gene on-target by matching the model's drug class to a
tuple of ARO drug-class keywords. A class missing from that map yields on_target=None,
which is not an error anywhere — the row is written, the column is simply blank, and
figure 04 (which filters `on_target == True`) silently drops the model.

That is exactly what happened: the map still held the pre-curation taxonomy, including
a `beta_lactams_carbapenems_others` key that no longer exists, and lacked carbapenems,
glycopeptides, macrolides, lincosamides, monobactams, phenicols, polymyxins and
glycylcyclines. 14 of 45 models could never be on-target, so the mechanism figure was
missing K. pneumoniae KPC/NDM and the entire E. faecium vanA cluster — two of the four
biology headlines — with nothing failing to signal it.

These tests fail loudly when the registry gains a class the map does not know, and when
the map keeps a key the registry has dropped.
"""
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

REGISTRY = PROJECT_ROOT / "config" / "registry" / "antibiotics.yaml"


def _registry_classes():
    yaml = pytest.importorskip("yaml")
    if not REGISTRY.exists():
        pytest.skip(f"antibiotic registry not found: {REGISTRY}")
    doc = yaml.safe_load(REGISTRY.read_text(encoding="utf-8")) or {}
    classes = doc.get("classes")
    if not isinstance(classes, dict) or not classes:
        pytest.skip("could not read `classes` from the antibiotic registry")
    return set(classes)


def test_class_keyword_map_covers_registry():
    from kb_tables import CLASS_TO_ARO_KEYWORD

    missing = sorted(_registry_classes() - set(CLASS_TO_ARO_KEYWORD))
    assert not missing, (
        "drug classes with no ARO keyword in kb_tables.CLASS_TO_ARO_KEYWORD: "
        f"{missing}. Models in these classes get on_target=None and vanish from "
        "figure 04 without any error being raised."
    )


def test_class_keyword_map_has_no_dead_keys():
    from kb_tables import CLASS_TO_ARO_KEYWORD

    dead = sorted(set(CLASS_TO_ARO_KEYWORD) - _registry_classes())
    assert not dead, (
        f"CLASS_TO_ARO_KEYWORD keys absent from the antibiotic registry: {dead}. "
        "A stale key is how the map drifted out of date the first time."
    )


def test_every_class_maps_to_nonempty_keywords():
    from kb_tables import CLASS_TO_ARO_KEYWORD

    empty = sorted(k for k, v in CLASS_TO_ARO_KEYWORD.items() if not v)
    assert not empty, f"classes mapped to an empty keyword tuple: {empty}"
