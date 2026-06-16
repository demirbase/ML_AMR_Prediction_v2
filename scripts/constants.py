#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Backward-compatibility shim.

The canonical antibiotic-class data now lives in config/registry/antibiotics.yaml
and is accessed through scripts/lib/registry.py (SCALE_MLOPS_PLAN.md §3). This
module is kept so existing ``from constants import ANTIBIOTIC_CLASSES`` imports
continue to work; it simply re-exports the registry-backed dictionary.

New code should prefer:
    from lib.registry import load_antibiotic_classes
"""

from lib.registry import load_antibiotic_classes

# {ClassDisplayName: [members]} — identical structure to the old hardcoded dict.
ANTIBIOTIC_CLASSES = load_antibiotic_classes()
