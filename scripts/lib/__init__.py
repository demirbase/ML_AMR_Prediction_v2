#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared library for the AMR Prediction pipeline (SCALE_MLOPS_PLAN.md §5).

Canonical home for code that was previously duplicated across the numbered
scripts. The top-level ``scripts/constants.py`` and ``scripts/utils.py`` remain
as thin backward-compatibility shims that re-export from here, so existing
``from utils import ...`` / ``from constants import ...`` imports keep working.

Submodules:
    registry      — organisms.yaml / antibiotics.yaml access (single source)
    config        — global config loader + {organism}/{antibiotic} path resolver
    chunking      — get_y_chunk (contiguous label slicing)
    io_utils      — run_command (shlex-based, never shell=True)
    run_metadata  — git hash / version capture, run_id generation (MLOps)
"""
