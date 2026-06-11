#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Centralised config loading + path resolution (SCALE_MLOPS_PLAN.md §4.2).

This module is the forward-looking, {organism}-aware path layer. It is
ADDITIVE: the existing numbered scripts keep reading the legacy ``paths:`` keys
directly and continue to work unchanged. New code (orchestrator, run metadata,
the migration script, future multi-organism runs) uses resolve_path() so that
adding an organism never requires touching path-construction code.

Public API:
    load_config()                                  -> dict   (global config.yaml)
    get_target(args=None)                          -> (organism, antibiotic)
    resolve_path(key, organism=, antibiotic=, run_id=) -> Path
"""

import os
from pathlib import Path

import yaml

# scripts/lib/config.py  ->  parents[2] == project root
PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONFIG_FILE = PROJECT_ROOT / "config" / "config.yaml"


def load_config(config_path=None):
    """Load and return the global config.yaml as a dict."""
    path = Path(config_path) if config_path else CONFIG_FILE
    if not path.exists():
        raise FileNotFoundError(f"Configuration file not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def get_target(args=None, config=None):
    """
    Resolve the (organism, antibiotic) target with the precedence:
        CLI args  >  environment variables  >  config.yaml defaults.

    This preserves the legacy "edit config.yaml, run the script" workflow
    (falls back to config), while enabling parameterised invocation.

    Args:
        args:   optional argparse.Namespace with .organism / .antibiotic.
        config: optional pre-loaded config dict (avoids re-reading the file).

    Returns:
        tuple(str, str): (organism, antibiotic)
    """
    cfg = config if config is not None else load_config()
    proj = cfg.get("project", {})

    organism = (
        getattr(args, "organism", None)
        or os.environ.get("AMR_ORGANISM")
        or proj.get("organism")
    )
    antibiotic = (
        getattr(args, "antibiotic", None)
        or os.environ.get("AMR_ANTIBIOTIC")
        or proj.get("target_antibiotic")
    )
    return organism, antibiotic


def resolve_path(key, organism=None, antibiotic=None, run_id=None, config=None):
    """
    Resolve a path template from config into an absolute Path.

    Looks the key up first in the new ``paths_organism:`` block (the
    {organism}-aware templates), then falls back to the legacy ``paths:`` block.
    Any ``{organism}`` / ``{antibiotic}`` / ``{run_id}`` placeholders present in
    the template are filled in.

    Args:
        key:        path key, e.g. "matrix_dir", "genomes_dir", "run_dir".
        organism:   organism slug (required if the template uses {organism}).
        antibiotic: antibiotic id (required if the template uses {antibiotic}).
        run_id:     run identifier (required if the template uses {run_id}).
        config:     optional pre-loaded config dict.

    Returns:
        Path: PROJECT_ROOT-anchored absolute path.
    """
    cfg = config if config is not None else load_config()
    paths_org = cfg.get("paths_organism", {}) or {}
    paths_legacy = cfg.get("paths", {}) or {}

    template = paths_org.get(key, paths_legacy.get(key))
    if template is None:
        raise KeyError(
            f"Path key '{key}' not found in config 'paths_organism:' or 'paths:'."
        )

    fmt = {}
    if organism is not None:
        fmt["organism"] = organism
    if antibiotic is not None:
        fmt["antibiotic"] = antibiotic
    if run_id is not None:
        fmt["run_id"] = run_id

    try:
        resolved = template.format(**fmt) if "{" in template else template
    except KeyError as missing:
        raise KeyError(
            f"Path template for '{key}' needs placeholder {missing} "
            f"but it was not provided (organism={organism}, "
            f"antibiotic={antibiotic}, run_id={run_id})."
        )

    return PROJECT_ROOT / resolved
