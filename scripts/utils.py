#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared utility functions for the AMR Prediction pipeline.

Centralises helpers that were previously duplicated (and allowed to diverge)
across 03/04/05/06:

    - get_y_chunk():  contiguous label slicing for a chunk id
    - run_command():  subprocess wrapper using shlex (NO shell=True)

Import as:
    from utils import get_y_chunk, run_command

This works because each pipeline script is launched directly
(``python scripts/0X_*.py``), which puts the ``scripts/`` directory on
``sys.path``.
"""

import shlex
import subprocess
import sys


def get_y_chunk(y_all, chunk_id, chunk_size, total_len):
    """
    Extract the label subset corresponding to a specific data chunk.

    Chunks are contiguous, fixed-size slices of the full label array, matching
    the row order written by 03_matrix_construction.py.

    Args:
        y_all:      Complete array of all labels (or any sliceable sequence).
        chunk_id:   Chunk identifier (0-indexed).
        chunk_size: Number of samples per chunk.
        total_len:  Total number of samples.

    Returns:
        The slice of ``y_all`` belonging to the requested chunk.
    """
    start = chunk_id * chunk_size
    end = min((chunk_id + 1) * chunk_size, total_len)
    return y_all[start:end]


def run_command(command, exit_on_error=True):
    """
    Execute an external command safely (NO shell interpretation).

    The command string is tokenised with ``shlex.split`` and executed with
    ``shell=False`` (the subprocess default). This prevents shell-injection:
    special characters in genome IDs / file paths (e.g. ``; rm -rf /``) are
    passed as literal argument tokens, never interpreted by a shell.

    Stdout is suppressed to keep console output clean; stderr is captured and
    printed on failure to aid debugging KMC / kmc_tools errors.

    Args:
        command (str):       Command line to execute.
        exit_on_error (bool): If True (default), call sys.exit(1) on failure
                              (preserves historical behaviour of 03). If False,
                              return False on failure instead.

    Returns:
        bool: True on success. False on failure when ``exit_on_error`` is False.
    """
    try:
        subprocess.run(
            shlex.split(command),
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,   # Capture stderr for diagnostics
            text=True
        )
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: Command failed: {command}")
        print(f"Return code: {e.returncode}")
        if e.stderr:
            # Print first 5 lines of stderr to avoid flooding the console
            for line in e.stderr.strip().splitlines()[:5]:
                print(f"  STDERR: {line}")
        if exit_on_error:
            sys.exit(1)
        return False
