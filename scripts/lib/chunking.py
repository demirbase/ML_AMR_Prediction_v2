#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Contiguous chunk label slicing — canonical single copy (SCALE_MLOPS_PLAN §5)."""


def get_y_chunk(y_all, chunk_id, chunk_size, total_len):
    """
    Extract the label subset corresponding to a specific data chunk.

    Chunks are contiguous, fixed-size slices of the full label array, matching
    the row order written by 03_matrix_construction.py (k-mers) / 03u_unitig_matrix.py
    (unitigs).

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
