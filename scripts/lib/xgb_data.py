#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Streaming XGBoost DMatrix construction from on-disk genome×feature chunks.

Feature-agnostic: the columns are k-mers (03) or unitigs (03u, the canonical
representation) — the streaming logic is identical. Replaces the legacy
incremental "1 tree per chunk" out-of-core regime. The full sparse matrix (e.g.
~109 GB for ~5k E. coli genomes × 50.8M k-mers, or ~8 GB for ~4.9M unitigs) is
never materialised in RAM: XGBoost pulls one chunk at a time through
``ChunkDMatrixIter`` and builds a single compact, quantised ``QuantileDMatrix``
(binary 0/1 data → ``max_bin=2`` → ~1 byte per non-zero). This enables standard
full-data gradient boosting with bounded memory (~one chunk + the histogram),
instead of training each tree on a single ~200-genome slice.

Used by 05_model_training (chunk-level train set) and 07b_feature_stability
(sample-level seed masks). Organism/antibiotic-agnostic: it simply streams
whatever chunk files it is handed, so the same code path serves every dataset.
"""

import numpy as np
import xgboost as xgb
from scipy.sparse import load_npz

from lib.chunking import get_y_chunk


class ChunkDMatrixIter(xgb.DataIter):
    """Yield genome×feature (k-mer/unitig) chunks one at a time for QuantileDMatrix construction.

    Args:
        files:       Ordered chunk file paths (``X_{ab}_part_{n}.npz``).
        y_all:       Full label array; row order matches 03's matrix output.
        chunk_size:  Rows per chunk (``preprocessing.chunk_size``). Chunk ``n``
                     occupies global rows ``[n*chunk_size : (n+1)*chunk_size]``
                     — the same contiguous convention as ``get_y_chunk``.
        row_mask:    Optional bool array over all rows; only ``True`` rows are
                     kept (sample-level splits, e.g. 07b seed holdout). ``None``
                     keeps every row.
        pos_weight:  Optional instance weight applied to positive-label rows
                     (typically the global neg/pos ratio); negatives get 1.0.
                     ``None`` trains unweighted.
        cache_prefix: Optional path prefix. When set, XGBoost spills the
                     quantised pages to disk here (external memory) instead of
                     holding the whole matrix in RAM — required when the full
                     train set is too large for the node (see build_quantile_dmatrix).
    """

    def __init__(self, files, y_all, chunk_size, *, row_mask=None, pos_weight=None,
                 cache_prefix=None):
        self._files = list(files)
        self._y_all = np.asarray(y_all)
        self._chunk_size = int(chunk_size)
        self._total = len(self._y_all)
        self._row_mask = row_mask
        self._pos_weight = pos_weight
        self._i = 0
        super().__init__(cache_prefix=cache_prefix)

    def reset(self):
        """Rewind to the first chunk (XGBoost calls this between passes)."""
        self._i = 0

    def next(self, input_data):
        """Push the next (masked) chunk; return 1 if data was supplied, else 0."""
        while self._i < len(self._files):
            f = self._files[self._i]
            self._i += 1
            chunk_id = int(f.stem.split('_')[-1])
            X = load_npz(f)
            y = np.asarray(get_y_chunk(self._y_all, chunk_id, self._chunk_size, self._total))
            if self._row_mask is not None:
                start = chunk_id * self._chunk_size
                local = self._row_mask[start:start + X.shape[0]]
                if not local.any():
                    continue  # this chunk contributes no rows to the split
                X = X[local]
                y = y[local]
            kwargs = {}
            if self._pos_weight is not None:
                kwargs["weight"] = np.where(y == 1, float(self._pos_weight), 1.0)
            input_data(data=X, label=y, **kwargs)
            return 1
        return 0


def global_pos_weight(files, y_all, chunk_size, row_mask=None):
    """Global neg/pos ratio over the given chunks (optionally row-masked).

    Returns 1.0 when there are no positives (degenerate split) so callers never
    divide by zero.
    """
    y_all = np.asarray(y_all)
    total = len(y_all)
    pos = neg = 0
    for f in files:
        cid = int(f.stem.split('_')[-1])
        y = get_y_chunk(y_all, cid, chunk_size, total)
        if row_mask is not None:
            start = cid * chunk_size
            y = y[row_mask[start:start + len(y)]]
        p = int(np.sum(y))
        pos += p
        neg += len(y) - p
    return (neg / pos) if pos > 0 else 1.0


def build_quantile_dmatrix(files, y_all, chunk_size, *, max_bin=2,
                           row_mask=None, pos_weight=None, cache_prefix=None, ref=None):
    """Stream the given chunks into a single quantised DMatrix.

    With ``cache_prefix`` set, builds an ``ExtMemQuantileDMatrix`` that spills
    pages to disk (external memory) — use this for the full train set, which can
    need far more RAM in-core than the node has (the ~109 GB matrix peaked >400 GB
    as a plain in-core QuantileDMatrix). Without it, builds an in-core
    ``QuantileDMatrix`` (fine for small subsets, e.g. HPO).

    ``ref`` must be the training DMatrix when building an *evaluation* matrix:
    XGBoost requires eval QuantileDMatrices to reuse the training quantiles.
    """
    it = ChunkDMatrixIter(files, y_all, chunk_size, row_mask=row_mask,
                          pos_weight=pos_weight, cache_prefix=cache_prefix)
    if cache_prefix is not None:
        return xgb.ExtMemQuantileDMatrix(it, max_bin=max_bin, ref=ref)
    return xgb.QuantileDMatrix(it, max_bin=max_bin, ref=ref)
