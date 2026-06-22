#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for 03u_unitig_matrix.py (rtab -> genome×unitig CSR transpose).

These exercise the core logic without invoking unitig-caller: a synthetic rtab
is written by hand and fed to rtab_to_chunks(). They verify (a) the unitig→column
order matches features.txt, (b) the rtab sample columns are mapped to the correct
output rows even when the rtab header order differs from the genome order, (c)
the absolute support filter drops singletons (below min_support) and zero-variance
core (present in every genome), and (d) the chunked CSR files reconstruct the
exact dense matrix.
"""

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import load_npz, vstack


def _write_rtab(path, header_samples, rows):
    """rows: list of (unitig_seq, [0/1 in header_samples order])."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("Unitig_sequence\t" + "\t".join(header_samples) + "\n")
        for seq, vals in rows:
            f.write(seq + "\t" + "\t".join(str(v) for v in vals) + "\n")


def _reconstruct(out_dir, antibiotic, n_chunks):
    parts = [load_npz(out_dir / f"X_{antibiotic}_part_{c}.npz") for c in range(n_chunks)]
    return vstack(parts).toarray()


@pytest.fixture
def mod(load_script):
    return load_script("03u_unitig_matrix.py")


def test_transpose_mapping_and_support_filter(mod, tmp_path):
    # Our canonical genome order:
    valid_genomes = ["g0", "g1", "g2", "g3"]
    valid_labels = [1, 0, 1, 0]
    # rtab header in a DIFFERENT order, to prove we map by name not position:
    header = ["g2", "g0", "g3", "g1"]
    # Values are written in header order.
    rows = [
        ("CORE", [1, 1, 1, 1]),   # present in all -> dropped (zero-variance core)
        ("SOLO", [0, 1, 0, 0]),   # only g0 -> support 1 -> dropped at min_support=2
        ("UAB",  [1, 1, 0, 0]),   # g2,g0 -> i.e. g0 & g2 present
        ("UCD",  [0, 0, 1, 1]),   # g3,g1 -> i.e. g1 & g3 present
    ]
    rtab = tmp_path / "unitigs.rtab"
    _write_rtab(rtab, header, rows)

    out_dir = tmp_path / "matrix_unitig"
    out_dir.mkdir()

    n_unitigs, n_chunks = mod.rtab_to_chunks(
        rtab, valid_genomes, valid_labels, out_dir,
        antibiotic="testdrug", chunk_size=2, min_support=2,
    )

    assert n_unitigs == 2          # CORE + SOLO dropped
    assert n_chunks == 2           # 4 genomes / chunk_size 2

    # features.txt == kept unitigs in column order (stream order: UAB, UCD)
    feats = (out_dir / "features.txt").read_text().splitlines()
    assert [ln.split("\t")[0] for ln in feats] == ["UAB", "UCD"]
    assert all(ln.endswith("\t1") for ln in feats)

    # Reconstructed dense matrix, rows in valid_genomes order, cols [UAB, UCD]
    dense = _reconstruct(out_dir, "testdrug", n_chunks)
    expected = np.array([
        [1, 0],   # g0: UAB present, UCD absent
        [0, 1],   # g1: UCD present
        [1, 0],   # g2: UAB present
        [0, 1],   # g3: UCD present
    ], dtype=np.int8)
    assert np.array_equal(dense, expected)
    assert dense.dtype == np.int8 and dense.max() <= 1

    # y / genomes csv match input order
    y = pd.read_csv(out_dir / "y_testdrug.csv")["label"].tolist()
    g = pd.read_csv(out_dir / "genomes_testdrug.csv")["Genome ID"].astype(str).tolist()
    assert y == valid_labels
    assert g == valid_genomes


def test_min_support_one_keeps_singletons(mod, tmp_path):
    valid_genomes = ["g0", "g1", "g2", "g3"]
    valid_labels = [1, 1, 0, 0]
    header = ["g0", "g1", "g2", "g3"]
    rows = [
        ("CORE", [1, 1, 1, 1]),   # still dropped (core) regardless of min_support
        ("SOLO", [1, 0, 0, 0]),   # kept at min_support=1
        ("PAIR", [0, 1, 1, 0]),
    ]
    rtab = tmp_path / "unitigs.rtab"
    _write_rtab(rtab, header, rows)
    out_dir = tmp_path / "m"
    out_dir.mkdir()

    n_unitigs, n_chunks = mod.rtab_to_chunks(
        rtab, valid_genomes, valid_labels, out_dir,
        antibiotic="testdrug", chunk_size=200, min_support=1,
    )
    assert n_unitigs == 2          # SOLO + PAIR (CORE dropped)
    feats = [ln.split("\t")[0] for ln in (out_dir / "features.txt").read_text().splitlines()]
    assert feats == ["SOLO", "PAIR"]
    dense = _reconstruct(out_dir, "testdrug", n_chunks)
    assert np.array_equal(dense, np.array([[1, 0], [0, 1], [0, 1], [0, 0]], dtype=np.int8))


def test_rejects_sample_not_in_genome_set(mod, tmp_path):
    rtab = tmp_path / "unitigs.rtab"
    _write_rtab(rtab, ["g0", "ZZZ"], [("UAB", [1, 1])])
    out_dir = tmp_path / "m"
    out_dir.mkdir()
    with pytest.raises(SystemExit):
        mod.rtab_to_chunks(rtab, ["g0", "g1"], [1, 0], out_dir,
                           antibiotic="testdrug", chunk_size=10, min_support=1)
