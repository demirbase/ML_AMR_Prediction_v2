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
from scipy.sparse import csr_matrix, load_npz, save_npz, vstack


def _write_store(tmp_path, store_genomes, unitig_seqs, presence):
    """presence[col] = iterable of store-row indices where unitig `col` is present."""
    dense = np.zeros((len(store_genomes), len(unitig_seqs)), dtype=np.int8)
    for col, rows in enumerate(presence):
        for r in rows:
            dense[r, col] = 1
    store = tmp_path / "unitig_all"
    store.mkdir()
    save_npz(store / "X_all_part_0.npz", csr_matrix(dense))
    pd.DataFrame({"Genome ID": store_genomes}).to_csv(store / "genomes_all.csv", index=False)
    with open(store / "features.txt", "w", encoding="utf-8") as f:
        for s in unitig_seqs:
            f.write(f"{s}\t1\n")
    return store


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


def test_subset_store_to_antibiotic(mod, tmp_path):
    # Organism store: 5 genomes × 4 unitigs.
    store_genomes = ["g0", "g1", "g2", "g3", "g4"]
    seqs = ["U0", "U1", "U2", "U3"]
    presence = [
        {0, 1, 2},        # U0 in g0,g1,g2
        {0, 3},           # U1 in g0,g3
        {1, 2, 3, 4},     # U2 in g1,g2,g3,g4
        {0, 2, 3},        # U3 in g0,g2,g3
    ]
    store = _write_store(tmp_path, store_genomes, seqs, presence)

    out_dir = tmp_path / "matrix_unitig"
    out_dir.mkdir()
    # Antibiotic subset (different order) g2,g0,g3 ; min_support=2, n_sel=3 -> max_support=2
    valid_genomes = ["g2", "g0", "g3"]
    valid_labels = [1, 0, 1]
    n_unitigs, n_chunks = mod.subset_store_to_antibiotic(
        store, out_dir, "testdrug", valid_genomes, valid_labels,
        chunk_size=2, min_support=2)

    # Support over the subset: U0=2, U1=2, U2=2, U3=3(==n_sel -> zero-variance core, dropped)
    assert n_unitigs == 3
    feats = [ln.split("\t")[0] for ln in (out_dir / "features.txt").read_text().splitlines()]
    assert feats == ["U0", "U1", "U2"]

    dense = _reconstruct(out_dir, "testdrug", n_chunks)
    expected = np.array([
        [1, 0, 1],   # g2: U0 yes, U1 no,  U2 yes
        [1, 1, 0],   # g0: U0 yes, U1 yes, U2 no
        [0, 1, 1],   # g3: U0 no,  U1 yes, U2 yes
    ], dtype=np.int8)
    assert np.array_equal(dense, expected)
    assert pd.read_csv(out_dir / "y_testdrug.csv")["label"].tolist() == valid_labels
    assert pd.read_csv(out_dir / "genomes_testdrug.csv")["Genome ID"].tolist() == valid_genomes
