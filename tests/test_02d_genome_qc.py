#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit tests for 02d_genome_qc.py (M15 CheckM2 + QUAST genome QC).

Exercise the pure logic without the tools: the CheckM2/QUAST TSV readers on
hand-written reports, and classify_row's threshold gating (incl. the "missing
metric is not a failure" rule)."""

import pytest

THR = {"completeness_min": 95.0, "contamination_max": 5.0,
       "n50_min": 50000, "max_contigs": 500}


@pytest.fixture
def mod(load_script):
    return load_script("02d_genome_qc.py")


def test_classify_row_pass_and_fail(mod):
    ok = mod.classify_row("g1", 99.0, 1.0, 120000, 80, 5_000_000, THR)
    assert ok["pass_overall"] is True
    assert all(ok[k] for k in ("pass_completeness", "pass_contamination",
                               "pass_n50", "pass_contigs"))
    # low completeness + high contamination + fragmented -> fail, per-check flags set
    bad = mod.classify_row("g2", 80.0, 9.0, 12000, 900, 5_200_000, THR)
    assert bad["pass_overall"] is False
    assert bad["pass_completeness"] is False
    assert bad["pass_contamination"] is False
    assert bad["pass_n50"] is False
    assert bad["pass_contigs"] is False


def test_classify_row_boundaries_inclusive(mod):
    # thresholds are inclusive (>=, <=)
    r = mod.classify_row("g", 95.0, 5.0, 50000, 500, 5_000_000, THR)
    assert r["pass_overall"] is True


def test_classify_row_missing_metric_not_a_failure(mod):
    # QUAST only (CheckM2 absent): completeness/contamination None, still passes
    # on the present QUAST checks.
    r = mod.classify_row("g", None, None, 120000, 50, 5_000_000, THR)
    assert r["pass_completeness"] is None and r["pass_contamination"] is None
    assert r["pass_overall"] is True
    # all metrics missing -> no present checks -> not a pass
    empty = mod.classify_row("g", None, None, None, None, None, THR)
    assert empty["pass_overall"] is False


def test_read_checkm2_and_quast(mod, tmp_path):
    cm = tmp_path / "checkm2"
    cm.mkdir()
    (cm / "quality_report.tsv").write_text(
        "Name\tCompleteness\tContamination\n"
        "562.100\t99.1\t0.5\n562.200\t72.0\t8.3\n", encoding="utf-8")
    got = mod._read_checkm2(cm)
    assert got["562.100"] == (99.1, 0.5)
    assert got["562.200"] == (72.0, 8.3)

    q = tmp_path / "quast"
    q.mkdir()
    (q / "transposed_report.tsv").write_text(
        "Assembly\t# contigs\tTotal length\tN50\n"
        "562.100\t80\t5000000\t120000\n562.200\t900\t5200000\t12000\n", encoding="utf-8")
    gq = mod._read_quast(q)
    assert gq["562.100"] == (120000.0, 80, 5000000)
    assert gq["562.200"] == (12000.0, 900, 5200000)


def test_readers_return_none_when_absent(mod, tmp_path):
    assert mod._read_checkm2(tmp_path / "nope") is None
    assert mod._read_quast(tmp_path / "nope") is None
