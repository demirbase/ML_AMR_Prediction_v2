#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Unit/smoke tests for 15_cross_antibiotic.py (S1 cross-antibiotic overlap).

Two layers, no external tools:
  * pure functions — ``_drug_family`` (β-lactam super-family collapse for H3)
    and ``hypergeom_sf`` (enrichment p / exact fallback);
  * a synthetic in-memory KB (real ``kb_schema``) exercising ``fill_drug_classes``,
    ``stable_sets`` and ``populate_overlap`` — the overlap-table write, the
    ``same_class`` (strict registry class) and ``same_drug_family`` (H3-ready)
    flags, and the Jaccard summary.
"""

import sqlite3

import pytest


@pytest.fixture
def mod(load_script):
    return load_script("15_cross_antibiotic.py")


# --- pure functions --------------------------------------------------------
def test_drug_family_collapses_betalactams(mod):
    # penicillins + cephalosporins + carbapenems + monobactams all collapse to
    # one family so an ampicillin/cefotaxime pair reads as "within β-lactam" (H3),
    # while a non-β-lactam class is left untouched.
    assert mod._drug_family("penicillins") == "beta_lactam"
    assert mod._drug_family("cephalosporins") == "beta_lactam"
    assert mod._drug_family("carbapenems") == "beta_lactam"
    assert mod._drug_family("monobactams") == "beta_lactam"
    assert mod._drug_family("quinolones") == "quinolones"
    assert mod._drug_family(None) is None


def test_hypergeom_sf_enrichment_and_degenerate(mod):
    # Universe N=10, two sets of 5, full overlap of 5 -> P(X>=5) = 1/C(10,5).
    p = mod.hypergeom_sf(k=5, N=10, K=5, n=5)
    assert p == pytest.approx(1 / 252, rel=1e-6)
    # More extreme overlap is never more probable than less extreme.
    assert mod.hypergeom_sf(5, 20, 5, 5) < mod.hypergeom_sf(2, 20, 5, 5)
    # Degenerate inputs -> None (no test, not a crash).
    assert mod.hypergeom_sf(0, 10, 5, 5) is None
    assert mod.hypergeom_sf(3, 0, 5, 5) is None


def _pair(a, b, fam, jac, ov, p=None):
    d = {"antibiotic_a": a, "antibiotic_b": b, "same_drug_family": fam,
         "jaccard": jac, "n_overlap": ov}
    if p is not None:
        d["hypergeom_p_enrichment"] = p
    return d


def test_h3_contrast_within_greater(mod):
    # ampicillin~cefotaxime (β-lactam, within) overlaps more than the cross-class
    # quinolone pairs -> H3 supported (descriptive).
    summaries = [
        _pair("ampicillin", "cefotaxime", True, 0.60, 3, p=1e-4),
        _pair("ampicillin", "ciprofloxacin", False, 0.14, 1),
        _pair("cefotaxime", "ciprofloxacin", False, 0.14, 1),
    ]
    h3 = mod.h3_contrast(summaries)
    assert h3["testable"] is True
    assert h3["verdict"] == "within_greater"
    assert h3["within_family"]["n_pairs"] == 1
    assert h3["cross_class"]["n_pairs"] == 2
    assert h3["within_family"]["mean_jaccard"] > h3["cross_class"]["mean_jaccard"]
    assert h3["within_family_min_p_enrichment"] == 1e-4


def test_h3_contrast_not_testable_single_cross_pair(mod):
    # The current canonical state: one cross-class pair, no within-family pair.
    h3 = mod.h3_contrast([_pair("ampicillin", "ciprofloxacin", False, 0.02, 4)])
    assert h3["testable"] is False
    assert h3["verdict"] is None
    assert h3["within_family"]["n_pairs"] == 0
    assert h3["within_family_min_p_enrichment"] is None


def test_h3_contrast_gene_family_level(mod):
    # The real cefotaxime finding: within-β-lactam (amp~cef) shares NO gene family
    # (TEM vs CTX-M/CMY), cross-class pairs overlap more -> H3 rejected at gene level.
    summaries = [
        {"antibiotic_a": "ampicillin", "antibiotic_b": "cefotaxime",
         "same_drug_family": True, "gene_family_jaccard": 0.0, "n_gene_family_overlap": 0},
        {"antibiotic_a": "cefotaxime", "antibiotic_b": "ciprofloxacin",
         "same_drug_family": False, "gene_family_jaccard": 0.33, "n_gene_family_overlap": 2},
    ]
    h3 = mod.h3_contrast(summaries, jaccard_key="gene_family_jaccard",
                         overlap_key="n_gene_family_overlap", p_key=None)
    assert h3["level"] == "gene_family"
    assert h3["testable"] is True
    assert h3["verdict"] == "within_not_greater"
    assert h3["within_family"]["mean_overlap"] == 0.0
    assert h3["within_family_min_p_enrichment"] is None


def test_gene_family_sets(mod, tmp_path):
    db = _kb(tmp_path, {"ampicillin": {1, 2}, "cefotaxime": {3, 4}})
    conn = sqlite3.connect(str(db))
    # ampicillin -> TEM; cefotaxime -> CTX-M + CMY (distinct β-lactamase families).
    rows = [(1, 1, "confirmed", "TEM beta-lactamase"),
            (3, 2, "confirmed", "CTX-M beta-lactamase"),
            (4, 2, "candidate", "CMY beta-lactamase"),
            (2, 1, "none", "should-be-ignored")]
    conn.executemany(
        """INSERT INTO blast_annotations(unitig_id, model_id, source_db, tier, aro_gene_family)
           VALUES (?,?,'card',?,?)""", rows)
    conn.commit()
    gf = mod.gene_family_sets(conn)
    assert gf["ampicillin"] == {"TEM beta-lactamase"}
    assert gf["cefotaxime"] == {"CTX-M beta-lactamase", "CMY beta-lactamase"}
    conn.close()


# --- synthetic KB ----------------------------------------------------------
def _kb(tmp_path, stable_by_ab):
    """Build a minimal populated KB. stable_by_ab: {antibiotic: set(unitig_id)}.

    Every referenced unitig is created; each antibiotic gets a run + model, and
    its stable unitigs are written as stable=1 gain_seed rows. Returns the path.
    """
    from lib.kb_schema import create_schema

    db = tmp_path / "amrk.db"
    conn = sqlite3.connect(str(db))
    create_schema(conn)
    all_uids = sorted({u for s in stable_by_ab.values() for u in s})
    for uid in all_uids:
        conn.execute("INSERT INTO unitigs(unitig_id, sequence, k) VALUES (?,?,?)",
                     (uid, f"SEQ{uid:04d}", 21))
    for mid, (ab, stable) in enumerate(stable_by_ab.items(), start=1):
        run_id = f"ecoli__{ab}__run"
        conn.execute("INSERT INTO pipeline_runs(run_id, organism, antibiotic) VALUES (?,?,?)",
                     (run_id, "ecoli", ab))
        conn.execute("INSERT INTO antibiotics(antibiotic, drug_class) VALUES (?,NULL)", (ab,))
        conn.execute("INSERT INTO models(model_id, run_id, antibiotic) VALUES (?,?,?)",
                     (mid, run_id, ab))
        for uid in sorted(stable):
            conn.execute(
                """INSERT INTO unitig_model_scores
                       (unitig_id, model_id, stable, selection_method)
                   VALUES (?,?,1,'gain_seed')""", (uid, mid))
    conn.commit()
    conn.close()
    return db


def test_stable_sets_and_overlap_cross_class(mod, tmp_path):
    db = _kb(tmp_path, {"ampicillin": {1, 2, 3, 4}, "ciprofloxacin": {3, 4, 5, 6, 7}})
    conn = sqlite3.connect(str(db))
    from lib.logging_utils import get_logger
    log = get_logger("test")

    classes = mod.fill_drug_classes(conn, log)
    assert classes["ampicillin"] == "penicillins"
    assert classes["ciprofloxacin"] == "quinolones"
    # drug_class actually written back to the antibiotics table.
    written = dict(conn.execute("SELECT antibiotic, drug_class FROM antibiotics").fetchall())
    assert written == {"ampicillin": "penicillins", "ciprofloxacin": "quinolones"}

    sets = mod.stable_sets(conn)
    assert sets == {"ampicillin": {1, 2, 3, 4}, "ciprofloxacin": {3, 4, 5, 6, 7}}

    summaries, union_all = mod.populate_overlap(conn, sets, classes, "ecoli", log)
    assert union_all == {1, 2, 3, 4, 5, 6, 7}
    assert len(summaries) == 1
    s = summaries[0]
    assert s["n_overlap"] == 2 and s["shared_unitig_ids"] == [3, 4]
    assert s["same_class"] is False and s["same_drug_family"] is False
    assert s["jaccard"] == pytest.approx(2 / 7)

    rows = conn.execute("SELECT unitig_id, antibiotic_a, antibiotic_b, same_class "
                        "FROM unitig_antibiotic_overlap ORDER BY unitig_id").fetchall()
    assert rows == [(3, "ampicillin", "ciprofloxacin", 0),
                    (4, "ampicillin", "ciprofloxacin", 0)]
    conn.close()


def test_same_drug_family_for_betalactam_pair(mod, tmp_path):
    # ampicillin (penicillins) + cefotaxime (cephalosporins): different registry
    # class (same_class=0) but the SAME β-lactam family (H3 within-class pair).
    db = _kb(tmp_path, {"ampicillin": {1, 2, 9}, "cefotaxime": {2, 9, 3}})
    conn = sqlite3.connect(str(db))
    from lib.logging_utils import get_logger
    log = get_logger("test")
    classes = mod.fill_drug_classes(conn, log)
    sets = mod.stable_sets(conn)
    summaries, _ = mod.populate_overlap(conn, sets, classes, "ecoli", log)
    s = summaries[0]
    assert s["same_class"] is False       # penicillins != cephalosporins
    assert s["same_drug_family"] is True   # both β-lactams -> H3-testable pair
    assert set(s["shared_unitig_ids"]) == {2, 9}
    conn.close()


def test_populate_overlap_is_idempotent(mod, tmp_path):
    db = _kb(tmp_path, {"ampicillin": {1, 2, 3}, "ciprofloxacin": {2, 3, 4}})
    conn = sqlite3.connect(str(db))
    from lib.logging_utils import get_logger
    log = get_logger("test")
    classes = mod.fill_drug_classes(conn, log)
    sets = mod.stable_sets(conn)
    mod.populate_overlap(conn, sets, classes, "ecoli", log)
    mod.populate_overlap(conn, sets, classes, "ecoli", log)  # full recompute, no dupes
    (n,) = conn.execute("SELECT COUNT(*) FROM unitig_antibiotic_overlap").fetchone()
    assert n == 2
    conn.close()
