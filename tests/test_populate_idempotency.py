#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audit Issue 3/24: populate must not duplicate blast_annotations /
validation_evidence (candidates+cpss double-insert or re-populate)."""

import sqlite3

from lib.kb_schema import create_schema, ensure_unique_indexes


def _kb(tmp_path):
    c = sqlite3.connect(str(tmp_path / "amrk.db"))
    create_schema(c)
    c.executescript("""
        INSERT INTO pipeline_runs(run_id, organism, antibiotic) VALUES ('R1','ecoli','ampicillin');
        INSERT INTO antibiotics(antibiotic) VALUES ('ampicillin');
        INSERT INTO models(model_id, run_id, antibiotic) VALUES (1,'R1','ampicillin');
        INSERT INTO unitigs(unitig_id, sequence, k) VALUES (1,'AAA',21),(2,'CCC',21);
    """)
    return c


def test_dedup_keeps_distinct_and_removes_true_dups(tmp_path):
    c = _kb(tmp_path)
    # unitig 1: two IDENTICAL card/TEM-1 rows (the candidates+cpss double-insert)
    # + one genuine multi-HSP (same gene, different identity/evalue) that must survive.
    c.executescript("""
        INSERT INTO blast_annotations(unitig_id,model_id,source_db,gene_symbol,identity_pct,coverage,evalue,tier)
            VALUES (1,1,'card','TEM-1',100.0,1.0,1e-40,'confirmed');
        INSERT INTO blast_annotations(unitig_id,model_id,source_db,gene_symbol,identity_pct,coverage,evalue,tier)
            VALUES (1,1,'card','TEM-1',100.0,1.0,1e-40,'confirmed');
        INSERT INTO blast_annotations(unitig_id,model_id,source_db,gene_symbol,identity_pct,coverage,evalue,tier)
            VALUES (1,1,'card','TEM-1',95.0,0.8,1e-20,'candidate');
        INSERT INTO validation_evidence(unitig_id,evidence_type,evidence_source,evidence_score,pipeline_run_id)
            VALUES (1,'blast','CARD 4.0.1',1e-40,'R1');
        INSERT INTO validation_evidence(unitig_id,evidence_type,evidence_source,evidence_score,pipeline_run_id)
            VALUES (1,'blast','CARD 4.0.1',1e-40,'R1');
    """)
    c.commit()
    assert c.execute("SELECT COUNT(*) FROM blast_annotations").fetchone()[0] == 3
    assert c.execute("SELECT COUNT(*) FROM validation_evidence").fetchone()[0] == 2

    ensure_unique_indexes(c)
    # true dups collapsed; the distinct multi-HSP row survives
    assert c.execute("SELECT COUNT(*) FROM blast_annotations").fetchone()[0] == 2
    assert c.execute("SELECT COUNT(*) FROM validation_evidence").fetchone()[0] == 1

    # idempotent on an already-clean DB
    ensure_unique_indexes(c)
    assert c.execute("SELECT COUNT(*) FROM blast_annotations").fetchone()[0] == 2
    c.close()


def test_unique_index_blocks_future_dups(tmp_path):
    c = _kb(tmp_path)
    ensure_unique_indexes(c)   # creates the UNIQUE indexes on an empty DB
    for _ in range(2):         # simulate candidates then cpss inserting the same row
        c.execute("INSERT OR IGNORE INTO blast_annotations"
                  "(unitig_id,model_id,source_db,gene_symbol,identity_pct,coverage,evalue,tier)"
                  " VALUES (1,1,'card','TEM-1',100.0,1.0,1e-40,'confirmed')")
        c.execute("INSERT OR IGNORE INTO validation_evidence"
                  "(unitig_id,evidence_type,evidence_source,evidence_score,pipeline_run_id)"
                  " VALUES (1,'blast','CARD 4.0.1',1e-40,'R1')")
    c.commit()
    assert c.execute("SELECT COUNT(*) FROM blast_annotations").fetchone()[0] == 1
    assert c.execute("SELECT COUNT(*) FROM validation_evidence").fetchone()[0] == 1
    c.close()


def test_null_unitig_evidence_not_collapsed(tmp_path):
    # concordance rows have unitig_id NULL + distinct evidence_type -> must all survive
    c = _kb(tmp_path)
    ensure_unique_indexes(c)
    for et in ("concordance_amrfinderplus", "concordance_resfinder", "head_to_head_model"):
        c.execute("INSERT OR IGNORE INTO validation_evidence"
                  "(unitig_id,evidence_type,evidence_source,evidence_score,pipeline_run_id)"
                  " VALUES (NULL,?,?,?,'R1')", (et, f"{et} src", 0.8))
    c.commit()
    assert c.execute("SELECT COUNT(*) FROM validation_evidence").fetchone()[0] == 3
    c.close()
