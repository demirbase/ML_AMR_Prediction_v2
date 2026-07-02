#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Smoke test for scripts/kb_report.py build_report against a synthetic KB."""

import sqlite3

import pytest

from lib import kb_queries as Q
from lib.kb_schema import create_schema


@pytest.fixture
def kb(tmp_path):
    db = tmp_path / "amrk.db"
    c = sqlite3.connect(str(db))
    create_schema(c)
    c.executescript("""
        INSERT INTO pipeline_runs(run_id, organism, antibiotic) VALUES ('R1','ecoli','ampicillin'),('R2','ecoli','cefotaxime');
        INSERT INTO antibiotics(antibiotic, drug_class) VALUES ('ampicillin','penicillins'),('cefotaxime','cephalosporins');
        INSERT INTO models(model_id, run_id, antibiotic, roc_auc, mcc, n_trees, auc_mean_seeds, auc_std_seeds)
            VALUES (1,'R1','ampicillin',0.924,0.72,146,0.951,0.011),(2,'R2','cefotaxime',0.969,0.844,45,0.955,0.020);
        INSERT INTO unitigs(unitig_id, sequence, k) VALUES (1,'AAA',21),(2,'CCC',21);
        INSERT INTO unitig_model_scores(unitig_id,model_id,stable,selection_method) VALUES (1,1,1,'cpss'),(1,2,1,'cpss');
        INSERT INTO blast_annotations(unitig_id,model_id,source_db,gene_symbol,tier,aro_gene_family)
            VALUES (1,1,'card','blaTEM-1','confirmed','TEM beta-lactamase');
        INSERT INTO validation_evidence(unitig_id,evidence_type,evidence_source,evidence_score,pipeline_run_id)
            VALUES (NULL,'concordance_resfinder','ResFinder 4.5.0 vs EUCAST/CLSI (bACC=0.934, kappa=0.859, n=4446)',0.859,'R1'),
                   (NULL,'head_to_head_model','unitig model vs tools on held-out test (bACC=0.873, kappa=0.707, n=800)',0.707,'R1'),
                   (1,'pyseer_lmm','pyseer LMM',1e-9,'R1');
        INSERT INTO kb_metadata(id,kb_schema_version,card_version,license,n_unitigs,n_models)
            VALUES (1,'0.4.0','4.0.1','CC-BY-4.0',2,2);
    """)
    c.commit(); c.close()
    return db


@pytest.fixture
def mod(load_script):
    return load_script("kb_report.py")


def test_build_report_contents(mod, kb):
    c = Q.connect(kb)
    md = mod.build_report(c)
    c.close()
    assert "# AMRK-DB — results summary" in md
    assert "schema** 0.4.0" in md and "CC-BY-4.0" in md
    # per-antibiotic performance
    assert "0.951±0.011" in md and "TEM beta-lactamase" in md
    assert "ampicillin" in md and "cefotaxime" in md
    # validation: resfinder kappa + pyseer count for ampicillin, head-to-head shown
    assert "0.859" in md and "bACC=0.873" in md
    # overlap: within-β-lactam amp~cef enumerated even with 0 shared
    assert "ampicillin ~ cefotaxime" in md
