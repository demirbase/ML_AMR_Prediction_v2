#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for lib/kb_queries.py (S8/S9 API backend) against a synthetic KB, plus
an optional FastAPI smoke test (skipped if fastapi is not installed)."""

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
        INSERT INTO pipeline_runs(run_id, organism, antibiotic) VALUES ('R1','ecoli','ampicillin');
        INSERT INTO pipeline_runs(run_id, organism, antibiotic) VALUES ('R2','ecoli','cefotaxime');
        INSERT INTO antibiotics(antibiotic, drug_class) VALUES ('ampicillin','penicillins'),('cefotaxime','cephalosporins');
        INSERT INTO models(model_id, run_id, antibiotic) VALUES (1,'R1','ampicillin'),(2,'R2','cefotaxime');
        INSERT INTO unitigs(unitig_id, sequence, k) VALUES (1,'AAA',21),(2,'CCC',21);
        INSERT INTO unitig_model_scores(unitig_id,model_id,gain,selection_frequency,stable,selection_method)
            VALUES (1,1,100.0,0.8,1,'cpss'),(2,1,10.0,0.3,0,'cpss'),(1,2,90.0,0.7,1,'cpss');
        INSERT INTO blast_annotations(unitig_id,model_id,source_db,gene_symbol,identity_pct,tier,aro_gene_family)
            VALUES (1,1,'card','blaTEM-1',100.0,'confirmed','TEM beta-lactamase');
        INSERT INTO unitig_background_frequency(unitig_id,model_id,prevalence_resistant,prevalence_susceptible,delta_prevalence,discriminative)
            VALUES (1,1,0.9,0.1,0.8,1);
        INSERT INTO variant_snp_check(unitig_id,model_id,card_model,snp,allele_class)
            VALUES (1,1,'x','S83L','wildtype');
        INSERT INTO unitig_antibiotic_overlap(unitig_id,organism,antibiotic_a,antibiotic_b,same_class)
            VALUES (1,'ecoli','ampicillin','cefotaxime',0);
        INSERT INTO validation_evidence(unitig_id,evidence_type,evidence_source,evidence_score,pipeline_run_id)
            VALUES (1,'blast','CARD 4.0.1',1e-40,'R1');
        INSERT INTO kb_metadata(id,kb_schema_version,card_version,zenodo_doi,license,n_unitigs,n_models)
            VALUES (1,'0.4.0','4.0.1',NULL,'CC-BY-4.0',2,2);
    """)
    c.commit(); c.close()
    return db


def test_metadata(kb):
    c = Q.connect(kb)
    m = Q.get_metadata(c)
    assert m["kb_schema_version"] == "0.4.0"
    assert m["license"] == "CC-BY-4.0"
    assert m["antibiotics"] == ["ampicillin", "cefotaxime"]
    c.close()


def test_stats(kb):
    c = Q.connect(kb)
    s = Q.get_stats(c)
    assert s["n_unitigs"] == 2 and s["n_models"] == 2
    amp = [x for x in s["per_antibiotic"] if x["antibiotic"] == "ampicillin"][0]
    assert amp["n_stable"] == 1 and amp["n_scored"] == 2
    assert {t["tier"]: t["n"] for t in s["blast_tiers"]}["confirmed"] == 1
    c.close()


def test_list_biomarkers_filters(kb):
    c = Q.connect(kb)
    assert len(Q.list_biomarkers(c, antibiotic="ampicillin")) == 2
    assert len(Q.list_biomarkers(c, antibiotic="ampicillin", stable_only=True)) == 1
    assert len(Q.list_biomarkers(c, min_stability=0.5)) == 2   # 0.8 (amp) + 0.7 (cef)
    # unitig 1 is stable in both models; its per-unitig confirmed BLAST attaches to
    # both antibiotic rows -> 2 confirmed rows, both blaTEM-1.
    conf = Q.list_biomarkers(c, tier="confirmed")
    assert len(conf) == 2 and all(r["gene_symbol"] == "blaTEM-1" for r in conf)
    # unitig 1 / ampicillin carries the confirmed gene via the best-hit join
    amp1 = [r for r in Q.list_biomarkers(c, antibiotic="ampicillin") if r["unitig_id"] == 1][0]
    assert amp1["tier"] == "confirmed" and amp1["gene_symbol"] == "blaTEM-1"
    c.close()


def test_get_unitig_full_chain(kb):
    c = Q.connect(kb)
    rec = Q.get_unitig(c, "AAA")
    assert rec["unitig"]["sequence"] == "AAA"
    assert len(rec["model_scores"]) == 2           # amp + cef
    assert rec["blast"][0]["gene_symbol"] == "blaTEM-1"
    assert rec["background_frequency"][0]["discriminative"] == 1
    assert rec["snp"][0]["snp"] == "S83L"
    assert rec["overlap"][0]["antibiotic_b"] == "cefotaxime"
    assert rec["evidence"][0]["evidence_type"] == "blast"
    assert Q.get_unitig(c, "ZZZ") is None
    c.close()


def test_overlap_order_independent(kb):
    c = Q.connect(kb)
    a = Q.get_overlap(c, "ampicillin", "cefotaxime")
    b = Q.get_overlap(c, "cefotaxime", "ampicillin")
    assert a["n_shared"] == 1 and b["n_shared"] == 1
    assert a["shared_unitigs"][0]["sequence"] == "AAA"
    c.close()


def test_fastapi_smoke(kb, load_script):
    pytest.importorskip("fastapi")
    from fastapi.testclient import TestClient
    api = load_script("kb_api.py")
    api.DB_PATH = kb
    client = TestClient(api._create_app())
    assert client.get("/api/v1/stats").json()["n_unitigs"] == 2
    assert client.get("/api/v1/metadata").json()["kb_schema_version"] == "0.4.0"
    r = client.get("/api/v1/kmers", params={"antibiotic": "ampicillin", "stable_only": True})
    assert r.json()["count"] == 1
    assert client.get("/api/v1/kmers/AAA").json()["blast"][0]["gene_symbol"] == "blaTEM-1"
    assert client.get("/api/v1/kmers/ZZZ").status_code == 404
