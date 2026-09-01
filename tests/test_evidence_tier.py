#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Composite evidence_tier (KB schema 0.7.0): the per-(unitig, model) grade that
folds the BLAST hit + the 5 statistical validation layers into one confidence
level, and — the point of the feature — surfaces `strong_novel` biomarkers the
BLAST-only tier hides as `none`. See populate_database.classify_evidence_tier."""

import sqlite3

import pandas as pd

from lib.kb_schema import create_schema
from populate_database import classify_evidence_tier, populate_evidence_tier


# ---- pure grading rule ----------------------------------------------------
def _tier(**kw):
    kw = {k: kw.get(k, 0) for k in ("blast_hit", "prevalence", "snp", "mda", "cpss", "pyseer")}
    return classify_evidence_tier(**kw)[0]


def test_confirmed_needs_blast_plus_two_stat():
    assert _tier(blast_hit=1, prevalence=1, cpss=1) == "confirmed"
    assert _tier(blast_hit=1, cpss=1, pyseer=1) == "confirmed"
    # a lone BLAST hit, or BLAST + only one stat layer, is NOT confirmed
    assert _tier(blast_hit=1) == "candidate"
    assert _tier(blast_hit=1, prevalence=1) == "candidate"


def test_strong_novel_is_no_gene_plus_cpss_pyseer_plus_one():
    # the flagship case: no CARD gene, lineage-aware backbone + a third layer
    assert _tier(cpss=1, pyseer=1, prevalence=1) == "strong_novel"
    assert _tier(cpss=1, pyseer=1, mda=1) == "strong_novel"
    assert _tier(cpss=1, pyseer=1, snp=1) == "strong_novel"
    # a known gene disqualifies it from *novel* (it becomes confirmed instead)
    assert _tier(blast_hit=1, cpss=1, pyseer=1, prevalence=1) == "confirmed"


def test_cpss_pyseer_pair_without_third_is_candidate_not_novel():
    t, n, layers, is_novel = classify_evidence_tier(0, 0, 0, 0, 1, 1)
    assert t == "candidate" and is_novel is False and n == 2


def test_candidate_weak_none_boundaries():
    assert _tier(prevalence=1, mda=1) == "candidate"      # >=2 stat, no gene, not the novel combo
    assert _tier(prevalence=1) == "weak"                  # exactly one stat layer
    assert _tier(cpss=1) == "weak"
    assert _tier() == "none"                              # no evidence at all


def test_layer_count_and_passed_names():
    t, n, layers, is_novel = classify_evidence_tier(1, 1, 0, 0, 1, 1)
    assert n == 4 and t == "confirmed"
    assert layers == "blast,prevalence,cpss,pyseer"       # blast first, then stat order


# ---- integration: populate_evidence_tier over a small KB ------------------
def _kb(tmp_path):
    c = sqlite3.connect(str(tmp_path / "amrk.db"))
    create_schema(c)
    c.executescript("""
        INSERT INTO pipeline_runs(run_id, organism, antibiotic) VALUES ('R1','ecoli','ampicillin');
        INSERT INTO antibiotics(antibiotic) VALUES ('ampicillin');
        INSERT INTO models(model_id, run_id, antibiotic) VALUES (1,'R1','ampicillin');
        INSERT INTO unitigs(unitig_id, sequence, k) VALUES
            (1,'AAA',21),(2,'CCC',21),(3,'GGG',21),(4,'TTT',21),(5,'AAG',21);
        -- universe: every candidate has a gain_seed score row
        INSERT INTO unitig_model_scores(unitig_id,model_id,selection_method,stable) VALUES
            (1,1,'gain_seed',0),(2,1,'gain_seed',0),(3,1,'gain_seed',0),
            (4,1,'gain_seed',0),(5,1,'gain_seed',0);
    """)
    return c


def test_populate_evidence_tier_grades_and_flags_novel(tmp_path):
    c = _kb(tmp_path)
    c.executescript("""
        -- u1: known TEM gene + prevalence + cpss  -> confirmed
        INSERT INTO blast_annotations(unitig_id,model_id,source_db,gene_symbol,tier)
            VALUES (1,1,'card','TEM-1','confirmed');
        INSERT INTO unitig_background_frequency(unitig_id,model_id,discriminative) VALUES (1,1,1);
        INSERT INTO unitig_model_scores(unitig_id,model_id,selection_method,stable) VALUES (1,1,'cpss',1);
        -- u2: NO gene, cpss + pyseer + prevalence  -> strong_novel
        INSERT INTO unitig_background_frequency(unitig_id,model_id,discriminative) VALUES (2,1,1);
        INSERT INTO unitig_model_scores(unitig_id,model_id,selection_method,stable) VALUES (2,1,'cpss',1);
        INSERT INTO validation_evidence(unitig_id,evidence_type,evidence_source,pipeline_run_id)
            VALUES (2,'pyseer_lmm','pyseer LMM','R1');
        -- u3: NO gene, cpss + pyseer only  -> candidate (no third layer)
        INSERT INTO unitig_model_scores(unitig_id,model_id,selection_method,stable) VALUES (3,1,'cpss',1);
        INSERT INTO validation_evidence(unitig_id,evidence_type,evidence_source,pipeline_run_id)
            VALUES (3,'pyseer_lmm','pyseer LMM','R1');
        -- u4: prevalence only  -> weak
        INSERT INTO unitig_background_frequency(unitig_id,model_id,discriminative) VALUES (4,1,1);
        -- u5: NO gene, cpss + pyseer + MDA (mda comes from perm_df)  -> strong_novel
        INSERT INTO unitig_model_scores(unitig_id,model_id,selection_method,stable) VALUES (5,1,'cpss',1);
        INSERT INTO validation_evidence(unitig_id,evidence_type,evidence_source,pipeline_run_id)
            VALUES (5,'pyseer_lmm','pyseer LMM','R1');
    """)
    c.commit()
    perm_df = pd.DataFrame({"kmer": ["AAG"], "permutation_significant": [1]})

    n = populate_evidence_tier(c, model_id=1, run_id="R1", perm_df=perm_df)
    assert n == 5
    got = dict(c.execute(
        "SELECT unitig_id, evidence_tier FROM unitig_evidence_tier ORDER BY unitig_id").fetchall())
    assert got == {1: "confirmed", 2: "strong_novel", 3: "candidate", 4: "weak", 5: "strong_novel"}
    novel = {r[0] for r in c.execute(
        "SELECT unitig_id FROM unitig_evidence_tier WHERE is_novel_candidate=1")}
    assert novel == {2, 5}
    # MDA layer was read from perm_df, not the KB
    layers5 = c.execute(
        "SELECT evidence_layers FROM unitig_evidence_tier WHERE unitig_id=5").fetchone()[0]
    assert "mda" in layers5 and "blast" not in layers5
    c.close()


def test_populate_evidence_tier_idempotent(tmp_path):
    c = _kb(tmp_path)
    c.execute("INSERT INTO unitig_background_frequency(unitig_id,model_id,discriminative) VALUES (1,1,1)")
    c.commit()
    populate_evidence_tier(c, 1, "R1", None)
    populate_evidence_tier(c, 1, "R1", None)   # re-populate must not duplicate (PK upsert)
    assert c.execute("SELECT COUNT(*) FROM unitig_evidence_tier WHERE unitig_id=1").fetchone()[0] == 1
    c.close()


# ---- step 11 SNP rows must join to real unitigs ---------------------------
# Regression: populate_snp used to fall back to step 11's `kmer_qseqid` column
# when no `kmer` column was present. That column holds the queried FASTA header
# (`Rank_n|Score_x|Feature_f...`), not a sequence, so every SNP row registered
# its own identifier string as a unitig. The layer then joined to nothing and
# graded 0 biomarkers while looking like it had run.
def test_attach_snp_sequences_maps_headers_to_kmers(tmp_path):
    from populate_database import _attach_snp_sequences

    d = tmp_path / "05_explainability"
    d.mkdir(parents=True)
    (d / "02_top_50_features_cipro.fasta").write_text(
        ">Rank_1|Score_9.0|Feature_f1\nACGTACGTAC\n"
        ">Rank_2|Score_8.0|Feature_f2\nTTTTGGGGCC\n", encoding="utf-8")
    df = pd.DataFrame({"kmer_qseqid": ["Rank_2|Score_8.0|Feature_f2",
                                       "Rank_9|Score_1.0|Feature_f9"],
                       "allele_class": ["resistant_allele", "wildtype"]})

    out = _attach_snp_sequences(df, tmp_path, "cipro")
    assert out.loc[0, "kmer"] == "TTTTGGGGCC"   # header resolved to its sequence
    assert pd.isna(out.loc[1, "kmer"])          # unmatched header stays unmapped


def test_populate_snp_never_registers_a_non_dna_unitig():
    from populate_database import populate_snp

    conn = sqlite3.connect(":memory:")
    create_schema(conn)
    conn.executescript("""
        INSERT INTO pipeline_runs(run_id, organism, antibiotic)
            VALUES ('R1','ecoli','ciprofloxacin');
        INSERT INTO antibiotics(antibiotic) VALUES ('ciprofloxacin');
        INSERT INTO models(model_id, run_id, antibiotic) VALUES (1,'R1','ciprofloxacin');
    """)
    snp = pd.DataFrame({"kmer": ["Rank_1|Score_9.0|Feature_f1", "ACGTACGTAC"],
                        "variant_gene": ["Ecol_gyrA_FLO"] * 2,
                        "snp": ["S83L"] * 2,
                        "allele_class": ["resistant_allele"] * 2})

    populate_snp(conn, 1, "R1", 31, snp)
    seqs = [r[0] for r in conn.execute("SELECT sequence FROM unitigs")]
    assert seqs == ["ACGTACGTAC"]   # the identifier string is skipped, not stored
