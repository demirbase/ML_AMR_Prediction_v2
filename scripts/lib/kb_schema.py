#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AMRK-DB knowledge-base schema (SQLite, stdlib only).

The schema follows docs/ROADMAP.md §1.1. It is intentionally plain SQL via the
stdlib ``sqlite3`` module (no SQLAlchemy/ORM) so populating the KB needs no extra
dependency or container rebuild; the DDL stays portable to PostgreSQL later
(the only SQLite-isms are ``INTEGER PRIMARY KEY`` autoincrement and the pragmas,
both trivially swapped).

Design notes
------------
* ``pipeline_runs`` is the provenance anchor — every model/score/evidence row
  links back to the exact run (git commit, CARD/KMC versions, config hash, seed)
  so any KB record is reproducible (ROADMAP §1.3, must-haves M6/M10).
* ``kmers`` is the deduplicated k-mer dictionary; everything else references it.
* ``validation_evidence`` is the generic evidence ledger (M11): one row per
  BLAST / background-frequency / SNP / permutation result, each tagged with its
  ``evidence_type``, ``evidence_source`` (incl. tool+version) and ``pipeline_run_id``.
* ``kb_metadata`` is a single-row table carrying ``kb_schema_version`` + FAIR
  fields (CARD version, Zenodo DOI, license) surfaced by the API ``/metadata``.

Bump ``KB_SCHEMA_VERSION`` (semantic versioning) on any schema change.
"""

KB_SCHEMA_VERSION = "0.3.0"

# Ordered DDL — parent tables before the children that reference them.
SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

-- Provenance: one row per pipeline execution that produced KB content. -------
CREATE TABLE IF NOT EXISTS pipeline_runs (
    run_id          TEXT PRIMARY KEY,        -- {org}__{ab}__{UTC}__{git7}
    organism        TEXT NOT NULL,
    antibiotic      TEXT NOT NULL,
    git_commit      TEXT,                    -- 40-char pipeline commit
    git_dirty       INTEGER,                 -- 0/1 working tree clean?
    card_version    TEXT,                    -- e.g. 4.0.1 (BLAST annotation source)
    kmc_version     TEXT,
    xgboost_version TEXT,
    random_seed     INTEGER,
    config_hash     TEXT,                    -- hash of config.yaml used
    min_support     INTEGER,                 -- effective (adaptive) feature filter
    n_genomes       INTEGER,
    created_at      TEXT                     -- ISO-8601 UTC
);

-- Antibiotic reference (class for cross-class vs within-class analysis). -----
CREATE TABLE IF NOT EXISTS antibiotics (
    antibiotic      TEXT PRIMARY KEY,        -- canonical id (registry spelling)
    drug_class      TEXT
);

-- One trained model per run, with held-out evaluation metrics. ---------------
CREATE TABLE IF NOT EXISTS models (
    model_id        INTEGER PRIMARY KEY,
    run_id          TEXT NOT NULL REFERENCES pipeline_runs(run_id),
    antibiotic      TEXT NOT NULL REFERENCES antibiotics(antibiotic),
    n_trees         INTEGER,
    operating_threshold REAL,
    roc_auc         REAL,
    roc_auc_ci_low  REAL,
    roc_auc_ci_high REAL,
    pr_auc          REAL,
    mcc             REAL,
    balanced_accuracy REAL,
    accuracy        REAL,
    auc_mean_seeds  REAL,                    -- 07b 5-seed mean
    auc_std_seeds   REAL,                    -- 07b 5-seed std
    UNIQUE(run_id, antibiotic)
);

-- Deduplicated k-mer dictionary. --------------------------------------------
CREATE TABLE IF NOT EXISTS kmers (
    kmer_id         INTEGER PRIMARY KEY,
    sequence        TEXT NOT NULL UNIQUE,    -- canonical k-mer
    k               INTEGER NOT NULL
);

-- Per-(kmer, model) importance + stability scores. --------------------------
CREATE TABLE IF NOT EXISTS kmer_model_scores (
    kmer_id              INTEGER NOT NULL REFERENCES kmers(kmer_id),
    model_id             INTEGER NOT NULL REFERENCES models(model_id),
    gain                 REAL,               -- XGBoost Gain importance
    in_gain_topn         INTEGER,            -- 0/1 in the single-model top-N
    selection_frequency  REAL,               -- selection frequency (method below)
    stable               INTEGER,            -- 0/1 selection_frequency >= threshold
    composite_score      REAL,               -- stability * log10(1/E) * identity
    mean_abs_shap        REAL,               -- mean |TreeSHAP| (CPSS rows; step 13)
    selection_method     TEXT,               -- 'gain_seed' (07/07b) | 'cpss' (step 13)
    PRIMARY KEY (kmer_id, model_id, selection_method)
);

-- BLAST hits (CARD local + NCBI remote), one row per (kmer, db, hit). --------
CREATE TABLE IF NOT EXISTS blast_annotations (
    annotation_id   INTEGER PRIMARY KEY,
    kmer_id         INTEGER NOT NULL REFERENCES kmers(kmer_id),
    model_id        INTEGER REFERENCES models(model_id),
    source_db       TEXT NOT NULL,           -- 'card' | 'ncbi'
    gene_symbol     TEXT,
    description     TEXT,
    identity_pct    REAL,
    coverage        REAL,                    -- alignment length / k
    evalue          REAL,
    tier            TEXT,                    -- confirmed | candidate | weak | none
    -- ARO/CARD ontology mapping (M16) — populated for CARD hits from 09's
    -- aro_index/card.json lookup; NULL for NCBI hits or unmapped CARD hits.
    aro_accession            TEXT,
    aro_gene_family          TEXT,
    aro_drug_class           TEXT,
    aro_resistance_mechanism TEXT
);

-- Resistant-vs-susceptible prevalence / discriminativeness (step 10). --------
CREATE TABLE IF NOT EXISTS kmer_background_frequency (
    kmer_id         INTEGER NOT NULL REFERENCES kmers(kmer_id),
    model_id        INTEGER NOT NULL REFERENCES models(model_id),
    prevalence_resistant   REAL,
    prevalence_susceptible REAL,
    prevalence_overall     REAL,
    delta_prevalence       REAL,
    odds_ratio             REAL,
    fisher_p               REAL,
    discriminative         INTEGER,          -- 0/1 |delta|>=min_delta AND p<alpha
    PRIMARY KEY (kmer_id, model_id)
);

-- CARD variant-model SNP allele check (step 11). ----------------------------
CREATE TABLE IF NOT EXISTS variant_snp_check (
    kmer_id         INTEGER NOT NULL REFERENCES kmers(kmer_id),
    model_id        INTEGER REFERENCES models(model_id),
    card_model      TEXT,
    snp             TEXT,                    -- e.g. S83L
    allele_class    TEXT,                    -- resistant_allele | wildtype | other | ambiguous
    PRIMARY KEY (kmer_id, model_id, card_model, snp)
);

-- Cross-antibiotic stable-k-mer overlap (step S1 / H3). ---------------------
CREATE TABLE IF NOT EXISTS kmer_antibiotic_overlap (
    kmer_id         INTEGER NOT NULL REFERENCES kmers(kmer_id),
    antibiotic_a    TEXT NOT NULL,
    antibiotic_b    TEXT NOT NULL,
    same_class      INTEGER,                 -- 0/1 within-class pair?
    PRIMARY KEY (kmer_id, antibiotic_a, antibiotic_b)
);

-- Generic evidence ledger — every validation result, fully attributed (M11). -
CREATE TABLE IF NOT EXISTS validation_evidence (
    evidence_id      INTEGER PRIMARY KEY,
    kmer_id          INTEGER REFERENCES kmers(kmer_id),
    evidence_type    TEXT NOT NULL,          -- blast | background_frequency | snp | permutation | temporal
    evidence_source  TEXT NOT NULL,          -- e.g. 'CARD 4.0.1', 'shuffled_labels', 'NCBI nt'
    evidence_score   REAL,                   -- E-value / delta-AUC / Fisher p ...
    pipeline_run_id  TEXT REFERENCES pipeline_runs(run_id)
);

-- Single-row KB metadata (FAIR; surfaced by API /metadata). -----------------
CREATE TABLE IF NOT EXISTS kb_metadata (
    id                INTEGER PRIMARY KEY CHECK (id = 1),
    kb_schema_version TEXT NOT NULL,
    card_version      TEXT,
    zenodo_doi        TEXT,
    license           TEXT DEFAULT 'CC-BY-4.0',
    created_at        TEXT,
    n_kmers           INTEGER,
    n_models          INTEGER
);

CREATE INDEX IF NOT EXISTS idx_kmers_sequence       ON kmers(sequence);
CREATE INDEX IF NOT EXISTS idx_blast_gene           ON blast_annotations(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_scores_stability     ON kmer_model_scores(selection_frequency);
CREATE INDEX IF NOT EXISTS idx_models_antibiotic    ON models(antibiotic);
"""


def create_schema(conn):
    """Create all tables/indexes on an open sqlite3 connection (idempotent)."""
    conn.executescript(SCHEMA_SQL)
    conn.commit()
