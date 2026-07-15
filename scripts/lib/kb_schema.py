#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AMRK-DB knowledge-base schema (SQLite, stdlib only).

The schema follows docs/ROADMAP.md §1.1. It is intentionally plain SQL via the
stdlib ``sqlite3`` module (no SQLAlchemy/ORM) so populating the KB needs no extra
dependency or container rebuild; the DDL stays portable to PostgreSQL later
(the only SQLite-isms are ``INTEGER PRIMARY KEY`` autoincrement and the pragmas,
both trivially swapped).

Feature unit = **unitig** (compacted de Bruijn graph paths from bcalm2/
unitig-caller; ROADMAP §0.1). The tables are named accordingly (``unitigs``,
``unitig_model_scores`` …); ``sequence`` is the unitig DNA and ``k`` the de
Bruijn k used to build it (21).

Design notes
------------
* ``pipeline_runs`` is the provenance anchor — every model/score/evidence row
  links back to the exact run (git commit, CARD/KMC versions, config hash, seed)
  so any KB record is reproducible (ROADMAP §1.3, must-haves M6/M10).
* ``unitigs`` is the deduplicated unitig dictionary; everything else references it.
* ``validation_evidence`` is the generic evidence ledger (M11): one row per
  BLAST / background-frequency / SNP / permutation / stability / pyseer result,
  each tagged with ``evidence_type``, ``evidence_source`` (incl. tool+version)
  and ``pipeline_run_id``.
* ``kb_metadata`` is a single-row table carrying ``kb_schema_version`` + FAIR
  fields (CARD version, Zenodo DOI, license) surfaced by the API ``/metadata``.

Bump ``KB_SCHEMA_VERSION`` (semantic versioning) on any schema change.
"""

KB_SCHEMA_VERSION = "0.7.1"

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
    -- 0.7.1: the tools the RESULTS depend on. Until now this table recorded kmc
    -- (a QC-only tool for the abandoned k-mer baseline) but not unitig-caller,
    -- which builds the features, nor PopPUNK, which defines the CV groups — so
    -- the KB could not say what produced its own lineage labels.
    unitig_caller_version TEXT,              -- builds the unitig features
    bcalm_version         TEXT,              -- compacted de Bruijn graph
    poppunk_version       TEXT,              -- defines the lineage-CV groups
    -- NOT redundant with poppunk_version: pinning PopPUNK does NOT pin its
    -- behaviour. Verified 2026-07-15 — a rebuild held poppunk at 2.7.8 while its
    -- network backend graph-tool went 2.98 -> 3.0, and E. coli re-clustered
    -- (324 -> 397 lineages, ARI 0.990). Different clusters = different folds =
    -- different AUCs, silently.
    graph_tool_version    TEXT,
    blast_version         TEXT,              -- BLAST+ used for the CARD/NCBI pass
    -- Reported by 14_pyseer_lmm into its own summary and read back by populate:
    -- pyseer ships in amr-tools.sif, populate runs in amr.sif, so only the step
    -- that invokes pyseer can honestly report its version.
    pyseer_version        TEXT,
    random_seed     INTEGER,
    config_hash     TEXT,                    -- data/config fingerprint
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
    auc_mean_seeds  REAL,                    -- lineage-CV / 5-seed mean
    auc_std_seeds   REAL,                    -- lineage-CV / 5-seed std
    cv_method       TEXT,                    -- lineage_group_kfold_Nfold (honest) | repeated_holdout_5seed (fallback)
    UNIQUE(run_id, antibiotic)
);

-- Deduplicated unitig dictionary. -------------------------------------------
CREATE TABLE IF NOT EXISTS unitigs (
    unitig_id       INTEGER PRIMARY KEY,
    sequence        TEXT NOT NULL UNIQUE,    -- unitig DNA sequence
    k               INTEGER NOT NULL         -- de Bruijn k used to build it
);

-- Per-(unitig, model) importance + stability scores. ------------------------
CREATE TABLE IF NOT EXISTS unitig_model_scores (
    unitig_id            INTEGER NOT NULL REFERENCES unitigs(unitig_id),
    model_id             INTEGER NOT NULL REFERENCES models(model_id),
    gain                 REAL,               -- XGBoost Gain importance
    in_gain_topn         INTEGER,            -- 0/1 in the single-model top-N
    selection_frequency  REAL,               -- selection frequency (method below)
    stable               INTEGER,            -- 0/1 selection_frequency >= threshold
    composite_score      REAL,               -- stability * log10(1/E) * identity
    mean_abs_shap        REAL,               -- mean |TreeSHAP| (CPSS rows; step 13)
    selection_method     TEXT,               -- 'gain_seed' (07/07b) | 'cpss' (step 13)
    PRIMARY KEY (unitig_id, model_id, selection_method)
);

-- BLAST hits (CARD local + NCBI remote), one row per (unitig, db, hit). ------
CREATE TABLE IF NOT EXISTS blast_annotations (
    annotation_id   INTEGER PRIMARY KEY,
    unitig_id       INTEGER NOT NULL REFERENCES unitigs(unitig_id),
    model_id        INTEGER REFERENCES models(model_id),
    source_db       TEXT NOT NULL,           -- 'card' | 'ncbi'
    gene_symbol     TEXT,
    description     TEXT,
    identity_pct    REAL,
    coverage        REAL,                    -- alignment length / unitig length
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
CREATE TABLE IF NOT EXISTS unitig_background_frequency (
    unitig_id       INTEGER NOT NULL REFERENCES unitigs(unitig_id),
    model_id        INTEGER NOT NULL REFERENCES models(model_id),
    prevalence_resistant   REAL,
    prevalence_susceptible REAL,
    prevalence_overall     REAL,
    delta_prevalence       REAL,
    odds_ratio             REAL,
    fisher_p               REAL,
    discriminative         INTEGER,          -- 0/1 |delta|>=min_delta AND p<alpha
    PRIMARY KEY (unitig_id, model_id)
);

-- CARD variant-model SNP allele check (step 11). ----------------------------
CREATE TABLE IF NOT EXISTS variant_snp_check (
    unitig_id       INTEGER NOT NULL REFERENCES unitigs(unitig_id),
    model_id        INTEGER REFERENCES models(model_id),
    card_model      TEXT,
    snp             TEXT,                    -- e.g. S83L
    allele_class    TEXT,                    -- resistant_allele | wildtype | other | ambiguous
    PRIMARY KEY (unitig_id, model_id, card_model, snp)
);

-- Cross-antibiotic stable-unitig overlap (step S1 / H3). --------------------
CREATE TABLE IF NOT EXISTS unitig_antibiotic_overlap (
    unitig_id       INTEGER NOT NULL REFERENCES unitigs(unitig_id),
    organism        TEXT NOT NULL,           -- 0.6.0: keep overlaps per-organism
    antibiotic_a    TEXT NOT NULL,
    antibiotic_b    TEXT NOT NULL,
    same_class      INTEGER,                 -- 0/1 within-class pair?
    PRIMARY KEY (unitig_id, organism, antibiotic_a, antibiotic_b)
);

-- Generic evidence ledger — every validation result, fully attributed (M11). -
CREATE TABLE IF NOT EXISTS validation_evidence (
    evidence_id      INTEGER PRIMARY KEY,
    unitig_id        INTEGER REFERENCES unitigs(unitig_id),
    evidence_type    TEXT NOT NULL,          -- blast | background_frequency | snp | permutation_mda | label_permutation | stability_selection | pyseer_lmm
    evidence_source  TEXT NOT NULL,          -- e.g. 'CARD 4.0.1', 'CPSS (B=100)', 'pyseer LMM'
    evidence_score   REAL,                   -- E-value / delta-AUC / Fisher p / selection freq ...
    pipeline_run_id  TEXT REFERENCES pipeline_runs(run_id)
);

-- Composite evidence tier (0.7.0) — one grade per (unitig, model) folding the
-- BLAST hit + the 5 statistical validation layers into a single confidence
-- level (see populate_database.classify_evidence_tier). This is ADDITIVE to and
-- independent of blast_annotations.tier (the BLAST-only layer-1 grade, kept for
-- backward-compat). Its purpose is to surface *novel* candidate biomarkers: a
-- unitig with no known CARD gene but CPSS-stable + pyseer-significant (the
-- lineage-aware backbone) grades `strong_novel` here while the BLAST-only tier
-- buries it as `none`.
CREATE TABLE IF NOT EXISTS unitig_evidence_tier (
    unitig_id          INTEGER NOT NULL REFERENCES unitigs(unitig_id),
    model_id           INTEGER NOT NULL REFERENCES models(model_id),
    evidence_tier      TEXT,     -- confirmed | strong_novel | candidate | weak | none
    n_evidence_layers  INTEGER,  -- 0..6 passed layers (blast + prevalence/snp/mda/cpss/pyseer)
    evidence_layers    TEXT,     -- csv of the passed layer names
    is_novel_candidate INTEGER,  -- 0/1: strong_novel (no known gene, CPSS+pyseer+>=1)
    PRIMARY KEY (unitig_id, model_id)
);

-- Single-row KB metadata (FAIR; surfaced by API /metadata). -----------------
CREATE TABLE IF NOT EXISTS kb_metadata (
    id                INTEGER PRIMARY KEY CHECK (id = 1),
    kb_schema_version TEXT NOT NULL,
    card_version      TEXT,
    zenodo_doi        TEXT,
    license           TEXT DEFAULT 'CC-BY-4.0',
    created_at        TEXT,
    n_unitigs         INTEGER,
    n_models          INTEGER
);

-- Organism reference (gram stain / phylum) for cross-phylum generalisation. --
-- Added 0.5.0. pipeline_runs.organism is the slug that keys here.
CREATE TABLE IF NOT EXISTS organisms (
    organism      TEXT PRIMARY KEY,       -- slug: ecoli, kpneumoniae, saureus
    display_name  TEXT,
    taxid         INTEGER,
    gram_stain    TEXT,                   -- 'negative' | 'positive'
    phylum        TEXT
);

-- External-validation concordance (M13): the model AND reference genotype tools
-- (AMRFinderPlus, ResFinder) scored vs EUCAST/CLSI phenotype on the model's
-- held-out TEST genomes (leakage-free). FDA ME/VME + Cohen's kappa + bACC.
-- Added 0.5.0; feeds the 'external validation' reviewer question directly.
CREATE TABLE IF NOT EXISTS external_concordance (
    model_id              INTEGER NOT NULL REFERENCES models(model_id),
    caller                TEXT NOT NULL,  -- 'model' | 'AMRFinderPlus' | 'ResFinder'
    reference             TEXT,           -- phenotype standard, e.g. 'EUCAST/CLSI'
    n_test                INTEGER,
    sensitivity           REAL,
    specificity           REAL,
    balanced_accuracy     REAL,
    cohen_kappa           REAL,
    major_error_rate      REAL,           -- FDA ME  (false-resistant)
    very_major_error_rate REAL,           -- FDA VME (false-susceptible)
    PRIMARY KEY (model_id, caller)
);

CREATE INDEX IF NOT EXISTS idx_unitigs_sequence    ON unitigs(sequence);
CREATE INDEX IF NOT EXISTS idx_blast_gene          ON blast_annotations(gene_symbol);
CREATE INDEX IF NOT EXISTS idx_scores_stability    ON unitig_model_scores(selection_frequency);
CREATE INDEX IF NOT EXISTS idx_models_antibiotic   ON models(antibiotic);
CREATE INDEX IF NOT EXISTS idx_evtier_model         ON unitig_evidence_tier(model_id, evidence_tier);
CREATE INDEX IF NOT EXISTS idx_evtier_novel         ON unitig_evidence_tier(is_novel_candidate);
"""


def _add_column(conn, table, col, decl):
    """Idempotent ALTER TABLE ADD COLUMN (SQLite has no ADD COLUMN IF NOT EXISTS)."""
    have = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
    if col not in have:
        conn.execute(f"ALTER TABLE {table} ADD COLUMN {col} {decl}")


def create_schema(conn):
    """Create all tables/indexes on an open sqlite3 connection (idempotent)."""
    conn.executescript(SCHEMA_SQL)
    # Additive migrations 0.4.0 -> 0.5.0: new columns on pre-existing tables
    # (CREATE TABLE IF NOT EXISTS won't alter an already-created table).
    _add_column(conn, "antibiotics", "mechanism_type", "TEXT")   # acquired | target_snp | mixed
    _add_column(conn, "antibiotics", "who_aware", "TEXT")         # Access | Watch | Reserve
    _add_column(conn, "models", "n_features", "INTEGER")          # # unitigs in the model's matrix
    # 0.6.0: unitig_antibiotic_overlap gains `organism` (the unified KB must not
    # merge e.g. gentamicin across organisms). It's a derived cache (rebuilt by
    # step 15), so if the old organism-less shape exists, recreate it empty.
    ov = {r[1] for r in conn.execute("PRAGMA table_info(unitig_antibiotic_overlap)")}
    if ov and "organism" not in ov:
        conn.execute("DROP TABLE unitig_antibiotic_overlap")
        conn.executescript(SCHEMA_SQL)   # re-creates only the dropped table (others IF NOT EXISTS)
    # 0.7.1: tool-version provenance for the tools that define the results. Added
    # here too so a pre-0.7.1 KB gains the columns instead of failing on INSERT.
    # They stay NULL for runs recorded before this landed — an honest "unknown",
    # which is the point: those rows genuinely cannot say what produced them.
    for col in ("unitig_caller_version", "bcalm_version", "poppunk_version",
                "graph_tool_version", "blast_version", "pyseer_version"):
        _add_column(conn, "pipeline_runs", col, "TEXT")
    conn.commit()


# Natural-key dedup + UNIQUE indexes (audit Issue 3/24). blast_annotations and
# validation_evidence have autoincrement PKs and were written with plain INSERTs,
# so populate_candidates + populate_cpss both writing the same unitig, or a
# re-populate, duplicated rows. The keys below distinguish genuine multi-HSP /
# multi-source rows (they differ in identity/coverage/evalue or evidence_score)
# from true duplicates (identical content), so deduping keeps every distinct hit.
_DEDUP_KEYS = [
    ("blast_annotations", "annotation_id",
     "unitig_id, model_id, source_db, gene_symbol, identity_pct, coverage, evalue"),
    ("validation_evidence", "evidence_id",
     "unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id"),
]


def ensure_unique_indexes(conn):
    """Deduplicate any legacy duplicate rows (keep the lowest id) then add UNIQUE
    indexes so future `INSERT OR IGNORE` writes cannot duplicate. Idempotent and
    safe on already-clean DBs. Must run before the populate inserts."""
    for tbl, pk, keys in _DEDUP_KEYS:
        conn.execute(f"DELETE FROM {tbl} WHERE {pk} NOT IN "
                     f"(SELECT MIN({pk}) FROM {tbl} GROUP BY {keys})")
    conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS ux_blast_natural ON "
                 "blast_annotations(unitig_id, model_id, source_db, gene_symbol, "
                 "identity_pct, coverage, evalue)")
    conn.execute("CREATE UNIQUE INDEX IF NOT EXISTS ux_evidence_natural ON "
                 "validation_evidence(unitig_id, evidence_type, evidence_source, "
                 "evidence_score, pipeline_run_id)")
    conn.commit()
