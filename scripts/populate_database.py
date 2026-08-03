#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Populate the AMRK-DB knowledge base from pipeline outputs (M8).

Reads the per-antibiotic artefacts produced by the pipeline and loads them into
a single SQLite knowledge base (schema: scripts/lib/kb_schema.py). Idempotent —
re-running an antibiotic refreshes its rows. One DB can hold many antibiotics, so
this is the substrate for the cross-antibiotic analysis (S1/H3) and the API (S8).

Inputs (globbed under results/models/runs; any missing input is skipped with a
note, so the KB can be built from a partial pipeline):
  runs/.../run_metadata.json               -> pipeline_runs
  models/.../manifest.json + 06 metrics    -> models
  results/.../10_repeated_holdout_summary  -> models (07b 5-seed AUC mean/std)
  results/.../07_kb_candidates_{ab}.csv     \\
  results/.../10_kmer_background_frequency   } -> unitigs, unitig_model_scores,
  results/.../11_variant_snp_check           /    blast_annotations,
                                                  unitig_background_frequency,
                                                  variant_snp_check, validation_evidence

Usage:
  python scripts/populate_database.py --organism ecoli --antibiotic ampicillin
  python scripts/populate_database.py            # organism/antibiotic from config
"""

import argparse
import datetime
import glob
import json
import os
import sqlite3
import sys
from pathlib import Path

import pandas as pd

from lib.config import load_config, resolve_path, get_target
from lib.kb_schema import KB_SCHEMA_VERSION, create_schema, ensure_unique_indexes
from lib.registry import (
    antibiotic_to_class,
    antibiotic_mechanism_type,
    antibiotic_who_aware,
    load_organisms,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _find(root, filename):
    """Return the most recent file named `filename` anywhere under `root`."""
    hits = sorted(Path(root).rglob(filename), key=lambda p: p.stat().st_mtime)
    return hits[-1] if hits else None


def _read_json(path):
    if path and Path(path).exists():
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return None


def _read_csv(path):
    if path and Path(path).exists():
        try:
            return pd.read_csv(path, encoding="utf-8")
        except Exception as e:
            print(f"  ⚠ could not read {path}: {e}")
    return None


def _f(x):
    """float-or-None (handles NaN / blanks)."""
    try:
        v = float(x)
        return None if v != v else v  # drop NaN
    except (TypeError, ValueError):
        return None


def _i(x):
    v = _f(x)
    return int(v) if v is not None else None


def _b(x):
    return 1 if str(x).strip().lower() in ("1", "true", "yes") else 0


def _s(x):
    """clean-string-or-None (drops NaN / blanks / literal 'nan')."""
    if x is None:
        return None
    s = str(x).strip()
    return None if s == "" or s.lower() == "nan" else s


# ---------------------------------------------------------------------------
# Per-table loaders
# ---------------------------------------------------------------------------
def unitig_id(conn, sequence, k):
    """INSERT-OR-IGNORE a k-mer and return its id (dedup on sequence)."""
    conn.execute("INSERT OR IGNORE INTO unitigs(sequence, k) VALUES (?,?)",
                 (sequence, k))
    row = conn.execute("SELECT unitig_id FROM unitigs WHERE sequence=?",
                       (sequence,)).fetchone()
    return row[0]


def populate_run(conn, organism, antibiotic, run_meta, card_version, min_support,
                 pyseer_version=None):
    """Insert/replace the pipeline_runs row; return run_id.

    ``pyseer_version`` comes from 14's summary rather than run_metadata: pyseer
    ships in amr-tools.sif while this runs in amr.sif, so collect_versions cannot
    see it and only the step that invoked it can report it honestly.
    """
    rm = run_meta or {}
    run_id = rm.get("run_id") or f"{organism}__{antibiotic}__unknown"
    versions = rm.get("versions", {}) if isinstance(rm.get("versions"), dict) else {}
    # Keys must match lib.run_metadata.build_run_metadata: git_commit_hash,
    # random_seed, started_at, data_fingerprint (a dict). config_hash column gets
    # the fingerprint's sha256 (or the whole dict JSON-encoded — SQLite can't bind
    # a dict directly, which crashed once run_metadata.json actually existed).
    fp = rm.get("data_fingerprint")
    if isinstance(fp, dict):
        config_hash = fp.get("sha256") if isinstance(fp.get("sha256"), str) \
            else json.dumps(fp, sort_keys=True)
    else:
        config_hash = rm.get("config_hash")
        if config_hash is not None and not isinstance(config_hash, str):
            config_hash = json.dumps(config_hash, sort_keys=True)
    # 0.7.1: versions of the tools the results actually depend on. `versions` comes
    # from lib.run_metadata.collect_versions, captured inside amr.sif — so pyseer
    # is absent there (it lives in amr-tools.sif) and is passed in separately from
    # 14's own summary, the only place that saw the binary that ran.
    conn.execute(
        """INSERT OR REPLACE INTO pipeline_runs
           (run_id, organism, antibiotic, git_commit, git_dirty, card_version,
            kmc_version, xgboost_version, unitig_caller_version, bcalm_version,
            poppunk_version, graph_tool_version, blast_version, pyseer_version,
            random_seed, config_hash, min_support, n_genomes, created_at)
           VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        (run_id, organism, antibiotic, rm.get("git_commit_hash"), _i(rm.get("git_dirty")),
         card_version, versions.get("kmc"), versions.get("xgboost"),
         versions.get("unitig_caller"), versions.get("bcalm"),
         versions.get("poppunk"), versions.get("graph_tool"), versions.get("blastn"),
         pyseer_version or versions.get("pyseer"),
         _i(rm.get("random_seed")), config_hash,
         _i(min_support), _i(rm.get("n_genomes")),
         rm.get("started_at") or datetime.datetime.now(datetime.timezone.utc).isoformat()),
    )
    return run_id


def populate_model(conn, run_id, antibiotic, drug_class, manifest, metrics, holdout):
    """Insert/replace the models row; return model_id."""
    conn.execute("INSERT OR IGNORE INTO antibiotics(antibiotic, drug_class) VALUES (?,?)",
                 (antibiotic, drug_class))
    m = (metrics or {}).get("metrics", metrics or {})
    ci = m.get("roc_auc_ci") or [None, None]
    auc_mean = auc_std = None
    cv_method = None
    if holdout is not None and "seed" in holdout.columns:
        idx = holdout.set_index("seed")["roc_auc"]
        auc_mean = _f(idx.get("MEAN"))
        auc_std = _f(idx.get("STD"))
        if "cv_method" in holdout.columns and len(holdout):
            cv_method = str(holdout["cv_method"].iloc[0])  # constant per model (07b)
    man = manifest or {}
    conn.execute(
        """INSERT OR REPLACE INTO models
           (run_id, antibiotic, n_trees, operating_threshold, roc_auc,
            roc_auc_ci_low, roc_auc_ci_high, pr_auc, mcc, balanced_accuracy,
            accuracy, auc_mean_seeds, auc_std_seeds, cv_method)
           VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
        (run_id, antibiotic, _i(man.get("n_trees")),
         _f((metrics or {}).get("operating_threshold") or man.get("threshold")),
         _f(m.get("roc_auc")), _f(ci[0]), _f(ci[1]), _f(m.get("pr_auc")),
         _f(m.get("mcc")), _f(m.get("balanced_accuracy")), _f(m.get("accuracy")),
         auc_mean, auc_std, cv_method),
    )
    return conn.execute("SELECT model_id FROM models WHERE run_id=? AND antibiotic=?",
                        (run_id, antibiotic)).fetchone()[0]


def populate_candidates(conn, model_id, run_id, k, cand_df, card_version):
    """Load per-k-mer scores + BLAST + background-frequency from the candidate
    table (10's output is a superset of 09's; falls back to 09)."""
    if cand_df is None or cand_df.empty:
        print("  ⚠ no candidate table (09/10) — skipping k-mer rows.")
        return 0
    has_bg = "discriminative" in cand_df.columns
    n = 0
    for _, r in cand_df.iterrows():
        seq = str(r.get("kmer", "")).strip()
        if not seq:
            continue
        kid = unitig_id(conn, seq, k)
        conn.execute(
            """INSERT OR REPLACE INTO unitig_model_scores
               (unitig_id, model_id, gain, in_gain_topn, selection_frequency,
                stable, composite_score, mean_abs_shap, selection_method)
               VALUES (?,?,?,?,?,?,?,?,?)""",
            (kid, model_id, _f(r.get("gain_score")), _b(r.get("in_gain_topN", 1)),
             _f(r.get("selection_frequency")), _b(r.get("stable")),
             _f(r.get("composite_score")), None, "gain_seed"),
        )
        # CARD BLAST annotation (best hit recorded in the candidate row)
        if str(r.get("card_gene", "")).strip():
            conn.execute(
                """INSERT OR IGNORE INTO blast_annotations
                   (unitig_id, model_id, source_db, gene_symbol, identity_pct,
                    evalue, tier, aro_accession, aro_gene_family, aro_drug_class,
                    aro_resistance_mechanism)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?)""",
                (kid, model_id, "card", str(r.get("card_gene")),
                 _f(r.get("card_identity")), _f(r.get("card_evalue")),
                 str(r.get("confidence_tier", "none")),
                 _s(r.get("aro_accession")), _s(r.get("aro_gene_family")),
                 _s(r.get("aro_drug_class")), _s(r.get("aro_resistance_mechanism"))),
            )
            conn.execute(
                """INSERT OR IGNORE INTO validation_evidence
                   (unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id)
                   VALUES (?,?,?,?,?)""",
                (kid, "blast", f"CARD {card_version}", _f(r.get("card_evalue")), run_id),
            )
        # Background frequency / discriminativeness (only present in 10's output)
        if has_bg:
            # Derive the prevalence gap when the upstream CSV does not carry it. Step 10
            # computed it internally but did not emit the column until 2026-08, so every
            # KB written before that has delta_prevalence NULL; deriving it here means a
            # re-populate fixes old runs without depending on which version of 10 ran.
            _dp = _f(r.get("delta_prevalence"))
            if _dp is None:
                _pr, _ps = _f(r.get("prevalence_resistant")), _f(r.get("prevalence_susceptible"))
                _dp = (_pr - _ps) if (_pr is not None and _ps is not None) else None
            conn.execute(
                """INSERT OR REPLACE INTO unitig_background_frequency
                   (unitig_id, model_id, prevalence_resistant, prevalence_susceptible,
                    prevalence_overall, delta_prevalence, odds_ratio, fisher_p,
                    discriminative) VALUES (?,?,?,?,?,?,?,?,?)""",
                (kid, model_id, _f(r.get("prevalence_resistant")),
                 _f(r.get("prevalence_susceptible")), _f(r.get("prevalence_overall")),
                 _dp, _f(r.get("odds_ratio")),
                 _f(r.get("fisher_p")), _b(r.get("discriminative"))),
            )
            conn.execute(
                """INSERT OR IGNORE INTO validation_evidence
                   (unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id)
                   VALUES (?,?,?,?,?)""",
                (kid, "background_frequency", "R-vs-S Fisher exact",
                 _f(r.get("fisher_p")), run_id),
            )
        n += 1
    return n


def populate_snp(conn, model_id, run_id, k, snp_df):
    """Load CARD variant-model SNP allele calls (step 11)."""
    if snp_df is None or snp_df.empty:
        return 0
    n = 0
    for _, r in snp_df.iterrows():
        seq = str(r.get("kmer", r.get("kmer_qseqid", ""))).strip()
        if not seq:
            continue
        kid = unitig_id(conn, seq, k)
        conn.execute(
            """INSERT OR REPLACE INTO variant_snp_check
               (unitig_id, model_id, card_model, snp, allele_class)
               VALUES (?,?,?,?,?)""",
            (kid, model_id, str(r.get("variant_gene", "")), str(r.get("snp", "")),
             str(r.get("allele_class", ""))),
        )
        conn.execute(
            """INSERT OR IGNORE INTO validation_evidence
               (unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id)
               VALUES (?,?,?,?,?)""",
            (kid, "snp", "CARD variant model", None, run_id),
        )
        n += 1
    return n


def populate_cpss(conn, model_id, run_id, k, cpss_df):
    """Load the CPSS-stable, CARD-annotated unitigs (steps 13/13b) into the KB:
    unitig_model_scores (selection_method='cpss', CPSS selection_frequency +
    mean|SHAP|), CARD blast_annotations (tier + coverage + ARO), and a
    stability_selection evidence row each."""
    if cpss_df is None or cpss_df.empty:
        return 0
    n = 0
    for _, r in cpss_df.iterrows():
        seq = str(r.get("kmer", "")).strip()
        if not seq:
            continue
        kid = unitig_id(conn, seq, k)
        conn.execute(
            """INSERT OR REPLACE INTO unitig_model_scores
               (unitig_id, model_id, gain, in_gain_topn, selection_frequency,
                stable, composite_score, mean_abs_shap, selection_method)
               VALUES (?,?,?,?,?,?,?,?,?)""",
            (kid, model_id, None, 0, _f(r.get("selection_frequency")),
             _b(r.get("stable")), _f(r.get("composite_score")),
             _f(r.get("mean_abs_shap")), "cpss"),
        )
        if str(r.get("card_gene", "")).strip():
            conn.execute(
                """INSERT OR IGNORE INTO blast_annotations
                   (unitig_id, model_id, source_db, gene_symbol, identity_pct,
                    coverage, evalue, tier, aro_accession, aro_gene_family,
                    aro_drug_class, aro_resistance_mechanism)
                   VALUES (?,?,?,?,?,?,?,?,?,?,?,?)""",
                (kid, model_id, "card", str(r.get("card_gene")),
                 _f(r.get("card_identity")), _f(r.get("coverage")),
                 _f(r.get("card_evalue")), str(r.get("confidence_tier", "none")),
                 _s(r.get("aro_accession")), _s(r.get("aro_gene_family")),
                 _s(r.get("aro_drug_class")), _s(r.get("aro_resistance_mechanism"))),
            )
        conn.execute(
            """INSERT OR IGNORE INTO validation_evidence
               (unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id)
               VALUES (?,?,?,?,?)""",
            (kid, "stability_selection", "CPSS (B=100, pi>=0.6, PFER-bounded)",
             _f(r.get("selection_frequency")), run_id),
        )
        n += 1
    return n


def populate_pyseer(conn, run_id, sig_df, threshold):
    """pyseer LMM (lineage-corrected) significance (step 14) -> validation_evidence,
    only for unitigs ALREADY in the KB (don't create bare kmers for genome-wide
    hits with no model/annotation)."""
    if sig_df is None or sig_df.empty:
        return 0
    pcol = "lrt-pvalue" if "lrt-pvalue" in sig_df.columns else "filter-pvalue"
    src = f"pyseer LMM lineage-corrected (Bonferroni {threshold:.2e})" if threshold else \
          "pyseer LMM lineage-corrected"
    n = 0
    for _, r in sig_df.iterrows():
        seq = str(r.get("variant", "")).strip()
        if not seq:
            continue
        row = conn.execute("SELECT unitig_id FROM unitigs WHERE sequence=?", (seq,)).fetchone()
        if not row:
            continue
        conn.execute(
            """INSERT OR IGNORE INTO validation_evidence
               (unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id)
               VALUES (?,?,?,?,?)""",
            (row[0], "pyseer_lmm", src, _f(r.get(pcol)), run_id),
        )
        n += 1
    return n


def populate_permutation(conn, run_id, k, perm_df, labelperm):
    """Permutation significance (step 12 / 12b) -> validation_evidence.

    Per-candidate MDA (test ROC-AUC drop) rows + one model-level label-permutation
    null row (unitig_id NULL, evidence_score = empirical p). Both are evidence, not
    per-kmer scores, so they live in the generic evidence ledger."""
    n = 0
    if perm_df is not None and not perm_df.empty:
        for _, r in perm_df.iterrows():
            seq = str(r.get("kmer", "")).strip()
            if not seq:
                continue
            kid = unitig_id(conn, seq, k)
            conn.execute(
                """INSERT OR IGNORE INTO validation_evidence
                   (unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id)
                   VALUES (?,?,?,?,?)""",
                (kid, "permutation_mda", "MDA test ROC-AUC drop (100 perms, BH-FDR)",
                 _f(r.get("mda_auc_drop")), run_id),
            )
            n += 1
    if labelperm:
        src = (f"label-shuffle null (N={labelperm.get('n_permutations')}, "
               f"real_auc={labelperm.get('real_roc_auc')}, "
               f"null_max={labelperm.get('null_auc_max')})")
        conn.execute(
            """INSERT OR IGNORE INTO validation_evidence
               (unitig_id, evidence_type, evidence_source, evidence_score, pipeline_run_id)
               VALUES (?,?,?,?,?)""",
            (None, "label_permutation", src, _f(labelperm.get("empirical_p")), run_id),
        )
        n += 1
    return n


# --- Composite evidence tier (0.7.0) ---------------------------------------
# Two of the five statistical layers — CPSS stability + pyseer LMM — are the
# lineage-aware, confounder-robust "novelty backbone" (literature E2/E3): a
# unitig passing BOTH with no known CARD gene is a strong *novel* candidate
# biomarker that the BLAST-only tier would hide as `none`.
def classify_evidence_tier(blast_hit, prevalence, snp, mda, cpss, pyseer):
    """Fold the BLAST hit + the 5 statistical validation layers into one
    per-(unitig, model) confidence grade. Returns
    ``(tier, n_layers, passed_layers_csv, is_novel)``.

    Grades (thesis Methods; HANDOFF §0.-5 item 0):
      confirmed     — known CARD gene (BLAST confirmed/candidate) + >=2 stat layers
      strong_novel  — NO known gene, but CPSS *and* pyseer + >=1 other stat layer
                      (surfaces novel biomarkers the BLAST-only tier hides)
      candidate     — a lone BLAST hit, OR >=2 stat layers (not the novel combo)
      weak          — exactly one statistical layer, no BLAST hit
      none          — no evidence
    """
    blast_hit = bool(blast_hit)
    stat = (("prevalence", prevalence), ("snp", snp), ("mda", mda),
            ("cpss", cpss), ("pyseer", pyseer))
    passed = (["blast"] if blast_hit else []) + [name for name, hit in stat if hit]
    n_stat = sum(1 for _, hit in stat if hit)
    n_layers = n_stat + (1 if blast_hit else 0)
    is_novel = bool((not blast_hit) and cpss and pyseer and (prevalence or mda or snp))
    if blast_hit and n_stat >= 2:
        tier = "confirmed"
    elif is_novel:
        tier = "strong_novel"
    elif blast_hit or n_stat >= 2:
        tier = "candidate"
    elif n_stat == 1:
        tier = "weak"
    else:
        tier = "none"
    return tier, n_layers, ",".join(passed), is_novel


def populate_evidence_tier(conn, model_id, run_id, perm_df):
    """Derive + load the composite ``evidence_tier`` for every unitig of ONE
    model. Reads the 5 KB-stored layers (blast / prevalence / snp / cpss /
    pyseer) and takes MDA significance from step 12's ``permutation_significant``
    column (the per-pass MDA flag is not otherwise stored in the KB —
    validation_evidence keeps the effect size ``mda_auc_drop`` for kb_tables)."""
    mda_sig = set()
    if perm_df is not None and not perm_df.empty and "permutation_significant" in perm_df.columns:
        for _, r in perm_df.iterrows():
            if _b(r.get("permutation_significant")):
                seq = str(r.get("kmer", "")).strip()
                if seq:
                    mda_sig.add(seq)
    rows = conn.execute(
        """SELECT s.unitig_id AS uid, u.sequence AS seq,
             MAX(CASE WHEN ba.tier IN ('confirmed','candidate') THEN 1 ELSE 0 END) AS blast_hit,
             MAX(COALESCE(bf.discriminative,0)) AS prevalence,
             MAX(CASE WHEN vs.allele_class='resistant_allele' THEN 1 ELSE 0 END) AS snp,
             MAX(CASE WHEN s.selection_method='cpss' AND s.stable=1 THEN 1 ELSE 0 END) AS cpss,
             MAX(CASE WHEN pe.unitig_id IS NOT NULL THEN 1 ELSE 0 END) AS pyseer
           FROM unitig_model_scores s
           JOIN unitigs u ON u.unitig_id = s.unitig_id
           LEFT JOIN blast_annotations ba
                  ON ba.unitig_id = s.unitig_id AND ba.model_id = s.model_id
           LEFT JOIN unitig_background_frequency bf
                  ON bf.unitig_id = s.unitig_id AND bf.model_id = s.model_id
           LEFT JOIN variant_snp_check vs
                  ON vs.unitig_id = s.unitig_id AND vs.model_id = s.model_id
           LEFT JOIN validation_evidence pe
                  ON pe.unitig_id = s.unitig_id AND pe.pipeline_run_id = ?
                     AND pe.evidence_type = 'pyseer_lmm'
           WHERE s.model_id = ?
           GROUP BY s.unitig_id, u.sequence""",
        (run_id, model_id),
    ).fetchall()
    n = 0
    for uid, seq, blast_hit, prevalence, snp, cpss, pyseer in rows:
        mda = 1 if seq in mda_sig else 0
        tier, n_layers, layers, is_novel = classify_evidence_tier(
            blast_hit, prevalence, snp, mda, cpss, pyseer)
        conn.execute(
            """INSERT OR REPLACE INTO unitig_evidence_tier
               (unitig_id, model_id, evidence_tier, n_evidence_layers,
                evidence_layers, is_novel_candidate)
               VALUES (?,?,?,?,?,?)""",
            (uid, model_id, tier, n_layers, layers, 1 if is_novel else 0),
        )
        n += 1
    return n


def update_metadata(conn, card_version):
    n_unitigs = conn.execute("SELECT COUNT(*) FROM unitigs").fetchone()[0]
    n_models = conn.execute("SELECT COUNT(*) FROM models").fetchone()[0]
    # Zenodo DOI (M10/FAIR): an AMR_ZENODO_DOI env override wins (mirrors
    # AMR_CARD_VERSION so a release can stamp the DOI without editing code);
    # otherwise PRESERVE any DOI already in the row — re-populating a new
    # antibiotic must not wipe the release DOI.
    prev = conn.execute("SELECT zenodo_doi FROM kb_metadata WHERE id = 1").fetchone()
    zenodo_doi = os.environ.get("AMR_ZENODO_DOI") or (prev[0] if prev else None)
    conn.execute(
        """INSERT OR REPLACE INTO kb_metadata
           (id, kb_schema_version, card_version, zenodo_doi, license, created_at,
            n_unitigs, n_models)
           VALUES (1,?,?,?,?,?,?,?)""",
        (KB_SCHEMA_VERSION, card_version, zenodo_doi, "CC-BY-4.0",
         datetime.datetime.now(datetime.timezone.utc).isoformat(), n_unitigs, n_models),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# 0.5.0 reference/meta populators.
# ---------------------------------------------------------------------------
def populate_organisms(conn):
    """Populate the organism reference table from the registry (single source of
    truth — organisms.yaml carries display_name / taxid / gram_stain / phylum)."""
    for slug, blk in load_organisms().items():
        conn.execute(
            "INSERT OR REPLACE INTO organisms(organism, display_name, taxid, gram_stain, phylum) "
            "VALUES (?,?,?,?,?)",
            (slug, blk.get("display_name"), blk.get("taxid"),
             blk.get("gram_stain"), blk.get("phylum")))


def populate_antibiotics_meta(conn):
    """Backfill WHO AWaRe + mechanism_type on whatever antibiotic rows exist so
    far, resolved from the registry (antibiotics.yaml) — no hard-coded tables."""
    for (ab,) in conn.execute("SELECT antibiotic FROM antibiotics").fetchall():
        conn.execute("UPDATE antibiotics SET mechanism_type=?, who_aware=? WHERE antibiotic=?",
                     (antibiotic_mechanism_type(ab), antibiotic_who_aware(ab), ab))


def _count_features(matrix_dir):
    f = Path(matrix_dir) / "features.txt"
    if f.exists():
        with open(f, encoding="utf-8") as fh:
            return sum(1 for _ in fh)
    return None


def populate_external_concordance(conn, model_id, organism, antibiotic):
    """M13 leakage-free head-to-head: our MODEL vs AMRFinderPlus vs ResFinder,
    all scored vs EUCAST/CLSI phenotype on the model's held-out TEST genomes
    (identical genome set for all three), from 16_concordance_summary_{org}.json.
    Replaces this model's rows so re-runs are idempotent."""
    conn.execute("DELETE FROM external_concordance WHERE model_id=?", (model_id,))
    summ = _read_json(_find(PROJECT_ROOT / "results" / organism,
                            f"16_concordance_summary_{organism}.json"))
    if not summ:
        return 0
    d = (summ.get("head_to_head_model_test_genomes") or {}).get(antibiotic)
    if not isinstance(d, dict):
        return 0
    n_test = d.get("n_common_test_genomes")
    cnt = 0
    for caller in ("model", "amrfinderplus", "resfinder"):
        s = d.get(caller)
        if not isinstance(s, dict):
            continue
        conn.execute(
            """INSERT OR REPLACE INTO external_concordance
               (model_id, caller, reference, n_test, sensitivity, specificity,
                balanced_accuracy, cohen_kappa, major_error_rate, very_major_error_rate)
               VALUES (?,?,?,?,?,?,?,?,?,?)""",
            (model_id, caller, "EUCAST/CLSI (held-out test)", n_test,
             _f(s.get("sensitivity")), _f(s.get("specificity")), _f(s.get("balanced_accuracy")),
             _f(s.get("cohen_kappa")), _f(s.get("major_error_rate")), _f(s.get("very_major_error_rate"))))
        cnt += 1
    return cnt


def _card_version_from_file(organism, antibiotic, config):
    """CARD DB provenance distilled from the file 08 writes (blastdbcmd -info) into
    one short string, e.g. 'CARD card.fna | 6,052 seqs | built Mar 21, 2026', so
    pipeline_runs.card_version records WHICH CARD snapshot annotated this model
    rather than a NULL. This is the auto-captured source; the env/config knobs stay
    as manual overrides. Best-effort — returns None if the file is absent."""
    try:
        d = resolve_path("dir_05_explainability", organism=organism,
                         antibiotic=antibiotic, config=config)
        txt = (d / "card_db_version.txt").read_text(encoding="utf-8")
    except Exception:
        return None
    name = date = seqs = None
    for raw in txt.splitlines():
        s = raw.strip()
        if s.startswith("Database:"):
            name = s.split(":", 1)[1].strip()
        elif s.startswith("Date:"):
            date = s.split(":", 1)[1].split("\t")[0].strip()
        elif "sequences;" in s:
            seqs = s.split("sequences;")[0].strip()
    if not name:
        return None
    parts = [f"CARD {name}"]
    if seqs:
        parts.append(f"{seqs} seqs")
    if date:
        parts.append(f"built {date}")
    return " | ".join(parts)


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="Populate the AMRK-DB knowledge base.")
    ap.add_argument("--organism", default=get_target(config=config)[0])
    ap.add_argument("--antibiotic", default=get_target(config=config)[1])
    ap.add_argument("--db", default=None, help="SQLite path (default: results/kb/amrk.db — "
                    "unified multi-organism KB; models.organism distinguishes rows)")
    args = ap.parse_args()
    organism, antibiotic = args.organism, args.antibiotic

    k_length = int(config["preprocessing"]["k_length"])
    # CARD version (M6) — AMR_CARD_VERSION env override wins over config.yaml so
    # HPC can record it without editing the (manually-tuned) config. Neither is
    # usually set, so fall back to the file 08 auto-writes (blastdbcmd -info): that
    # keeps card_version a real value instead of the NULL it was.
    card_version = (os.environ.get("AMR_CARD_VERSION")
                    or (config.get("blast", {}) or {}).get("card_version")
                    or _card_version_from_file(organism, antibiotic, config))
    drug_class = antibiotic_to_class(antibiotic)  # registry class_id (e.g. 'penicillins')

    # Resolve roots (organism/antibiotic-scoped) and glob the artefacts.
    models_dir = resolve_path("models_dir", organism=organism, antibiotic=antibiotic, config=config)
    results_root = PROJECT_ROOT / "results" / organism / antibiotic
    runs_root = PROJECT_ROOT / "runs" / organism / antibiotic

    run_meta = _read_json(_find(runs_root, "run_metadata.json"))
    manifest = _read_json(models_dir / "manifest.json")
    metrics = _read_json(_find(results_root, f"09_metrics_{antibiotic}.json"))
    holdout = _read_csv(_find(results_root, f"10_repeated_holdout_summary_{antibiotic}.csv"))
    # 10's output is the richest per-k-mer table; fall back to 09's candidates.
    cand = _read_csv(_find(results_root, f"10_kmer_background_frequency_{antibiotic}.csv"))
    if cand is None:
        cand = _read_csv(_find(results_root, f"07_kb_candidates_{antibiotic}.csv"))
    snp = _read_csv(_find(results_root, f"11_variant_snp_check_{antibiotic}.csv"))
    # Permutation significance (step 12 MDA + step 12b label-permutation null).
    perm_df = _read_csv(_find(results_root, f"12_permutation_test_{antibiotic}.csv"))
    labelperm = _read_json(_find(results_root, f"12b_label_permutation_summary_{antibiotic}.json"))
    # CPSS-stable, CARD-annotated unitigs (steps 13/13b).
    cpss = _read_csv(_find(results_root, f"13_stable_kb_candidates_{antibiotic}.csv"))
    # pyseer LMM lineage-corrected significance (step 14).
    pyseer_sig = _read_csv(_find(results_root, f"14_pyseer_significant_{antibiotic}.csv"))
    pyseer_sum = _read_json(_find(results_root, f"14_pyseer_summary_{antibiotic}.json"))

    # Adaptive min_support actually used (from pipeline_runs if present, else config)
    min_support = (run_meta or {}).get("min_support")

    # Unified multi-organism KB (schema tags each model with organism). Per-organism
    # KBs are deprecated; pass --db to override. results_root above stays
    # organism/antibiotic-scoped because the per-run artefacts live there.
    db_path = Path(args.db) if args.db else (PROJECT_ROOT / "results" / "kb" / "amrk.db")
    db_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print(f"POPULATE AMRK-DB  ({organism} / {antibiotic})  ->  {db_path}")
    print("=" * 70)
    print(f"  run_metadata : {'yes' if run_meta else 'MISSING'}")
    print(f"  manifest/06  : {'yes' if manifest else 'MISSING'} / {'yes' if metrics else 'MISSING'}")
    print(f"  07b holdout  : {'yes' if holdout is not None else 'MISSING'}")
    print(f"  candidates   : {'10 (with background)' if (cand is not None and 'discriminative' in cand.columns) else ('09' if cand is not None else 'MISSING')}")
    print(f"  11 SNP       : {'yes' if snp is not None else 'absent (ok)'}")
    print(f"  permutation  : MDA {'yes' if perm_df is not None else 'absent'} / "
          f"label-perm {'yes' if labelperm else 'absent'}")
    print(f"  CPSS stable  : {'yes' if cpss is not None else 'absent'}")
    print(f"  pyseer LMM   : {'yes' if pyseer_sig is not None else 'absent'}")

    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    create_schema(conn)
    # Dedup any legacy duplicates + add UNIQUE indexes so the blast/evidence
    # INSERT OR IGNOREs below are idempotent (audit Issue 3/24).
    ensure_unique_indexes(conn)
    matrix_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic, config=config)
    try:
        populate_organisms(conn)                                  # 0.5.0 reference table
        run_id = populate_run(conn, organism, antibiotic, run_meta, card_version,
                              min_support,
                              pyseer_version=(pyseer_sum or {}).get("pyseer_version"))
        model_id = populate_model(conn, run_id, antibiotic, drug_class, manifest, metrics, holdout)
        n_k = populate_candidates(conn, model_id, run_id, k_length, cand, card_version)
        n_s = populate_snp(conn, model_id, run_id, k_length, snp)
        n_p = populate_permutation(conn, run_id, k_length, perm_df, labelperm)
        n_c = populate_cpss(conn, model_id, run_id, k_length, cpss)
        n_l = populate_pyseer(conn, run_id, pyseer_sig,
                              (pyseer_sum or {}).get("bonferroni_threshold"))
        # 0.7.0: composite evidence tier (runs last — needs every other layer
        # already written for this model, esp. pyseer/cpss/blast).
        n_et = populate_evidence_tier(conn, model_id, run_id, perm_df)
        n_novel = conn.execute(
            "SELECT COUNT(*) FROM unitig_evidence_tier WHERE model_id=? AND is_novel_candidate=1",
            (model_id,)).fetchone()[0]
        # 0.5.0: model feature count + antibiotic meta + external concordance
        nf = _count_features(matrix_dir)
        conn.execute("UPDATE models SET n_features=? WHERE model_id=?", (nf, model_id))
        populate_antibiotics_meta(conn)                           # after the ab row exists
        n_ext = populate_external_concordance(conn, model_id, organism, antibiotic)
        update_metadata(conn, card_version)
        conn.commit()
    finally:
        conn.close()

    print(f"\n  ✓ run_id={run_id}  model_id={model_id}  | unitigs loaded: {n_k} "
          f"| SNP rows: {n_s} | permutation evidence: {n_p} | CPSS stable: {n_c} "
          f"| pyseer LMM evidence: {n_l} | n_features: {nf} | external concordance: {n_ext}")
    print(f"  ✓ evidence_tier graded: {n_et} unitigs | novel candidates (strong_novel): {n_novel}")
    print(f"  ✓ KB written: {db_path}  (schema {KB_SCHEMA_VERSION})")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
