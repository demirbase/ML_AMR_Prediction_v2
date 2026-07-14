#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Read-only query layer over the AMRK-DB knowledge base (S8/S9 API backend).

Pure ``sqlite3`` functions returning plain dicts/lists — no web framework, so
they are unit-testable without FastAPI. ``kb_api.py`` is a thin FastAPI wrapper
that just exposes these over HTTP. Schema: ``lib/kb_schema.py`` (unitigs,
unitig_model_scores, blast_annotations, unitig_background_frequency,
variant_snp_check, unitig_antibiotic_overlap, validation_evidence, models,
pipeline_runs, kb_metadata).
"""

import sqlite3


def connect(db_path):
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    return conn


def _rows(conn, sql, params=()):
    return [dict(r) for r in conn.execute(sql, params).fetchall()]


def get_metadata(conn):
    """FAIR machine-readable KB metadata (S9): schema version, DOI, license, sizes."""
    r = conn.execute("SELECT * FROM kb_metadata WHERE id = 1").fetchone()
    meta = dict(r) if r else {}
    meta["antibiotics"] = [x["antibiotic"] for x in
                           conn.execute("SELECT antibiotic FROM models ORDER BY antibiotic")]
    return meta


def get_stats(conn):
    """Aggregate counts for the /stats endpoint."""
    one = lambda sql: conn.execute(sql).fetchone()[0]
    per_ab = _rows(conn,
        """SELECT m.antibiotic,
                  COUNT(DISTINCT s.unitig_id) AS n_scored,
                  COUNT(DISTINCT CASE WHEN s.stable=1 THEN s.unitig_id END) AS n_stable
             FROM models m LEFT JOIN unitig_model_scores s ON s.model_id = m.model_id
            GROUP BY m.antibiotic ORDER BY m.antibiotic""")
    tiers = _rows(conn,
        "SELECT tier, COUNT(*) AS n FROM blast_annotations WHERE tier IS NOT NULL "
        "GROUP BY tier ORDER BY n DESC")
    return {
        "n_unitigs": one("SELECT COUNT(*) FROM unitigs"),
        "n_models": one("SELECT COUNT(*) FROM models"),
        "n_evidence": one("SELECT COUNT(*) FROM validation_evidence"),
        "per_antibiotic": per_ab,
        "blast_tiers": tiers,
    }


def list_biomarkers(conn, antibiotic=None, min_stability=None, tier=None,
                    stable_only=False, limit=200, offset=0):
    """Filterable biomarker list (unitig × model), with best BLAST gene/tier.

    Joins the per-(unitig,model) scores to the model's antibiotic and the unitig's
    best confirmed/candidate BLAST hit. Filters are all optional."""
    where, params = ["1=1"], []
    if antibiotic:
        where.append("m.antibiotic = ?"); params.append(antibiotic)
    if min_stability is not None:
        where.append("s.selection_frequency >= ?"); params.append(float(min_stability))
    if stable_only:
        where.append("s.stable = 1")
    if tier:
        where.append("ba.tier = ?"); params.append(tier)
    sql = f"""
        SELECT u.unitig_id, u.sequence, m.antibiotic, s.selection_method,
               s.gain, s.selection_frequency, s.stable, s.composite_score,
               s.mean_abs_shap, ba.gene_symbol, ba.tier, ba.identity_pct,
               ba.aro_accession, ba.aro_gene_family, ba.aro_drug_class
          FROM unitig_model_scores s
          JOIN models m   ON m.model_id = s.model_id
          JOIN unitigs u  ON u.unitig_id = s.unitig_id
          LEFT JOIN blast_annotations ba
                 ON ba.unitig_id = s.unitig_id
                AND ba.tier IN ('confirmed','candidate')
                AND ba.annotation_id = (
                    SELECT annotation_id FROM blast_annotations b2
                     WHERE b2.unitig_id = s.unitig_id
                       AND b2.tier IN ('confirmed','candidate')
                     ORDER BY b2.identity_pct DESC LIMIT 1)
         WHERE {' AND '.join(where)}
         GROUP BY u.unitig_id, m.model_id, s.selection_method
         ORDER BY s.selection_frequency DESC, s.gain DESC
         LIMIT ? OFFSET ?"""
    return _rows(conn, sql, (*params, int(limit), int(offset)))


def get_unitig(conn, sequence):
    """Full evidence chain for one unitig (by exact sequence), or None."""
    u = conn.execute("SELECT * FROM unitigs WHERE sequence = ?", (sequence,)).fetchone()
    if not u:
        return None
    uid = u["unitig_id"]
    return {
        "unitig": dict(u),
        "model_scores": _rows(conn,
            """SELECT m.antibiotic, s.selection_method, s.gain, s.selection_frequency,
                      s.stable, s.composite_score, s.mean_abs_shap
                 FROM unitig_model_scores s JOIN models m ON m.model_id = s.model_id
                WHERE s.unitig_id = ?""", (uid,)),
        "blast": _rows(conn,
            "SELECT source_db, gene_symbol, identity_pct, coverage, evalue, tier, "
            "aro_accession, aro_gene_family, aro_drug_class, aro_resistance_mechanism "
            "FROM blast_annotations WHERE unitig_id = ?", (uid,)),
        "background_frequency": _rows(conn,
            "SELECT m.antibiotic, bf.prevalence_resistant, bf.prevalence_susceptible, "
            "bf.delta_prevalence, bf.odds_ratio, bf.fisher_p, bf.discriminative "
            "FROM unitig_background_frequency bf JOIN models m ON m.model_id=bf.model_id "
            "WHERE bf.unitig_id = ?", (uid,)),
        "snp": _rows(conn,
            "SELECT card_model, snp, allele_class FROM variant_snp_check "
            "WHERE unitig_id = ?", (uid,)),
        "overlap": _rows(conn,
            "SELECT antibiotic_a, antibiotic_b, same_class FROM unitig_antibiotic_overlap "
            "WHERE unitig_id = ?", (uid,)),
        "evidence": _rows(conn,
            "SELECT evidence_type, evidence_source, evidence_score, pipeline_run_id "
            "FROM validation_evidence WHERE unitig_id = ?", (uid,)),
    }


def get_overlap(conn, ab1, ab2, organism=None):
    """Cross-antibiotic shared stable unitigs for a pair (order-independent).

    The overlap table is organism-aware (schema 0.6.0): pass ``organism`` to keep
    a same-drug pair from being merged across species; None returns all organisms.
    """
    org_clause = " AND o.organism = ?" if organism else ""
    params = (ab1, ab2, ab2, ab1) + ((organism,) if organism else ())
    rows = _rows(conn,
        f"""SELECT o.unitig_id, o.organism, u.sequence, o.same_class,
                  (SELECT group_concat(DISTINCT gene_symbol) FROM blast_annotations b
                    WHERE b.unitig_id = o.unitig_id AND b.tier IN ('confirmed','candidate')
                      AND gene_symbol IS NOT NULL) AS genes
             FROM unitig_antibiotic_overlap o JOIN unitigs u ON u.unitig_id = o.unitig_id
            WHERE ((antibiotic_a = ? AND antibiotic_b = ?)
               OR (antibiotic_a = ? AND antibiotic_b = ?)){org_clause}""",
        params)
    return {"antibiotic_a": ab1, "antibiotic_b": ab2, "organism": organism,
            "n_shared": len(rows), "shared_unitigs": rows}
