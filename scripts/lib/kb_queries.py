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


def _table_exists(conn, name):
    """Guard for pre-0.7.0 KBs that lack unitig_evidence_tier — keep the API
    working (evidence_tier fields just come back NULL) until a re-populate adds
    the table, instead of failing with 'no such table'."""
    return conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (name,)
    ).fetchone() is not None


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
                    stable_only=False, evidence_tier=None, novel_only=False,
                    limit=200, offset=0):
    """Filterable biomarker list (unitig × model), with best BLAST gene/tier AND
    the composite evidence_tier (0.7.0).

    Joins the per-(unitig,model) scores to the model's antibiotic, the unitig's
    best confirmed/candidate BLAST hit, and its composite evidence tier. Filters
    are all optional: ``tier`` is the BLAST-only layer-1 grade; ``evidence_tier``
    is the composite grade; ``novel_only`` keeps only strong_novel candidates."""
    has_et = _table_exists(conn, "unitig_evidence_tier")
    if (evidence_tier or novel_only) and not has_et:
        return []   # pre-0.7.0 KB: no composite tier to filter on
    where, params = ["1=1"], []
    if antibiotic:
        where.append("m.antibiotic = ?"); params.append(antibiotic)
    if min_stability is not None:
        where.append("s.selection_frequency >= ?"); params.append(float(min_stability))
    if stable_only:
        where.append("s.stable = 1")
    if tier:
        where.append("ba.tier = ?"); params.append(tier)
    if evidence_tier:
        where.append("et.evidence_tier = ?"); params.append(evidence_tier)
    if novel_only:
        where.append("et.is_novel_candidate = 1")
    et_select = ("et.evidence_tier, et.n_evidence_layers, et.evidence_layers, "
                 "et.is_novel_candidate" if has_et else
                 "NULL AS evidence_tier, NULL AS n_evidence_layers, "
                 "NULL AS evidence_layers, NULL AS is_novel_candidate")
    et_join = ("LEFT JOIN unitig_evidence_tier et "
               "ON et.unitig_id = s.unitig_id AND et.model_id = s.model_id"
               if has_et else "")
    sql = f"""
        SELECT u.unitig_id, u.sequence, m.antibiotic, s.selection_method,
               s.gain, s.selection_frequency, s.stable, s.composite_score,
               s.mean_abs_shap, ba.gene_symbol, ba.tier, ba.identity_pct,
               ba.aro_accession, ba.aro_gene_family, ba.aro_drug_class,
               {et_select}
          FROM unitig_model_scores s
          JOIN models m   ON m.model_id = s.model_id
          JOIN unitigs u  ON u.unitig_id = s.unitig_id
          {et_join}
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


def list_novel_candidates(conn, antibiotic=None, organism=None, limit=200, offset=0):
    """The flagship 0.7.0 query: `strong_novel` biomarkers — CPSS-stable +
    pyseer-significant unitigs with NO known CARD gene, which the BLAST-only tier
    hides as `none`. Returns unitig × model rows ordered by evidence breadth."""
    if not _table_exists(conn, "unitig_evidence_tier"):
        return []   # pre-0.7.0 KB
    where, params = ["et.is_novel_candidate = 1"], []
    if antibiotic:
        where.append("m.antibiotic = ?"); params.append(antibiotic)
    if organism:
        where.append("r.organism = ?"); params.append(organism)
    sql = f"""
        SELECT u.unitig_id, u.sequence, r.organism, m.antibiotic,
               et.evidence_tier, et.n_evidence_layers, et.evidence_layers
          FROM unitig_evidence_tier et
          JOIN models m   ON m.model_id = et.model_id
          JOIN pipeline_runs r ON r.run_id = m.run_id
          JOIN unitigs u  ON u.unitig_id = et.unitig_id
         WHERE {' AND '.join(where)}
         GROUP BY et.unitig_id, et.model_id
         ORDER BY et.n_evidence_layers DESC, u.unitig_id
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
        "evidence_tier": _rows(conn,
            "SELECT m.antibiotic, et.evidence_tier, et.n_evidence_layers, "
            "et.evidence_layers, et.is_novel_candidate "
            "FROM unitig_evidence_tier et JOIN models m ON m.model_id = et.model_id "
            "WHERE et.unitig_id = ?", (uid,)) if _table_exists(conn, "unitig_evidence_tier") else [],
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
