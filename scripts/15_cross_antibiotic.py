#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cross-antibiotic stable-unitig overlap (Should-have S1; ROADMAP §1.6 / H3).

Runs *after* populate_database.py — it operates directly on the populated KB
(``results/{org}/kb/amrk.db``), which already holds every antibiotic's stable
unitig set in ``unitig_model_scores`` (``stable = 1``). For every pair of
antibiotics it computes the overlap of the stable sets, records each shared
unitig in ``unitig_antibiotic_overlap`` (the schema's S1/H3 substrate), and
writes a descriptive report (CSV + JSON).

This is the substrate the API ``/overlap?ab1=&ab2=`` (S8) reads, and the input
to H3 ("within-class β-lactam overlap > cross-class overlap").

**Hypergeometric significance is DEFERRED by design** (``--with-test``, off by
default). Rationale: with the current canonical pair (ampicillin = penicillin,
ciprofloxacin = fluoroquinolone) there is only *one* cross-class pair and *no*
within-class pair, so H3 cannot yet be concluded and a single p-value would be
descriptive at best. The test machinery is implemented and wired so that adding
a within-class β-lactam (e.g. cefotaxime) makes ``--with-test`` immediately
meaningful — no code change. See ROADMAP §1.6 for the universe definition.

Reification-safe wording (S10): unitigs shared across antibiotics are
"co-selected" / carry a "shared statistical signal" — not "cause" resistance.

Usage:
  python scripts/15_cross_antibiotic.py                      # organism from config
  python scripts/15_cross_antibiotic.py --organism ecoli
  python scripts/15_cross_antibiotic.py --with-test          # opt into the (preliminary) hypergeometric p
"""

import argparse
import csv
import datetime
import itertools
import json
import sqlite3
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config  # noqa: E402
from lib.logging_utils import get_logger  # noqa: E402
from lib.registry import antibiotic_to_class  # noqa: E402

# β-lactam super-family: penicillins, cephalosporins, carbapenems and the other
# β-lactam registry classes all share the β-lactamase / PBP mechanism, so an
# "is this a within-β-lactam pair?" question (H3) spans several *registry*
# classes. ``same_class`` in the DB stays the strict registry-class identity
# (schema semantics); this coarser grouping is surfaced only in the summary so
# H3 is ready the moment a second β-lactam (cefotaxime) enters the KB.
_BETA_LACTAM_CLASSES = {
    "penicillins",
    "cephalosporins",
    "beta_lactams_carbapenems_others",
    "monobactams",
}


def _drug_family(class_id):
    """Coarse family for the H3 within-/cross- contrast (β-lactams collapse)."""
    if class_id in _BETA_LACTAM_CLASSES:
        return "beta_lactam"
    return class_id


def fill_drug_classes(conn, logger):
    """Backfill ``antibiotics.drug_class`` from the registry (populate leaves it
    NULL). Returns {antibiotic: class_id}."""
    classes = {}
    for (ab,) in conn.execute("SELECT antibiotic FROM antibiotics").fetchall():
        cid = antibiotic_to_class(ab)
        classes[ab] = cid
        conn.execute("UPDATE antibiotics SET drug_class = ? WHERE antibiotic = ?",
                     (cid, ab))
        if cid is None:
            logger.warning("antibiotic '%s' not found in registry -> drug_class NULL", ab)
    return classes


def stable_sets(conn):
    """Return {antibiotic: set(unitig_id)} of stable unitigs (any selection
    method) per antibiotic, joining through the model. Only antibiotics that
    actually have a model + stable unitigs appear."""
    rows = conn.execute(
        """SELECT m.antibiotic, s.unitig_id
             FROM unitig_model_scores s
             JOIN models m ON m.model_id = s.model_id
            WHERE s.stable = 1"""
    ).fetchall()
    sets = {}
    for ab, uid in rows:
        sets.setdefault(ab, set()).add(uid)
    return sets


def gene_for_unitig(conn, unitig_id):
    """Best confirmed/candidate gene symbol(s) for a unitig, '' if none."""
    rows = conn.execute(
        """SELECT DISTINCT gene_symbol FROM blast_annotations
            WHERE unitig_id = ? AND tier IN ('confirmed','candidate')
              AND gene_symbol IS NOT NULL""",
        (unitig_id,),
    ).fetchall()
    return ";".join(sorted({r[0] for r in rows if r[0]}))


def hypergeom_sf(k, N, K, n):
    """P(X >= k) for the overlap of two stable sets of sizes K and n drawn from
    a universe of N unitigs (enrichment). Uses scipy if available, else an exact
    math.comb fallback. Returns None if the universe is degenerate."""
    if N <= 0 or K <= 0 or n <= 0 or k <= 0:
        return None
    try:
        from scipy.stats import hypergeom
        return float(hypergeom.sf(k - 1, N, K, n))
    except Exception:
        from math import comb
        denom = comb(N, n)
        if denom == 0:
            return None
        kmax = min(K, n)
        p = sum(comb(K, i) * comb(N - K, n - i) for i in range(k, kmax + 1)) / denom
        return float(p)


def populate_overlap(conn, sets, classes, logger):
    """Recompute ``unitig_antibiotic_overlap`` for every antibiotic pair.

    Returns a list of per-pair summary dicts (descriptive)."""
    conn.execute("DELETE FROM unitig_antibiotic_overlap")  # full recompute (idempotent)
    abs_sorted = sorted(sets)
    pairs = list(itertools.combinations(abs_sorted, 2))
    union_all = set().union(*sets.values()) if sets else set()
    summaries = []
    for a, b in pairs:
        sa, sb = sets[a], sets[b]
        shared = sorted(sa & sb)
        union = sa | sb
        same_class = int(classes.get(a) is not None and classes.get(a) == classes.get(b))
        same_family = int(_drug_family(classes.get(a)) is not None
                          and _drug_family(classes.get(a)) == _drug_family(classes.get(b)))
        for uid in shared:
            conn.execute(
                """INSERT OR REPLACE INTO unitig_antibiotic_overlap
                       (unitig_id, antibiotic_a, antibiotic_b, same_class)
                   VALUES (?,?,?,?)""",
                (uid, a, b, same_class),
            )
        jaccard = len(shared) / len(union) if union else 0.0
        summaries.append({
            "antibiotic_a": a,
            "antibiotic_b": b,
            "class_a": classes.get(a),
            "class_b": classes.get(b),
            "same_class": bool(same_class),
            "same_drug_family": bool(same_family),
            "n_stable_a": len(sa),
            "n_stable_b": len(sb),
            "n_overlap": len(shared),
            "jaccard": round(jaccard, 6),
            "shared_unitig_ids": shared,
        })
        logger.info("%s vs %s | stable %d/%d | overlap %d | jaccard %.4f | same_class=%d",
                    a, b, len(sa), len(sb), len(shared), jaccard, same_class)
    return summaries, union_all


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="Cross-antibiotic stable-unitig overlap (S1/H3).")
    ap.add_argument("--organism", default=config.get("project", {}).get("organism", "ecoli"))
    ap.add_argument("--db", default=None, help="SQLite path (default: results/{org}/kb/amrk.db)")
    ap.add_argument("--with-test", action="store_true",
                    help="Also compute the (preliminary) hypergeometric p per pair. "
                         "DEFERRED by default — see module docstring / ROADMAP §1.6.")
    args = ap.parse_args()
    organism = args.organism

    db_path = Path(args.db) if args.db else (PROJECT_ROOT / "results" / organism / "kb" / "amrk.db")
    out_dir = db_path.parent
    logger = get_logger("s1-cross-antibiotic")

    if not db_path.exists():
        logger.error("KB not found: %s (run populate_database.py first)", db_path)
        sys.exit(1)

    logger.info("=" * 64)
    logger.info("S1 cross-antibiotic overlap  (%s)  ->  %s", organism, db_path)
    logger.info("=" * 64)

    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    try:
        classes = fill_drug_classes(conn, logger)
        sets = stable_sets(conn)
        if len(sets) < 2:
            logger.warning("only %d antibiotic(s) with stable unitigs in the KB — "
                           "need >=2 for any overlap. Nothing to do.", len(sets))
            conn.commit()
            return
        summaries, union_all = populate_overlap(conn, sets, classes, logger)

        # Attach shared-unitig gene annotations for the report.
        for s in summaries:
            s["shared_genes"] = {uid: gene_for_unitig(conn, uid) for uid in s["shared_unitig_ids"]}

        # Hypergeometric significance — DEFERRED unless explicitly requested.
        N = len(union_all)  # ROADMAP §1.6 universe: union of all stable unitigs
        has_within = any(s["same_drug_family"] for s in summaries)
        if args.with_test:
            for s in summaries:
                s["hypergeom_universe_N"] = N
                s["hypergeom_p_enrichment"] = hypergeom_sf(
                    s["n_overlap"], N, s["n_stable_a"], s["n_stable_b"])
            logger.warning("PRELIMINARY hypergeometric test computed over union-universe "
                           "N=%d (ROADMAP §1.6). Caveat: per-antibiotic unitig "
                           "vocabularies differ; H3 needs a within-class pair.", N)
        else:
            logger.info("Hypergeometric test DEFERRED (use --with-test to compute a "
                        "preliminary p). within-class/β-lactam pair present: %s.",
                        has_within)

        conn.commit()
    finally:
        conn.close()

    # ---- write descriptive outputs --------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / f"15_cross_antibiotic_overlap_{organism}.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["antibiotic_a", "antibiotic_b", "same_class", "same_drug_family",
                    "unitig_id", "shared_gene"])
        for s in summaries:
            for uid in s["shared_unitig_ids"]:
                w.writerow([s["antibiotic_a"], s["antibiotic_b"], int(s["same_class"]),
                            int(s["same_drug_family"]), uid, s["shared_genes"].get(uid, "")])

    summary_path = out_dir / f"15_cross_antibiotic_summary_{organism}.json"
    summary_doc = {
        "organism": organism,
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "n_antibiotics": len(sets),
        "antibiotics": {ab: {"class": classes.get(ab), "n_stable": len(s)}
                        for ab, s in sets.items()},
        "hypergeometric_test": "computed" if args.with_test else "deferred",
        "pairs": [{k: v for k, v in s.items() if k != "shared_genes"} for s in summaries],
    }
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary_doc, fh, indent=2)

    logger.info("✓ overlap table populated; %d pair(s) summarised", len(summaries))
    logger.info("✓ %s", csv_path)
    logger.info("✓ %s", summary_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
