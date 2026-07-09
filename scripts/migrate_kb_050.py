#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Backfill the AMRK-DB from schema 0.4.0 -> 0.5.0.

Additive, in-place, and SAFE on an already-populated KB: it never re-inserts a
model/evidence row (so it cannot hit populate_database's model INSERT FK issue).
It only:
  * creates the new tables/columns (via create_schema)         — organisms,
    external_concordance, antibiotics.{mechanism_type,who_aware}, models.n_features
  * fills the `organisms` reference table (gram stain / phylum)
  * backfills antibiotic meta (WHO AWaRe + acquired/target-SNP)
  * sets models.n_features (unitig count from each model's matrix features.txt)
  * loads external_concordance from 16_concordance_{organism}.csv (M13), if present
  * stamps card_version where it was left NULL

Run inside the container with the unitig representation resolved, e.g.:
    AMR_FEATURE_REPR=unitig AMR_CARD_VERSION=4.0.1 \
      apptainer exec --no-home $SIF python scripts/migrate_kb_050.py --db results/kb/amrk.db
"""
import argparse
import sqlite3
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.kb_schema import create_schema, ensure_unique_indexes, KB_SCHEMA_VERSION  # noqa: E402
from lib.config import load_config, resolve_path  # noqa: E402
from populate_database import (  # noqa: E402
    populate_organisms, populate_antibiotics_meta, _count_features,
    populate_external_concordance, update_metadata,
)


def main():
    ap = argparse.ArgumentParser(description="Backfill AMRK-DB to schema 0.5.0 (additive, safe).")
    ap.add_argument("--db", default="results/kb/amrk.db")
    ap.add_argument("--card-version", default="4.0.1")
    args = ap.parse_args()
    config = load_config()

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys = ON")
    create_schema(conn)           # new tables + columns (idempotent)
    ensure_unique_indexes(conn)
    populate_organisms(conn)
    populate_antibiotics_meta(conn)

    rows = conn.execute(
        "SELECT m.model_id, p.organism, m.antibiotic FROM models m "
        "JOIN pipeline_runs p USING(run_id) ORDER BY m.model_id").fetchall()
    n_ext = 0
    for mid, org, ab in rows:
        md = resolve_path("matrix_dir", organism=org, antibiotic=ab, config=config)
        conn.execute("UPDATE models SET n_features=? WHERE model_id=?", (_count_features(md), mid))
        n_ext += populate_external_concordance(conn, mid, org, ab)

    conn.execute("UPDATE pipeline_runs SET card_version=? WHERE card_version IS NULL OR card_version=''",
                 (args.card_version,))
    update_metadata(conn, args.card_version)
    conn.commit()

    print(f"KB backfilled to schema {KB_SCHEMA_VERSION}")
    q = lambda s: conn.execute(s).fetchone()[0]  # noqa: E731
    print(f"  organisms:                 {q('SELECT COUNT(*) FROM organisms')}")
    print(f"  antibiotics w/ AWaRe:      {q('SELECT COUNT(*) FROM antibiotics WHERE who_aware IS NOT NULL')}")
    print(f"  antibiotics w/ mechanism:  {q('SELECT COUNT(*) FROM antibiotics WHERE mechanism_type IS NOT NULL')}")
    print(f"  models w/ n_features:      {q('SELECT COUNT(*) FROM models WHERE n_features IS NOT NULL')} / {len(rows)}")
    print(f"  external_concordance rows: {q('SELECT COUNT(*) FROM external_concordance')}")
    print(f"  card_version (metadata):   {q('SELECT card_version FROM kb_metadata')}")
    conn.close()


if __name__ == "__main__":
    main()
