#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""AMRK-DB REST API (Should-have S8/S9; ROADMAP §1.8).

A minimal, read-only FastAPI over the SQLite knowledge base — enough for the
publication "the KB is queryable via an open API" criterion (Database Oxford /
Briefings). The data-access logic lives in ``lib/kb_queries`` (pure sqlite3,
unit-tested without a web server); this module is only the thin HTTP layer:
versioned ``/api/v1`` routes, CORS, and auto OpenAPI docs at ``/docs``.

Endpoints (ROADMAP §1.8):
    GET /api/v1/metadata                     FAIR metadata (schema ver, DOI, license) — S9
    GET /api/v1/stats                        aggregate counts
    GET /api/v1/kmers?antibiotic=&tier=&evidence_tier=&novel_only=&min_stability=&stable_only=&limit=&offset=
    GET /api/v1/kmers/{sequence}             one unitig's full evidence chain
    GET /api/v1/novel?antibiotic=&organism=  strong_novel biomarkers (0.7.0) — no known gene, CPSS+pyseer
    GET /api/v1/overlap?ab1=&ab2=            cross-antibiotic shared stable unitigs

Run:
    pip install fastapi uvicorn
    AMR_KB_DB=results/kb/amrk.db uvicorn scripts.kb_api:app --reload
    # or: python scripts/kb_api.py   (serves on :8000)
"""

import os
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import kb_queries as Q  # noqa: E402
from lib.kb_schema import KB_SCHEMA_VERSION  # noqa: E402

DB_PATH = Path(os.environ.get("AMR_KB_DB", PROJECT_ROOT / "results" / "kb" / "amrk.db"))


def _create_app():
    from fastapi import FastAPI, HTTPException, Query
    from fastapi.middleware.cors import CORSMiddleware

    app = FastAPI(
        title="AMRK-DB API",
        version=KB_SCHEMA_VERSION,   # single source of truth — never hardcode a copy
        description="Read-only API over the stability-filtered, lineage-validated, "
                    "unitig-resolution AMR biomarker knowledge base for ESKAPEE pathogens.",
    )
    app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["GET"],
                       allow_headers=["*"])

    def _conn():
        if not DB_PATH.exists():
            raise HTTPException(503, f"KB not available at {DB_PATH}")
        return Q.connect(DB_PATH)

    @app.get("/api/v1/metadata")
    def metadata():
        c = _conn()
        try:
            return Q.get_metadata(c)
        finally:
            c.close()

    @app.get("/api/v1/stats")
    def stats():
        c = _conn()
        try:
            return Q.get_stats(c)
        finally:
            c.close()

    @app.get("/api/v1/kmers")
    def kmers(antibiotic: str | None = None, tier: str | None = None,
              min_stability: float | None = None, stable_only: bool = False,
              evidence_tier: str | None = None, novel_only: bool = False,
              limit: int = Query(200, le=2000), offset: int = 0):
        c = _conn()
        try:
            rows = Q.list_biomarkers(c, antibiotic=antibiotic, tier=tier,
                                     min_stability=min_stability, stable_only=stable_only,
                                     evidence_tier=evidence_tier, novel_only=novel_only,
                                     limit=limit, offset=offset)
            return {"count": len(rows), "results": rows}
        finally:
            c.close()

    @app.get("/api/v1/novel")
    def novel(antibiotic: str | None = None, organism: str | None = None,
              limit: int = Query(200, le=2000), offset: int = 0):
        """strong_novel biomarkers (0.7.0): CPSS-stable + pyseer-significant with
        no known CARD gene — the candidates the BLAST-only tier hides as `none`."""
        c = _conn()
        try:
            rows = Q.list_novel_candidates(c, antibiotic=antibiotic, organism=organism,
                                           limit=limit, offset=offset)
            return {"count": len(rows), "results": rows}
        finally:
            c.close()

    @app.get("/api/v1/kmers/{sequence}")
    def kmer(sequence: str):
        c = _conn()
        try:
            rec = Q.get_unitig(c, sequence)
            if rec is None:
                raise HTTPException(404, f"unitig not found: {sequence}")
            return rec
        finally:
            c.close()

    @app.get("/api/v1/overlap")
    def overlap(ab1: str, ab2: str, organism: str = None):
        c = _conn()
        try:
            return Q.get_overlap(c, ab1, ab2, organism=organism)
        finally:
            c.close()

    @app.get("/")
    def root():
        return {"name": "AMRK-DB API", "version": KB_SCHEMA_VERSION, "docs": "/docs",
                "endpoints": ["/api/v1/metadata", "/api/v1/stats", "/api/v1/kmers",
                              "/api/v1/kmers/{sequence}", "/api/v1/novel",
                              "/api/v1/overlap?ab1=&ab2="]}

    return app


# Module-level `app` for uvicorn (`uvicorn scripts.kb_api:app`). If FastAPI is not
# installed, importing lib.kb_queries still works (tested standalone); `app` is None.
try:
    app = _create_app()
except ImportError:
    app = None


if __name__ == "__main__":
    if app is None:
        sys.exit("FastAPI not installed. `pip install fastapi uvicorn` first.")
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=int(os.environ.get("AMR_API_PORT", 8000)))
