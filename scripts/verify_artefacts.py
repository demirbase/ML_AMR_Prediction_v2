#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Verify a delivered artefact tree: KB + tidy tables + thesis figures.

Answers one question — *is this tree the thing we think it is?* — for whichever
copy you point it at (laptop, HPC scratch, an unpacked Zenodo tarball). Run it
before every release, and after every sync, because the failures this project
actually hit were all silent:

  * a figure regenerated into an empty canvas over a good one (data root pruned)
  * a table that differed between machines because a sort tie fell to filesystem
    order
  * `on_target` blank for 14 of 45 models, so a figure dropped them without a word
  * a corrected figure that never re-rendered because argparse rejected a flag
    while the calling loop read the exit as success

Nothing here re-derives the science; it checks that the artefacts are present,
paired, non-empty, mutually consistent, and identical to a reference when one is
given.

Usage:
    python scripts/verify_artefacts.py --root results
    python scripts/verify_artefacts.py --root $AMR_WORK/results --expect-md5 md5.txt

`--expect-md5` takes `MD5␠␠FILENAME` lines (the output of `md5sum results/tables/*.csv`,
basenames only are compared) and turns a cross-machine comparison into a check
you can run on one machine instead of eyeballing two lists.

Exit code is 0 only when every check passes.
"""
import argparse
import csv
import hashlib
import json
import sqlite3
import sys
from pathlib import Path

EXPECTED_TABLES = {
    "models_summary.csv": 45,
    "kb_overview.csv": 45,
    "biomarkers.csv": 3571,
    "mechanisms.csv": None,          # row count varies with the CARD snapshot
    "cv_comparison.csv": 45,
    "h3_gene_family_overlap.csv": 138,
    "novel_ncbi_context.csv": 23,
    "novel_ncbi_hits.csv": None,
    "fair_mapping.csv": 15,          # the 15 FAIR principles, one row each
    # kb_tables_thesis.py. The counts below are structural, not sample-dependent:
    # 6 organisms, 7 evidence layers, 9+2+4 provenance items, 45 models, 9 limitations.
    # headline_biomarkers is left unpinned for the same reason as mechanisms.csv --
    # it moves with the CARD snapshot.
    "lineage_summary.csv": 6,
    "evidence_accounting.csv": 7,
    "provenance_tools.csv": 15,
    "hyperparameters.csv": 45,
    "limitations.csv": 9,
    "headline_biomarkers.csv": None,
}

# KB table -> expected row count. None = "must exist, count not pinned";
# 0 = "expected to be EMPTY" (see METHODOLOGY §5.2 — these are documented gaps,
# so a sudden non-zero count is as much a surprise as a sudden zero).
EXPECTED_KB = {
    "models": 45, "pipeline_runs": 45, "organisms": 7, "antibiotics": 22,
    "unitigs": 3844, "unitig_model_scores": 3613, "unitig_evidence_tier": 3571,
    "blast_annotations": 3611, "unitig_background_frequency": 2409,
    "variant_snp_check": 953, "validation_evidence": 10530, "kb_metadata": 1,
    "external_concordance": 0, "unitig_antibiotic_overlap": 0,
}

EXPECTED_TIERS = {"confirmed": 349, "strong_novel": 23, "candidate": 942,
                  "weak": 1920, "none": 337}

N_FIGURES = 41                  # 37 data + 2 schematics (37 pipeline, 38 KB schema)
                                # + 39 evidence combinations, 40 structure vs inflation
MIN_PNG_BYTES = 20_000          # a blank matplotlib canvas lands far below this


class Report:
    def __init__(self):
        self.fails, self.warns, self.n = [], [], 0

    def check(self, ok, label, detail=""):
        self.n += 1
        mark = "PASS" if ok else "FAIL"
        print(f"  [{mark}] {label}" + (f" — {detail}" if detail else ""))
        if not ok:
            self.fails.append(label)

    def warn(self, label, detail=""):
        self.warns.append(label)
        print(f"  [WARN] {label}" + (f" — {detail}" if detail else ""))


def verify_kb(root, rep):
    db = root / "kb" / "amrk.db"
    print("\nKB")
    if not db.exists():
        rep.check(False, "amrk.db present", str(db))
        return
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        have = {r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        for tbl, want in EXPECTED_KB.items():
            if tbl not in have:
                rep.check(False, f"table {tbl}", "missing")
                continue
            got = con.execute(f'SELECT COUNT(*) FROM "{tbl}"').fetchone()[0]
            rep.check(got == want, f"{tbl} rows", f"{got} (expected {want})")

        meta = con.execute(
            "SELECT kb_schema_version, zenodo_doi FROM kb_metadata").fetchone()
        rep.check(meta[0] == "0.7.1", "schema version", meta[0])
        rep.check(bool(meta[1]), "Zenodo DOI stamped", meta[1] or "(EMPTY)")

        tiers = dict(con.execute(
            "SELECT evidence_tier, COUNT(*) FROM unitig_evidence_tier "
            "GROUP BY evidence_tier"))
        rep.check(tiers == EXPECTED_TIERS, "evidence tiers",
                  " · ".join(f"{k} {v}" for k, v in sorted(tiers.items())))

        cv = dict(con.execute(
            "SELECT cv_method, COUNT(*) FROM models GROUP BY cv_method"))
        rep.check(cv == {"lineage_group_kfold_5fold": 45},
                  "all models lineage-CV (no fallback)", str(cv))

        dirty, commits, seeds = con.execute(
            "SELECT SUM(git_dirty), COUNT(DISTINCT git_commit), "
            "COUNT(DISTINCT random_seed) FROM pipeline_runs").fetchone()
        rep.check(seeds == 1, "single random seed", f"{seeds} distinct")
        if dirty:
            rep.warn("git_dirty on every run",
                     f"{dirty}/45 — known, documented in METHODOLOGY §5.3")
        # The two empty tables are expected; say so out loud so nobody "fixes" it.
        print("       (external_concordance and unitig_antibiotic_overlap are "
              "expected to be empty — METHODOLOGY §5.2)")
    finally:
        con.close()


def verify_tables(root, rep, expect_md5=None):
    print("\nTABLES")
    tdir = root / "tables"
    for name, want_rows in EXPECTED_TABLES.items():
        f = tdir / name
        if not f.exists():
            rep.check(False, name, "missing")
            continue
        with open(f, newline="", encoding="utf-8") as fh:
            rows = list(csv.reader(fh))
        n = len(rows) - 1
        ok = (want_rows is None) or (n == want_rows)
        rep.check(ok, f"{name} rows",
                  f"{n}" + (f" (expected {want_rows})" if want_rows else ""))

    js = tdir / "h3_summary.json"
    if js.exists():
        h = json.loads(js.read_text(encoding="utf-8"))
        rep.check(h.get("n_pairs") == 138 and h.get("H3_supported") is True,
                  "h3_summary", f"{h.get('n_pairs')} pairs · "
                                f"supported={h.get('H3_supported')}")
    else:
        rep.check(False, "h3_summary.json", "missing")

    # on_target must never be blank: a blank drops the model from figure 04.
    mech = tdir / "mechanisms.csv"
    if mech.exists():
        rows = list(csv.DictReader(open(mech, encoding="utf-8")))
        blank = sum(1 for r in rows if r.get("on_target") not in ("True", "False"))
        models = len({(r["organism"], r["antibiotic"]) for r in rows
                      if r.get("on_target") == "True"})
        rep.check(blank == 0, "mechanisms on_target complete",
                  f"{blank} blank · {models} models on-target")

    # A column copied from one tidy table into another must survive the trip as the
    # SAME TEXT. pandas' default CSV float parser is accurate to an ULP, not exact, and
    # which value it rounds differs by version: headline_biomarkers.csv differed between
    # laptop and HPC because the laptop mangled one delta_prevalence and the container
    # mangled one composite_score. Every other table matched, so nothing else caught it.
    hb, bm = tdir / "headline_biomarkers.csv", tdir / "biomarkers.csv"
    if hb.exists() and bm.exists():
        src = {}
        for r in csv.DictReader(open(bm, encoding="utf-8")):
            src[(r["unitig_id"], r["model_id"])] = r
        cols = ["delta_prevalence", "odds_ratio", "pyseer_lrt_p",
                "selection_frequency", "composite_score", "identity_pct", "coverage"]
        drift, checked = [], 0
        for r in csv.DictReader(open(hb, encoding="utf-8")):
            ref = src.get((r["unitig_id"], r["model_id"]))
            if not ref:
                continue
            for c in cols:
                if c in r and c in ref:
                    checked += 1
                    if r[c] != ref[c]:
                        drift.append(f"{c} {ref[c]!r}->{r[c]!r}")
        rep.check(not drift, "headline floats byte-identical to biomarkers.csv",
                  f"{checked} values checked" if not drift
                  else f"{len(drift)} differ, e.g. {drift[0]} "
                       f"(read_csv needs float_precision='round_trip')")

    if expect_md5:
        ref = {}
        for line in Path(expect_md5).read_text(encoding="utf-8").splitlines():
            parts = line.split()
            if len(parts) >= 2:
                ref[Path(parts[-1]).name] = parts[0]
        for name, want in sorted(ref.items()):
            f = tdir / name
            got = hashlib.md5(f.read_bytes()).hexdigest() if f.exists() else None
            rep.check(got == want, f"md5 {name}",
                      "identical" if got == want else f"{got} != {want}")


def verify_figures(root, rep):
    print("\nFIGURES")
    fdir = root / "figures"
    pngs = sorted(p for p in fdir.glob("*.png"))
    rep.check(len(pngs) == N_FIGURES, "PNG count", f"{len(pngs)} (expected {N_FIGURES})")
    rep.check(len(list(fdir.glob("*.pdf"))) == N_FIGURES, "PDF count",
              str(len(list(fdir.glob("*.pdf")))))

    missing_pdf = [p.name for p in pngs if not p.with_suffix(".pdf").exists()]
    rep.check(not missing_pdf, "every PNG has a PDF", ", ".join(missing_pdf) or "ok")

    tiny = [f"{p.name} ({p.stat().st_size // 1024} KB)"
            for p in pngs if p.stat().st_size < MIN_PNG_BYTES]
    rep.check(not tiny, "no suspiciously small PNG", ", ".join(tiny) or "ok")

    # Ink coverage is the check that actually catches a blank canvas; it needs
    # Pillow, which ships with matplotlib, so it is present wherever figures were
    # generated. Degrade to the size check rather than failing when it is not.
    try:
        import numpy as np
        from PIL import Image
        blank = []
        for p in pngs:
            a = np.asarray(Image.open(p).convert("L"), dtype=np.uint8)
            if float((a < 250).mean()) < 0.010:
                blank.append(p.name)
        rep.check(not blank, "no near-empty figure (ink >= 1%)",
                  ", ".join(blank) or f"{len(pngs)} checked")
    except ImportError:
        rep.warn("ink check skipped", "Pillow/numpy unavailable — size check only")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default="results", help="results/ root to verify")
    ap.add_argument("--expect-md5", help="file of 'MD5  path' lines for the tables")
    args = ap.parse_args()

    root = Path(args.root)
    print(f"VERIFYING {root.resolve()}")
    rep = Report()
    verify_kb(root, rep)
    verify_tables(root, rep, args.expect_md5)
    verify_figures(root, rep)

    print("\n" + "=" * 60)
    print(f"{rep.n - len(rep.fails)}/{rep.n} checks passed"
          + (f" · {len(rep.warns)} warning(s)" if rep.warns else ""))
    if rep.fails:
        print("FAILED:")
        for f in rep.fails:
            print(f"  - {f}")
    print("=" * 60)
    return 1 if rep.fails else 0


if __name__ == "__main__":
    sys.exit(main())
