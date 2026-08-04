#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 18 — genomic context for the KB's `strong_novel` biomarkers (NCBI nt pass).

Why this step exists
--------------------
`strong_novel` means "CPSS-stable + pyseer-significant + discriminative, but NO
curated CARD determinant". That is an explicit *knowledge gap*, not a discovery
claim (METHODOLOGY §4.4). The gap is answerable: step 08 already BLASTed the
same query FASTA against NCBI `nt` (organism-restricted via entrez_query) and
wrote `04_ncbi_blast_results_{ab}.tsv` — ~231k alignments across the 45 models.
Those hits never enter the KB (`populate_database.py` only ever writes
`source_db='card'`), so the novel biomarkers currently have no stated locus.

This step joins the two and answers, per novel unitig: *where in the genome does
this sequence sit, and what is annotated there?*

IMPORTANT — this is CONTEXT, not evidence
-----------------------------------------
The NCBI pass is restricted to the study organism and searches `nt`, a sequence
archive rather than a curated resistance catalogue. An `nt` hit therefore says
"this sequence occurs here in this species", NOT "this is a resistance gene".
Nothing here changes an `evidence_tier`; the KB is not modified. Read the output
as annotation for the Results/Discussion narrative.

Two limits to quote alongside any number from this step:
  * Step 08 ran the remote pass with `max_target_seqs=50` (the default in
    `08_blast_annotation.py`; absent from config.yaml). The retained alignments
    are therefore a BOUNDED sample of `nt` — plasmid/chromosome proportions are
    proportions of what BLAST kept, not a census of the database. (`max_target_seqs`
    is also not a "top 50 best hits" filter; it bounds the search, so borderline
    subjects can be missed.)
  * Queries are 31–106 bp. Full-length 100%-identity matches at these lengths are
    unambiguous placements, but they identify a LOCUS, never a mechanism.

Method
------
  1. Read the `strong_novel` (unitig, model) pairs from the KB. `unitig_model_scores`
     holds up to two rows per pair (gain + CPSS paths), so the CPSS frequency is
     aggregated with GROUP BY before joining — a plain JOIN double-counts.
  2. Per model, parse `02_top_{N}_features_{ab}.fasta` to map sequence -> qseqid
     (header `Rank_R|Score_S|Feature_fI`). Matching is done on the sequence, so
     the header scheme can change without breaking this step. The reverse
     complement is tried as a fallback.
  3. Stream `04_ncbi_blast_results_{ab}.tsv` (outfmt 6: qseqid sseqid pident
     length mismatch gapopen qstart qend sstart send evalue bitscore qlen stitle),
     keep the rows for those qseqids, compute query coverage = length / qlen.
  4. Rank hits by bitscore and report the best ones per unitig.

Output
------
    results/tables/novel_ncbi_context.csv    one row per novel unitig (best hit
                                             + the top-N subject titles)
    results/tables/novel_ncbi_hits.csv       every retained alignment (long form)
    stdout                                   human-readable summary

Usage (TRUBA, login node is fine — stdlib only, no pandas)
    python3 scripts/18_novel_ncbi_context.py \
        --kb $AMR_WORK/results/kb/amrk.db \
        --results-root $AMR_WORK/results \
        --out $AMR_WORK/results/tables
"""

import argparse
import csv
import sqlite3
import sys
from pathlib import Path

# outfmt 6 column order fixed by 08_blast_annotation.py (OUTFMT).
TSV_COLS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "stitle"]

_COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq):
    return seq.translate(_COMP)[::-1]


def novel_biomarkers(kb_path):
    """The strong_novel (organism, antibiotic, sequence) rows, CPSS-deduplicated."""
    con = sqlite3.connect(kb_path)
    con.row_factory = sqlite3.Row
    try:
        return con.execute("""
            WITH skor AS (
                SELECT unitig_id, model_id, MAX(selection_frequency) AS cpss_freq
                FROM unitig_model_scores GROUP BY unitig_id, model_id)
            SELECT p.organism, m.antibiotic, u.sequence,
                   LENGTH(u.sequence) AS unitig_len,
                   e.evidence_layers, s.cpss_freq
            FROM unitig_evidence_tier e
            JOIN unitigs u       ON u.unitig_id = e.unitig_id
            JOIN models m        ON m.model_id  = e.model_id
            JOIN pipeline_runs p ON p.run_id    = m.run_id
            LEFT JOIN skor s     ON s.unitig_id = e.unitig_id
                                AND s.model_id  = e.model_id
            WHERE e.evidence_tier = 'strong_novel'
            ORDER BY p.organism, m.antibiotic, u.sequence""").fetchall()
    finally:
        con.close()


def read_fasta_index(fasta_path):
    """sequence -> qseqid for one model's step-07 query FASTA."""
    index, header = {}, None
    with open(fasta_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].split()[0]  # BLAST truncates qseqid at whitespace
            elif header:
                index.setdefault(line.upper(), header)
                header = None
    return index


def replicon_stats(hits, majority=0.8):
    """Replicon context over ALL retained alignments, not just the best hit.

    `nt` subject titles name the replicon ("… chromosome, complete genome",
    "… plasmid pXYZ"), so the hit distribution answers *plasmid or chromosome?*
    A single best hit silently picks one side: 8 of the 23 novel biomarkers here
    land on both replicon types at comparable rates, which is a mobile-element
    signature (IS/transposon/integron backbone), not a bad alignment. Hence the
    call is made on the whole distribution and stays `mixed` unless one side
    holds a `majority` share.
    """
    n = len(hits)
    if not n:
        return 0, 0, None, "no_hit"
    titles = [h["stitle"].lower() for h in hits]
    n_pl = sum(1 for t in titles if "plasmid" in t)
    n_ch = sum(1 for t in titles
               if "plasmid" not in t and ("chromosome" in t or "genome" in t))
    frac = max(n_pl, n_ch) / n
    call = ("plasmid" if n_pl > n_ch else "chromosome") if frac >= majority else "mixed"
    return n_pl, n_ch, round(frac, 3), call


def hits_for(tsv_path, wanted_qseqids):
    """Alignments whose qseqid is in `wanted_qseqids`, streamed."""
    out = []
    with open(tsv_path, encoding="utf-8", errors="replace") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < len(TSV_COLS) or row[0] not in wanted_qseqids:
                continue
            r = dict(zip(TSV_COLS, row[:len(TSV_COLS) - 1] + ["\t".join(row[len(TSV_COLS) - 1:])]))
            try:
                r["bitscore"] = float(r["bitscore"])
                r["pident"] = float(r["pident"])
                qlen = float(r["qlen"])
                r["qcov"] = round(float(r["length"]) / qlen, 3) if qlen else None
            except (TypeError, ValueError):
                continue
            out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kb", required=True, help="path to amrk.db")
    ap.add_argument("--results-root", required=True,
                    help="results/ root holding {organism}/{antibiotic}/05_explainability")
    ap.add_argument("--out", required=True, help="directory for the two CSVs")
    ap.add_argument("--top", type=int, default=3,
                    help="subject titles to keep per unitig (default 3)")
    ap.add_argument("--majority", type=float, default=0.8,
                    help="share of alignments one replicon type needs before the "
                         "call is plasmid/chromosome rather than mixed (default 0.8)")
    args = ap.parse_args()

    results_root, out_dir = Path(args.results_root), Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = novel_biomarkers(args.kb)
    if not rows:
        print("No strong_novel rows in the KB — nothing to annotate.")
        return
    print(f"strong_novel biomarkers in KB : {len(rows)}")

    # Group by model so each FASTA/TSV pair is opened once.
    by_model = {}
    for r in rows:
        by_model.setdefault((r["organism"], r["antibiotic"]), []).append(r)
    print(f"models involved               : {len(by_model)}\n")

    summary, long_form, missing = [], [], []
    for (organism, antibiotic), group in sorted(by_model.items()):
        expl = results_root / organism / antibiotic / "05_explainability"
        fastas = sorted(expl.glob(f"02_top_*_features_{antibiotic}.fasta"))
        tsvs = sorted(expl.glob(f"04_ncbi_blast_results_{antibiotic}.tsv"))
        if not fastas or not tsvs:
            print(f"[!] {organism}/{antibiotic}: "
                  f"{'FASTA' if not fastas else 'NCBI TSV'} missing — skipped")
            missing.extend({**dict(r), "reason": "input file missing"} for r in group)
            continue

        seq2id = read_fasta_index(fastas[0])
        wanted = {}
        for r in group:
            seq = r["sequence"].upper()
            qid = seq2id.get(seq) or seq2id.get(revcomp(seq))
            if qid is None:
                missing.append({**dict(r), "reason": "sequence not in step-07 query FASTA"})
                continue
            wanted.setdefault(qid, []).append(r)

        hits = hits_for(tsvs[0], set(wanted)) if wanted else []
        per_q = {}
        for h in hits:
            per_q.setdefault(h["qseqid"], []).append(h)

        print(f"{organism}/{antibiotic}: {len(group)} novel · "
              f"{len(wanted)} mapped to FASTA · {len(hits)} nt alignments")

        for qid, members in wanted.items():
            hs = sorted(per_q.get(qid, []), key=lambda x: x["bitscore"], reverse=True)
            for r in members:
                base = {"organism": organism, "antibiotic": antibiotic,
                        "unitig_len": r["unitig_len"], "cpss_freq": r["cpss_freq"],
                        "evidence_layers": r["evidence_layers"],
                        "qseqid": qid, "sequence": r["sequence"], "n_nt_hits": len(hs)}
                n_pl, n_ch, frac, call = replicon_stats(hs, args.majority)
                base.update({"n_plasmid_hits": n_pl, "n_chromosome_hits": n_ch,
                             "dominant_fraction": frac, "replicon_call": call})
                if hs:
                    b = hs[0]
                    base.update({
                        "best_pident": b["pident"], "best_qcov": b["qcov"],
                        "best_evalue": b["evalue"], "best_bitscore": b["bitscore"],
                        "best_subject": b["stitle"][:300],
                        "top_subjects": " || ".join(h["stitle"][:120] for h in hs[:args.top]),
                    })
                else:
                    base.update({"best_pident": None, "best_qcov": None,
                                 "best_evalue": None, "best_bitscore": None,
                                 "best_subject": "", "top_subjects": ""})
                summary.append(base)
                for h in hs[:args.top]:
                    long_form.append({"organism": organism, "antibiotic": antibiotic,
                                      "qseqid": qid, "sequence": r["sequence"],
                                      **{c: h[c] for c in TSV_COLS if c != "qseqid"},
                                      "qcov": h["qcov"]})

    def write_csv(path, records):
        if not records:
            return
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(records[0].keys()))
            w.writeheader()
            w.writerows(records)
        print(f"  -> {path}  ({len(records)} rows)")

    print("\n" + "=" * 78)
    write_csv(out_dir / "novel_ncbi_context.csv", summary)
    write_csv(out_dir / "novel_ncbi_hits.csv", long_form)

    annotated = [s for s in summary if s["n_nt_hits"]]
    print(f"\nnovel biomarkers annotated    : {len(annotated)}/{len(rows)}")
    if missing:
        print(f"unmapped (no FASTA entry/file): {len(missing)}")
        for m in missing[:10]:
            print(f"    {m['organism']}/{m['antibiotic']} — {m['reason']}")

    if annotated:
        calls = {}
        for s in annotated:
            calls[s["replicon_call"]] = calls.get(s["replicon_call"], 0) + 1
        print("replicon calls                : "
              + " · ".join(f"{k} {v}" for k, v in sorted(calls.items())))
        print("\nBest nt hit per novel biomarker "
              f"(replicon call uses all hits, >={args.majority:.0%} majority)")
        print("-" * 78)
        for s in sorted(annotated, key=lambda x: (x["organism"], x["antibiotic"])):
            print(f"{s['organism'][:18]:18s} {s['antibiotic'][:22]:22s} "
                  f"{s['unitig_len']:>4} bp  id={s['best_pident']:.1f}% "
                  f"cov={s['best_qcov']}  E={s['best_evalue']}  "
                  f"[{s['replicon_call']} {s['dominant_fraction']:.0%} "
                  f"of {s['n_nt_hits']}]")
            print(f"    {s['best_subject'][:150]}")
    print()


if __name__ == "__main__":
    sys.exit(main())
