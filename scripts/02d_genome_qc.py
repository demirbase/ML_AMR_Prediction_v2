#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Assembly-level genome QC — CheckM2 + QUAST (must-have M15; ROADMAP §0.2).

Hybrid genome QC on the raw assemblies, complementing 02b's k-mer complexity
IQR advisory: **CheckM2** (completeness / contamination via a machine-learning
model over DIAMOND hits) + **QUAST** (N50, contig count, total length). This is
an *organism-level* step — it QCs every ``{genome_id}.fna`` once and the result
covers all antibiotics.

CheckM2 lives in ``amr-checkm2.sif`` (pins python<3.9) and QUAST in
``amr-tools.sif``; the config-driven Python prep/post needs ``amr.sif``
(yaml+pandas). So — like 14_pyseer_lmm.py — this script runs only the Python
halves and the SLURM job chains the tool CLIs between them:

    amr.sif:       02d_genome_qc.py --mode prep   # genome dir + out dirs -> paths.sh
    amr-checkm2.sif: checkm2 predict --input $GENOMES_DIR -x fna \
                        --output-directory $CHECKM2_OUT --threads N
    amr-tools.sif: quast.py $GENOMES_DIR/*.fna -o $QUAST_OUT --threads N --no-plots
    amr.sif:       02d_genome_qc.py --mode post   # merge + threshold + report

``prep`` writes ``02d_qc_paths_{organism}.sh`` (GENOMES_DIR / CHECKM2_OUT /
QUAST_OUT / N_GENOMES) so the SLURM script sources the config-resolved paths
instead of hardcoding them.

``post`` merges CheckM2 ``quality_report.tsv`` + QUAST ``transposed_report.tsv``,
applies the thresholds, and writes a per-genome table, a summary JSON (pass rate
+ distributions — the Methods "data quality" statement), and an exclusion list
of failing genome IDs. The exclusion list is *advisory*: QC does NOT rebuild the
matrix or retrain — only re-run the pipeline (03u -> training) with these
excluded if a material fraction fails and you choose stricter data.

Usage:
  # on amr.sif
  python scripts/02d_genome_qc.py --mode prep   [--organism ecoli]
  # ... run CheckM2 + QUAST (SLURM, their containers) ...
  python scripts/02d_genome_qc.py --mode post   [--completeness-min 95 --contamination-max 5 \
                                                  --n50-min 50000 --max-contigs 500]
"""

import argparse
import datetime
import json
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config, resolve_path  # noqa: E402
from lib.logging_utils import get_logger  # noqa: E402

logger = get_logger("m15-genome-qc")


def _qc_dirs(organism, config):
    """Resolve the (genomes_dir, qc_out_dir) pair, both organism-scoped."""
    genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    base = resolve_path("dir_global_exploration", organism=organism, config=config)
    qc_out = Path(base) / "genome_qc"
    return Path(genomes_dir), qc_out


def do_prep(organism, config):
    genomes_dir, qc_out = _qc_dirs(organism, config)
    fna = sorted(genomes_dir.glob("*.fna"))
    if not fna:
        logger.error("no .fna assemblies under %s", genomes_dir)
        sys.exit(1)
    qc_out.mkdir(parents=True, exist_ok=True)
    checkm2_out = qc_out / "checkm2"
    quast_out = qc_out / "quast"

    paths_sh = qc_out / f"02d_qc_paths_{organism}.sh"
    with open(paths_sh, "w", encoding="utf-8") as f:
        f.write(f'GENOMES_DIR="{genomes_dir}"\n')
        f.write(f'CHECKM2_OUT="{checkm2_out}"\n')
        f.write(f'QUAST_OUT="{quast_out}"\n')
        f.write(f'N_GENOMES={len(fna)}\n')
    logger.info("prep: %d assemblies in %s", len(fna), genomes_dir)
    logger.info("prep: CheckM2 -> %s | QUAST -> %s", checkm2_out, quast_out)
    logger.info("prep: SLURM sources -> %s", paths_sh)


def _read_checkm2(checkm2_out):
    """Return {genome_id: (completeness, contamination)} from quality_report.tsv."""
    f = Path(checkm2_out) / "quality_report.tsv"
    if not f.exists():
        return None
    # Read 'Name' as str — genome ids like "562.100" must not be parsed as
    # floats (562.1) or the merge with genome ids silently breaks.
    df = pd.read_csv(f, sep="\t", dtype={"Name": str})
    # CheckM2 'Name' is the file stem (== genome_id); columns Completeness/Contamination.
    out = {}
    for _, r in df.iterrows():
        out[str(r["Name"])] = (float(r["Completeness"]), float(r["Contamination"]))
    return out


def _read_quast(quast_out):
    """Return {genome_id: (n50, n_contigs, total_len)} from transposed_report.tsv."""
    f = Path(quast_out) / "transposed_report.tsv"
    if not f.exists():
        return None
    # 'Assembly' as str — same genome-id float-parsing hazard as CheckM2.
    df = pd.read_csv(f, sep="\t", dtype={"Assembly": str})
    # QUAST 'Assembly' is the file stem; robust column lookup (names vary slightly).
    def col(*cands):
        for c in cands:
            if c in df.columns:
                return c
        return None
    c_n50 = col("N50")
    c_ctg = col("# contigs", "# contigs (>= 0 bp)")
    c_len = col("Total length", "Total length (>= 0 bp)")
    out = {}
    for _, r in df.iterrows():
        out[str(r["Assembly"])] = (float(r[c_n50]), int(r[c_ctg]), int(r[c_len]))
    return out


def classify_row(gid, comp, cont, n50, nctg, tlen, thr):
    """Apply thresholds to one genome's metrics -> row dict.

    A check passes if its metric is present AND within threshold; a MISSING
    metric is ``None`` (its tool didn't run) and is not counted as a failure.

    ``pass_overall`` (the EXCLUSION gate) is CheckM2 completeness+contamination ONLY.
    QUAST N50/contigs are reported for information but are ADVISORY, not exclusion
    criteria: unitig features are ~30-60 bp and survive fragmented assemblies intact,
    so contiguity says little about presence/absence AMR content. Gating on N50 (>=50 kb)
    dropped ~63% of E. faecium (N50<50 kb draft assemblies whose AMR genes are fine),
    which would gut the VRE panel. CheckM2 completeness (missing content) and
    contamination (chimeric/mixed) are the quality signals that actually matter."""
    pass_comp = None if comp is None else comp >= thr["completeness_min"]
    pass_cont = None if cont is None else cont <= thr["contamination_max"]
    pass_n50 = None if n50 is None else n50 >= thr["n50_min"]   # advisory only
    pass_ctg = None if nctg is None else nctg <= thr["max_contigs"]  # advisory only
    # Fail ONLY if a CheckM2 gate is affirmatively violated. A missing metric (None,
    # e.g. CheckM2 didn't assess this genome) does not fail — "missing metric is not a
    # failure" — and N50/contigs never gate (advisory). In the real run CheckM2 is
    # present for every genome, so this is exactly "completeness>=95 AND contamination<=5".
    overall = (pass_comp is not False) and (pass_cont is not False)
    return {
        "genome_id": gid, "completeness": comp, "contamination": cont,
        "n50": n50, "n_contigs": nctg, "total_length": tlen,
        "pass_completeness": pass_comp, "pass_contamination": pass_cont,
        "pass_n50": pass_n50, "pass_contigs": pass_ctg, "pass_overall": overall,
    }


def do_post(organism, config, thr):
    genomes_dir, qc_out = _qc_dirs(organism, config)
    checkm2 = _read_checkm2(qc_out / "checkm2")
    quast = _read_quast(qc_out / "quast")
    if checkm2 is None and quast is None:
        logger.error("neither CheckM2 nor QUAST reports found under %s — run the "
                     "tools first (see module docstring / SLURM).", qc_out)
        sys.exit(1)
    if checkm2 is None:
        logger.warning("CheckM2 report missing — completeness/contamination skipped.")
    if quast is None:
        logger.warning("QUAST report missing — N50/contig checks skipped.")

    ids = sorted(set(checkm2 or {}) | set(quast or {}))
    rows = []
    for gid in ids:
        comp, cont = (checkm2 or {}).get(gid, (None, None))
        n50, nctg, tlen = (quast or {}).get(gid, (None, None, None))
        rows.append(classify_row(gid, comp, cont, n50, nctg, tlen, thr))
    df = pd.DataFrame(rows)

    qc_out.mkdir(parents=True, exist_ok=True)
    csv_path = qc_out / f"02d_genome_qc_{organism}.csv"
    df.to_csv(csv_path, index=False, encoding="utf-8")

    n = len(df)
    n_pass = int(df["pass_overall"].sum())
    fails = df.loc[~df["pass_overall"], "genome_id"].tolist()
    exclude_path = qc_out / f"02d_genome_qc_exclude_{organism}.txt"
    with open(exclude_path, "w", encoding="utf-8") as f:
        f.write("\n".join(fails) + ("\n" if fails else ""))

    # The canonical outlier file 03u actually consumes:
    # dir_global_exploration/global_qc_outliers.csv, a CSV with a 'Genome' column.
    # 02d previously wrote ONLY the exclude .txt above (different name, in the qc/
    # subdir, headerless) — 03u looks for global_qc_outliers.csv one level up, so it
    # never found it and the CheckM2/QUAST QC was computed but SILENTLY NOT applied
    # to the models. Write the exact artifact 03u reads (an empty-but-headed file
    # when nothing fails, so 03u reads it cleanly and excludes zero).
    outliers_path = qc_out.parent / "global_qc_outliers.csv"
    pd.DataFrame({"Genome": fails}).to_csv(outliers_path, index=False, encoding="utf-8")
    logger.info("post: wrote %d QC outliers -> %s", len(fails), outliers_path)

    def _dist(colname):
        s = df[colname].dropna()
        if s.empty:
            return None
        return {"min": float(s.min()), "median": float(s.median()),
                "mean": round(float(s.mean()), 3), "max": float(s.max())}

    summary = {
        "organism": organism,
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "thresholds": thr,
        "n_genomes": n,
        "n_pass": n_pass,
        "n_fail": n - n_pass,
        "pass_rate": round(n_pass / n, 4) if n else None,
        "n_fail_completeness": int((df["pass_completeness"] == False).sum()),  # noqa: E712
        "n_fail_contamination": int((df["pass_contamination"] == False).sum()),  # noqa: E712
        "n_fail_n50": int((df["pass_n50"] == False).sum()),  # noqa: E712
        "n_fail_contigs": int((df["pass_contigs"] == False).sum()),  # noqa: E712
        "distributions": {k: _dist(k) for k in
                          ("completeness", "contamination", "n50", "n_contigs")},
        "tools": {"checkm2": checkm2 is not None, "quast": quast is not None},
    }
    summary_path = qc_out / f"02d_genome_qc_summary_{organism}.json"
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    logger.info("post: %d genomes | pass %d (%.1f%%) | fail %d",
                n, n_pass, 100 * n_pass / n if n else 0, n - n_pass)
    logger.info("post: fail breakdown — completeness %d, contamination %d, N50 %d, contigs %d",
                summary["n_fail_completeness"], summary["n_fail_contamination"],
                summary["n_fail_n50"], summary["n_fail_contigs"])
    logger.info("✓ %s", csv_path)
    logger.info("✓ %s", summary_path)
    logger.info("✓ %s (%d IDs — advisory exclusion list)", exclude_path, len(fails))


def main():
    config = load_config()
    ap = argparse.ArgumentParser(description="Genome assembly QC — CheckM2 + QUAST (M15).")
    ap.add_argument("--mode", choices=["prep", "post"], required=True)
    ap.add_argument("--organism", default=config.get("project", {}).get("organism", "ecoli"))
    # Thresholds (ROADMAP §0.2): completeness >=95-99%, contamination <=5% (strict
    # 2-3%), N50 >=50kb, contig upper bound. Defaults are the lenient end; tighten
    # per the write-up.
    ap.add_argument("--completeness-min", type=float, default=95.0)
    ap.add_argument("--contamination-max", type=float, default=5.0)
    ap.add_argument("--n50-min", type=int, default=50000)
    ap.add_argument("--max-contigs", type=int, default=500)
    args = ap.parse_args()

    if args.mode == "prep":
        do_prep(args.organism, config)
    else:
        thr = {"completeness_min": args.completeness_min,
               "contamination_max": args.contamination_max,
               "n50_min": args.n50_min, "max_contigs": args.max_contigs}
        do_post(args.organism, config, thr)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        print(f"\nFATAL ERROR: {e}")
        traceback.print_exc()
        sys.exit(1)
