#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 00a — Download + clean BV-BRC AMR data, then fetch genome assemblies.

Pipeline:
    1. Fetch the genome_amr table from the BV-BRC Data API (paginated), filtered
       to taxon_id / public / evidence="Laboratory Method".            (or --raw-csv)
    2. Clean it (EUCAST/CLSI only, Resistant/Susceptible only, antibiotic name
       normalisation, majority -> newest-year conflict resolution) via lib.bvbrc.
    3. Download each surviving genome's contigs FASTA as {genome_id}.fna from the
       BV-BRC FTP, in parallel, with retries and resume.

Outputs (organism resolved from the registry via taxon_id):
    data/external/{org}/metadata/BVBRC_genome_amr.csv     raw API/web table
    data/external/{org}/metadata/amr_cleaned_long.csv     cleaned (genome,antibiotic,label)
    data/external/{org}/metadata/download_manifest.json   provenance (query, date, counts)
    data/raw/{org}/genomes/{genome_id}.fna                assemblies
    logs/{org}/00_download.log                            full log
    logs/{org}/cleaning_report.json                       per-step row/pair counts
    logs/{org}/download_report.csv                        per-genome status

Examples:
    python scripts/00a_download_bvbrc.py --organism ecoli --max-genomes 20   # dry run
    python scripts/00a_download_bvbrc.py --organism ecoli                    # full
    python scripts/00a_download_bvbrc.py --organism ecoli --retry-failed     # retry failures
    python scripts/00a_download_bvbrc.py --organism ecoli --raw-csv path.csv # skip API fetch

NOTE: the BV-BRC Data API caps deep pagination. If the automatic fetch is
truncated for a very large table, download the filtered table from the website
("DOWNLOAD" on the AMR Phenotypes tab, including the Testing Standard column) and
pass it with --raw-csv.

CLI BACKEND — IMPORTANT: sourcing the BV-BRC environment
(`source /Applications/BV-BRC.app/user-env.sh`) puts the p3-* tools on PATH but
ALSO replaces `python` with BV-BRC's bundled (old) interpreter. Run this script
with your project's Python explicitly, e.g.:

    source /Applications/BV-BRC.app/user-env.sh        # p3-* on PATH
    /path/to/conda/envs/<env>/bin/python scripts/00a_download_bvbrc.py --backend cli ...

The script runs under your Python (with pandas); the p3-* subprocess calls
inherit the sourced PATH / Perl env and work normally.
"""

import argparse
import csv
import datetime
import io
import json
import logging
import shutil
import subprocess
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import registry                                  # noqa: E402
from lib.bvbrc import clean_amr_table, standardise_columns  # noqa: E402
from lib.config import load_config, resolve_path           # noqa: E402

API_URL = "https://www.bv-brc.org/api/genome_amr/"
FTP_FASTA = "https://ftp.bv-brc.org/genomes/{gid}/{gid}.fna"
API_BATCH = 25000
SELECT_FIELDS = ["genome_id", "genome_name", "antibiotic", "resistant_phenotype",
                 "testing_standard", "testing_standard_year"]

log = logging.getLogger("bvbrc")


# ---------------------------------------------------------------------------
def setup_logging(log_path):
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log.setLevel(logging.INFO)
    log.handlers.clear()
    fmt = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    fh = logging.FileHandler(log_path, encoding="utf-8"); fh.setFormatter(fmt)
    sh = logging.StreamHandler(); sh.setFormatter(fmt)
    log.addHandler(fh); log.addHandler(sh)


def _http_get(url, timeout=120):
    req = urllib.request.Request(url, headers={"Accept": "text/csv",
                                               "User-Agent": "amr-pipeline/00a"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def fetch_amr_table(taxid):
    """Fetch the full filtered genome_amr table via paginated API requests."""
    base_rql = (
        f"and(eq(taxon_id,{taxid}),eq(public,true),"
        f'eq(evidence,"Laboratory Method"))'
        f"&select({','.join(SELECT_FIELDS)})&sort(+genome_id)"
    )
    frames, start = [], 0
    while True:
        rql = f"{base_rql}&limit({API_BATCH},{start})"
        url = API_URL + "?" + rql.replace(" ", "%20")
        log.info(f"API fetch: rows {start}..{start + API_BATCH}")
        try:
            text = _http_get(url)
        except Exception as e:
            log.error(f"API request failed at start={start}: {e}")
            log.error("If the table is large, use the website DOWNLOAD and pass --raw-csv.")
            break
        batch = pd.read_csv(io.StringIO(text)) if text.strip() else pd.DataFrame()
        if batch.empty:
            break
        frames.append(batch)
        if len(batch) < API_BATCH:
            break
        start += API_BATCH
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---- CLI backend (BV-BRC p3-* tools; handles large result sets) -----------
def cli_available():
    return bool(shutil.which("p3-get-genome-drugs") and shutil.which("p3-genome-fasta"))


def fetch_amr_table_cli(taxid, limit_genomes=0):
    """
    Fetch genome_amr via the BV-BRC CLI in two stages.

    Stage 1 lists the taxon's public genome_ids (fast — ids only). Stage 2 looks
    up the AMR drug records for those genomes. When ``limit_genomes`` > 0 only the
    first N genomes are queried, so a dry run finishes in seconds instead of
    pulling drug records for all ~85k E. coli genomes.
    """
    # Stage 1: genome id list
    log.info(f"CLI stage 1/2: listing public genome_ids for taxon {taxid}...")
    r1 = subprocess.run(
        ["p3-all-genomes", "--eq", f"taxon_id,{taxid}", "--eq", "public,true", "--attr", "genome_id"],
        capture_output=True, text=True)
    if r1.returncode != 0 or not r1.stdout.strip():
        log.error(f"p3-all-genomes failed: {(r1.stderr or '')[:300]}")
        return pd.DataFrame()
    lines = r1.stdout.splitlines()
    header, ids = lines[0], lines[1:]
    log.info(f"  found {len(ids)} genomes")
    if limit_genomes and limit_genomes < len(ids):
        ids = ids[:limit_genomes]
        log.info(f"  [dry-run] querying drug records for first {len(ids)} genomes only")

    # Stage 2: drug records for those genomes (genome_ids fed via stdin)
    log.info(f"CLI stage 2/2: p3-get-genome-drugs for {len(ids)} genomes "
             f"(this is the slow step)...")
    cmd2 = ["p3-get-genome-drugs", "--eq", "evidence,Laboratory Method",
            "--attr", "genome_id", "--attr", "genome_name", "--attr", "antibiotic",
            "--attr", "resistant_phenotype", "--attr", "testing_standard",
            "--attr", "testing_standard_year", "--attr", "evidence"]
    stdin_text = "\n".join([header] + ids) + "\n"
    r2 = subprocess.run(cmd2, input=stdin_text, capture_output=True, text=True)
    if r2.returncode != 0:
        log.error(f"p3-get-genome-drugs failed: {(r2.stderr or '')[:500]}")
        return pd.DataFrame()
    if not r2.stdout.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(r2.stdout), sep="\t")  # p3 tools emit TSV


def download_one_cli(gid, dest_dir, retries=3):
    """Download {gid}.fna via `p3-genome-fasta`. Returns (gid, status, size, error)."""
    out = dest_dir / f"{gid}.fna"
    if out.exists() and out.stat().st_size > 0:
        return gid, "skipped", out.stat().st_size, ""
    last = ""
    for _ in range(retries):
        try:
            r = subprocess.run(["p3-genome-fasta", str(gid)],
                               capture_output=True, text=True, timeout=180)
            if r.returncode == 0 and r.stdout.lstrip().startswith(">"):
                out.write_text(r.stdout, encoding="utf-8")
                return gid, "downloaded", len(r.stdout), ""
            last = (r.stderr or "empty/non-FASTA").strip()[:200]
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
    return gid, "failed", 0, last


# ---------------------------------------------------------------------------
def download_one(gid, dest_dir, retries=3, timeout=120):
    """Download {gid}.fna. Returns (gid, status, size, error)."""
    out = dest_dir / f"{gid}.fna"
    if out.exists() and out.stat().st_size > 0:
        return gid, "skipped", out.stat().st_size, ""
    url = FTP_FASTA.format(gid=gid)
    last = ""
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "amr-pipeline/00a"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                data = r.read()
            if not data or not data.lstrip().startswith(b">"):
                last = "empty or non-FASTA response"
                continue
            out.write_bytes(data)
            return gid, "downloaded", len(data), ""
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
    return gid, "failed", 0, last


def download_all(genome_ids, dest_dir, workers, report_path, download_fn=download_one):
    dest_dir.mkdir(parents=True, exist_ok=True)
    results = []
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(download_fn, gid, dest_dir): gid for gid in genome_ids}
        done = 0
        for fut in as_completed(futs):
            gid, status, size, err = fut.result()
            results.append({"genome_id": gid, "status": status, "size": size, "error": err})
            done += 1
            if status == "failed":
                log.warning(f"  download failed: {gid} ({err})")
            if done % 100 == 0 or done == len(genome_ids):
                log.info(f"  downloaded {done}/{len(genome_ids)}")
    _write_report(results, report_path)
    counts = pd.DataFrame(results)["status"].value_counts().to_dict()
    log.info(f"Download status: {counts}")
    return results


def _write_report(results, report_path):
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["genome_id", "status", "size", "error"])
        w.writeheader()
        w.writerows(results)


def retry_failed(report_path, dest_dir, workers, download_fn=download_one):
    if not report_path.exists():
        log.error(f"No download report at {report_path}; nothing to retry.")
        return
    prev = pd.read_csv(report_path, dtype={"genome_id": str})
    failed = prev[prev["status"] == "failed"]["genome_id"].tolist()
    if not failed:
        log.info("No failed downloads to retry.")
        return
    log.info(f"Retrying {len(failed)} failed downloads...")
    new = download_all(failed, dest_dir, workers, dest_dir.parent / "_retry_tmp.csv", download_fn)
    # merge: replace prior rows for retried ids
    new_df = pd.DataFrame(new)
    merged = pd.concat([prev[~prev["genome_id"].isin(failed)], new_df], ignore_index=True)
    _write_report(merged.to_dict("records"), report_path)
    (dest_dir.parent / "_retry_tmp.csv").unlink(missing_ok=True)


# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description="Download + clean BV-BRC AMR data and assemblies.")
    p.add_argument("--organism", default=None, help="registry slug (default: config project.organism)")
    p.add_argument("--raw-csv", default=None, help="use an existing raw AMR CSV instead of the API")
    p.add_argument("--max-genomes", type=int, default=0,
                   help="cap candidate genomes (0=all). With the CLI backend this also "
                        "limits the (slow) drug-record fetch, so dry runs finish fast.")
    p.add_argument("--workers", type=int, default=12, help="parallel download threads")
    p.add_argument("--skip-download", action="store_true", help="fetch+clean only; no FASTA download")
    p.add_argument("--retry-failed", action="store_true", help="re-attempt only previously failed downloads")
    p.add_argument("--backend", choices=["auto", "cli", "api"], default="auto",
                   help="data transport: BV-BRC CLI (p3-*), HTTP API, or auto-detect (default)")
    args = p.parse_args()

    config = load_config()
    organism = args.organism or config.get("project", {}).get("organism", "ecoli")
    try:
        org = registry.get_organism(organism)
    except KeyError as e:
        print(e); sys.exit(1)
    taxid = org.get("taxid")

    # Choose transport. CLI is preferred (no deep-pagination limit; official
    # p3-genome-fasta for assemblies); falls back to the HTTP API if the p3 tools
    # are not on PATH.
    use_cli = args.backend == "cli" or (args.backend == "auto" and cli_available())
    if args.backend == "cli" and not cli_available():
        print("ERROR: --backend cli requested but p3-* tools not found on PATH.")
        sys.exit(1)
    download_fn = download_one_cli if use_cli else download_one

    genomes_dir = resolve_path("raw_genomes_dir", organism=organism, config=config)
    meta_dir = resolve_path("metadata_file", organism=organism, config=config).parent
    logs_dir = resolve_path("logs_dir", organism=organism, antibiotic="_global", config=config).parent
    meta_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(logs_dir / "00_download.log")

    raw_csv = meta_dir / "BVBRC_genome_amr.csv"
    cleaned_csv = meta_dir / "amr_cleaned_long.csv"
    dl_report = logs_dir / "download_report.csv"

    log.info("=" * 70)
    log.info(f"BV-BRC download — organism={organism} taxid={taxid} "
             f"backend={'cli' if use_cli else 'api'}")
    log.info("=" * 70)

    if args.retry_failed:
        retry_failed(dl_report, genomes_dir, args.workers, download_fn)
        log.info("Retry complete. Re-run 00_prepare_metadata.py to refresh the matrix.")
        return

    # 1) obtain raw table
    if args.raw_csv:
        log.info(f"Loading raw table from {args.raw_csv}")
        raw = pd.read_csv(args.raw_csv, sep=None, engine="python", dtype=str)
    elif use_cli:
        log.info("Fetching genome_amr via BV-BRC CLI...")
        raw = fetch_amr_table_cli(taxid, limit_genomes=args.max_genomes)
        if raw.empty:
            log.error("CLI fetch returned nothing. Check p3-login/connectivity or use --backend api / --raw-csv.")
            sys.exit(1)
        standardise_columns(raw).to_csv(raw_csv, index=False)
        log.info(f"Raw table saved: {raw_csv} ({len(raw)} rows)")
    else:
        log.info("Fetching genome_amr from BV-BRC HTTP API...")
        raw = fetch_amr_table(taxid)
        if raw.empty:
            log.error("No rows fetched. Provide --raw-csv from the website DOWNLOAD instead.")
            sys.exit(1)
        standardise_columns(raw).to_csv(raw_csv, index=False)
        log.info(f"Raw table saved: {raw_csv} ({len(raw)} rows)")

    # 2) clean
    cleaned, report = clean_amr_table(raw)
    cleaned.to_csv(cleaned_csv, index=False)
    (logs_dir / "cleaning_report.json").write_text(json.dumps(report, indent=2))
    log.info(f"Cleaning report: {report}")
    log.info(f"Cleaned long table: {cleaned_csv} "
             f"({report['n_genomes']} genomes × {report['n_antibiotics']} antibiotics)")

    genome_ids = sorted(cleaned["genome_id"].astype(str).unique())
    if args.max_genomes and args.max_genomes < len(genome_ids):
        genome_ids = genome_ids[:args.max_genomes]
        log.info(f"[dry-run] limiting to first {len(genome_ids)} genomes")

    # provenance manifest
    manifest = {
        "organism": organism, "taxon_id": taxid,
        "source": "BV-BRC genome_amr",
        "backend": "raw-csv" if args.raw_csv else ("cli" if use_cli else "api"),
        "filter": 'eq(public,true) & eq(evidence,"Laboratory Method")',
        "cleaning": "EUCAST/CLSI only; Resistant/Susceptible only; name-normalised; "
                    "majority->newest-year conflict resolution",
        "downloaded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "n_candidate_genomes": len(genome_ids),
        "cleaning_report": report,
    }
    (meta_dir / "download_manifest.json").write_text(json.dumps(manifest, indent=2))

    # 3) download assemblies
    if args.skip_download:
        log.info("--skip-download set; not fetching FASTA.")
    else:
        log.info(f"Downloading {len(genome_ids)} assemblies -> {genomes_dir}")
        download_all(genome_ids, genomes_dir, args.workers, dl_report, download_fn)

    log.info("Done. Next: python scripts/00_prepare_metadata.py "
             "(builds amr_phenotypes.csv from the present .fna files).")


if __name__ == "__main__":
    main()
