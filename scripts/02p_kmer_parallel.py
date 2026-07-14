#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 02 (parallel) — genome-level parallel KMC k-mer counting.

A drop-in alternative to 02_kmer_extraction.py for HPC / many-core nodes. The
serial step 02 runs one KMC process at a time, so on a 20–40 core scheduler
allocation the CPUs sit idle and the job is flagged as low-efficiency (and may
be auto-cancelled, e.g. on TRUBA). This version runs many single-threaded KMC
processes concurrently (one per genome), saturating the allocation and finishing
far faster, while writing the exact same per-genome `{genome_id}.kmc_pre/.kmc_suf`
databases that steps 02b/03 consume.

Concurrency = preprocessing.threads (set this equal to the SLURM --cpus-per-task).
Each worker runs `kmc -t1 -m4` in its own temp subdirectory (no collisions) and
is resumable (genomes already counted are skipped).

Run identically to step 02 (organism comes from config.yaml):
    python scripts/02p_kmer_parallel.py
"""

import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import get_target, load_config, resolve_path, resolve_tool  # noqa: E402


def main() -> None:
    cfg = load_config()
    # AMR_ORGANISM env overrides config (parallel per-organism KMC, like 03u).
    organism = get_target(config=cfg)[0]
    k_length = cfg["preprocessing"]["k_length"]
    workers = int(cfg["preprocessing"]["threads"])

    raw_dir = resolve_path("raw_genomes_dir", organism=organism, config=cfg)
    out_dir = resolve_path("kmc_outputs_dir", organism=organism, config=cfg)
    tmp_dir = out_dir / "tmp"
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    kmc = resolve_tool("kmc_bin", "kmc", config=cfg)
    if not kmc:
        sys.exit("ERROR: KMC not found. Install KMC (conda install -c bioconda kmc) "
                 "so `kmc` is on PATH, or set AMR_KMC_BIN.")

    genomes = sorted(raw_dir.glob("*.fna"))
    if not genomes:
        sys.exit(f"ERROR: no .fna genomes in {raw_dir}")
    print(f"Parallel KMC: {len(genomes)} genomes | workers={workers} | per-KMC -t1 -m4",
          flush=True)

    def work(genome: Path) -> str:
        gid = genome.stem
        if (out_dir / f"{gid}.kmc_pre").exists() and (out_dir / f"{gid}.kmc_suf").exists():
            return "skip"
        wtmp = tmp_dir / f"{gid}_t"
        wtmp.mkdir(parents=True, exist_ok=True)
        lst = tmp_dir / f"{gid}.lst"
        lst.write_text(str(genome) + "\n")
        try:
            r = subprocess.run(
                [kmc, f"-k{k_length}", "-m4", "-t1", "-ci1", "-fm",
                 f"@{lst}", str(out_dir / gid), str(wtmp)],
                capture_output=True, text=True)
            ok = r.returncode == 0
        finally:
            lst.unlink(missing_ok=True)
            shutil.rmtree(wtmp, ignore_errors=True)
        return "ok" if ok else "fail"

    counts = {"ok": 0, "skip": 0, "fail": 0}
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(work, g): g for g in genomes}
        for fut in as_completed(futures):
            counts[fut.result()] += 1
            done += 1
            if done % 200 == 0 or done == len(genomes):
                print(f"  {done}/{len(genomes)}  "
                      f"ok={counts['ok']} skip={counts['skip']} fail={counts['fail']}",
                      flush=True)

    print(f"DONE: ok={counts['ok']} skip={counts['skip']} fail={counts['fail']} "
          f"total={len(genomes)}", flush=True)
    if counts["fail"]:
        print("WARNING: some genomes failed KMC counting (see above).")


if __name__ == "__main__":
    main()
