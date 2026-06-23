#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Annotation Orchestrator — Step 08

This script coordinates the biological validation of the top k-mer features
identified in Step 07 (07_explainability.py) by delegating to a Nextflow
pipeline (08_blast_pipeline.nf) that runs two parallel BLAST searches:

  1. CARD Local BLAST:
     Queries the Comprehensive Antibiotic Resistance Database (CARD).
     Directly tests whether the top k-mers overlap with documented
     resistance genes (e.g., gyrA, parC for fluoroquinolones).
     Requires a pre-built local blastn database.

  2. NCBI Remote BLAST:
     Queries the full NCBI nucleotide (nt) database over the internet.
     Captures novel or uncharacterised resistance determinants not yet
     in CARD. Uses -remote flag — no local database needed.

Why Nextflow?
    Both BLAST searches are embarrassingly parallel and independent.
    Nextflow manages the parallel execution, retries, and output staging
    automatically, while this Python script provides the project-standard
    CLI experience (config loading, step printing, ✓ checkmarks).

Output Files (inside analysis_results/{antibiotic}/05_explainability/):
    03_card_blast_results_{antibiotic}.tsv   — CARD local hits
    04_ncbi_blast_results_{antibiotic}.tsv   — NCBI remote hits

Prerequisite Setup (CARD database):
    Download CARD nucleotide FASTA:
        wget https://card.mcmaster.ca/latest/data/nucleotide_fasta_protein_homolog_model.fasta
    Build blastn database:
        makeblastdb -in <file>.fasta -dbtype nucl -out data/blast_db/card_nt/card
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import subprocess
import shutil
import sys
import os
import yaml
from pathlib import Path

# Ensure the conda environment's bin is on PATH so shutil.which() finds
# blastn even when this script is launched via the full python interpreter
# path (which bypasses the activated environment's PATH export).
_conda_bin = Path(sys.executable).parent
os.environ['PATH'] = str(_conda_bin) + os.pathsep + os.environ.get('PATH', '')


# ============================================================================
# LOAD CONFIGURATION FROM YAML
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH  = PROJECT_ROOT / "config" / "config.yaml"

if not CONFIG_PATH.exists():
    raise FileNotFoundError(
        f"Configuration file not found: {CONFIG_PATH}\n"
        f"Please ensure config.yaml exists in the config/ directory."
    )

with open(CONFIG_PATH, 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)

# Extract project-level identifiers
TARGET_ANTIBIOTIC = config['project']['target_antibiotic']
ORGANISM          = config.get('project', {}).get('organism', 'ecoli')
TOP_N             = config['analysis']['top_n_features']

# Organism-aware path resolution (SCALE_MLOPS_PLAN §4.2)
from lib.config import resolve_path

# Resolve BLAST parameters from config
blast_cfg   = config.get('blast', {})
CARD_DB_DIR = PROJECT_ROOT / blast_cfg.get('card_db_dir', 'data/blast_db/card_nt')
CARD_DB     = CARD_DB_DIR / blast_cfg.get('card_db_name', 'card')
EVALUE      = blast_cfg.get('evalue',    10)
WORD_SIZE   = blast_cfg.get('word_size', 11)
THREADS     = blast_cfg.get('threads',   8)
# BLAST task is chosen from the ACTUAL query length (see choose_blast_task): the
# 'blastn-short' params are tuned for queries <50 bp (k-mers AND short unitigs),
# 'blastn' for longer. Picking by feature type was wrong — unitigs can be short
# (~30-50 bp), where 'blastn' finds nothing. blast.task overrides the auto choice.
BLAST_TASK_OVERRIDE = blast_cfg.get('task')


def choose_blast_task(fasta_path, override=None, short_max=50):
    """Pick the BLAST task from the MEDIAN query length: 'blastn-short' when the
    bulk of queries are short (median < short_max bp), else 'blastn'. Median (not
    max) because a few long queries shouldn't force 'blastn' on a short-dominated
    set — 'blastn-short' (with word_size 7) finds full-length hits even for the
    longer ones, whereas 'blastn' misses the short ones. ``override`` wins."""
    if override:
        return override
    lens = []
    try:
        for line in open(fasta_path, encoding='utf-8'):
            s = line.strip()
            if s and not s.startswith('>'):
                lens.append(len(s))
    except OSError:
        return 'blastn-short'
    if not lens:
        return 'blastn-short'
    lens.sort()
    median = lens[len(lens) // 2]
    return 'blastn-short' if median < short_max else 'blastn'

# Resolve I/O paths (organism-aware)
EXPLAINABILITY_DIR = resolve_path('dir_05_explainability', organism=ORGANISM,
                                  antibiotic=TARGET_ANTIBIOTIC, config=config)
# Filename must track top_n_features from config — 07 writes 02_top_{TOP_N}_features.
# Hardcoding 50 silently broke this step whenever top_n_features != 50.
FASTA_INPUT = EXPLAINABILITY_DIR / f"02_top_{TOP_N}_features_{TARGET_ANTIBIOTIC}.fasta"

# Nextflow pipeline path
PIPELINE_PATH = PROJECT_ROOT / "scripts" / "08_blast_pipeline.nf"

# Expected output files (for final confirmation print)
CARD_OUT  = EXPLAINABILITY_DIR / f"03_card_blast_results_{TARGET_ANTIBIOTIC}.tsv"
NCBI_OUT  = EXPLAINABILITY_DIR / f"04_ncbi_blast_results_{TARGET_ANTIBIOTIC}.tsv"


# ============================================================================
# MAIN ORCHESTRATION FUNCTION
# ============================================================================
def main() -> None:
    """
    Orchestrate the BLAST annotation pipeline for AMR k-mer features.

    Workflow:
        1. Validate tool availability (nextflow, blastn)
        2. Validate input FASTA from Step 07
        3. Validate CARD local database
        4. Execute Nextflow pipeline (CARD + NCBI in parallel)
        5. Confirm output files were created
    """
    print("=" * 80)
    print(f"BLAST ANNOTATION: {TARGET_ANTIBIOTIC.upper()} — K-MER BIOLOGICAL VALIDATION")
    print("=" * 80)
    print(f"  Target antibiotic : {TARGET_ANTIBIOTIC}")
    print(f"  Top-N features    : {TOP_N}")
    print(f"  E-value threshold : {EVALUE}")
    print(f"  Word size         : {WORD_SIZE}")
    print(f"  BLAST threads     : {THREADS}")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # STEP 1: Validate required tools
    # -------------------------------------------------------------------------
    print("\n[STEP 1/4] Checking required tool availability...")

    missing_tools = []
    for tool in ("nextflow", "blastn"):
        path = shutil.which(tool)
        if path:
            print(f"  ✓ {tool:12s} found: {path}")
        else:
            print(f"  ✗ {tool:12s} NOT FOUND")
            missing_tools.append(tool)

    if missing_tools:
        print("\nERROR: The following required tools are not installed or not on PATH:")
        for tool in missing_tools:
            if tool == "nextflow":
                print("  • nextflow  → Install: https://www.nextflow.io/docs/latest/getstarted.html")
                print("                Quick:   curl -s https://get.nextflow.io | bash && mv nextflow /usr/local/bin/")
            elif tool == "blastn":
                print("  • blastn    → Install BLAST+: https://www.ncbi.nlm.nih.gov/books/NBK569861/")
                print("                macOS:   brew install blast")
                print("                conda:   conda install -c bioconda blast")
        sys.exit(1)

    # -------------------------------------------------------------------------
    # STEP 2: Validate input FASTA from Step 07
    # -------------------------------------------------------------------------
    print("\n[STEP 2/4] Validating input FASTA from Step 07...")

    if not FASTA_INPUT.exists():
        print(f"  ✗ FASTA not found: {FASTA_INPUT}")
        print(f"\n  Run feature extraction first:")
        print(f"    python scripts/07_explainability.py")
        sys.exit(1)

    fasta_lines = FASTA_INPUT.read_text(encoding='utf-8').strip().splitlines()
    seq_count   = sum(1 for l in fasta_lines if l.startswith('>'))
    print(f"  ✓ FASTA input     : {FASTA_INPUT.name}")
    print(f"  ✓ Sequences       : {seq_count}")

    # -------------------------------------------------------------------------
    # STEP 3: Validate CARD local database
    # -------------------------------------------------------------------------
    print("\n[STEP 3/4] Validating CARD local database...")

    # A blastn nucleotide DB may be single-volume (<db>.nhr) or split into
    # multiple volumes (<db>.00.nhr, <db>.01.nhr, ...) described by an alias
    # file (<db>.nal). Checking only ".nhr" gave a false "missing" verdict on
    # multi-volume CARD databases. Accept any of these layouts.
    card_present = (
        (CARD_DB.parent / (CARD_DB.name + ".nhr")).exists()
        or (CARD_DB.parent / (CARD_DB.name + ".nal")).exists()
        or any(CARD_DB.parent.glob(CARD_DB.name + ".*.nhr"))
    )
    if not card_present:
        print(f"  ⚠ CARD database not found at: {CARD_DB}")
        print(f"    CARD local BLAST will fail. To build the database:")
        print(f"      1. Download: https://card.mcmaster.ca/download")
        print(f"      2. makeblastdb -in <card.fna> -dbtype nucl -out {CARD_DB}")
        print(f"    Continuing — NCBI remote BLAST will still run.\n")
    else:
        print(f"  ✓ CARD database   : {CARD_DB}")

        # Record the CARD database provenance for reproducibility (was P-14:
        # "CARD version not recorded"). blastdbcmd -info reports the build date
        # and sequence counts; we persist it next to the BLAST outputs so the
        # exact database snapshot used can be cited in Methods.
        try:
            info = subprocess.run(
                ["blastdbcmd", "-db", str(CARD_DB), "-info"],
                capture_output=True, text=True, check=False
            )
            if info.returncode == 0 and info.stdout.strip():
                EXPLAINABILITY_DIR.mkdir(parents=True, exist_ok=True)
                version_file = EXPLAINABILITY_DIR / "card_db_version.txt"
                version_file.write_text(info.stdout, encoding='utf-8')
                first_line = info.stdout.strip().splitlines()[0]
                print(f"    ↳ CARD DB info recorded: {version_file.name} ({first_line})")
        except Exception as e:
            print(f"    ⚠ Could not record CARD DB version: {e}")

    # -------------------------------------------------------------------------
    # STEP 4: Execute Nextflow pipeline
    # -------------------------------------------------------------------------
    print("\n[STEP 4/4] Launching Nextflow pipeline (CARD + NCBI in parallel)...")
    print("=" * 80)

    blast_task = choose_blast_task(FASTA_INPUT, BLAST_TASK_OVERRIDE)
    # blastn-short needs a small word_size (7) to seed short queries; the config
    # word_size (11+) truncated/missed full-length hits on ~30-50 bp unitigs.
    word_size = 7 if blast_task == "blastn-short" else WORD_SIZE
    cmd = [
        "nextflow", "run", str(PIPELINE_PATH),
        "--fasta",      str(FASTA_INPUT),
        "--card_db",    str(CARD_DB),
        "--outdir",     str(EXPLAINABILITY_DIR),
        "--antibiotic", TARGET_ANTIBIOTIC,
        "--threads",    str(THREADS),
        "--evalue",     str(EVALUE),
        "--word_size",  str(word_size),
        "--task",       str(blast_task),
    ]
    print(f"  BLAST task: {blast_task} | word_size: {word_size} "
          f"(auto from median query length; override=blast.task)")

    print(f"  Command: {' '.join(cmd)}\n")

    # Force English locale for the JVM that Nextflow spawns.
    # On systems with a Turkish locale, Java's String.toLowerCase() converts
    # 'I' → 'ı' (dotless-i), which breaks Nextflow's errorStrategy keyword
    # matching ('ignore' fails because the JVM sees 'ıgnore').
    # (os is already imported at module level.)
    nxf_env = os.environ.copy()
    nxf_env['NXF_OPTS'] = nxf_env.get('NXF_OPTS', '') + ' -Duser.language=en -Duser.country=US'

    result = subprocess.run(cmd, cwd=str(PROJECT_ROOT), env=nxf_env)

    if result.returncode != 0:
        print("\nERROR: Nextflow pipeline exited with a non-zero status.")
        print(f"  Return code: {result.returncode}")
        print("  Check Nextflow logs in the .nextflow.log file for details.")
        sys.exit(result.returncode)

    # -------------------------------------------------------------------------
    # COMPLETION: Confirm output files
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("BLAST ANNOTATION COMPLETE")
    print("=" * 80)
    print("\nOutput files:")

    for out_path in (CARD_OUT, NCBI_OUT):
        if out_path.exists():
            size_kb = out_path.stat().st_size / 1024
            print(f"  ✓ {out_path.name}  ({size_kb:.1f} KB)")
        else:
            print(f"  ⚠ Not found: {out_path.name}  (check Nextflow logs)")

    print(f"\nAll outputs in: {EXPLAINABILITY_DIR}")
    print("\nNext step:")
    print("  Run 09_biological_summary.py — it grades every hit into")
    print("  confirmed / candidate / weak tiers (thresholds in config.yaml →")
    print("  analysis.confidence_tiers), joins 07b stability, and computes the")
    print("  known-mechanism recovery rate, composite score and novel fraction.")
    print("=" * 80)


# ============================================================================
# SCRIPT ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    main()
