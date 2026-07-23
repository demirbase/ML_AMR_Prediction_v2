#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BLAST Annotation Orchestrator — Step 08

This script coordinates the biological validation of the top unitig/k-mer
features identified in Step 07 (07_explainability.py) by running blastn directly
via subprocess (no Nextflow) in two BLAST searches:

  1. CARD Local BLAST:
     Queries the Comprehensive Antibiotic Resistance Database (CARD).
     Directly tests whether the top k-mers overlap with documented
     resistance genes (e.g., gyrA, parC for fluoroquinolones).
     Requires a pre-built local blastn database.

  2. NCBI Remote BLAST:
     Queries the full NCBI nucleotide (nt) database over the internet.
     Captures novel or uncharacterised resistance determinants not yet
     in CARD. Uses -remote flag — no local database needed.

Both BLAST searches are embarrassingly parallel and independent. They are driven
from this script via subprocess — an earlier Nextflow orchestration (and the
section here explaining why it was used) was removed in the M9 review; the
pipeline is deliberately pure Python.

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

from lib.config import get_target  # early: env>config target before module globals

# Extract project-level identifiers
TARGET_ANTIBIOTIC = get_target(config=config)[1]
ORGANISM          = get_target(config=config)[0]
TOP_N             = config['analysis']['top_n_features']

# Organism-aware path resolution (SCALE_MLOPS_PLAN §4.2)
from lib.config import env_bool, resolve_path
from lib.registry import get_organism

# Resolve BLAST parameters from config
blast_cfg   = config.get('blast', {})
CARD_DB_DIR = PROJECT_ROOT / blast_cfg.get('card_db_dir', 'data/blast_db/card_nt')
CARD_DB     = CARD_DB_DIR / blast_cfg.get('card_db_name', 'card')
# Escape hatch for a deliberate NCBI-remote-only run. Env rather than a CLI flag
# because 08 has no argparse — the target itself comes from AMR_ORGANISM/
# AMR_ANTIBIOTIC — so this matches how every other knob reaches a SLURM job.
ALLOW_MISSING_CARD_DB = env_bool('AMR_ALLOW_MISSING_CARD_DB', False)
EVALUE      = blast_cfg.get('evalue',    10)
WORD_SIZE   = blast_cfg.get('word_size', 11)
THREADS     = blast_cfg.get('threads',   8)
# BLAST task is chosen from the ACTUAL query length (see choose_blast_task): the
# 'blastn-short' params are tuned for queries <50 bp (k-mers AND short unitigs),
# 'blastn' for longer. Picking by feature type was wrong — unitigs can be short
# (~30-50 bp), where 'blastn' finds nothing. blast.task overrides the auto choice.
BLAST_TASK_OVERRIDE = blast_cfg.get('task')

# NCBI remote pass parameters — DECOUPLED from the local CARD pass. The public
# NCBI BLAST server kills 'blastn-short' + word_size 7 over nt with SIGXCPU, so
# the remote pass uses 'blastn' + word_size 11 (sufficient for the high-identity
# genomic-context hits we want) and restricts the search to the study organism.
NCBI_TASK       = blast_cfg.get('ncbi_task',       'blastn')
NCBI_WORD_SIZE  = blast_cfg.get('ncbi_word_size',  11)
MAX_TARGET_SEQS = blast_cfg.get('max_target_seqs', 50)


def organism_entrez_query(organism_id):
    """Build an NCBI entrez_query that restricts the remote search to the study
    organism, derived from the registry (never hardcoded, so it auto-adjusts per
    organism). Prefer a TAXID filter ``txid<N>[Organism:exp]`` because it has NO
    spaces and passes as a single CLI token; a scientific-name filter like
    ``Escherichia coli[organism]`` contains a space that would be word-split into a
    broken ``-entrez_query`` argument. ':exp' explodes the taxon to include all
    descendant strains. Returns '' (no restriction) if the organism is unknown or
    has neither taxid nor name."""
    try:
        block = get_organism(organism_id)
    except Exception:
        return ""
    taxid = block.get('taxid')
    if taxid:
        return f"txid{taxid}[Organism:exp]"
    name = (block.get('display_name') or "").strip()
    return f"{name}[organism]" if name else ""


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
        1. Validate tool availability (blastn)
        2. Validate input FASTA from Step 07
        3. Validate CARD local database
        4. Run the two BLAST passes via subprocess (CARD local, NCBI remote)
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
    for tool in ("blastn",):
        path = shutil.which(tool)
        if path:
            print(f"  ✓ {tool:12s} found: {path}")
        else:
            print(f"  ✗ {tool:12s} NOT FOUND")
            missing_tools.append(tool)

    if missing_tools:
        print("\nERROR: The following required tools are not installed or not on PATH:")
        for tool in missing_tools:
            if tool == "blastn":
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
        # HARD FAIL by default. The local CARD pass is where tier, gene_symbol and
        # the ARO mapping come from — the KB's entire biology layer. Skipping it
        # does not produce a smaller KB, it produces a HOLLOW one: the pipeline
        # still exits 0, populate still writes rows, and nothing downstream
        # notices that every annotation is empty. That is the failure mode you
        # only discover after a full re-populate, which is exactly when it costs
        # the most. A warning is not enough for something this silent.
        msg = (
            f"CARD database not found at: {CARD_DB}\n"
            f"  The local CARD BLAST pass supplies tier / gene_symbol / ARO — without\n"
            f"  it the KB's biology layer would be silently empty, so this is fatal.\n"
            f"  To build the database:\n"
            f"    1. Download: https://card.mcmaster.ca/download\n"
            f"    2. makeblastdb -in <card.fna> -dbtype nucl -out {CARD_DB}\n"
            f"  If you really do want an NCBI-remote-only run, pass --allow-missing-card-db."
        )
        if not ALLOW_MISSING_CARD_DB:
            sys.exit(f"ERROR: {msg}")
        print(f"  ⚠ {msg}")
        print(f"    --allow-missing-card-db given — continuing with NCBI remote only.\n")
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
    # STEP 4: Run the two BLAST passes directly (CARD local + NCBI remote)
    # -------------------------------------------------------------------------
    # blastn is called via subprocess — NO Nextflow. The earlier .nf orchestration
    # was deleted but main() still shelled out to `nextflow run <the missing .nf>`,
    # so 08 could not run at all (nextflow is not even in the pinned amr.sif). The
    # pipeline only ever ran two blastn passes; doing them here directly is simpler
    # and drops the whole JVM/Nextflow dependency.
    print("\n[STEP 4/4] Running BLAST (CARD local + NCBI remote)...")
    print("=" * 80)

    # 09 reads outfmt-6 with exactly these columns (09_biological_summary.TSV_COLS);
    # 'qlen' before 'stitle' lets 09 compute query coverage for the tier cutoffs.
    OUTFMT = ("6 qseqid sseqid pident length mismatch gapopen qstart qend "
              "sstart send evalue bitscore qlen stitle")

    blast_task = choose_blast_task(FASTA_INPUT, BLAST_TASK_OVERRIDE)
    # blastn-short needs a small word_size (7) to seed short queries; the config
    # word_size (11+) truncated/missed full-length hits on ~30-50 bp unitigs.
    word_size = 7 if blast_task == "blastn-short" else WORD_SIZE
    # NCBI remote pass: organism-restricted entrez_query (registry taxid/name).
    entrez_query = organism_entrez_query(ORGANISM)
    print(f"  CARD task: {blast_task} | word_size: {word_size} "
          f"(auto from median query length; override=blast.task)")
    print(f"  NCBI task: {NCBI_TASK} | word_size: {NCBI_WORD_SIZE} | "
          f"max_target_seqs: {MAX_TARGET_SEQS} | "
          f"entrez_query: {entrez_query or '(none)'}")

    def run_blast(cmd, label, out_path):
        print(f"\n  → {label}: {' '.join(cmd)}")
        r = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
        if r.returncode != 0:
            sys.exit(f"ERROR: {label} blastn exited with status {r.returncode}.")
        # blastn writes an empty output file when there are no hits — a valid,
        # meaningful result (a novel k-mer with no DB match), NOT an error; 09
        # already tolerates empty TSVs. Only a truly missing file is fatal.
        if not out_path.exists():
            sys.exit(f"ERROR: {label} produced no output file at {out_path}.")
        print(f"    ✓ {out_path.name} ({out_path.stat().st_size / 1024:.1f} KB)")

    # CARD local pass (skipped only when the DB is absent AND --allow-missing-card-db;
    # otherwise STEP 3 already hard-failed).
    if card_present:
        run_blast(
            ["blastn", "-query", str(FASTA_INPUT), "-db", str(CARD_DB),
             "-out", str(CARD_OUT), "-outfmt", OUTFMT,
             "-task", blast_task, "-word_size", str(word_size),
             "-evalue", str(EVALUE), "-max_target_seqs", str(MAX_TARGET_SEQS),
             "-num_threads", str(THREADS)],
            "CARD local", CARD_OUT)
    else:
        CARD_OUT.write_text("", encoding="utf-8")  # empty file so 09 finds the path
        print("\n  → CARD local: SKIPPED (no DB, --allow-missing-card-db).")

    # NCBI remote pass. -remote runs server-side, so -num_threads is NOT allowed
    # (blastn errors if both are given). This is the slow pass (~10-20 min).
    ncbi_cmd = ["blastn", "-query", str(FASTA_INPUT), "-db", "nt", "-remote",
                "-out", str(NCBI_OUT), "-outfmt", OUTFMT,
                "-task", NCBI_TASK, "-word_size", str(NCBI_WORD_SIZE),
                "-evalue", str(EVALUE), "-max_target_seqs", str(MAX_TARGET_SEQS)]
    if entrez_query:
        ncbi_cmd += ["-entrez_query", entrez_query]
    print("\n  (NCBI remote BLAST over nt can take ~10-20 min — not a hang.)")
    run_blast(ncbi_cmd, "NCBI remote", NCBI_OUT)

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
            print(f"  ⚠ Not found: {out_path.name}  (unexpected — blastn should have written it)")

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
