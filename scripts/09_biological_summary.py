#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Biological Summary Report Generator — Step 09

Generates a Markdown report (05_final_biological_report.md) that maps the
top k-mer features (from Step 07) to their biological meaning via:

  1. CARD local BLAST results  → acquired resistance gene names
  2. NCBI remote BLAST results → core-genome / SNP context

For NCBI hits, instead of using the generic 'stitle' (which typically says
"complete genome"), this script queries NCBI Entrez efetch in real-time to
retrieve the specific gene/product name that overlaps the matched coordinates.

API behaviour is throttled (0.3 s between calls) and fully wrapped in
try/except so the script never crashes from network errors.
"""

# ============================================================================
# LIBRARY IMPORTS
# ============================================================================
import sys
import re
import time
import yaml
import pandas as pd
from pathlib import Path
from Bio import Entrez, SeqIO

# ============================================================================
# GLOBAL CONSTANTS
# ============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# NCBI Entrez identification (email / optional api_key) is configured at runtime
# from config.yaml in configure_entrez(); see main(). A fake/placeholder email
# (e.g. user@example.com) violates NCBI's Terms of Use and risks an IP ban, so
# we never hardcode one here.

# ----------------------------------------------------------------------------
# BLAST confidence tiers for short (k-mer) alignments
# ----------------------------------------------------------------------------
# A 21-mer producing an E-value of 1.5 is NOT "confirmed" homology. We grade
# every hit so the report cannot present weak matches as findings (was P-07).
CONFIRMED_MAX_EVALUE = 1e-3
CONFIRMED_MIN_IDENT  = 95.0
CANDIDATE_MAX_EVALUE = 1.0
CANDIDATE_MIN_IDENT  = 90.0
# Hits weaker than the candidate tier are dropped from the report entirely.
REPORT_MAX_EVALUE    = CANDIDATE_MAX_EVALUE


def classify_confidence(pident, evalue):
    """Grade a BLAST hit into confirmed / candidate / weak by E-value + identity."""
    try:
        pident = float(pident)
        evalue = float(evalue)
    except (TypeError, ValueError):
        return "weak"
    if evalue <= CONFIRMED_MAX_EVALUE and pident >= CONFIRMED_MIN_IDENT:
        return "confirmed"
    if evalue <= CANDIDATE_MAX_EVALUE and pident >= CANDIDATE_MIN_IDENT:
        return "candidate"
    return "weak"


def configure_entrez(config):
    """Configure NCBI Entrez identity from config; warn (don't crash) if unset."""
    ncbi_cfg = config.get('ncbi', {}) or {}
    email = (ncbi_cfg.get('entrez_email') or "").strip()
    api_key = (ncbi_cfg.get('api_key') or "").strip()

    if not email:
        print("WARNING: ncbi.entrez_email is not set in config.yaml.")
        print("         NCBI may rate-limit or ban requests without a valid e-mail.")
        print("         Set config['ncbi']['entrez_email'] before running Step 09.")
    else:
        Entrez.email = email

    # An API key raises the NCBI rate limit from 3 to 10 requests/sec and is the
    # recommended way to avoid throttling/bans for many sequential efetch calls.
    if api_key:
        Entrez.api_key = api_key
        print("  ✓ NCBI api_key configured (higher rate limit enabled).")


# ============================================================================
# CONFIGURATION LOADER
# ============================================================================
def load_config(config_path=None):
    if config_path is None:
        config_path = PROJECT_ROOT / "config" / "config.yaml"
    try:
        with open(config_path, "r") as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Error loading config: {e}")
        sys.exit(1)


# ============================================================================
# CARD HELPER
# ============================================================================
def extract_card_gene(sseqid):
    """Extract AMR gene symbol from CARD sseqid.

    Example:
        gb|NG_068181.1|+|100-925|ARO:3006096|OXA-909  →  OXA-909
    """
    sseqid = str(sseqid)
    if '|' in sseqid:
        return sseqid.split('|')[-1].strip()
    return sseqid


# ============================================================================
# NCBI STITLE CLEANER  (used as fallback when Entrez lookup fails)
# ============================================================================
def clean_ncbi_stitle(stitle):
    """Strip generic genome-level metadata from an NCBI stitle string."""
    stitle = str(stitle)
    patterns = [
        r",\s*complete genome",
        r"\s*complete genome",
        r",\s*complete sequence",
        r"\s*complete sequence",
        r"\s*genome assembly,\s*chromosome:\s*\w+",
        r"\s*genome assembly,\s*chromosome",
        r"\s*genome assembly.*",
        r"\s*chromosome,.*",
        r"\s*chromosome.*",
        r"\s*plasmid.*",
        r"\s+DNA,.*",
        r"\s+DNA.*",
        r",\s*partial cds",
        r"\s*partial cds",
        r"\s*gene for.*",
    ]
    cleaned = stitle
    for pattern in patterns:
        cleaned = re.sub(pattern, '', cleaned, flags=re.IGNORECASE)
    return cleaned.strip()


# ============================================================================
# ENTREZ COORDINATE-BASED GENE NAME LOOKUP
# ============================================================================
def _extract_accession(sseqid: str) -> str:
    """Return the bare accession number from a raw BLAST sseqid field.

    Handles formats such as:
        gi|123456|gb|NZ_CP012345.1|    →  NZ_CP012345.1
        ref|NZ_CP012345.1|             →  NZ_CP012345.1
        NZ_CP012345.1                  →  NZ_CP012345.1
    """
    sseqid = str(sseqid).strip()
    if '|' in sseqid:
        parts = [p for p in sseqid.split('|') if p.strip()]
        # The accession is the last non-empty token
        return parts[-1].strip()
    return sseqid


def fetch_gene_name_at_coords(sseqid: str, sstart: int, send: int,
                               stitle: str) -> str:
    """Query NCBI Entrez to find the gene/product overlapping [sstart, send].

    Parameters
    ----------
    sseqid  : raw BLAST subject sequence ID (accession, possibly pipe-delimited)
    sstart  : subject alignment start (1-based)
    send    : subject alignment end   (1-based)
    stitle  : original BLAST stitle column (used for fallback organism label)

    Returns
    -------
    A human-readable string in one of three formats:
        "GeneName (OrganismName)"          – successful Entrez lookup
        "Intergenic Region (OrganismName)" – no CDS/gene feature at coords
        "API Error (OrganismName)"         – network / parse error
    """
    organism_label = clean_ncbi_stitle(stitle)
    accession = _extract_accession(sseqid)

    # Reverse-strand hits have sstart > send; normalise for efetch
    seq_start = min(sstart, send)
    seq_stop  = max(sstart, send)

    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=accession,
            rettype="gb",
            retmode="text",
            seq_start=seq_start,
            seq_stop=seq_stop,
        )
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Walk features looking for an annotated gene or product qualifier
        for feature in record.features:
            if feature.type not in ("CDS", "gene", "rRNA", "tRNA", "ncRNA",
                                    "misc_RNA", "misc_feature"):
                continue

            qualifiers = feature.qualifiers

            # Prefer /gene over /product (shorter, canonical symbol)
            gene_name = (
                qualifiers.get("gene",    [None])[0]
                or qualifiers.get("product", [None])[0]
            )
            if gene_name:
                return f"{gene_name} ({organism_label})"

        # No annotated feature found at these coordinates
        return f"Intergenic Region ({organism_label})"

    except Exception:
        return f"API Error ({organism_label})"


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("Loading configuration...")
    config = load_config()
    antibiotic = config['project']['target_antibiotic']
    top_n = config.get('analysis', {}).get('top_n_features', 50)

    # Configure NCBI Entrez identity (email / api_key) from config — never a
    # hardcoded placeholder e-mail (NCBI ToS / ban risk).
    configure_entrez(config)

    # ------------------------------------------------------------------
    # Resolve paths (organism-aware — SCALE_MLOPS_PLAN §4.2)
    # ------------------------------------------------------------------
    from lib.config import resolve_path
    organism = config.get('project', {}).get('organism', 'ecoli')
    explain_dir = resolve_path('dir_05_explainability', organism=organism,
                               antibiotic=antibiotic, config=config)
    if not explain_dir.exists():
        print(f"Error: Directory {explain_dir} does not exist.")
        sys.exit(1)

    # Track top_n_features from config (07 writes 01_top_{top_n}_features).
    csv_file  = explain_dir / f"01_top_{top_n}_features_{antibiotic}.csv"
    card_file = explain_dir / f"03_card_blast_results_{antibiotic}.tsv"
    ncbi_file = explain_dir / f"04_ncbi_blast_results_{antibiotic}.tsv"
    out_file  = explain_dir / "05_final_biological_report.md"

    if not csv_file.exists():
        print(f"Error: Cannot find top features CSV at {csv_file}")
        sys.exit(1)

    print(f"Reading {csv_file}...")
    df_features = pd.read_csv(csv_file)

    # ------------------------------------------------------------------
    # TSV column schema shared by both BLAST result files
    # ------------------------------------------------------------------
    tsv_cols = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'stitle',
    ]

    # ------------------------------------------------------------------
    # Load & filter CARD results  (gene names extracted offline from sseqid)
    # ------------------------------------------------------------------
    if card_file.exists() and card_file.stat().st_size > 0:
        print(f"Reading {card_file}...")
        df_card = pd.read_csv(card_file, sep='\t', header=None, names=tsv_cols)
        df_card['pident'] = pd.to_numeric(df_card['pident'], errors='coerce')
        df_card['evalue'] = pd.to_numeric(df_card['evalue'], errors='coerce')
        # Drop clearly-insignificant hits (E > 1.0 for a 21-mer is noise), then
        # grade the survivors into confirmed/candidate tiers. The previous
        # E ≤ 50 filter let weak/meaningless alignments through as "findings".
        df_card = df_card[
            (df_card['pident'] >= CANDIDATE_MIN_IDENT) & (df_card['evalue'] <= REPORT_MAX_EVALUE)
        ].copy()
        df_card['Gene_Match'] = df_card['sseqid'].apply(extract_card_gene)
        df_card['Confidence'] = df_card.apply(
            lambda r: classify_confidence(r['pident'], r['evalue']), axis=1)
    else:
        print(f"Warning: {card_file} is missing or empty.")
        df_card = pd.DataFrame(columns=tsv_cols + ['Gene_Match'])

    # ------------------------------------------------------------------
    # Load & filter NCBI results  (gene names resolved at report-write time)
    # Gene_Match column is intentionally left blank here; it will be
    # populated row-by-row inside the report loop below.
    # ------------------------------------------------------------------
    if ncbi_file.exists() and ncbi_file.stat().st_size > 0:
        print(f"Reading {ncbi_file}...")
        df_ncbi = pd.read_csv(ncbi_file, sep='\t', header=None, names=tsv_cols)
        df_ncbi['pident'] = pd.to_numeric(df_ncbi['pident'], errors='coerce')
        df_ncbi['evalue'] = pd.to_numeric(df_ncbi['evalue'], errors='coerce')
        df_ncbi['sstart'] = pd.to_numeric(df_ncbi['sstart'], errors='coerce').fillna(0).astype(int)
        df_ncbi['send']   = pd.to_numeric(df_ncbi['send'],   errors='coerce').fillna(0).astype(int)
        df_ncbi = df_ncbi[
            (df_ncbi['pident'] >= CANDIDATE_MIN_IDENT) & (df_ncbi['evalue'] <= REPORT_MAX_EVALUE)
        ].copy()
        df_ncbi['Confidence'] = df_ncbi.apply(
            lambda r: classify_confidence(r['pident'], r['evalue']), axis=1)
    else:
        print(f"Warning: {ncbi_file} is missing or empty.")
        df_ncbi = pd.DataFrame(columns=tsv_cols + ['Confidence'])

    # ------------------------------------------------------------------
    # Generate Markdown report
    # ------------------------------------------------------------------
    print(f"Generating markdown report at {out_file}...")

    with open(out_file, "w") as f:
        f.write("# Final Biological Report\n")
        f.write(f"**Target Antibiotic:** {antibiotic.capitalize()}\n\n")
        f.write("**Confidence tiers** (short-alignment BLAST grading):\n")
        f.write(f"- `confirmed` — E ≤ {CONFIRMED_MAX_EVALUE:g} and identity ≥ {CONFIRMED_MIN_IDENT:g}%\n")
        f.write(f"- `candidate` — E ≤ {CANDIDATE_MAX_EVALUE:g} and identity ≥ {CANDIDATE_MIN_IDENT:g}%\n")
        f.write(f"- Weaker hits (E > {REPORT_MAX_EVALUE:g}) are excluded from this report.\n\n")
        f.write("---\n\n")

        for _, row in df_features.iterrows():
            rank     = int(row['Rank'])
            score    = float(row['Gain_Score'])
            feat_id  = str(row['Feature_ID'])
            sequence = str(row['Kmer_Sequence'])

            # Reconstruct the query ID to match BLAST qseqid column
            q_id = f"Rank_{rank}|Score_{score:.4f}|Feature_{feat_id}"

            f.write(f"### Rank {rank}: {sequence} (Gain: {score:.4f})\n")

            # --------------------------------------------------------------
            # CARD hits — no Entrez call needed, offline gene symbol lookup
            # --------------------------------------------------------------
            f.write("**CARD Hits (Acquired Resistance / Plasmids):**\n")
            card_hits = df_card[df_card['qseqid'] == q_id].head(10)
            if not card_hits.empty:
                for _, hit in card_hits.iterrows():
                    f.write(
                        f"- {hit['Gene_Match']}, "
                        f"Identity: {hit['pident']}%, "
                        f"E-value: {hit['evalue']} "
                        f"[{hit.get('Confidence', 'n/a')}]\n"
                    )
            else:
                f.write("*No high-confidence hits*\n")

            # --------------------------------------------------------------
            # NCBI hits — real-time Entrez coordinate lookup (top 10 only)
            # --------------------------------------------------------------
            f.write("**NCBI Hits (Core Genome / SNPs):**\n")
            ncbi_hits = df_ncbi[df_ncbi['qseqid'] == q_id].head(10)
            if not ncbi_hits.empty:
                for _, hit in ncbi_hits.iterrows():
                    gene_label = fetch_gene_name_at_coords(
                        sseqid=hit['sseqid'],
                        sstart=int(hit['sstart']),
                        send=int(hit['send']),
                        stitle=hit['stitle'],
                    )
                    f.write(
                        f"- {gene_label}, "
                        f"Identity: {hit['pident']}%, "
                        f"E-value: {hit['evalue']} "
                        f"[{hit.get('Confidence', 'n/a')}]\n"
                    )
                    # Be polite to the NCBI API. With an api_key the allowance is
                    # 10 req/s (0.1s); without one it is 3 req/s, so 0.34s is the
                    # safe floor. Use the looser delay when no key is configured.
                    time.sleep(0.1 if getattr(Entrez, 'api_key', None) else 0.34)
            else:
                f.write("*No high-confidence hits*\n")

            f.write("\n")

    print("Generation complete!")


# ============================================================================
# ENTRY POINT
# ============================================================================
if __name__ == "__main__":
    main()
