#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 11 — k-mer-centric CARD variant-model SNP allele check.

A homolog-model CARD hit (step 08/09) only proves a k-mer lies in a known ARG
*region*; for SNP-mediated resistance (e.g. gyrA/parC fluoroquinolone mutations,
rpoB rifampicin) the gene is present in every isolate and only a specific point
mutation confers resistance. This step asks the precise question:

    Does the candidate k-mer span a known CARD resistance SNP position, and if
    so, does it carry the RESISTANT allele or the wildtype base?

Method:
  1. BLAST the candidate k-mers (07 FASTA) against CARD's *protein variant
     model* nucleotide sequences (blastn-short), capturing the aligned query
     and subject strings + subject coordinates.
  2. Parse card.json for each variant model's resistance SNPs (protein
     positions, e.g. "S83L"); map protein position p -> CDS nt codon (3p-2..3p).
  3. For every hit covering a SNP codon, read the k-mer's bases aligned to that
     codon (strand-aware), translate, and classify:
         resistant_allele | wildtype | other_variant | partial/ambiguous.

This turns "k-mer is in gyrA" into "k-mer carries gyrA S83L (resistant)" — the
difference between a lineage/wildtype signal and a true resistance determinant.

Requires the full CARD data download (not shipped):
    mkdir -p data/external/card && cd data/external/card
    curl -L https://card.mcmaster.ca/latest/data -o card-data.tar.bz2
    tar xjf card-data.tar.bz2
If the variant FASTA / card.json are absent the step exits cleanly with a note.

Output:
    results/{org}/{ab}/05_explainability/11_variant_snp_check_{ab}.csv
"""

import io
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd

from lib.config import load_config, resolve_path, resolve_tool, get_target

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Ensure the conda environment's bin is on PATH so blastn/makeblastdb resolve
# even when launched via the full interpreter path (matches step 08).
os.environ['PATH'] = str(Path(sys.executable).parent) + os.pathsep + os.environ.get('PATH', '')

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
_CODON_TABLE = {  # standard genetic code (DNA codons -> 1-letter AA; '*' = stop)
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S',
    'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A',
    'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAT': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
    'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

_SNP_RE = re.compile(r'^([A-Z])(\d+)([A-Z])$')

# BLAST tabular columns we request (aligned strings let us read the SNP codon
# directly, strand-aware, without manual coordinate gymnastics).
OUTFMT_COLS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
               "sstart", "send", "sstrand", "qseq", "sseq", "evalue", "bitscore"]


# --------------------------------------------------------------------------
# Pure helpers (unit-tested)
# --------------------------------------------------------------------------
def parse_snp_token(token):
    """'S83L' -> ('S', 83, 'L'); returns None for non-substitution tokens."""
    m = _SNP_RE.match(str(token).strip())
    if not m:
        return None
    return m.group(1), int(m.group(2)), m.group(3)


def protein_pos_to_codon_nt(protein_pos):
    """1-based protein position -> (start, end) 1-based CDS nt codon range."""
    start = (int(protein_pos) - 1) * 3 + 1
    return start, start + 2


def translate_codon(codon):
    """Translate a 3-nt DNA codon to a 1-letter amino acid ('X' if unknown)."""
    return _CODON_TABLE.get(str(codon).upper(), 'X')


def query_codon_from_alignment(sstart, send, sstrand, qseq_aln, sseq_aln,
                               codon_positions):
    """
    Read the QUERY (k-mer) bases aligned to the given subject CDS positions,
    returned in CDS (plus) sense. None if the codon is not fully covered or a
    gap falls in it.

    BLAST reports qseq/sseq in the alignment's orientation. For a plus-strand
    subject hit the query base IS the CDS-sense base; for a minus-strand hit the
    CDS-sense base is the complement of the aligned query char (and subject
    coordinates run downward). We walk the alignment columns tracking the
    subject coordinate and collect the query bases at the requested positions.
    """
    wanted = set(int(p) for p in codon_positions)
    plus = (str(sstrand).lower() != "minus") and not (int(sstart) > int(send))
    subj = int(sstart)
    step = 1 if plus else -1
    collected = {}
    for qc, sc in zip(str(qseq_aln), str(sseq_aln)):
        if sc == '-':                      # insertion in query vs subject: no subj advance
            continue
        if subj in wanted:
            if qc == '-':                  # deletion in query at a SNP position -> ambiguous
                return None
            base = qc if plus else qc.translate(_COMPLEMENT)
            collected[subj] = base.upper()
        subj += step
    if len(collected) != len(wanted):
        return None                        # codon not fully covered by the alignment
    return ''.join(collected[p] for p in sorted(wanted))


def classify_allele(query_codon, wt_aa, mut_aa):
    """Classify the k-mer's codon against a CARD SNP (wt -> mut)."""
    if not query_codon or len(query_codon) != 3:
        return "ambiguous"
    aa = translate_codon(query_codon)
    if aa == mut_aa:
        return "resistant_allele"
    if aa == wt_aa:
        return "wildtype"
    return "other_variant"


# --------------------------------------------------------------------------
# CARD parsing
# --------------------------------------------------------------------------
def parse_card_variant_snps(card_json_path):
    """Return {ARO_accession: [(wt, pos, mut), ...]} for protein variant models."""
    with open(card_json_path, encoding="utf-8") as f:
        data = json.load(f)
    out = {}
    for v in data.values():
        if not isinstance(v, dict):
            continue
        if v.get("model_type") != "protein variant model":
            continue
        aro = str(v.get("ARO_accession") or "").replace("ARO:", "")
        if not aro:
            continue
        snp_param = (v.get("model_param", {}) or {}).get("snp", {}) or {}
        values = snp_param.get("param_value", {}) or {}
        snps = []
        for entry in values.values():
            # param_value maps an id -> the SNP token directly ("S83L"), though
            # some CARD releases nest it as {"param_value": "S83L"}.
            tok = entry.get("param_value") if isinstance(entry, dict) else entry
            parsed = parse_snp_token(tok) if tok else None
            if parsed:
                snps.append(parsed)
        if snps:
            out[aro] = snps
    return out


def aro_from_sseqid(sseqid):
    """Extract bare ARO accession from a CARD sseqid (…|ARO:3003297|gyrA…)."""
    for tok in str(sseqid).split('|'):
        if tok.startswith("ARO:"):
            return tok.split("ARO:")[1].strip()
    return None


def gene_from_sseqid(sseqid):
    parts = [p for p in str(sseqid).split('|') if p.strip()]
    return parts[-1].split(' ')[0] if parts else str(sseqid)


# --------------------------------------------------------------------------
def main():
    config = load_config()
    antibiotic = get_target(config=config)[1]
    organism = get_target(config=config)[0]
    top_n = config.get('analysis', {}).get('top_n_features', 50)
    blast_cfg = config.get('blast', {})

    explain_dir = resolve_path('dir_05_explainability', organism=organism,
                               antibiotic=antibiotic, config=config)
    fasta = explain_dir / f"02_top_{top_n}_features_{antibiotic}.fasta"
    variant_fasta = PROJECT_ROOT / blast_cfg.get(
        'card_variant_fasta', 'data/external/card/nucleotide_fasta_protein_variant_model.fasta')
    card_json = PROJECT_ROOT / blast_cfg.get('card_json', 'data/external/card/card.json')
    out_path = explain_dir / f"11_variant_snp_check_{antibiotic}.csv"

    print("=" * 80)
    print(f"VARIANT-MODEL SNP CHECK: {antibiotic.upper()} ({organism})")
    print("=" * 80)

    # Graceful skip if CARD variant data is not downloaded
    missing = [str(p) for p in (variant_fasta, card_json) if not p.exists()]
    if missing:
        print("  CARD variant data not found — skipping (this step is optional):")
        for m in missing:
            print(f"    missing: {m}")
        print("  Download: curl -L https://card.mcmaster.ca/latest/data -o card-data.tar.bz2 "
              "&& tar xjf it under data/external/card/")
        return
    if not fasta.exists() or fasta.stat().st_size == 0:
        print(f"  No candidate FASTA at {fasta}; run 07 first. Nothing to check.")
        return

    blastn = resolve_tool('blast_bin', 'blastn', config=config)
    makeblastdb = resolve_tool('makeblastdb_bin', 'makeblastdb', config=config)
    if not blastn or not makeblastdb:
        print("  blastn/makeblastdb not on PATH (conda install -c bioconda blast). Skipping.")
        return

    # 1) Build the variant-model BLAST DB once (cached)
    db_dir = variant_fasta.parent / "variant_db"
    db_dir.mkdir(parents=True, exist_ok=True)
    db = db_dir / "variant"
    if not (db.parent / (db.name + ".nhr")).exists() and not any(db.parent.glob(db.name + ".*.nhr")):
        print("  Building CARD variant-model BLAST DB...")
        subprocess.run([makeblastdb, "-in", str(variant_fasta), "-dbtype", "nucl",
                        "-out", str(db)], check=True, capture_output=True, text=True)

    # 2) BLAST candidate k-mers vs variant models (short-query mode)
    print("  BLAST candidate k-mers vs CARD variant models...")
    r = subprocess.run(
        [blastn, "-query", str(fasta), "-db", str(db), "-task", "blastn-short",
         "-dust", "no", "-evalue", str(blast_cfg.get('evalue', 10)),
         "-word_size", str(blast_cfg.get('word_size', 11)),
         "-outfmt", "6 " + " ".join(OUTFMT_COLS)],
        capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  blastn failed: {(r.stderr or '')[:300]}")
        return
    if not r.stdout.strip():
        print("  No variant-model BLAST hits for the candidate k-mers.")
        pd.DataFrame(columns=["kmer_qseqid", "variant_gene", "aro", "snp",
                              "observed_codon", "observed_aa", "allele_class"]).to_csv(out_path, index=False)
        print(f"  Saved (empty): {out_path.name}")
        return
    hits = pd.read_csv(io.StringIO(r.stdout), sep="\t", names=OUTFMT_COLS)

    # 3) CARD SNP annotations
    print("  Parsing CARD variant SNP annotations...")
    aro_snps = parse_card_variant_snps(card_json)

    # 4) For each hit, test every SNP codon it covers
    rows = []
    for _, h in hits.iterrows():
        aro = aro_from_sseqid(h['sseqid'])
        snps = aro_snps.get(aro, [])
        if not snps:
            continue
        s0, s1 = int(h['sstart']), int(h['send'])
        lo, hi = min(s0, s1), max(s0, s1)
        for wt, pos, mut in snps:
            c_start, c_end = protein_pos_to_codon_nt(pos)
            if c_start < lo or c_end > hi:
                continue                                  # codon not within the alignment span
            qcodon = query_codon_from_alignment(
                h['sstart'], h['send'], h['sstrand'], h['qseq'], h['sseq'],
                (c_start, c_start + 1, c_start + 2))
            allele = classify_allele(qcodon, wt, mut)
            rows.append({
                "kmer_qseqid": h['qseqid'],
                "variant_gene": gene_from_sseqid(h['sseqid']),
                "aro": aro,
                "snp": f"{wt}{pos}{mut}",
                "pident": h['pident'], "evalue": h['evalue'],
                "observed_codon": qcodon or "",
                "observed_aa": translate_codon(qcodon) if qcodon else "",
                "allele_class": allele,
            })

    out = pd.DataFrame(rows, columns=["kmer_qseqid", "variant_gene", "aro", "snp",
                                      "pident", "evalue", "observed_codon",
                                      "observed_aa", "allele_class"])
    out.to_csv(out_path, index=False)

    n_res = int((out['allele_class'] == 'resistant_allele').sum()) if len(out) else 0
    n_wt = int((out['allele_class'] == 'wildtype').sum()) if len(out) else 0
    print(f"  Variant-model SNP codons spanned: {len(out)}  "
          f"(resistant_allele={n_res}, wildtype={n_wt})")
    if n_res:
        for _, r2 in out[out['allele_class'] == 'resistant_allele'].iterrows():
            print(f"    ✓ {r2['variant_gene']} {r2['snp']} — k-mer carries RESISTANT allele "
                  f"({r2['kmer_qseqid']})")
    print(f"  Saved: {out_path.name}")
    print("=" * 80)


if __name__ == "__main__":
    main()
