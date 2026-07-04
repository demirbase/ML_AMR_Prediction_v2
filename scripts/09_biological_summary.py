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
import os
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
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
from lib.config import load_config, get_target  # noqa: E402  (canonical loader; audit Issue 23)

# NCBI Entrez identification (email / optional api_key) is configured at runtime
# from config.yaml in configure_entrez(); see main(). A fake/placeholder email
# (e.g. user@example.com) violates NCBI's Terms of Use and risks an IP ban, so
# we never hardcode one here.

# ----------------------------------------------------------------------------
# BLAST confidence tiers for short (k-mer) alignments
# ----------------------------------------------------------------------------
# A 21-mer producing an E-value of 1.5 is NOT "confirmed" homology, so every hit
# is graded (was P-07). Thresholds are read from config.yaml (analysis.
# confidence_tiers) so they are citeable in Methods and easy to defend; a hard-
# coded fallback keeps the script runnable if the keys are absent.
DEFAULT_TIERS = {
    "confirmed": {"min_identity": 95.0, "min_coverage": 0.95, "max_evalue": 1.0},
    "candidate": {"min_identity": 90.0, "min_coverage": 0.80, "max_evalue": 10.0},
    "weak":      {"min_identity": 80.0, "min_coverage": 0.60, "max_evalue": 50.0},
}
DEFAULT_REPORT_MAX_EVALUE = 50.0
DEFAULT_KMER_LENGTH = 21


def load_tiers(config):
    """Return (tiers, report_max_evalue, weak_min_identity, weak_min_cov, k_length, stability_threshold)."""
    analysis = config.get("analysis", {}) or {}
    tiers = analysis.get("confidence_tiers", DEFAULT_TIERS) or DEFAULT_TIERS
    report_max = float(analysis.get("report_max_evalue", DEFAULT_REPORT_MAX_EVALUE))
    weak = tiers.get("weak", DEFAULT_TIERS["weak"])
    weak_min_ident = float(weak["min_identity"])
    weak_min_cov = float(weak.get("min_coverage", DEFAULT_TIERS["weak"]["min_coverage"]))
    k_length = int(config.get("preprocessing", {}).get("k_length", DEFAULT_KMER_LENGTH))
    stability_threshold = float(analysis.get("stability_threshold", 0.6))
    return tiers, report_max, weak_min_ident, weak_min_cov, k_length, stability_threshold


def classify_confidence(pident, evalue, length, query_length, tiers):
    """
    Grade a BLAST hit into confirmed / candidate / weak (best-first) or 'none'.

    The primary, database-size-independent criteria are IDENTITY and COVERAGE
    (alignment length / QUERY length). ``query_length`` is k for a k-mer query
    and the unitig length for a unitig query (08 emits `qlen`); using a fixed k
    for unitigs would make coverage meaningless (always >1). E-value is only a
    loose secondary gate (scales with DB size; not comparable CARD vs NCBI nt).
    """
    try:
        pident = float(pident)
        evalue = float(evalue)
        coverage = float(length) / float(query_length) if query_length else 0.0
    except (TypeError, ValueError, ZeroDivisionError):
        return "none"
    for tier in ("confirmed", "candidate", "weak"):
        t = tiers.get(tier)
        if (t and pident >= float(t["min_identity"])
                and coverage >= float(t.get("min_coverage", 0.0))
                and evalue <= float(t["max_evalue"])):
            return tier
    return "none"


# BLAST outfmt-6 columns. 08 now emits `qlen` (query length) before `stitle`;
# older outputs lack it (handled by read_blast_tsv -> qlen = NaN -> k_length).
TSV_COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'stitle']


def read_blast_tsv(path):
    """Read a BLAST outfmt-6 TSV, tolerating presence/absence of the `qlen`
    column (forward/backward compatible). Always returns a frame with all
    TSV_COLS; a missing `qlen` is filled NaN so callers fall back to k_length.
    """
    df = pd.read_csv(path, sep='\t', header=None)
    ncol = df.shape[1]
    if ncol == len(TSV_COLS):
        df.columns = TSV_COLS
    elif ncol == len(TSV_COLS) - 1:                     # legacy: no qlen
        df.columns = [c for c in TSV_COLS if c != 'qlen']
        df['qlen'] = float('nan')
    else:                                               # best effort
        df.columns = (TSV_COLS + [f'extra{i}' for i in range(ncol)])[:ncol]
        if 'qlen' not in df.columns:
            df['qlen'] = float('nan')
    return df


def configure_entrez(config):
    """Configure NCBI Entrez identity from config; warn (don't crash) if unset."""
    ncbi_cfg = config.get('ncbi', {}) or {}
    # AMR_ENTREZ_EMAIL / AMR_ENTREZ_API_KEY env overrides let HPC set the Entrez
    # identity without editing config.yaml (which carries manual HPC tuning that
    # must not be overwritten). Env wins over config; both fall back gracefully.
    email = (os.environ.get('AMR_ENTREZ_EMAIL') or ncbi_cfg.get('entrez_email') or "").strip()
    api_key = (os.environ.get('AMR_ENTREZ_API_KEY') or ncbi_cfg.get('api_key') or "").strip()

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


def aro_from_sseqid(sseqid):
    """Extract the ARO accession digits from a CARD sseqid (ROADMAP §0.2 M16).

    'gb|NG_068181.1|+|100-925|ARO:3006096|OXA-909' -> '3006096'; '' if absent.
    """
    mobj = re.search(r'ARO:(\d+)', str(sseqid))
    return mobj.group(1) if mobj else ''


def load_aro_index(path):
    """Load CARD's aro_index.tsv -> {aro_accession_digits: {gene_family,
    drug_class, resistance_mechanism, short_name}} for the ARO ontology mapping
    (ROADMAP §0.2 M16). Returns {} if the file is absent (full CARD download only),
    so the KB simply gets empty ARO fields rather than failing. Column names are
    resolved case-insensitively because they drift slightly between CARD versions.
    """
    path = Path(path)
    if not path.exists():
        return {}
    df = pd.read_csv(path, sep='\t', dtype=str).fillna('')

    def col(*names):
        for n in names:
            for c in df.columns:
                if c.strip().lower() == n.lower():
                    return c
        return None

    aro_c = col('ARO Accession')
    if aro_c is None:
        return {}
    fam_c, drug_c = col('AMR Gene Family'), col('Drug Class')
    mech_c, name_c = col('Resistance Mechanism'), col('CARD Short Name', 'ARO Name')
    out = {}
    for _, r in df.iterrows():
        acc = str(r[aro_c]).replace('ARO:', '').strip()
        if acc:
            out[acc] = {
                'gene_family': r[fam_c] if fam_c else '',
                'drug_class': r[drug_c] if drug_c else '',
                'resistance_mechanism': r[mech_c] if mech_c else '',
                'short_name': r[name_c] if name_c else '',
            }
    return out


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
# KB-CANDIDATE TABLE + QUANTITATIVE VALIDATION (M7 / composite / H4)
# ============================================================================
import json
import math


def best_card_hit(df_card, q_id):
    """Lowest-E-value CARD hit for a feature, or None. Returns a dict."""
    if df_card.empty:
        return None
    hits = df_card[df_card['qseqid'] == q_id]
    if hits.empty:
        return None
    row = hits.loc[hits['evalue'].idxmin()]
    return {
        'gene': row.get('Gene_Match', ''),
        'identity': float(row['pident']),
        'evalue': float(row['evalue']),
        'tier': row.get('Confidence', 'none'),
        'sseqid': row.get('sseqid', ''),
    }


def composite_score(selection_frequency, identity_pct, evalue):
    """ROADMAP §1.4: stability × log10(1/E) × (identity/100). NaN if inputs missing."""
    try:
        sf = float(selection_frequency)
        idp = float(identity_pct)
        ev = float(evalue)
    except (TypeError, ValueError):
        return float('nan')
    if not (sf == sf) or ev <= 0:          # sf NaN or non-positive E
        return float('nan')
    # E-value can exceed 1 for short k-mers -> log10(1/E) goes negative; clamp the
    # contribution at 0 so a weak hit cannot produce a negative composite score.
    return sf * max(0.0, math.log10(1.0 / ev)) * (idp / 100.0)


def build_kb_candidates(df_features, df_card, stability_threshold, aro_index=None):
    """
    Join the candidate k-mers (07: gain top-N ∪ stable) with their best CARD hit
    and compute the composite score. Returns (DataFrame, metrics dict).

    When ``aro_index`` is provided (CARD aro_index.tsv), each CARD hit is mapped to
    its ARO ontology (ROADMAP §0.2 M16): accession + gene family + drug class +
    resistance mechanism — the 5-field schema for the KB.
    """
    aro_index = aro_index or {}
    rows = []
    for _, feat in df_features.iterrows():
        rank = int(feat['Rank'])
        score = float(feat['Gain_Score'])
        feat_id = str(feat['Feature_ID'])
        q_id = f"Rank_{rank}|Score_{score:.4f}|Feature_{feat_id}"
        sel_freq = feat.get('selection_frequency', float('nan'))
        try:
            sel_freq = float(sel_freq)
        except (TypeError, ValueError):
            sel_freq = float('nan')
        is_stable = bool(feat.get('stable', False))

        hit = best_card_hit(df_card, q_id)
        aro_acc = aro_from_sseqid(hit['sseqid']) if hit else ''
        aro = aro_index.get(aro_acc, {})
        rows.append({
            'rank': rank,
            'kmer': feat.get('Kmer_Sequence', ''),
            'feature_id': feat_id,
            'gain_score': score,
            'in_gain_topN': bool(feat.get('in_gain_topN', True)),
            'selection_frequency': sel_freq,
            'stable': is_stable,
            'card_gene': hit['gene'] if hit else '',
            'card_identity': hit['identity'] if hit else float('nan'),
            'card_evalue': hit['evalue'] if hit else float('nan'),
            'confidence_tier': hit['tier'] if hit else 'none',
            'has_card_hit': hit is not None,
            'composite_score': composite_score(sel_freq, hit['identity'], hit['evalue']) if hit else float('nan'),
            # ARO ontology (M16): empty when no hit / no aro_index.tsv present.
            'aro_accession': f"ARO:{aro_acc}" if aro_acc else '',
            'aro_gene_family': aro.get('gene_family', ''),
            'aro_drug_class': aro.get('drug_class', ''),
            'aro_resistance_mechanism': aro.get('resistance_mechanism', ''),
        })
    cols = ['rank', 'kmer', 'feature_id', 'gain_score', 'in_gain_topN',
            'selection_frequency', 'stable', 'card_gene', 'card_identity',
            'card_evalue', 'confidence_tier', 'has_card_hit', 'composite_score',
            'aro_accession', 'aro_gene_family', 'aro_drug_class',
            'aro_resistance_mechanism']
    kb = pd.DataFrame(rows, columns=cols)
    if not kb.empty:
        kb = kb.sort_values(['composite_score', 'gain_score'],
                            ascending=False, na_position='last')

    # ---- quantitative metrics (M7 recovery rate, H2, H4 novel fraction) -----
    n_features = len(kb)
    stable = kb[kb['stable']] if n_features else kb
    n_stable = len(stable)
    n_stable_confirmed = int((stable['confidence_tier'] == 'confirmed').sum()) if n_stable else 0
    n_stable_with_hit = int(stable['has_card_hit'].sum()) if n_stable else 0
    n_stable_novel = n_stable - n_stable_with_hit
    metrics = {
        'n_candidate_features': n_features,
        'n_in_gain_topN': int(kb['in_gain_topN'].sum()) if n_features else 0,
        'n_stable': n_stable,
        'stability_threshold': stability_threshold,
        'tier_counts_all': (kb['confidence_tier'].value_counts().to_dict() if n_features else {}),
        # M7 / H2: of the reproducible (stable) k-mers, fraction that map to a
        # known ARG at confirmed confidence. H2 accepts >= 0.40.
        'known_mechanism_recovery_rate': (n_stable_confirmed / n_stable) if n_stable else None,
        'H2_pass': ((n_stable_confirmed / n_stable) >= 0.40) if n_stable else None,
        # H4: fraction of stable k-mers with NO CARD hit (novel candidates).
        'novel_candidate_fraction': (n_stable_novel / n_stable) if n_stable else None,
        'n_stable_confirmed': n_stable_confirmed,
        'n_stable_novel': n_stable_novel,
    }
    return kb, metrics


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("Loading configuration...")
    config = load_config()
    antibiotic = get_target(config=config)[1]
    top_n = config.get('analysis', {}).get('top_n_features', 50)
    tiers, report_max_evalue, weak_min_ident, weak_min_cov, k_length, stability_threshold = load_tiers(config)

    # Configure NCBI Entrez identity (email / api_key) from config — never a
    # hardcoded placeholder e-mail (NCBI ToS / ban risk).
    configure_entrez(config)

    # ------------------------------------------------------------------
    # Resolve paths (organism-aware — SCALE_MLOPS_PLAN §4.2)
    # ------------------------------------------------------------------
    from lib.config import resolve_path
    organism = get_target(config=config)[0]
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
    tsv_cols = TSV_COLS

    # ------------------------------------------------------------------
    # Load & filter CARD results  (gene names extracted offline from sseqid)
    # ------------------------------------------------------------------
    if card_file.exists() and card_file.stat().st_size > 0:
        print(f"Reading {card_file}...")
        df_card = read_blast_tsv(card_file)
        df_card['pident'] = pd.to_numeric(df_card['pident'], errors='coerce')
        df_card['evalue'] = pd.to_numeric(df_card['evalue'], errors='coerce')
        df_card['length'] = pd.to_numeric(df_card['length'], errors='coerce')
        df_card['qlen']   = pd.to_numeric(df_card['qlen'], errors='coerce')
        # Effective query length: the real qlen for unitigs, else k_length (k-mers).
        df_card['qlen_eff'] = df_card['qlen'].where(df_card['qlen'] > 0, k_length)
        # Keep everything down to the weak tier (identity + coverage floors,
        # E ≤ report_max_evalue) and grade each hit. Weak hits are kept and
        # FLAGGED rather than dropped, for transparency (ROADMAP Risk-4 / §1.4).
        df_card = df_card[
            (df_card['pident'] >= weak_min_ident)
            & (df_card['length'] >= weak_min_cov * df_card['qlen_eff'])
            & (df_card['evalue'] <= report_max_evalue)
        ].copy()
        df_card['Gene_Match'] = df_card['sseqid'].apply(extract_card_gene)
        df_card['Confidence'] = df_card.apply(
            lambda r: classify_confidence(r['pident'], r['evalue'], r['length'], r['qlen_eff'], tiers), axis=1)
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
        df_ncbi = read_blast_tsv(ncbi_file)
        df_ncbi['pident'] = pd.to_numeric(df_ncbi['pident'], errors='coerce')
        df_ncbi['evalue'] = pd.to_numeric(df_ncbi['evalue'], errors='coerce')
        df_ncbi['length'] = pd.to_numeric(df_ncbi['length'], errors='coerce')
        df_ncbi['qlen']   = pd.to_numeric(df_ncbi['qlen'], errors='coerce')
        df_ncbi['qlen_eff'] = df_ncbi['qlen'].where(df_ncbi['qlen'] > 0, k_length)
        df_ncbi['sstart'] = pd.to_numeric(df_ncbi['sstart'], errors='coerce').fillna(0).astype(int)
        df_ncbi['send']   = pd.to_numeric(df_ncbi['send'],   errors='coerce').fillna(0).astype(int)
        df_ncbi = df_ncbi[
            (df_ncbi['pident'] >= weak_min_ident)
            & (df_ncbi['length'] >= weak_min_cov * df_ncbi['qlen_eff'])
            & (df_ncbi['evalue'] <= report_max_evalue)
        ].copy()
        df_ncbi['Confidence'] = df_ncbi.apply(
            lambda r: classify_confidence(r['pident'], r['evalue'], r['length'], r['qlen_eff'], tiers), axis=1)
    else:
        print(f"Warning: {ncbi_file} is missing or empty.")
        df_ncbi = pd.DataFrame(columns=tsv_cols + ['Confidence'])

    # ------------------------------------------------------------------
    # KB-candidate table + quantitative validation (M7 / composite / H4)
    # ------------------------------------------------------------------
    # ARO ontology mapping (M16): aro_index.tsv ships in the full CARD download
    # next to card.json; absent -> empty mapping (KB just gets blank ARO fields).
    aro_index_path = config.get('blast', {}).get('aro_index')
    if aro_index_path:
        aro_index_path = PROJECT_ROOT / aro_index_path
    else:
        card_json = config.get('blast', {}).get('card_json', 'data/external/card/card.json')
        aro_index_path = PROJECT_ROOT / Path(card_json).parent / 'aro_index.tsv'
    aro_index = load_aro_index(aro_index_path)
    if aro_index:
        print(f"  ✓ ARO ontology: {len(aro_index)} accessions ({Path(aro_index_path).name})")
    else:
        print(f"  ⚠ ARO index not found ({aro_index_path}); KB ARO fields left blank "
              f"(needs the full CARD download).")

    kb, metrics = build_kb_candidates(df_features, df_card, stability_threshold,
                                      aro_index=aro_index)
    kb_path = explain_dir / f"07_kb_candidates_{antibiotic}.csv"
    kb.to_csv(kb_path, index=False)
    metrics_path = explain_dir / f"08_validation_metrics_{antibiotic}.json"
    metrics_path.write_text(json.dumps(metrics, indent=2))
    print(f"  ✓ KB candidates: {kb_path.name}  |  metrics: {metrics_path.name}")

    def _pct(x):
        return "n/a" if x is None else f"{100 * x:.1f}%"

    # ------------------------------------------------------------------
    # Generate Markdown report
    # ------------------------------------------------------------------
    print(f"Generating markdown report at {out_file}...")

    with open(out_file, "w") as f:
        f.write("# Final Biological Report\n")
        f.write(f"**Target Antibiotic:** {antibiotic.capitalize()}\n\n")

        # --- Quantitative validation summary (top of report) ---------------
        ct = metrics['tier_counts_all']
        f.write("## Quantitative validation summary\n\n")
        f.write(f"- Candidate k-mers analysed: **{metrics['n_candidate_features']}** "
                f"(gain top-N: {metrics['n_in_gain_topN']}, "
                f"stable ≥ {stability_threshold:g}: {metrics['n_stable']})\n")
        f.write(f"- Known-mechanism recovery rate (M7 / H2): "
                f"**{_pct(metrics['known_mechanism_recovery_rate'])}** of stable k-mers "
                f"are confirmed CARD ARGs "
                f"(H2 accept ≥ 40% → **{'PASS' if metrics['H2_pass'] else 'fail' if metrics['H2_pass'] is not None else 'n/a'}**)\n")
        f.write(f"- Novel candidate fraction (H4): "
                f"**{_pct(metrics['novel_candidate_fraction'])}** of stable k-mers have no CARD hit\n")
        f.write(f"- CARD tier distribution (best hit, all candidates): "
                f"confirmed={ct.get('confirmed', 0)}, candidate={ct.get('candidate', 0)}, "
                f"weak={ct.get('weak', 0)}, none={ct.get('none', 0)}\n\n")
        f.write("> Composite score = stability × log10(1/E-value) × (identity/100); "
                "see `07_kb_candidates` for the ranked table.\n\n")
        f.write("**Methodological notes:**\n")
        # Reification-safe framing (ROADMAP §0.2 S10 / Takefuji 2025): the model
        # reports statistical association, NOT causation.
        f.write("- **Reification caveat:** a high Gain/stability score means a feature is "
                "*statistically associated with / predictive of* resistance — it does **not** "
                "*cause / determine / confer* resistance. This catalogue is a stability-filtered "
                "signal list, not a causal claim; confidence tiers explicitly separate the "
                "**statistical** layer from the **homology-based biological** layer, and weak hits "
                "are reported transparently, not claimed as findings.\n")
        f.write("- CARD hits use the **homolog model** (gene presence). A hit indicates the "
                "feature lies in a known ARG region, **not** that a resistance-conferring SNP is "
                "present — SNP/variant mechanisms (e.g. gyrA/parC fluoroquinolone mutations) are "
                "not confirmed here and need CARD's variant model / RGI.\n")
        if aro_index:
            f.write("- CARD hits are mapped to the **ARO ontology** (accession + gene family + "
                    "drug class + resistance mechanism) in `07_kb_candidates` (M16).\n")
        f.write("- Tiers grade on **identity + coverage** (alignment length / query length); "
                "E-value is a loose secondary gate only, as it is not comparable between CARD and "
                "NCBI nt (different database sizes).\n\n")

        f.write("**Confidence tiers** (coverage = alignment length / query length; checked best-first):\n")
        for tier in ("confirmed", "candidate", "weak"):
            t = tiers.get(tier, {})
            f.write(f"- `{tier}` — identity ≥ {float(t.get('min_identity', 0)):g}%, "
                    f"coverage ≥ {float(t.get('min_coverage', 0)) * 100:g}%, "
                    f"E ≤ {float(t.get('max_evalue', 0)):g}\n")
        f.write(f"- Hits below the weak floor or with E > {report_max_evalue:g} are excluded.\n\n")
        f.write("---\n\n")

        for _, row in df_features.iterrows():
            rank     = int(row['Rank'])
            score    = float(row['Gain_Score'])
            feat_id  = str(row['Feature_ID'])
            sequence = str(row['Kmer_Sequence'])

            # Reconstruct the query ID to match BLAST qseqid column
            q_id = f"Rank_{rank}|Score_{score:.4f}|Feature_{feat_id}"

            # Provenance flags: gain top-N membership + 07b stability (if present)
            sel_freq = row.get('selection_frequency', float('nan'))
            try:
                sel_freq = float(sel_freq)
            except (TypeError, ValueError):
                sel_freq = float('nan')
            flags = []
            if bool(row.get('in_gain_topN', True)):
                flags.append("gain-topN")
            if bool(row.get('stable', False)):
                flags.append(f"stable (freq={sel_freq:.2f})")
            elif sel_freq == sel_freq:      # not NaN
                flags.append(f"freq={sel_freq:.2f}")
            flag_str = f" — {', '.join(flags)}" if flags else ""

            f.write(f"### Rank {rank}: {sequence} (Gain: {score:.4f}){flag_str}\n")

            # --------------------------------------------------------------
            # CARD hits — no Entrez call needed, offline gene symbol lookup
            # --------------------------------------------------------------
            f.write("**CARD Hits (Acquired Resistance / Plasmids):**\n")
            card_hits = df_card[df_card['qseqid'] == q_id].sort_values('evalue').head(10)
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
            ncbi_hits = df_ncbi[df_ncbi['qseqid'] == q_id].sort_values('evalue').head(10)
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
