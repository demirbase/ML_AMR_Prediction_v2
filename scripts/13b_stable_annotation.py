#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Annotate the CPSS-stable unitig set — Step 13b

Bridges step 13 (CPSS stability selection + SHAP) to the knowledge base. The
stable unitigs are a *different, statistically-grounded* set from the Gain-top
candidates, so they need their own CARD annotation. This reuses step 09's exact
tiering (identity + COVERAGE, so the spurious ~14 bp hits a short query throws
against the huge CARD DB are graded 'none', not over-claimed) and ARO ontology
mapping — no logic is duplicated.

Flow
----
1. Read the stable set from ``13_stability_selection_{ab}.csv``.
2. BLAST it against the local CARD homolog DB (``blastn-short``, word_size 7 —
   correct for short unitig queries; same as step 08's CARD pass).
3. For each stable unitig: best CARD hit -> identity/coverage tier (09's
   ``classify_confidence``) + gene symbol + ARO accession/family/class/mechanism.
4. Merge in the CPSS ``selection_frequency`` and mean |SHAP|, compute the
   composite score, and write a KB-ready table.

Output (results/{org}/{ab}/05_explainability/)
    13_stable_kb_candidates_{ab}.csv  — KB-ready: kmer, CPSS freq, SHAP, CARD gene,
        identity, coverage, tier, ARO fields, composite_score
"""

import argparse
import importlib
import subprocess
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config, resolve_path, resolve_tool, get_target  # noqa: E402

_b = importlib.import_module("09_biological_summary")  # reuse tiering + ARO


def run_card_blast(fasta, card_db, out_tsv):
    """blastn-short vs the local CARD DB, emitting 09's expected outfmt-6 cols."""
    blastn = resolve_tool("blastn", "blastn")
    cmd = [blastn, "-query", str(fasta), "-db", str(card_db),
           "-task", "blastn-short", "-word_size", "7", "-dust", "no",
           "-outfmt", "6 " + " ".join(_b.TSV_COLS),
           "-max_target_seqs", "5", "-evalue", "10", "-out", str(out_tsv)]
    subprocess.run(cmd, check=True)


def annotate(stable_df, blast_df, tiers, aro_index):
    """Best CARD hit per stable unitig -> tiered, ARO-mapped, CPSS-scored record."""
    rows = []
    for _, s in stable_df.iterrows():
        fidx = int(s["feature_index"])
        kmer = str(s.get("kmer", ""))
        rec = {
            "feature_index": fidx, "kmer": kmer,
            "selection_frequency": s.get("selection_frequency"),
            "mean_abs_shap": s.get("mean_abs_shap"),
            "stable": 1,
            "card_gene": "", "card_identity": None, "coverage": None,
            "card_evalue": None, "confidence_tier": "none", "has_card_hit": 0,
            "composite_score": None,
            "aro_accession": "", "aro_gene_family": "",
            "aro_drug_class": "", "aro_resistance_mechanism": "",
        }
        # qseqid is a bare integer here, so read_blast_tsv infers it as int64 —
        # compare as strings on both sides.
        hits = (blast_df[blast_df["qseqid"].astype(str) == str(fidx)]
                if not blast_df.empty else blast_df)
        if hits is not None and not hits.empty:
            best = hits.loc[hits["evalue"].idxmin()]
            qlen = float(best["qlen"]) if best["qlen"] == best["qlen"] else (len(kmer) or 1)
            tier = _b.classify_confidence(best["pident"], best["evalue"],
                                          best["length"], qlen, tiers)
            aro_acc = _b.aro_from_sseqid(best["sseqid"])
            aro = aro_index.get(aro_acc, {})
            rec.update(
                card_gene=_b.extract_card_gene(best["sseqid"]),
                card_identity=float(best["pident"]),
                coverage=(float(best["length"]) / qlen if qlen else 0.0),
                card_evalue=float(best["evalue"]),
                confidence_tier=tier, has_card_hit=1,
                composite_score=_b.composite_score(s.get("selection_frequency"),
                                                   best["pident"], best["evalue"]),
                aro_accession=(f"ARO:{aro_acc}" if aro_acc else ""),
                aro_gene_family=aro.get("gene_family", ""),
                aro_drug_class=aro.get("drug_class", ""),
                aro_resistance_mechanism=aro.get("resistance_mechanism", ""),
            )
        rows.append(rec)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description="Annotate CPSS-stable unitigs (CARD tier + ARO).")
    ap.add_argument("--organism", default=None)
    ap.add_argument("--antibiotic", default=None)
    args = ap.parse_args()

    config = load_config()
    organism, antibiotic = get_target(args, config=config)

    out_dir = resolve_path("dir_05_explainability", organism=organism,
                           antibiotic=antibiotic, config=config)
    stable_csv = out_dir / f"13_stability_selection_{antibiotic}.csv"
    if not stable_csv.exists():
        print(f"ERROR: {stable_csv} not found — run 13_stability_selection.py first.")
        sys.exit(1)

    df = pd.read_csv(stable_csv, encoding="utf-8")
    stable_df = df[(df["stable"] == 1) & (df["kmer"].astype(str).str.len() > 0)].copy()
    print(f"  stable unitigs to annotate: {len(stable_df)}")

    # FASTA keyed by feature_index (clean join key).
    fasta = out_dir / f"13_stable_features_{antibiotic}.fasta"
    with open(fasta, "w", encoding="utf-8") as fh:
        for _, r in stable_df.iterrows():
            fh.write(f">{int(r['feature_index'])}\n{r['kmer']}\n")

    blast_cfg = config.get("blast", {})
    card_db = (PROJECT_ROOT / blast_cfg.get("card_db_dir", "data/blast_db/card_nt")
               / blast_cfg.get("card_db_name", "card"))
    blast_tsv = out_dir / f"13_stable_card_hits_{antibiotic}.tsv"
    print("  BLAST vs local CARD (blastn-short)...", flush=True)
    run_card_blast(fasta, card_db, blast_tsv)

    blast_df = _b.read_blast_tsv(blast_tsv)
    tiers, *_ = _b.load_tiers(config)
    aro_index = _b.load_aro_index(
        PROJECT_ROOT / "data" / "external" / "card" / "aro_index.tsv")

    out = annotate(stable_df, blast_df, tiers, aro_index)
    out = out.sort_values(["confidence_tier", "selection_frequency", "mean_abs_shap"],
                          ascending=[True, False, False])
    out_csv = out_dir / f"13_stable_kb_candidates_{antibiotic}.csv"
    out.to_csv(out_csv, index=False, encoding="utf-8")

    tier_counts = out["confidence_tier"].value_counts().to_dict()
    n_card = int((out["has_card_hit"] == 1).sum())
    print("\n" + "=" * 70)
    print(f"  annotated {len(out)} stable unitigs | CARD-hit (any) {n_card} | tiers {tier_counts}")
    print(f"  ✓ {out_csv.name}")
    print("=" * 70)


if __name__ == "__main__":
    main()
