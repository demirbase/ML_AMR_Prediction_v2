#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""One-command thesis-ready results summary from the AMRK-DB knowledge base.

Reads a populated ``amrk.db`` and renders a Markdown snapshot of everything the
thesis Results section needs — per-antibiotic model performance (lineage-CV +
single-split AUC, MCC, tree count), CPSS-stable counts, confirmed CARD gene
families, statistical validation (pyseer LMM), external concordance (M13
AMRFinderPlus/ResFinder + model head-to-head), and the cross-antibiotic overlap
(S1/H3). Pure read-only over ``lib/kb_queries`` + a few report aggregates, so it
is reproducible and stays in sync with whatever is actually in the KB.

Usage:
  python scripts/kb_report.py                       # -> stdout + results/{org}/kb/KB_REPORT_{org}.md
  python scripts/kb_report.py --db path/to/amrk.db --organism ecoli
"""

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib import kb_queries as Q  # noqa: E402


def _f(x, nd=3):
    return "—" if x is None else f"{x:.{nd}f}"


def build_report(conn):
    """Return the KB results summary as a Markdown string."""
    m = Q.get_metadata(conn)
    s = Q.get_stats(conn)
    L = []
    L.append("# AMRK-DB — results summary")
    L.append("")
    L.append(f"- **schema** {m.get('kb_schema_version','?')} · **CARD** "
             f"{m.get('card_version','?')} · **license** {m.get('license','?')} · "
             f"**Zenodo DOI** {m.get('zenodo_doi') or '_reserved_'}")
    L.append(f"- **{s['n_unitigs']}** unitigs · **{s['n_models']}** models · "
             f"**{s['n_evidence']}** evidence rows · antibiotics: "
             f"{', '.join(m.get('antibiotics', []))}")

    # --- per-antibiotic model performance -----------------------------------
    models = conn.execute(
        """SELECT m.model_id, m.antibiotic, a.drug_class, m.roc_auc, m.mcc, m.n_trees,
                  m.auc_mean_seeds, m.auc_std_seeds
             FROM models m LEFT JOIN antibiotics a ON a.antibiotic = m.antibiotic
            ORDER BY m.antibiotic""").fetchall()
    L.append("\n## Model performance (per antibiotic)")
    L.append("| antibiotic | class | lineage-CV ROC-AUC | single-split AUC | MCC | trees | "
             "CPSS stable | confirmed gene families |")
    L.append("|---|---|---|---|---|---|---|---|")
    for md in models:
        mid = md["model_id"]
        n_stable = conn.execute(
            "SELECT COUNT(DISTINCT unitig_id) FROM unitig_model_scores "
            "WHERE model_id=? AND stable=1 AND selection_method='cpss'", (mid,)).fetchone()[0]
        fams = [r[0] for r in conn.execute(
            "SELECT DISTINCT aro_gene_family FROM blast_annotations "
            "WHERE model_id=? AND tier='confirmed' AND aro_gene_family IS NOT NULL",
            (mid,)).fetchall()]
        cv = (f"{_f(md['auc_mean_seeds'])}±{_f(md['auc_std_seeds'])}"
              if md["auc_mean_seeds"] is not None else "—")
        L.append(f"| {md['antibiotic']} | {md['drug_class'] or '—'} | {cv} | "
                 f"{_f(md['roc_auc'])} | {_f(md['mcc'])} | {md['n_trees'] or '—'} | "
                 f"{n_stable} | {', '.join(sorted(fams)) or '—'} |")

    # --- statistical + external validation (validation_evidence) ------------
    run_ab = dict(conn.execute("SELECT run_id, antibiotic FROM models").fetchall())

    def ev_for(ab, etype, agg="count"):
        rows = conn.execute(
            "SELECT evidence_source, evidence_score, pipeline_run_id FROM validation_evidence "
            "WHERE evidence_type=?", (etype,)).fetchall()
        rows = [r for r in rows if run_ab.get(r[2]) == ab]
        if agg == "count":
            return len(rows)
        return rows

    L.append("\n## Statistical & external validation")
    L.append("| antibiotic | pyseer-LMM sig | AMRFinderPlus (κ) | ResFinder (κ) | "
             "model head-to-head |")
    L.append("|---|---|---|---|---|")
    for ab in m.get("antibiotics", []):
        pys = ev_for(ab, "pyseer_lmm")
        afp = ev_for(ab, "concordance_amrfinderplus", "rows")
        rf = ev_for(ab, "concordance_resfinder", "rows")
        h2h = ev_for(ab, "head_to_head_model", "rows")
        L.append(f"| {ab} | {pys} | {_f(afp[0][1]) if afp else '—'} | "
                 f"{_f(rf[0][1]) if rf else '—'} | "
                 f"{(h2h[0][0].split('(')[-1].rstrip(') ')) if h2h else '—'} |")

    # --- cross-antibiotic overlap (S1/H3) ----------------------------------
    # Enumerate ALL antibiotic pairs (the overlap table only stores pairs with a
    # shared unitig, so a 0-overlap pair — e.g. the within-β-lactam ampicillin~
    # cefotaxime, the H3 crux — must be shown explicitly as 0).
    import itertools
    counts, same = {}, {}
    for r in conn.execute("SELECT antibiotic_a, antibiotic_b, COUNT(*) n, MAX(same_class) sc "
                          "FROM unitig_antibiotic_overlap GROUP BY antibiotic_a, antibiotic_b"):
        counts[frozenset((r[0], r[1]))] = r["n"]
        same[frozenset((r[0], r[1]))] = r["sc"]
    abs_list = m.get("antibiotics", [])
    L.append("\n## Cross-antibiotic stable-unitig overlap (S1/H3)")
    if len(abs_list) >= 2:
        L.append("| pair | same registry class | shared stable unitigs |")
        L.append("|---|---|---|")
        for a, b in itertools.combinations(abs_list, 2):
            key = frozenset((a, b))
            n = counts.get(key, 0)
            sc = "yes" if same.get(key) else "no"
            L.append(f"| {a} ~ {b} | {sc} | {n} |")
        L.append("\n_H3 (within-class β-lactam overlap > cross-class) is evaluated in "
                 "`15_cross_antibiotic_summary_{organism}.json` at unitig and ARO "
                 "gene-family level._".replace("{organism}", m.get("organism", "ecoli")))
    else:
        L.append("_need ≥2 antibiotics for overlap (run 15_cross_antibiotic.py)_")
    L.append("")
    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description="Thesis-ready KB results summary.")
    ap.add_argument("--organism", default="ecoli")
    ap.add_argument("--db", default=None, help="KB path (default: results/{org}/kb/amrk.db)")
    args = ap.parse_args()
    db_path = Path(args.db) if args.db else (
        PROJECT_ROOT / "results" / args.organism / "kb" / "amrk.db")
    if not db_path.exists():
        sys.exit(f"KB not found: {db_path}")
    conn = Q.connect(db_path)
    try:
        report = build_report(conn)
    finally:
        conn.close()
    print(report)
    out = db_path.parent / f"KB_REPORT_{args.organism}.md"
    out.write_text(report, encoding="utf-8")
    print(f"\n[written] {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
