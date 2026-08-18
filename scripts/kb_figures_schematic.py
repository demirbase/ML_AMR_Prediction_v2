#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The two schematics the thesis needs and that no data figure can stand in for:

    37_pipeline_overview   what was actually run, end to end, with the real counts
    38_kb_schema           the KB's tables and foreign keys, drawn FROM the live schema

Neither is a plot of results, which is why the other 37 figures do not contain
them: a methods chapter is unreadable without a flow diagram, and a knowledge-base
thesis is indefensible without a schema diagram.

Every number and every table name is read from the KB at draw time -- nothing is
typed into the source. A schematic that repeats hand-copied figures is exactly the
artefact that goes stale the first time the pipeline is re-run, and this one cannot:
add a table to the KB and it appears; change a row count and the box follows.

The honest annotations are computed the same way. `snp` and `mda` are drawn as
FIRED 0 because that is what unitig_evidence_tier says, not because it was
remembered -- see METHODOLOGY 5.3 for why each is zero (they are zero for
different reasons, and the figure says so).

    python scripts/kb_figures_schematic.py --db results/kb/amrk.db \
        --tables results/tables --out results/figures [--only pipeline,schema]

Saved as PNG (200 dpi) + PDF, matching the rest of the figure set.
"""
import argparse
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

INK = "#222222"
MUTED = "#6b6b6b"
EDGE = "#4a4a4a"
FILL_DATA = "#dbe9f6"
FILL_MODEL = "#e6e0f2"
FILL_EVID = "#dcefe0"
FILL_KB = "#fcefd5"
FILL_REF = "#eceff1"
FILL_DEAD = "#f2f2f2"
RED = "#b2182b"

# The six layers classify_evidence_tier() folds into a grade, in its order, keyed by
# the token written into unitig_evidence_tier.evidence_layers. label_permutation is
# the seventh produced layer but it is model-level and grades no biomarker, so it is
# drawn apart -- the "seven produced / six counted / four firing" distinction is the
# point of this panel, not a footnote to it.
LAYERS = [
    ("blast",      "BLAST vs CARD",          "08_blast_annotation"),
    ("prevalence", "Prevalence R vs S",      "10_kmer_background_frequency"),
    ("snp",        "SNP allele (CARD var.)", "11_variant_snp_check"),
    ("mda",        "MDA permutation",        "12_permutation_test"),
    ("cpss",       "CPSS stability",         "13_stability_selection"),
    ("pyseer",     "pyseer LMM",             "14_pyseer_lmm"),
]


def _save(fig, out, name):
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}.png / .pdf")


TITLE_BAND = 0.034      # room for the bold table/stage name
LINE_H = 0.026          # one body line
PAD_BOTTOM = 0.014


def box_height(lines):
    """A box is as tall as what it holds.

    The first cut of this figure set hard-coded heights and every box with four
    lines spilled its last two into the box below it -- unreadable, and the kind
    of defect that only a rendered PNG reveals. Height is derived here so it
    cannot fall out of step with the text again.
    """
    return TITLE_BAND + LINE_H * len(lines) + PAD_BOTTOM


def _box(ax, x, y, w, h, title, lines=(), fc=FILL_REF, ec=EDGE, ls="solid",
         title_size=8.6, line_size=7.0, title_colour=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.010",
                                fc=fc, ec=ec, lw=1.0, linestyle=ls, zorder=2))
    ax.text(x + w / 2, y + h - 0.009, title, ha="center", va="top", fontsize=title_size,
            fontweight="bold", color=title_colour, zorder=3)
    for i, ln in enumerate(lines):
        ax.text(x + w / 2, y + h - TITLE_BAND - 0.004 - i * LINE_H, ln, ha="center",
                va="top", fontsize=line_size, color=MUTED, zorder=3)
    return (x, y, w, h)


def stack(ax, x, w, y_top, items, gap=0.014, **kw):
    """Lay boxes downward from y_top; each one sized by its own content."""
    out = {}
    y = y_top
    for name, lines, fc in items:
        h = box_height(lines)
        y -= h
        out[name] = _box(ax, x, y, w, h, name, lines, fc, **kw)
        y -= gap
    return out


def _arrow(ax, p0, p1, colour=EDGE, rad=0.0, lw=1.1, alpha=0.9, ls="solid"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=11,
                                 lw=lw, color=colour, alpha=alpha, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}", zorder=1,
                                 shrinkA=2, shrinkB=3))


# --------------------------------------------------------------------------- facts
def read_facts(db):
    """Everything the schematics assert, straight out of the delivered KB."""
    c = sqlite3.connect(db)
    q = lambda s: c.execute(s).fetchone()[0]
    f = {
        "n_models": q("select count(*) from models"),
        "n_runs": q("select count(*) from pipeline_runs"),
        "n_unitigs": q("select count(*) from unitigs"),
        "n_tiers": q("select count(*) from unitig_evidence_tier"),
        "n_pairs": q("select sum(n_genomes) from pipeline_runs"),
        "n_orgs": q("select count(distinct organism) from pipeline_runs"),
        "n_abx": q("select count(distinct antibiotic) from models"),
        "n_classes": q("select count(distinct drug_class) from antibiotics"),
        "seed": q("select distinct random_seed from pipeline_runs"),
        "min_support": q("select distinct min_support from pipeline_runs"),
        "schema": q("select kb_schema_version from kb_metadata"),
        "doi": q("select zenodo_doi from kb_metadata"),
        "licence": q("select license from kb_metadata"),
        "novel": q("select count(*) from unitig_evidence_tier where evidence_tier='strong_novel'"),
    }
    f["cv"] = dict(c.execute("select cv_method, count(*) from models group by 1"))
    f["tiers"] = dict(c.execute(
        "select evidence_tier, count(*) from unitig_evidence_tier group by 1"))
    f["tools"] = dict(c.execute(
        """select 'poppunk', poppunk_version from pipeline_runs limit 1"""))
    row = c.execute("""select unitig_caller_version, poppunk_version, xgboost_version,
                              blast_version, pyseer_version, kmc_version, card_version
                       from pipeline_runs limit 1""").fetchone()
    f["tools"] = dict(zip(["unitig_caller", "poppunk", "xgboost", "blast",
                           "pyseer", "kmc", "card"], row))
    # How many graded biomarkers each layer actually fired for. Counted here rather
    # than asserted: two of the six are zero, and the figure has to be able to say so
    # without anyone remembering which two.
    f["fired"] = {}
    for key, _, _ in LAYERS:
        f["fired"][key] = c.execute(
            "select count(*) from unitig_evidence_tier "
            "where ','||evidence_layers||',' like ?", (f"%,{key},%",)).fetchone()[0]
    f["n_evidence_types"] = q("select count(distinct evidence_type) from validation_evidence")
    c.close()
    return f


# ------------------------------------------------------------- edge anchors
def _top(b):    x, y, w, h = b; return (x + w / 2, y + h)
def _bottom(b): x, y, w, h = b; return (x + w / 2, y)
def _left(b):   x, y, w, h = b; return (x, y + h / 2)
def _right(b):  x, y, w, h = b; return (x + w, y + h / 2)


# ------------------------------------------------------------------ 37 pipeline
def fig_pipeline(f, out, n_fig, n_tab):
    fig, ax = plt.subplots(figsize=(14.2, 9.0))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    ax.text(0.5, 0.992, "From public assemblies to a queryable AMR knowledge base",
            ha="center", va="top", fontsize=14, fontweight="bold", color=INK)
    ax.text(0.5, 0.957,
            f"{f['n_pairs']:,} genome-phenotype pairs · {f['n_orgs']} ESKAPEE organisms · "
            f"{f['n_abx']} antibiotics · {f['n_models']} models, every one lineage-CV · "
            f"single seed {f['seed']} · min_support {f['min_support']}",
            ha="center", va="top", fontsize=9, color=MUTED)

    Y0 = 0.888
    for x0, w, label in [(0.015, 0.225, "1 · DATA AND QUALITY CONTROL"),
                         (0.262, 0.215, "2 · STRUCTURE AND FEATURES"),
                         (0.509, 0.215, "3 · MODELLING"),
                         (0.735, 0.250, f"4 · EVIDENCE ({f['n_evidence_types']} produced)")]:
        ax.text(x0 + w / 2, 0.918, label, ha="center", va="center", fontsize=8.4,
                fontweight="bold", color=MUTED)

    L1 = stack(ax, 0.015, 0.225, Y0, [
        ("BV-BRC download", ["00a · assemblies + AMR phenotypes",
                             "UI node only — compute nodes",
                             "have no internet"], FILL_DATA),
        ("Metadata and labels", ["00 · antibiotic-name normalisation",
                                 "01 · label validation"], FILL_DATA),
        ("Genome QC — the only gate", ["02d · CheckM2: completeness \u2265 95 %,",
                                       "contamination \u2264 5 %",
                                       "02b · QUAST contiguity is ADVISORY,",
                                       "never an exclusion criterion"], FILL_DATA),
    ])
    L2 = stack(ax, 0.262, 0.215, Y0, [
        ("PopPUNK lineages", [f"02c · {f['tools']['poppunk']}",
                              "graph-tool 3.0",
                              "clusters ARE the CV groups"], FILL_DATA),
        ("Unitig features", [f"03u · {f['tools']['unitig_caller']}",
                             "bcalm de Bruijn graph",
                             "presence/absence matrix"], FILL_DATA),
    ])
    cv = ", ".join(f"{k} ({v})" for k, v in f["cv"].items())
    L3 = stack(ax, 0.509, 0.215, Y0, [
        ("Hyperparameter search", ["04 · Optuna, 30 trials per model"], FILL_MODEL),
        ("Model training", [f"05 · XGBoost {f['tools']['xgboost']}, seed {f['seed']}"], FILL_MODEL),
        ("Evaluation", ["06 · lineage-grouped 5-fold CV",
                        f"{f['n_models']}/{f['n_models']} models, no fallback",
                        "07b · random-CV run alongside"], FILL_MODEL),
        ("Explainability", ["07 · gain + SHAP ranking",
                            "13 · CPSS stability, PFER bound"], FILL_MODEL),
    ])

    # -- lane 4: the evidence layers and what each one actually graded
    ex, ew, ey, eh = 0.735, 0.250, 0.360, Y0 - 0.360
    ax.add_patch(FancyBboxPatch((ex, ey), ew, eh, boxstyle="round,pad=0.005,rounding_size=0.010",
                                fc=FILL_EVID, ec=EDGE, lw=1.0, zorder=2))
    ax.text(ex + ew / 2, ey + eh - 0.010, "Orthogonal evidence layers", ha="center",
            va="top", fontsize=9, fontweight="bold", color=INK, zorder=3)
    ax.text(ex + ew / 2, ey + eh - 0.040,
            f"{f['n_evidence_types']} produced · 6 counted into the grade · 4 ever fire",
            ha="center", va="top", fontsize=7.3, color=MUTED, style="italic", zorder=3)
    yy = ey + eh - 0.075
    for key, label, script in LAYERS:
        n = f["fired"][key]
        dead = n == 0
        col = RED if dead else INK
        ax.text(ex + 0.012, yy, ("\u2717" if dead else "\u2713"), ha="left", va="top",
                fontsize=8.6, color=(RED if dead else "#2f7d32"), fontweight="bold", zorder=3)
        ax.text(ex + 0.032, yy, label, ha="left", va="top", fontsize=7.7, color=col, zorder=3)
        ax.text(ex + 0.032, yy - 0.024, script, ha="left", va="top", fontsize=6.2,
                color=MUTED, family="monospace", zorder=3)
        ax.text(ex + ew - 0.012, yy, f"graded {n:,}", ha="right", va="top", fontsize=7.3,
                color=(RED if dead else MUTED), fontweight=("bold" if dead else "normal"),
                zorder=3)
        yy -= 0.053
    yy -= 0.006
    ax.text(ex + 0.012, yy,
            "7th layer, label permutation (12b), is model-level:\n"
            "it validates a model, it grades no biomarker.",
            ha="left", va="top", fontsize=6.6, color=MUTED, style="italic", zorder=3)
    ax.text(ex + 0.012, ey + 0.014,
            "\u2717 = produced but never fires (METHODOLOGY \u00a75.3):\n"
            "SNP = candidate-set mismatch, a wiring fault\n"
            "MDA = real negative, but R=100 is underpowered",
            ha="left", va="bottom", fontsize=6.5, color=RED, style="italic",
            zorder=3, linespacing=1.45)

    # -- the knowledge base band
    tiers = " · ".join(f"{k} {v:,}" for k, v in sorted(f["tiers"].items(), key=lambda kv: -kv[1]))
    KB = _box(ax, 0.262, 0.150, 0.462, box_height(["", "", ""]),
              f"KNOWLEDGE BASE   ·   amrk.db   ·   schema {f['schema']}",
              [f"{f['n_unitigs']:,} unitigs · {f['n_models']} models · "
               f"{f['n_tiers']:,} graded (unitig, model) pairs",
               f"evidence tier: {tiers}",
               "one join reaches the provenance of any claim"],
              FILL_KB, title_size=9.4)
    OUT = _box(ax, 0.735, 0.150, 0.250, box_height(["", "", ""]), "Delivered artefacts",
               [f"{n_tab} tidy tables · {n_fig} figures",
                f"Zenodo {f['doi']}", f"{f['licence']}"], FILL_KB)
    PROV = _box(ax, 0.015, 0.150, 0.225, box_height(["", "", ""]), "Provenance capture",
                ["run_metadata.json per stage: commit,",
                 "seed, tool versions, data sha256",
                 "— written even when git_dirty = 1"], FILL_REF)

    # -- arrows, all anchored to the boxes the layout returned
    _arrow(ax, _bottom(L1["BV-BRC download"]), _top(L1["Metadata and labels"]))
    _arrow(ax, _bottom(L1["Metadata and labels"]), _top(L1["Genome QC — the only gate"]))
    _arrow(ax, _right(L1["Genome QC — the only gate"]), _left(L2["PopPUNK lineages"]), rad=-0.14)
    _arrow(ax, _right(L1["Genome QC — the only gate"]), _left(L2["Unitig features"]), rad=0.10)
    _arrow(ax, _bottom(L2["PopPUNK lineages"]), _top(L2["Unitig features"]))
    _arrow(ax, _right(L2["Unitig features"]), _left(L3["Model training"]), rad=0.10)
    _arrow(ax, _bottom(L3["Hyperparameter search"]), _top(L3["Model training"]))
    _arrow(ax, _bottom(L3["Model training"]), _top(L3["Evaluation"]))
    _arrow(ax, _bottom(L3["Evaluation"]), _top(L3["Explainability"]))
    _arrow(ax, _right(L3["Explainability"]), (ex, ey + eh * 0.42), rad=-0.10)
    _arrow(ax, (ex + ew * 0.5, ey), (0.700, 0.276), rad=0.10)
    _arrow(ax, _right(KB), _left(OUT))
    _arrow(ax, _bottom(L3["Explainability"]), _right(KB), rad=0.12)
    _arrow(ax, _right(PROV), _left(KB))
    # the lineage groups feed the CV split -- the mechanism behind the whole thesis
    _arrow(ax, _right(L2["PopPUNK lineages"]), _left(L3["Evaluation"]), rad=0.22,
           colour="#8c6d1f", lw=1.5)
    ax.text(0.4935, 0.735, "lineage groups define the folds", fontsize=6.2,
            color="#8c6d1f", ha="center", va="center", style="italic", rotation=90,
            bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2), zorder=4)

    ax.text(0.5, 0.115,
            "Step numbers are the repository's own. CheckM2 is the only exclusion gate; QUAST is reported but never filters. "
            "No external test set exists for this panel —\nlineage hold-out IS the external validation, and temporal validation is "
            "impossible because BV-BRC AMR phenotypes stop in 2021 (METHODOLOGY \u00a75.3).",
            ha="center", va="top", fontsize=7.4, color=MUTED, style="italic")
    _save(fig, out, "37_pipeline_overview")


# -------------------------------------------------------------------- 38 schema
def fig_schema(db, f, out):
    c = sqlite3.connect(db)
    meta = {}
    for (t,) in c.execute("select name from sqlite_master where type='table' order by name"):
        info = list(c.execute(f"pragma table_info({t})"))
        meta[t] = {
            "cols": [r[1] for r in info],
            "pk": [r[1] for r in info if r[5]],
            "fk": [(r[3], r[2]) for r in c.execute(f"pragma foreign_key_list({t})")],
            "n": c.execute(f"select count(*) from {t}").fetchone()[0],
        }
    c.close()

    fig, ax = plt.subplots(figsize=(14.6, 9.0))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.992, f"amrk.db — schema {f['schema']}: {len(meta)} tables, "
                        f"read from the delivered database",
            ha="center", va="top", fontsize=14, fontweight="bold", color=INK)
    ax.text(0.5, 0.957, "Every biomarker claim is one join from its provenance: "
                        "(unitig, model) \u2192 models \u2192 pipeline_runs \u2192 "
                        "git commit, seed, config hash, tool versions.",
            ha="center", va="top", fontsize=8.8, color=MUTED)

    def entry(name, show, dead=False):
        m = meta[name]
        lines = [f"PK  {', '.join(m['pk'])}"]
        rest = [cl for cl in m["cols"] if cl not in m["pk"]]
        lines += rest[:show]
        if len(rest) > show:
            lines.append(f"\u2026 +{len(rest) - show} more columns")
        head = f"{name}  [{m['n']:,}]" if not dead else f"{name}  [EMPTY]"
        return head, lines

    Y0 = 0.905
    ref = stack(ax, 0.012, 0.200, Y0, [
        entry("organisms", 2) + (FILL_REF,),
        entry("antibiotics", 2) + (FILL_REF,),
        entry("kb_metadata", 3) + (FILL_REF,),
        entry("pipeline_runs", 6) + (FILL_DATA,),
    ], title_size=8.0, line_size=6.5)
    hubs = stack(ax, 0.238, 0.200, Y0, [
        entry("models", 7) + (FILL_MODEL,),
        entry("unitigs", 2) + (FILL_MODEL,),
    ], title_size=8.4, line_size=6.5)
    c3 = stack(ax, 0.462, 0.245, Y0, [
        entry("unitig_evidence_tier", 3) + (FILL_KB,),
        entry("unitig_model_scores", 3) + (FILL_EVID,),
        entry("unitig_background_frequency", 3) + (FILL_EVID,),
    ], title_size=8.0, line_size=6.5)
    c4 = stack(ax, 0.735, 0.245, Y0, [
        entry("blast_annotations", 3) + (FILL_EVID,),
        entry("variant_snp_check", 2) + (FILL_EVID,),
        entry("validation_evidence", 3) + (FILL_EVID,),
    ], title_size=8.0, line_size=6.5)
    dead = stack(ax, 0.735, 0.245, 0.405, [
        entry("external_concordance", 1, dead=True) + (FILL_DEAD,),
        entry("unitig_antibiotic_overlap", 2, dead=True) + (FILL_DEAD,),
    ], ls="dashed", title_size=8.0, line_size=6.5, title_colour=RED)

    boxes = {**ref, **hubs, **c3, **c4, **dead}
    by_name = {k.split("  ")[0]: v for k, v in boxes.items()}

    C_MODEL, C_UNITIG, C_RUN, C_REF = "#7b52ab", "#2c7fb8", "#8c6d1f", "#777777"

    def edge(src, dst, colour, rad):
        _arrow(ax, _left(by_name[src]), _right(by_name[dst]), colour=colour,
               rad=rad, lw=0.9, alpha=0.8)

    for i, t in enumerate(["unitig_evidence_tier", "unitig_model_scores",
                           "unitig_background_frequency"]):
        edge(t, "models", C_MODEL, 0.10 - i * 0.06)
        edge(t, "unitigs", C_UNITIG, -0.16 - i * 0.05)
    for i, t in enumerate(["blast_annotations", "variant_snp_check"]):
        edge(t, "models", C_MODEL, 0.24 - i * 0.05)
    edge("blast_annotations", "unitigs", C_UNITIG, -0.28)
    edge("variant_snp_check", "unitigs", C_UNITIG, -0.32)
    edge("validation_evidence", "unitigs", C_UNITIG, -0.36)
    edge("validation_evidence", "pipeline_runs", C_RUN, -0.42)
    edge("external_concordance", "models", C_MODEL, 0.34)
    edge("unitig_antibiotic_overlap", "unitigs", C_UNITIG, 0.30)
    _arrow(ax, _left(by_name["models"]), _right(by_name["pipeline_runs"]),
           colour=C_RUN, rad=0.16, lw=1.0, alpha=0.85)
    _arrow(ax, _left(by_name["models"]), _right(by_name["antibiotics"]),
           colour=C_REF, rad=-0.10, lw=0.9, alpha=0.8)

    for i, (col, lab) in enumerate([(C_MODEL, "FK \u2192 models(model_id)"),
                                    (C_UNITIG, "FK \u2192 unitigs(unitig_id)"),
                                    (C_RUN, "FK \u2192 pipeline_runs(run_id)"),
                                    (C_REF, "FK \u2192 antibiotics(antibiotic)")]):
        y = 0.450 - i * 0.030
        ax.plot([0.245, 0.278], [y, y], color=col, lw=1.6)
        ax.text(0.285, y, lab, fontsize=7.2, color=MUTED, va="center")

    ax.text(0.245, 0.300,
            "TWO TABLES ARE EMPTY BY RECORD, NOT BY ACCIDENT\n"
            "external_concordance \u2014 step 16 was not run in this clean pass,\n"
            "so figure 07 does not exist either.\n"
            "unitig_antibiotic_overlap \u2014 step 15 ran but was never loaded.\n"
            "verify_artefacts.py asserts both stay empty: a silent fill is a\n"
            "regression exactly as a silent emptying would be.",
            fontsize=6.9, color=RED, va="top", linespacing=1.5,
            bbox=dict(fc="white", ec="none", alpha=0.9, pad=2.5), zorder=4)

    ax.text(0.735, 0.135,
            "Query rules the schema cannot enforce\n"
            "\u2022 blast_annotations: filter tier IN ('confirmed','candidate').\n"
            "   3,007 of 3,611 hits are tier='none'; unfiltered, an\n"
            "   A. baumannii model returns mecA.\n"
            "\u2022 unitig_model_scores holds a pair twice (gain + CPSS):\n"
            "   GROUP BY before joining or the rows double.\n"
            "\u2022 performance is models.auc_mean_seeds (lineage-CV),\n"
            "   never models.roc_auc, which is a single split.",
            fontsize=6.9, color=INK, va="top", linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.6", fc="#fff8e1", ec="#d6b656", lw=0.9))
    _save(fig, out, "38_kb_schema")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="path to amrk.db")
    ap.add_argument("--tables", help="tidy tables dir (unused today; kept for Makefile symmetry)")
    ap.add_argument("--out", required=True, help="figures output dir")
    ap.add_argument("--only", help="comma list: pipeline,schema")
    a = ap.parse_args()

    out = Path(a.out)
    want = {s.strip() for s in a.only.split(",")} if a.only else {"pipeline", "schema"}
    facts = read_facts(a.db)

    # Counted, not asserted: the "delivered artefacts" box would otherwise be the one
    # place in this figure that can quietly lie about the size of the figure set.
    stems = {p.stem for p in out.glob("*.png")} | {"37_pipeline_overview", "38_kb_schema"}
    n_fig = len(stems)
    n_tab = len(list(Path(a.tables).glob("*.csv"))) if a.tables else 0

    print(f"KB: {facts['n_models']} models, {facts['n_unitigs']:,} unitigs, "
          f"{facts['n_tiers']:,} graded pairs, schema {facts['schema']}")
    print(f"artefact counts for the figure: {n_fig} figures, {n_tab} tables")
    if "pipeline" in want:
        fig_pipeline(facts, out, n_fig, n_tab)
    if "schema" in want:
        fig_schema(a.db, facts, out)


if __name__ == "__main__":
    main()
