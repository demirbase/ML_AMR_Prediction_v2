#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AMRK-DB explorer — a local Streamlit UI over the knowledge-base SQLite file.

This is the queryable interface for the unitig AMR biomarker knowledge base
(ROADMAP S8/N1): browse the stable/confirmed biomarkers, inspect each unitig's
multi-layer evidence chain (BLAST/CARD + ARO, R-vs-S discriminativeness, CPSS
stability, permutation, pyseer LMM), and see the run provenance.

Run locally (not part of the HPC pipeline / container):
    pip install streamlit pandas
    streamlit run scripts/kb_app.py
Then point the sidebar at your amrk.db (default: results/kb/amrk.db — the
unified multi-organism KB).
"""

import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="AMRK-DB — AMR Unitig Knowledge Base",
                   page_icon="🧬", layout="wide")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = PROJECT_ROOT / "results" / "kb" / "amrk.db"  # unified multi-organism KB


@st.cache_data(show_spinner=False)
def load_tables(db_path: str, mtime: float):
    """Load every KB table into a dict of DataFrames (mtime busts the cache)."""
    con = sqlite3.connect(db_path)
    try:
        names = [r[0] for r in con.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")]
        return {n: pd.read_sql_query(f"SELECT * FROM '{n}'", con) for n in names}
    finally:
        con.close()


def best_annotation(blast: pd.DataFrame) -> pd.DataFrame:
    """One row per unitig: the best CARD hit (confirmed > candidate > weak > none,
    then lowest E-value)."""
    if blast.empty:
        return pd.DataFrame(columns=["unitig_id"])
    tier_rank = {"confirmed": 0, "candidate": 1, "weak": 2, "none": 3}
    b = blast.copy()
    b["_t"] = b["tier"].map(tier_rank).fillna(9)
    b["_e"] = pd.to_numeric(b["evalue"], errors="coerce").fillna(1e9)
    b = b.sort_values(["unitig_id", "_t", "_e"])
    return b.groupby("unitig_id", as_index=False).first()


# --- sidebar: DB selection -------------------------------------------------
st.sidebar.title("🧬 AMRK-DB")
db_path = st.sidebar.text_input("Veritabanı yolu (amrk.db)", str(DEFAULT_DB))
if not Path(db_path).exists():
    st.warning(f"Veritabanı bulunamadı: `{db_path}`\n\n"
               "Kenar çubuğundan `amrk.db` yolunu gir (Drive yedeğinden indirdiğin dosya).")
    st.stop()

T = load_tables(db_path, Path(db_path).stat().st_mtime)
meta = T.get("kb_metadata", pd.DataFrame())
scores = T.get("unitig_model_scores", pd.DataFrame())
unitigs = T.get("unitigs", pd.DataFrame())
blast = T.get("blast_annotations", pd.DataFrame())
bg = T.get("unitig_background_frequency", pd.DataFrame())
evidence = T.get("validation_evidence", pd.DataFrame())
models = T.get("models", pd.DataFrame())
runs = T.get("pipeline_runs", pd.DataFrame())
overlap = T.get("unitig_antibiotic_overlap", pd.DataFrame())

# --- header / FAIR metadata ------------------------------------------------
st.title("AMR Unitig Biyobelirteç Bilgi Tabanı")
if not meta.empty:
    m = meta.iloc[0]
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Şema", m.get("kb_schema_version", "—"))
    c2.metric("CARD", m.get("card_version", "—") or "—")
    c3.metric("Unitig", int(m.get("n_unitigs", 0) or 0))
    c4.metric("Model", int(m.get("n_models", 0) or 0))
    c5.metric("Lisans", m.get("license", "—"))

# --- build the biomarker view ----------------------------------------------
ann = best_annotation(blast)
view = scores.merge(unitigs, on="unitig_id", how="left")
if not ann.empty:
    view = view.merge(
        ann[["unitig_id", "source_db", "gene_symbol", "tier", "identity_pct",
             "coverage", "aro_gene_family", "aro_drug_class"]],
        on="unitig_id", how="left")
if not bg.empty:
    view = view.merge(
        bg[["unitig_id", "discriminative", "prevalence_resistant",
            "prevalence_susceptible", "fisher_p"]],
        on="unitig_id", how="left")

# --- sidebar: filters ------------------------------------------------------
st.sidebar.header("Filtreler")
methods = sorted(view["selection_method"].dropna().unique()) if "selection_method" in view else []
sel_method = st.sidebar.multiselect("Seçim yöntemi", methods, default=methods)
tiers = ["confirmed", "candidate", "weak", "none"]
sel_tier = st.sidebar.multiselect("Güven seviyesi (CARD)", tiers,
                                  default=["confirmed", "candidate", "weak"])
stable_only = st.sidebar.checkbox("Sadece kararlı (stable)", value=False)
search = st.sidebar.text_input("Gen / dizi ara").strip().lower()

f = view.copy()
if sel_method:
    f = f[f["selection_method"].isin(sel_method)]
if "tier" in f.columns and sel_tier:
    f = f[f["tier"].isin(sel_tier) | f["tier"].isna() & ("none" in sel_tier)]
if stable_only and "stable" in f.columns:
    f = f[f["stable"] == 1]
if search:
    mask = f.get("gene_symbol", pd.Series(dtype=str)).fillna("").str.lower().str.contains(search)
    mask = mask | f["sequence"].fillna("").str.lower().str.contains(search)
    f = f[mask]

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs(
    ["🔬 Biyobelirteçler", "🧩 Kanıt zinciri", "📊 Model & Provenance",
     "🔗 Çapraz-antibiyotik (H3)", "✅ Dış doğrulama (M13)", "🗄️ Ham tablolar"])

with tab1:
    st.caption(f"{len(f)} unitig (filtreli). Güven seviyesi CARD identity+coverage'a dayanır.")
    cols = [c for c in ["sequence", "selection_method", "selection_frequency", "stable",
                        "mean_abs_shap", "gain", "composite_score", "gene_symbol", "tier",
                        "identity_pct", "coverage", "aro_gene_family", "aro_drug_class",
                        "discriminative", "prevalence_resistant", "prevalence_susceptible",
                        "fisher_p"] if c in f.columns]
    show = f[cols].sort_values(
        [c for c in ["stable", "selection_frequency", "mean_abs_shap"] if c in cols],
        ascending=False)
    st.dataframe(show, use_container_width=True, height=520,
                 column_config={"sequence": st.column_config.TextColumn("unitig", width="medium")})
    if "tier" in f.columns:
        st.write("**Güven dağılımı:**",
                 f["tier"].fillna("none").value_counts().to_dict())

with tab2:
    st.caption("Bir unitig seç → tüm doğrulama kanıtları (BLAST/CARD, ayırt edicilik, "
               "CPSS kararlılık, permütasyon, pyseer LMM).")
    confirmed = view[view.get("tier", pd.Series(dtype=str)) == "confirmed"] if "tier" in view else view
    opts = confirmed if not confirmed.empty else view
    label = opts.apply(lambda r: f"{r.get('gene_symbol') or '—'}  |  {str(r['sequence'])[:40]}…", axis=1)
    pick = st.selectbox("Unitig", options=list(opts["unitig_id"]),
                        format_func=lambda uid: label[opts.index[opts["unitig_id"] == uid][0]])
    if pick is not None and not evidence.empty:
        ev = evidence[evidence["unitig_id"] == pick][
            ["evidence_type", "evidence_source", "evidence_score"]]
        seq = unitigs.loc[unitigs["unitig_id"] == pick, "sequence"].iloc[0]
        st.code(seq, language="text")
        st.dataframe(ev, use_container_width=True, height=320)
        st.write(f"**{len(ev)} kanıt satırı / {ev['evidence_type'].nunique()} tür**")

with tab3:
    orgs = T.get("organisms", pd.DataFrame())
    if not orgs.empty:
        used = set(runs["organism"].dropna()) if not runs.empty else set()
        st.subheader("Organizmalar")
        st.dataframe(orgs[orgs["organism"].isin(used)] if used else orgs,
                     use_container_width=True)
    st.subheader("Modeller")
    if not models.empty:
        st.dataframe(models, use_container_width=True)
    st.subheader("Provenance (pipeline_runs)")
    if not runs.empty:
        st.dataframe(runs, use_container_width=True)
    st.caption("git_commit + config_hash + seed → her KB kaydı tam tekrarlanabilir.")

with tab4:
    st.caption("Antibiyotikler arası paylaşılan kararlı unitig'ler (S1). H3: "
               "sınıf-içi (β-laktam) overlap > sınıf-arası? Tüm çiftler listelenir "
               "(0 paylaşım = within-β-laktam ampicillin~cefotaxime, H3'ün özü).")
    abx = sorted(models["antibiotic"].dropna().unique()) if not models.empty else []
    if len(abx) >= 2 and not overlap.empty:
        import itertools as _it
        cnt = (overlap.groupby(["antibiotic_a", "antibiotic_b"])
               .size().to_dict())
        sc = (overlap.groupby(["antibiotic_a", "antibiotic_b"])["same_class"]
              .max().to_dict())
        rows = []
        for a, b in _it.combinations(abx, 2):
            n = cnt.get((a, b), cnt.get((b, a), 0))
            same = sc.get((a, b), sc.get((b, a), 0))
            rows.append({"çift": f"{a} ~ {b}",
                         "aynı registry sınıfı": "evet" if same else "hayır",
                         "paylaşılan kararlı unitig": int(n)})
        st.dataframe(pd.DataFrame(rows), use_container_width=True)
        shared = overlap.merge(unitigs, on="unitig_id", how="left")
        if not ann.empty:
            shared = shared.merge(ann[["unitig_id", "gene_symbol", "tier"]],
                                  on="unitig_id", how="left")
        st.write("**Paylaşılan unitig'ler:**")
        st.dataframe(shared[[c for c in ["antibiotic_a", "antibiotic_b", "sequence",
                     "gene_symbol", "tier", "same_class"] if c in shared.columns]],
                     use_container_width=True)
    else:
        st.info("Overlap tablosu boş / <2 antibiyotik. `15_cross_antibiotic.py` çalıştır.")

with tab5:
    st.caption("Dış doğrulama (M13, şema 0.5.0 `external_concordance` tablosu): "
               "AMRFinderPlus / ResFinder (ve varsa model) vs EUCAST/CLSI fenotip — "
               "held-out test genomlarında dengeli doğruluk (bACC), Cohen κ, FDA ME/VME.")
    ext = T.get("external_concordance", pd.DataFrame())
    if not ext.empty:
        e = ext.copy()
        if not models.empty:
            e["antibiyotik"] = e["model_id"].map(dict(zip(models["model_id"], models["antibiotic"])))
        cols = [c for c in ["antibiyotik", "caller", "reference", "n_test", "sensitivity",
                            "specificity", "balanced_accuracy", "cohen_kappa",
                            "major_error_rate", "very_major_error_rate"] if c in e.columns]
        st.dataframe(e[cols].sort_values([c for c in ["antibiyotik", "caller"] if c in cols]),
                     use_container_width=True, height=320)
        st.caption(f"{e['model_id'].nunique()} model × {e['caller'].nunique()} araç. "
                   "(Model bACC + K. pneu için `16_external_concordance.py`'yi o modellerde çalıştır.)")
    else:
        st.info("`external_concordance` tablosu boş. `16_external_concordance.py` çalıştır + "
                "`migrate_kb_050.py` ile KB'ye yükle.")

with tab6:
    # Count/version read from the KB itself — a hardcoded "13 tables (schema 0.6.0)"
    # here was wrong on both counts once 0.7.0 added unitig_evidence_tier.
    st.caption(f"KB'nin {len(T)} ham tablosu. Her tablonun/kolonun anlamı: "
               "`docs/KB_ACIKLAMA.md`.")
    _order = ["kb_metadata", "organisms", "antibiotics", "pipeline_runs", "models",
              "unitigs", "unitig_model_scores", "blast_annotations",
              "unitig_background_frequency", "variant_snp_check",
              "unitig_evidence_tier",
              "unitig_antibiotic_overlap", "validation_evidence", "external_concordance"]
    names = [t for t in _order if t in T] + [t for t in sorted(T) if t not in _order]
    tname = st.selectbox("Tablo", names)
    df = T.get(tname, pd.DataFrame())
    st.write(f"**{tname}** — {len(df)} satır × {len(df.columns)} kolon")
    st.dataframe(df, use_container_width=True, height=500)
    st.download_button(f"{tname}.csv indir", df.to_csv(index=False).encode("utf-8"),
                       file_name=f"{tname}.csv", mime="text/csv")

st.sidebar.caption("ROADMAP S8/N1 · CC-BY-4.0")
