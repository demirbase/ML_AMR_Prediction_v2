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
Then point the sidebar at your amrk.db (default: results/<org>/kb/amrk.db).
"""

import sqlite3
from pathlib import Path

import pandas as pd
import streamlit as st

st.set_page_config(page_title="AMRK-DB — AMR Unitig Knowledge Base",
                   page_icon="🧬", layout="wide")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = PROJECT_ROOT / "results" / "ecoli" / "kb" / "amrk.db"


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

tab1, tab2, tab3 = st.tabs(["🔬 Biyobelirteçler", "🧩 Kanıt zinciri", "📊 Model & Provenance"])

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
    st.subheader("Modeller")
    if not models.empty:
        st.dataframe(models, use_container_width=True)
    st.subheader("Provenance (pipeline_runs)")
    if not runs.empty:
        st.dataframe(runs, use_container_width=True)
    st.caption("git_commit + config_hash + seed → her KB kaydı tam tekrarlanabilir.")

st.sidebar.caption("ROADMAP S8/N1 · CC-BY-4.0")
