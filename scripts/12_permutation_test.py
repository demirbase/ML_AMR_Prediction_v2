#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Permutation feature-importance test (MDA) — Step 12  (ROADMAP §0.2 / must-have M9)

Why
---
BLAST says *which gene* a candidate unitig maps to and step 10 says *whether it
discriminates* R vs S — but neither asks whether the **model actually relies on
that feature**. M9 answers that with a **permutation importance / Mean Decrease
in Accuracy (MDA)** test: keep the trained model fixed (NO retraining — ROADMAP
§0.2 explicitly rules out the hundreds-of-refits null), permute one candidate
feature's column in the held-out test set, and measure how much the test ROC-AUC
drops. A feature the model genuinely uses → consistent AUC drop → low empirical
p; an irrelevant feature → ~no change.

Method (defensible + cheap)
---------------------------
* Reuse the **exact held-out test set** from the experiment config (same set 06
  reported AUC on) via 06's loaders — so the baseline AUC matches the headline.
* Features are **binary** (presence/absence, max_bin=2), so permuting a column =
  redistributing its 1s among the test genomes at random, keeping the count.
* XGBoost prediction only changes for the genomes whose feature value flipped, so
  each permutation re-predicts **only the changed rows** (fast) and overlays them
  on the baseline prediction vector before recomputing AUC.
* A candidate the model never splits on (not in ``get_score``) cannot change the
  prediction → MDA = 0 by construction; we skip its permutations.
* Significance: empirical ``p = (1 + #{perm_auc >= baseline}) / (R + 1)`` over
  ``R`` permutations, then **Benjamini-Hochberg FDR** across candidates (ROADMAP
  §0.2: Q < 0.05).

Output (results/{org}/{ab}/05_explainability/)
    12_permutation_test_{ab}.csv      — per-candidate mda_auc_drop, perm_p, perm_q, significant
    12_permutation_summary_{ab}.json  — baseline AUC, R, #significant, params
"""

import argparse
import importlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xgboost as xgb
from scipy.sparse import csr_matrix
from sklearn.metrics import roc_auc_score

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

from lib.config import load_config, resolve_path, get_target  # noqa: E402

# Reuse 06's EXACT held-out test-set loaders (single source of truth for the
# split) — importing by a digit-leading module name needs importlib.
_ev = importlib.import_module("06_evaluation")


def benjamini_hochberg(pvals):
    """BH-FDR q-values for a 1-D array of p-values (NaN-safe; NaN stays NaN)."""
    p = np.asarray(pvals, dtype=float)
    q = np.full_like(p, np.nan)
    mask = ~np.isnan(p)
    idx = np.where(mask)[0]
    if idx.size == 0:
        return q
    pm = p[idx]
    order = np.argsort(pm)
    ranked = pm[order]
    m = ranked.size
    q_sorted = ranked * m / (np.arange(1, m + 1))
    # enforce monotonicity from the largest p downwards
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_sorted = np.clip(q_sorted, 0, 1)
    q_back = np.empty_like(q_sorted)
    q_back[order] = q_sorted
    q[idx] = q_back
    return q


def load_candidates(results_root, antibiotic):
    """Per-candidate (kmer, matrix column index). Prefer step-10's table (it
    carries feature_index); fall back to 07's KB candidates."""
    def _find(name):
        hits = sorted(Path(results_root).rglob(name))
        return hits[-1] if hits else None

    path = _find(f"10_kmer_background_frequency_{antibiotic}.csv")
    if path is None:
        path = _find(f"07_kb_candidates_{antibiotic}.csv")
    if path is None:
        raise FileNotFoundError(
            f"No candidate table (10_kmer_background_frequency / 07_kb_candidates) "
            f"for {antibiotic} under {results_root}"
        )
    df = pd.read_csv(path, encoding="utf-8")
    if "feature_index" not in df.columns:
        raise KeyError(
            f"{path.name} has no 'feature_index' column — cannot map candidates to "
            f"matrix columns. Re-run step 10 (it emits feature_index)."
        )
    return df, path


def main():
    ap = argparse.ArgumentParser(description="MDA permutation importance (M9).")
    ap.add_argument("--organism", default=None)
    ap.add_argument("--antibiotic", default=None)
    ap.add_argument("--n-permutations", type=int, default=100,
                    help="permutations per candidate (ROADMAP: 100)")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    config = load_config()
    organism, antibiotic = get_target(args, config=config)
    R = args.n_permutations
    rng = np.random.default_rng(args.seed)

    matrix_dir = resolve_path("matrix_dir", organism=organism, antibiotic=antibiotic, config=config)
    models_dir = resolve_path("models_dir", organism=organism, antibiotic=antibiotic, config=config)
    out_dir = resolve_path("dir_05_explainability", organism=organism, antibiotic=antibiotic, config=config)
    out_dir.mkdir(parents=True, exist_ok=True)
    results_root = PROJECT_ROOT / "results" / organism / antibiotic

    print("=" * 78)
    print(f"PERMUTATION IMPORTANCE (MDA)  —  {organism} / {antibiotic}  |  R={R}")
    print("=" * 78)

    # --- model -------------------------------------------------------------
    model_path = next(
        (p for p in [
            models_dir / f"xgboost_{antibiotic}_final_v2.json",
            models_dir / f"xgboost_{antibiotic}.json",
        ] if p.exists()), None)
    if model_path is None:
        print(f"ERROR: no model in {models_dir}"); sys.exit(1)
    model = xgb.Booster()
    model.load_model(str(model_path))
    used_cols = {int(k[1:]) for k in model.get_score(importance_type="gain")}
    print(f"  model: {model_path.name}  | trees={model.num_boosted_rounds()} "
          f"| features used by model={len(used_cols)}")

    # --- exact held-out test set (06's loaders) ----------------------------
    y_all = pd.read_csv(matrix_dir / f"y_{antibiotic}.csv", encoding="utf-8")["label"].values
    chunk_files = sorted(matrix_dir.glob(f"X_{antibiotic}_part_*.npz"),
                         key=lambda x: int(x.stem.split("_")[-1]))
    test_filenames, _thr = _ev.load_test_files_from_config()
    X_test, y_test, _ids = _ev.load_test_data(y_all, chunk_files, test_filenames)
    X_test = X_test.tocsr()
    X_csc = X_test.tocsc()
    y_test = np.asarray(y_test).astype(int)
    n_test = X_test.shape[0]

    baseline_pred = model.inplace_predict(X_test)
    baseline_auc = roc_auc_score(y_test, baseline_pred)
    print(f"  test genomes={n_test}  | baseline ROC-AUC={baseline_auc:.4f}")

    cand, cand_path = load_candidates(results_root, antibiotic)
    print(f"  candidates: {len(cand)} (from {cand_path.name})\n")

    rows = []
    for _, c in cand.iterrows():
        fidx = int(c["feature_index"])
        kmer = str(c.get("kmer", ""))
        rec = {"rank": c.get("rank"), "kmer": kmer, "feature_index": fidx,
               "card_gene": c.get("card_gene"), "stable": c.get("stable"),
               "used_by_model": int(fidx in used_cols)}

        col = np.asarray(X_csc.getcol(fidx).todense()).ravel()
        ones = np.where(col != 0)[0]
        c_ones = ones.size
        # Unused feature, or no variation in the test set → MDA undefined/0.
        if fidx not in used_cols or c_ones == 0 or c_ones == n_test:
            rec.update(mda_auc_drop=0.0, perm_p=1.0, n_perm=0)
            rows.append(rec); continue

        # Only column `fidx` varies between permutations; every other feature is
        # fixed. So a genome's prediction depends solely on its feature-`fidx`
        # value — precompute BOTH states once (pred1 = has the feature, pred0 =
        # does not), then each permutation is a per-row pick. 2 model calls per
        # candidate instead of R, and exact (not an approximation).
        # Set column `fidx` with pure sparse arithmetic (no LIL — tolil on a
        # ~n_test x 4.9M matrix is Python-slow and hit the node CPU ulimit).
        # Binary/missing semantics preserved: present = stored 1; absent =
        # *unstored* (matches training, where absence is an implicit/missing 0) —
        # so force-absent removes the entry (eliminate_zeros), not a stored 0.
        rows0 = np.where(col == 0)[0]   # currently absent -> force present for pred1
        rows1 = np.where(col != 0)[0]   # currently present -> force absent for pred0
        pred1 = baseline_pred.copy()
        if rows0.size:
            sub = X_test[rows0]
            e = csr_matrix((np.ones(rows0.size),
                            (np.arange(rows0.size), np.full(rows0.size, fidx))),
                           shape=sub.shape)
            pred1[rows0] = model.inplace_predict(sub.maximum(e).tocsr())
        pred0 = baseline_pred.copy()
        if rows1.size:
            sub = X_test[rows1]
            e = csr_matrix((np.ones(rows1.size),
                            (np.arange(rows1.size), np.full(rows1.size, fidx))),
                           shape=sub.shape)
            m = sub - e
            m.eliminate_zeros()
            pred0[rows1] = model.inplace_predict(m.tocsr())

        perm_aucs = np.empty(R)
        for r in range(R):
            assign = np.zeros(n_test, dtype=bool)
            assign[rng.choice(n_test, size=c_ones, replace=False)] = True
            perm_aucs[r] = roc_auc_score(y_test, np.where(assign, pred1, pred0))

        mda = baseline_auc - float(perm_aucs.mean())
        p = (1 + int(np.sum(perm_aucs >= baseline_auc))) / (R + 1)
        rec.update(mda_auc_drop=mda, perm_p=p, n_perm=R)
        rows.append(rec)
        print(f"  f{fidx:<8} {kmer[:24]:<24} MDA={mda:+.4f}  p={p:.3f}  "
              f"({c.get('card_gene') or '—'})")

    res = pd.DataFrame(rows)
    res["perm_q"] = benjamini_hochberg(res["perm_p"].where(res["n_perm"] > 0).values)
    res["permutation_significant"] = (res["perm_q"] < 0.05).fillna(False).astype(int)
    res = res.sort_values(["permutation_significant", "mda_auc_drop"],
                          ascending=[False, False])

    csv_path = out_dir / f"12_permutation_test_{antibiotic}.csv"
    res.to_csv(csv_path, index=False, encoding="utf-8")

    n_tested = int((res["n_perm"] > 0).sum())
    n_sig = int(res["permutation_significant"].sum())
    summary = {
        "antibiotic": antibiotic, "organism": organism,
        "n_permutations": R, "seed": args.seed,
        "baseline_roc_auc": baseline_auc,
        "n_candidates": int(len(res)),
        "n_tested": n_tested,                      # used-by-model & variable
        "n_significant_fdr": n_sig,                # BH-FDR Q<0.05
        "model": model_path.name,
    }
    json_path = out_dir / f"12_permutation_summary_{antibiotic}.json"
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print("\n" + "=" * 78)
    print(f"  baseline AUC {baseline_auc:.4f} | tested {n_tested}/{len(res)} "
          f"candidates | significant (Q<0.05): {n_sig}")
    print(f"  ✓ {csv_path.name}\n  ✓ {json_path.name}")
    print("=" * 78)


if __name__ == "__main__":
    main()
