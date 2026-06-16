# Technical Review & Remediation Record

**Project:** ML AMR Prediction Framework v2
**Scope:** Consolidated record of the code review / audit findings and their
resolution. This single document supersedes and merges three earlier working
documents (now archived, see end of file):

- `technical_review_report.md` — static analysis (findings B01–B22)
- `AUDIT_ISSUES.md` — second-pass audit (CRITICAL/HIGH/MEDIUM/LOW)
- `AMR_Project_Analysis.md` — scientific + reproducibility analysis (P-01–P-15)

All concrete code/config/documentation findings below are **resolved** on the
`fix/amr-audit-remediation` branch. The scientific gaps and large MLOps items
are tracked as deliberate future work (see `docs/SCALE_MLOPS_PLAN.md` and
`docs/ROADMAP.md`).

---

## 1. Summary

The pipeline (alignment-free k-mer features → out-of-core XGBoost → BLAST
annotation) was architecturally sound but carried a number of correctness,
security, reproducibility and consistency defects. The remediation:

- fixed **1 data-leakage bug**, **1 security bug**, and ~30 correctness/quality
  defects across all 12 scripts;
- removed duplicated code into a shared `scripts/lib/` package + registry;
- introduced an organism-aware path layer, run metadata/manifests, and a
  layered test suite (smoke / unit / synthetic end-to-end);
- the end-to-end test then surfaced and fixed 3 further runtime robustness bugs
  (`base_score`, `n_estimators` off-by-one, empty feature importances).

---

## 2. Critical findings (resolved)

| ID | Finding | File | Resolution |
|----|---------|------|------------|
| P-01 / B02 | **Threshold data leakage** — Youden's J fit on the *test* set and written back to config | `06_evaluation.py` | Removed the test-set fit+write; the unbiased, training-derived threshold is applied; Youden's J is logged for information only. A test asserts 06 never mutates the threshold. |
| P-02 / B01 | **Shell injection** — `subprocess.run(..., shell=True)` | `03_matrix_construction.py` | Replaced with the shared `run_command()` (`shlex.split`, `shell=False`). |
| P-03 / B03 | `eval_metric` hardcoded `aucpr`, overriding config & mislabelling the score | `04_optimization.py` | Read from `config['xgboost_params']['eval_metric']`. |
| P-04 / B06 | `colsample_bytree` search space `[0.05, 0.3]` ~100× off the stated `1/√p` heuristic | `04_optimization.py`, `METHODOLOGY.md` | Search space derived dynamically as a log-scale window around `1/√p` (`compute_colsample_range`). |
| AUDIT | **R/S double-count** — `counts.get(1.0)+counts.get(1)` on a float index | `01_data_validation.py` | Use the float keys only. |
| AUDIT | Hardcoded `top_50` filenames break when `top_n_features ≠ 50` | `08`, `09` | Use `TOP_N` from config. |

## 3. Runtime bugs found by the end-to-end test (resolved)

| Finding | File | Resolution |
|---------|------|------------|
| Config/data path mismatch after the organism migration | `config.yaml` + all scripts | Scripts resolve paths via `lib.config.resolve_path(organism, antibiotic)`. |
| `base_score must be in (0,1)` on a pure/weighted chunk (XGBoost ≥2.0 auto base_score) | `05_model_training.py` | Pin `base_score=0.5`. |
| `n_estimators=0` — `best_iteration` (0-indexed) stored without `+1` → 0-tree model | `04`, `05` | Store `max(1, best_iteration+1)`; `05` guards with `max(1, …)`. |
| Empty feature importances → `list index out of range` | `07_explainability.py` | Detect empty importance, write empty outputs, exit cleanly. |

## 4. High / medium findings (resolved)

- **O(N²) Gram matrix** in SVD QC → `TruncatedSVD` on the sparse stack (`03b`).
- **BLAST confidence tiers** — replaced the over-permissive `evalue ≤ 50` with
  confirmed/candidate/weak tiers (E-value + identity) (`09`).
- **NCBI Entrez** — email/api_key read from config (`ncbi:` section); no
  hardcoded placeholder e-mail (`09`).
- **PR-AUC** unified on `average_precision_score` (`06`).
- **KMC resume** checks both `.kmc_pre` and `.kmc_suf` (`02`).
- **NCBI Nextflow** process gets `errorStrategy 'ignore'` (`08_blast_pipeline.nf`).
- **CARD DB version** recorded via `blastdbcmd -info` (`08`).
- Duplicated `get_y_chunk`/`run_command`/`ANTIBIOTIC_CLASSES` → `scripts/lib/`
  + `config/registry/`.
- Numerous quality fixes: `mkdir` dup, `scale_pos_weight` loop pop, global-RNG
  reseed, f-string header, docstrings (`k=21`), unused imports, `print("\n="*60)`,
  feature-index parsing, step-label numbering, config comment, etc.

## 5. Reproducibility / robustness additions

- **Bootstrap 95% CIs** for ROC-AUC / PR-AUC (`06`) — interim variance estimate.
- **Run metadata** (`runs/.../run_metadata.json`), **model manifest**
  (`models/.../manifest.json`), **metrics** (`metrics.json`).
- **Layered test suite** (`tests/`): smoke + unit (seconds) and an opt-in
  synthetic end-to-end run (minutes) — see `tests/README.md`.
- `environment.yml` (conda incl. KMC/BLAST/Nextflow) + comprehensive
  `requirements.txt`.

---

## 6. Open scientific gaps (future work — not code bugs)

These are deferred by design (see `docs/ROADMAP.md`, `docs/SCALE_MLOPS_PLAN.md`):

- **Cross-validation** — currently a single stratified split + bootstrap CIs;
  full k-fold / repeated holdout is future work (P-05).
- **Feature stability** (multi-seed selection frequency) — `07b` not yet
  implemented (P-06).
- **Cross-antibiotic overlap** analysis — `10` not yet implemented.
- **Phylogenetic / MLST bias control**, external/temporal validation.
- **Knowledge Base** (SQLite→Postgres + API) — not yet implemented.

## 7. Known data gaps (operational, not code)

- `features.txt` (the k-mer vocabulary) is gitignored (>100 MB) and may need
  regeneration via step 03 before re-running 07–09 (P-08).
- Only **cefotaxime** and **gentamicin** matrices are present on disk;
  **ampicillin** and **ciprofloxacin** matrices must be regenerated (02→03).

---

*Archived source documents (kept locally under `archive/`, excluded from git):
`technical_review_report.md`, `AUDIT_ISSUES.md`, `AMR_Project_Analysis.md`.*
