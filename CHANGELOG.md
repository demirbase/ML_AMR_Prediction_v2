# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- Knowledge-base layer (M8, the thesis contribution): `lib/kb_schema.py` (SQLite
  DDL, 11 tables per ROADMAP §1.1, stdlib only — no new dependency) +
  `populate_database.py` (loads pipeline outputs — run_metadata, manifest, 06
  metrics, 07b holdout, 09/10 candidates+background, 11 SNP — into one queryable
  `results/{org}/kb/amrk.db`; idempotent, multi-antibiotic, graceful on missing
  inputs). `KB_SCHEMA_VERSION=0.1.0`. FastAPI (S8) to follow.
- `lib/xgb_data.py` — `ChunkDMatrixIter` (streaming `xgb.DataIter`) +
  `build_quantile_dmatrix` / `global_pos_weight`: build a single in-core
  `QuantileDMatrix` from on-disk chunks without materialising the full sparse
  matrix (binary data + `max_bin=2` → ~1 byte/non-zero). Supports sample-level
  row masks and global class weighting; shared by steps 05 and 07b.
- Research-software-engineering scaffolding: `LICENSE` (MIT), `CITATION.cff`,
  `pyproject.toml` (PEP 621 metadata + ruff/mypy/pytest config), GitHub Actions
  CI (ruff + unit/smoke tests on Python 3.10–3.12), `.pre-commit-config.yaml`,
  `CONTRIBUTING.md`, this changelog, a `Makefile`, and a `run_pipeline.py`
  orchestrator.
- `lib/logging_utils.py` — standard logger factory for the orchestrator and new code.
- Step 10 `kmer_background_frequency.py` — resistant-vs-susceptible prevalence,
  Fisher's exact test, and a discriminativeness flag (ROADMAP §1.1).
- Step 11 `variant_snp_check.py` — k-mer-centric CARD variant-model SNP allele
  check (resistant allele vs wildtype).
- Step 09: KB-candidate table, composite score, known-mechanism recovery rate
  (M7), novel-candidate fraction (H4); confidence tiers moved to `config.yaml`.

### Changed
- Feature filter (step 03): `min_support` is now **data-adaptive** —
  `max(min_support_floor=5, ceil(min_prevalence=0.01 * n_genomes))` — so it scales
  with dataset size across antibiotics/organisms (small sets fall back to the floor
  and keep all markers; large sets get de-confounding + faster training). An
  explicit integer `preprocessing.min_support` still overrides. (config knobs:
  `min_support`, `min_support_floor`, `min_prevalence`.)
- Training regime (steps 05 and 07b): replaced the epoch-based 1-tree-per-chunk
  incremental warm-start with **standard full-data gradient boosting** over a
  streaming **`ExtMemQuantileDMatrix`** (external memory). Every tree now sees
  the whole training set (stronger fit; saturates HPC cores, fixing the
  low-CPU-efficiency warning). Quantised pages spill to fast scratch
  (`cache_prefix`) so the matrix never has to fit in RAM — an in-core
  `QuantileDMatrix` of the full train set peaked >400 GB and OOM-killed a 384 GB
  node. Class imbalance handled once via a global `neg/pos` instance weight.
  Resolves the documented "04 vs 05 training regimes differ" caveat. On-disk
  chunking and the chunk-level train/test split are unchanged, so no re-run of
  03/04 is required.
- HPO (step 04): runs trials concurrently (`training.optuna_threads_per_trial`)
  over a `QuantileDMatrix` HPO subset, to use all allocated cores without OOM.
- BLAST: CARD search uses `blastn-short -dust no`; confidence tiers grade on
  identity + coverage (database-size-independent), E-value secondary.
- Tool discovery is PATH-aware (`lib.config.resolve_tool`): conda/module on
  Linux/HPC, bundled macOS binary only as a Darwin fallback.
- BV-BRC step 00a: certifi-based SSL, URL-encoded API queries, API-sampled dry
  runs, batched CLI fetch.
- Pipeline order is `07b → 07` so the candidate set includes stable k-mers.

### Fixed
- Out-of-core HPO/training `base_score=0.5` (pure-class chunk error).
- Removed hardcoded macOS KMC paths in steps 02b/03.

### Repository hygiene
- `.gitignore` hardened; generated data, matrices, models, results and the full
  CARD bundle are no longer tracked (only the small CARD homolog BLAST DB is).
