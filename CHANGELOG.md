# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
- **M15 genome QC executed (CheckM2 + QUAST).** `02d_genome_qc.py` run on all 5470
  assemblies: **97.1% pass** (5312/5470) at completeness≥95 / contamination≤5 /
  N50≥50kb / contigs≤500; 158 fails (mostly low N50). Per-genome table + summary
  JSON + advisory exclusion list. Fails <3% → not retrained (data-quality
  statement for Methods). Also: `kb_report.py` one-command thesis results summary;
  `kb_app.py` H3-overlap + M13-concordance tabs; S. aureus taxid fix (927→1280).
- **KB REST API (S8) + FAIR metadata endpoint (S9).** `scripts/kb_api.py` (FastAPI,
  CORS, auto OpenAPI at `/docs`) over `scripts/lib/kb_queries.py` (pure sqlite3,
  unit-tested without a web server): `/api/v1/kmers` (filter by antibiotic / tier /
  min_stability / stable_only), `/kmers/{sequence}` (full evidence chain),
  `/overlap`, `/stats`, and `/metadata` (schema version + Zenodo DOI + CC-BY-4.0
  license + counts — the FAIR machine-readable access point). Run with
  `uvicorn scripts.kb_api:app`.
- **Reification safeguard (S10).** `METHODOLOGY.md §4.4` — associational-not-causal
  wording policy + three structural safeguards (layered orthogonal evidence,
  measured confounding incl. the H3 negative finding, provenance-over-assertion).
- **M13 external-validation concordance.** `scripts/16_external_concordance.py`
  + `scripts/lib/concordance.py` (balanced accuracy, sensitivity, specificity,
  Cohen's κ, McNemar, FDA major/very-major error bands). AMRFinderPlus
  2026-05-15.1 + ResFinder 4.5.0 run on all 5468 genomes; per-antibiotic
  genotype-vs-phenotype concordance **and a leakage-free model-vs-tool
  head-to-head** on the held-out test split (`06_evaluation.py` now saves
  `16_model_preds_{ab}.csv`). The unitig model's balanced accuracy (amp 0.873 /
  cef 0.925 / cip 0.928) matches ResFinder for cefotaxime & ciprofloxacin and
  beats AMRFinderPlus for ciprofloxacin & ampicillin — external evidence the
  signal is mechanism-driven, not lineage memorisation.
- **Third antibiotic (cefotaxime) + H3 result.** Canonical 03u→populate for
  cefotaxime (model_id 3, lineage-CV 0.9546±0.020, CTX-M-276/278 confirmed).
  `scripts/15_cross_antibiotic.py` (S1) computes cross-antibiotic stable-unitig
  overlap + the H3 within/cross contrast at both unitig and ARO gene-family
  level. **H3 rejected**: within-β-lactam (ampicillin=TEM vs cefotaxime=CTX-M/CMY)
  shares no gene family — same class, distinct enzymes — a biologically
  substantive negative finding.
- **M10 FAIR/Zenodo prep.** `populate_database.py` preserves `zenodo_doi` across
  re-populates and honours the `AMR_ZENODO_DOI` env override; added `.zenodo.json`,
  refreshed `CITATION.cff`, a README "Data availability & KB versioning" section,
  and `docs/RELEASE_ZENODO.md` deposit checklist. (Deposit itself pending.)
- **Second antibiotic + multi-antibiotic KB (ciprofloxacin).** Canonical 04→populate
  run for ciprofloxacin appended to the same `amrk.db` (model_id 2). lineage-CV
  0.9496±0.007. **SNP showcase:** step 11 detects `gyrA S83L` + `parC S80I` as
  resistant alleles — the fluoroquinolone positive control — while CARD-homolog
  recovery is ~0 (resistance is a target-gene SNP, not an acquired gene), the exact
  complement of ampicillin's acquired-gene mechanism.
- **`13_stability_selection.py --base-trees`** — decouple the CPSS base selector
  (sparse, default 10 trees) from the final model so the PFER bound stays tight
  (a 66/146-tree model over-selects → PFER ~100; 10 trees → ~2.7).
- **`scripts/kb_app.py`** — local Streamlit explorer over `amrk.db`: filterable
  biomarker table, per-unitig multi-layer evidence chain, model + provenance
  (ROADMAP S8/N1). `pip install streamlit pandas`.
- **Canonical reproducible re-run** of the whole pipeline from step 04 (HPO
  included) so every artefact carries provenance: `pipeline_runs` now stamps
  git_commit + random_seed + config_hash(sha256) + CARD version (run_id is no
  longer `…__unknown`).

### Changed (KB schema 0.3.0 → 0.4.0 — renamed k-mer → unitig)
- Tables/columns renamed to the actual feature unit: `kmers`→`unitigs`,
  `kmer_id`→`unitig_id`, `kmer_model_scores`→`unitig_model_scores`,
  `kmer_background_frequency`→`unitig_background_frequency`,
  `kmer_antibiotic_overlap`→`unitig_antibiotic_overlap`, `kb_metadata.n_kmers`→
  `n_unitigs`. Candidate-CSV column reads stay `kmer` (on-disk name); output
  filenames unchanged.

### Fixed
- `populate_run` crashed once `run_metadata.json` actually existed: wrong keys
  (`git_commit`/`seed`/`created_at`) and a dict (`data_fingerprint`) bound to a
  TEXT column. Aligned to `git_commit_hash`/`random_seed`/`started_at` and store
  `data_fingerprint.sha256` in `config_hash`.
- `12b_label_permutation_test.py` test/label misalignment when the experiment
  config's `test_files` are non-ascending (REAL AUC collapsed to ~0.49) — build
  the test matrix in ascending chunk order.

- **M9 permutation significance tests.** `12_permutation_test.py` — MDA
  (per-feature permutation importance): model fixed, permute each candidate
  unitig's column in the held-out test set, measure ROC-AUC drop + BH-FDR
  (ROADMAP §0.2). `12b_label_permutation_test.py` — label-permutation null
  (ROADMAP §1.7): shuffle all labels, retrain frozen-HP 8 trees over a built-once
  streamed `QuantileDMatrix` (swap labels via `set_label`/`set_weight`), build the
  null ROC-AUC distribution + empirical p. Ampicillin: baseline AUC 0.9534; MDA
  0/51 significant under Q<0.05 (expected — unitig redundancy); label-perm real
  ≫ null (~0.50) → model highly significant.
- **M16 ARO/CARD ontology in the KB.** `blast_annotations` gains
  `aro_accession`, `aro_gene_family`, `aro_drug_class`, `aro_resistance_mechanism`
  (kb_schema `0.1.0`→`0.2.0`); `populate_database.py` writes them for CARD hits.
- **M6 CARD version recorded** — `AMR_CARD_VERSION` env override →
  `kb_metadata.card_version` (e.g. `4.0.1`), no config.yaml edit on HPC.
- `AMR_ENTREZ_EMAIL` / `AMR_ENTREZ_API_KEY` env overrides for step 09 Entrez.

### Changed
- **Step 08 NCBI remote BLAST decoupled from CARD.** The public NCBI server
  SIGXCPU-kills `blastn-short`/word7 over `nt`; the remote pass now uses
  `blastn`/word11 + a taxid `-entrez_query` (`txid<N>[Organism:exp]`, from the
  registry — a scientific-name value breaks the Nextflow CLI launcher) +
  `-max_target_seqs`; CARD local pass unchanged. Step 08 also sets
  `NXF_ANSI_LOG=false` so backgrounded Nextflow doesn't stall (SIGTTOU).

### KB schema
- `KB_SCHEMA_VERSION` `0.1.0` → `0.2.0` (added ARO columns to `blast_annotations`).

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
