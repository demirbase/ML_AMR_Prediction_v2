# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Fixed
- **`kb_overview.cpss_n_stable` counted more than CPSS.** `kb_tables.py` counted
  every `unitig_model_scores` row with `stable=1`, without filtering
  `selection_method` and without deduplicating, so the column — and the axis of
  figure 02, labelled "CPSS stable unitigs (π≥0.6)" — reported **2,060** rows
  over 2,045 distinct pairs, of which **856 came from the `gain_seed` path**.
  The CPSS figure is **1,204**, which is what all 45 `13_stability_summary`
  JSONs and the `cpss` row of `evidence_accounting.csv` report. Now filtered to
  `selection_method='cpss'` and deduplicated; `kb_overview.csv` and figure 02
  regenerated, and the column agrees with the per-model summaries in 45/45.
- **`limitations.csv` row 9 carried a pre-rebuild claim** that the `snp` layer
  grades 0 biomarkers. It now reports the graded count from the KB (18).
- **`snp` evidence layer reached no biomarker — a load-time defect, not a
  set-membership one.** Step 11 reports its hits by the FASTA header it queried
  (`Rank_n|Score_x|Feature_f...`), not by sequence.
  `populate_database.populate_snp` fell back to that header column whenever no
  `kmer` column was present and handed it to `unitig_id()`, which registered the
  identifier string itself as a unitig. The layer therefore joined to nothing
  while appearing to have run, and 335 non-DNA rows accumulated in `unitigs`.
  The previous diagnosis in METHODOLOGY/HANDOFF — that step 11's candidate set
  and the graded universe are disjoint — was **wrong**; the two sets overlap.
  `_attach_snp_sequences()` now resolves the header back to its k-mer through
  the same FASTA, and `populate_snp` refuses any value that is not DNA.
  Step 11 was **not** re-run: its outputs were correct throughout.
  Effect on the delivered KB: `snp` fires on **18** of 3571 (unitig, model)
  pairs · all 21 `resistant_allele` calls now reach the graded universe ·
  5 biomarkers move `weak` → `candidate` (`candidate` 942→947, `weak`
  1920→1915; `confirmed`, `strong_novel` and `none` unchanged) · `unitigs`
  3844→**3509** · the accounting becomes **7 produced / 6 counted / 5 firing**,
  leaving `mda` as the only dead layer (a genuine underpowered null).
  The 18 are the canonical target-site substitutions a homolog-model BLAST
  cannot distinguish: *gyrA* S83L (*E. coli*), *gyrA* T83I (*P. aeruginosa*),
  *gyrA* S84L and *parC* S80F (*S. aureus*), *parC* S84L (*A. baumannii*).
  Regression tests added in `tests/test_evidence_tier.py`.
  ⚠️ The KB now differs from the archived v0.7.1 and requires a new Zenodo
  version; `amrk.db`, `results/tables/` and figures 06/35/39 were regenerated.
  `models.n_trees`/`n_features` were carried over from the pre-rebuild database
  because the production `manifest.json` files are not available locally.
  `kb_metadata.zenodo_doi` is deliberately left NULL: the previous value
  (`10.5281/zenodo.21789464`) identifies an archive that no longer holds this
  content, so it must be repopulated when the new version is minted.

### Added
- **Step 18 — genomic context for `strong_novel` biomarkers.**
  `scripts/18_novel_ncbi_context.py` joins the KB's novel set to the
  organism-restricted NCBI `nt` alignments step 08 already produced but never
  loaded (`populate_database.py` writes `source_db='card'` only). All **23/23**
  novel biomarkers align at 100% identity over full query length, so none is an
  assembly artefact. Replicon call over *all* retained alignments (≥80% majority,
  not the single best hit): 10 chromosomal, 5 plasmid, 8 mixed — the mixed class
  being a mobile-element signature. Read-only: the KB and its DOI are untouched.

### Fixed (documentation — no code/KB change)
- **`METHODOLOGY.md` reconciled with the delivered run.** Removed the stale
  "Nextflow BLAST pipeline" description (step 08 has been pure-Python
  `subprocess` since the M9 review); corrected the H3 result from *negative* to
  the delivered positive contrast; recorded that the canonical unitig path used a
  **fixed `min_support = 10`** on all 45 runs rather than the adaptive formula in
  §2.4; separated CARD (enters the KB) from NCBI (context only); bumped schema
  0.4.0 → 0.7.1. Added **§5 "What the delivered run actually did"** — panel and
  parameters as executed, the evidence ladder as executed, and the limitations
  that must be stated.
- **Evidence-layer count pinned.** The project produces **7** orthogonal
  analyses, `classify_evidence_tier()` counts **6**, and **4** actually fired in
  the delivered KB (max observed `n_evidence_layers` = 4, reached by 7
  biomarkers). Two designed layers contributed nothing, for different reasons:
  `snp` **0/3571** because step 11 scans the step-07 FASTA while the KB's scored
  universe comes from steps 10+13 — the two sets are disjoint (953 rows / 335
  unitigs, zero overlap), silently stranding 21 `resistant_allele` calls; `mda`
  **0/3571** as a genuine null (sets overlap 2409/2409, nothing survives BH-FDR
  at q<0.05, with R=100 permutations underpowered at ~2400 candidates).
  Corrected in `METHODOLOGY.md`, `docs/KB_ACIKLAMA.md`, `docs/KB_KAVRAMLAR.md`.
- **`docs/KB_ACIKLAMA.md` / `docs/KB_KAVRAMLAR.md` brought to 0.7.1.** Both still
  described the 21-model / 2-organism / 2363-unitig 0.6.0 KB. All row counts,
  the CARD version string, the organism and drug-class lists, and the schema
  version were corrected against the shipped database. The advisor-question
  section quoted M13 concordance results (bACC 0.926 vs AMRFinderPlus 0.538) as
  findings — `external_concordance` is **empty** in 0.7.1, so those are now
  marked as superseded 0.6.0-era numbers that must not be used in the thesis.
- **CARD tier filter documented.** 3007 of 3611 `blast_annotations` rows sit at
  `tier='none'` (mean coverage 0.38, E-values to 9.3), including all 2035 rows
  whose `gene_symbol` is the literal `"nan"`. Any biological claim must filter
  `tier IN ('confirmed','candidate')`; unfiltered joins return, for example,
  staphylococcal *mecA* under an *A. baumannii* model.

---

## [0.7.1] - 2026-08-04

Archived on Zenodo: **[10.5281/zenodo.21789464](https://doi.org/10.5281/zenodo.21789464)**
(Dataset, CC-BY-4.0). 45 models · 6 ESKAPEE organisms · 14 antibiotic classes ·
78,556 genome–phenotype pairs · lineage-aware CV on 45/45 models (mean ROC-AUC
0.842) · 3571 tiered biomarkers.

> **How to read the entries below.** They accumulated across the project's whole
> development, from the 1-antibiotic pilot through the 2- and 3-antibiotic KBs to
> the final 45-model rebuild, and record *when each capability landed* — not that
> every one of them contributed an artefact to the shipped database. Three
> capabilities listed here produced results in earlier KB builds that are **absent
> from the 0.7.1 database**: **M13 external concordance** (`external_concordance`
> ships with 0 rows — the AMRFinderPlus/ResFinder numbers quoted below are from
> the 2-organism build and must not be cited as 0.7.1 results), **step 15
> cross-antibiotic overlap** (`unitig_antibiotic_overlap` ships with 0 rows), and
> **step 11's SNP allele check** (953 rows land in `variant_snp_check` but share no
> unitig with any scored model, so the `snp` evidence layer is never applied). See
> `METHODOLOGY.md §5.2`.

### Added
- **M15 genome QC executed (CheckM2 + QUAST).** `02d_genome_qc.py` run on all
  **17,742** assemblies of the final 6-organism panel: **98.7% pass**
  (17,516/17,742). Per-genome table + summary JSON + advisory exclusion list.
  Fails <2% → not retrained (data-quality statement for Methods).
  **The enforced gate is completeness≥95 / contamination≤5 only.** N50≥50 kb and
  contigs≤500 are computed and reported but deliberately NOT enforced: an N50 gate
  removes **1,305 of 2,078 E. faecium genomes (63%)**, i.e. it selects on assembly
  provenance rather than genome quality for a species whose BV-BRC entries are
  routinely short-contig drafts. Figure 12 now plots all four criteria, marking
  which two are enforced. (The earlier "97.1% pass (5312/5470), 158 fails mostly
  low N50" line described the 2-organism era run and no longer holds.) Also: `kb_report.py` one-command thesis results summary;
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
  measured confounding, provenance-over-assertion). **Note:** the H3 example in
  §4.4 was written against the 2-organism run and read as a negative finding; the
  45-model re-analysis supersedes it (see below).
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
  level. **H3 rejected at this scale**: within-β-lactam (ampicillin=TEM vs
  cefotaxime=CTX-M/CMY) shares no gene family — same class, distinct enzymes.
  **SUPERSEDED by the 45-model analysis** (`scripts/17_h3_gene_family_overlap.py`):
  over 138 pairs, same-class pairs share more ARO gene families than cross-class
  pairs (mean overlap 0.84 vs 0.29, Mann–Whitney p=0.0015, `H3_supported: true`).
  The claim rests on that contrast — 0/138 individual pairs survive
  Benjamini–Yekutieli, and only 5 within-class pairs exist by panel construction.
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
