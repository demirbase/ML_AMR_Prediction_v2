# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project follows
[Semantic Versioning](https://semver.org/).

## [Unreleased]

### Added
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
