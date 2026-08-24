# ============================================================================
# Makefile — ML_AMR_Prediction_v2 developer & pipeline shortcuts
# ============================================================================
# `make help` lists targets. PYTHON can be overridden, e.g.:
#     make test PYTHON=/opt/anaconda3/envs/bitirme_vol2/bin/python
# ============================================================================
PYTHON ?= python
ORG    ?= ecoli
AB     ?= ampicillin

.DEFAULT_GOAL := help
.PHONY: help setup dev-install lint format typecheck test test-all \
        pipeline data features train biology tables figures clean-pyc

# Thesis artefact paths. Override on the command line if the KB lives elsewhere,
# e.g. `make figures KB=$$AMR_WORK/results/kb/amrk.db`.
KB      ?= results/kb/amrk.db
RESULTS ?= results
TABLES  ?= results/tables
FIGURES ?= results/figures
# kb_tables_thesis reads these two as well: the PopPUNK cluster CSVs for the lineage
# summary, and run_metadata.json for the per-model hyperparameters. Both tables are
# skipped (never written empty) when the inputs are absent.
DATA    ?= data/processed
RUNS    ?= runs

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
	  awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-14s\033[0m %s\n", $$1, $$2}'

setup:  ## Create the conda environment (KMC/BLAST/Nextflow + Python deps)
	conda env create -f environment.yml

dev-install:  ## Install the dev/QA toolchain (ruff, mypy, pre-commit, pytest)
	$(PYTHON) -m pip install -e ".[dev]"
	pre-commit install

lint:  ## Lint with ruff
	ruff check .

format:  ## Auto-format with ruff
	ruff format .

typecheck:  ## Type-check the lib/ core with mypy
	mypy

test:  ## Run unit + smoke tests (fast)
	$(PYTHON) -m pytest -ra

test-all:  ## Run everything incl. the synthetic integration test (needs KMC/xgboost)
	$(PYTHON) -m pytest -m "unit or smoke or integration" -ra

pipeline:  ## Run the analysis core (01->10) for ORG/AB from config
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --antibiotic $(AB)

data:  ## Acquire BV-BRC data + build the phenotype matrix (00a, 00)
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --only 00a 00

features:  ## k-mer counting + matrix construction (02, 02b, 03)
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --only 02 02b 03

train:  ## HPO + train + evaluate (04, 05, 06)
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --only 04 05 06

biology:  ## Stability + explainability + BLAST + reports (07b..11)
	$(PYTHON) scripts/run_pipeline.py --organism $(ORG) --only 07b 07 08 09 10 11

tables:  ## Rebuild the tidy thesis tables from the KB (kb_tables + H3 + CV comparison)
	$(PYTHON) scripts/kb_tables.py --db $(KB) --results $(RESULTS) --out $(TABLES)
	$(PYTHON) scripts/17_h3_gene_family_overlap.py --db $(KB) --tables $(TABLES) --figures $(FIGURES)
	$(PYTHON) scripts/kb_cv_comparison.py --tables $(TABLES) --results $(RESULTS) --out $(FIGURES)
	$(PYTHON) scripts/18_novel_ncbi_context.py --kb $(KB) --results-root $(RESULTS) --out $(TABLES)
	$(PYTHON) scripts/kb_fair_mapping.py --db $(KB) --out $(TABLES)
	$(PYTHON) scripts/kb_tables_thesis.py --db $(KB) --tables $(TABLES) --data $(DATA) --runs $(RUNS)

figures: tables  ## Rebuild every thesis figure (run `make tables` implicitly)
	# Each script takes a DIFFERENT set of flags -- kb_figures_model has no --db, and
	# passing one makes argparse exit non-zero while a loop that only greps stdout for
	# a checkmark reports success. That is how a corrected figure silently kept its old
	# rendering for a whole review round. Spell the flags out, once, here.
	$(PYTHON) scripts/kb_figures.py         --tables $(TABLES) --results $(RESULTS) --db $(KB) --out $(FIGURES)
	$(PYTHON) scripts/kb_figures_data.py    --tables $(TABLES) --results $(RESULTS) --db $(KB) --out $(FIGURES)
	$(PYTHON) scripts/kb_figures_model.py   --tables $(TABLES) --results $(RESULTS) --out $(FIGURES)
	$(PYTHON) scripts/kb_figures_biology.py --tables $(TABLES) --db $(KB) --out $(FIGURES)
	$(PYTHON) scripts/kb_figures_schematic.py --db $(KB) --tables $(TABLES) --out $(FIGURES)

clean-pyc:  ## Remove Python caches
	find . -type d -name __pycache__ -prune -exec rm -rf {} + ; \
	find . -type f -name '*.py[co]' -delete
