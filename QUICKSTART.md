# ML AMR Prediction Framework — Quickstart

Spin up the `ML_AMR_Prediction_v2` pipeline: validate → extract → construct →
train → evaluate → biologically annotate, for a chosen organism + antibiotic.

## 1. Prerequisites & setup

### Environment

- **Python** 3.10+
- **conda (recommended):** `conda env create -f environment.yml && conda activate amr-prediction`
  (installs KMC, BLAST+ and Nextflow too), **or** `pip install -r requirements.txt`
  and install the external tools yourself.

### External tools

- [ ] **KMC** ≥3.2 — k-mer counting (steps 02/03). Binaries expected under `bin/bin/`.
- [ ] **BLAST+** ≥2.12 — homology search (step 08). `brew install blast` / `conda install -c bioconda blast`.
- [ ] **Nextflow** ≥22.10 — runs `08_blast_pipeline.nf`. `curl -s https://get.nextflow.io | bash`.

### Inputs (organism-scoped layout)

For organism slug `ecoli` (set in `config/config.yaml → project.organism`):

- [ ] **Raw genomes:** `.fna` files in `data/raw/ecoli/genomes/`
- [ ] **Metadata:** phenotype labels in `data/external/ecoli/metadata/amr_phenotypes.csv`
      (`Genome ID` column + one column per antibiotic, values 0/1)
- [ ] **CARD DB** (step 08): `data/external/blast_db/card_nt/card.*`
- [ ] **NCBI e-mail** (step 09): set `config/config.yaml → ncbi.entrez_email`

> Adding a new organism: add a block to `config/registry/organisms.yaml`
> (`enabled: true`), drop its data in `data/raw/{organism}/` and
> `data/external/{organism}/`, then run the pipeline — no code changes.

---

## 2. Configure the target

Edit `config/config.yaml`:

```yaml
project:
  organism: "ecoli"            # registry slug (config/registry/organisms.yaml)
  target_antibiotic: "gentamicin"
```

Key parameters also live here: `preprocessing` (k_length, min_support,
chunk_size), `training` (n_trials, test/validation fractions), `analysis`
(top_n_features), `blast`, `ncbi`.

---

## 3. Run the pipeline

```bash
python scripts/01_data_validation.py     # metadata validation + EDA
python scripts/02_kmer_extraction.py     # KMC k-mer counting
python scripts/03_matrix_construction.py # sparse .npz matrix chunks
python scripts/04_optimization.py        # Optuna HPO -> config_{antibiotic}.yaml
python scripts/05_model_training.py      # out-of-core XGBoost -> model + manifest
python scripts/06_evaluation.py          # metrics, ROC/PR, bootstrap CIs
python scripts/07_explainability.py      # top k-mers -> CSV + FASTA
python scripts/08_blast_annotation.py    # BLAST vs CARD + NCBI (Nextflow)
python scripts/09_biological_summary.py  # confidence-tiered biological report
```

---

## 4. Expected outputs

Per-run provenance lands in `runs/{organism}/{antibiotic}/{run_id}/`
(`run_metadata.json`, `metrics.json`); the trained model + `manifest.json` in
`models/{organism}/{antibiotic}/`.

Analysis artifacts are routed to `results/{organism}/{antibiotic}/`:

- `01_data_exploration/` — class distribution, missingness, co-occurrence
- `02_matrix_qc/` — sparsity, prevalence, SVD separability
- `03_model_optimization/` — Optuna history & importance
- `04_evaluation/` — confusion matrix, ROC/PR, calibration, metrics CSV, bootstrap CIs
- `05_explainability/` — top-feature CSV/FASTA + CARD/NCBI BLAST TSVs + final report

Each plot is saved with its underlying `.csv` for figure reproduction.

> `results/`, `logs/` and `runs/` are generated (not version-controlled).

---

## 5. Validate without a full run

```bash
pytest                 # smoke + unit (seconds)
pytest -m integration  # tiny synthetic end-to-end (minutes); needs xgboost/KMC
```

See `tests/README.md` for details.
