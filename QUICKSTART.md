# ML AMR Prediction Framework — Quickstart

End-to-end guide: **install → acquire data → run the pipeline → read the
results**, for a chosen organism + antibiotic. Works on macOS, Linux and HPC
clusters. Commands are copy-pasteable.

> **Pipeline at a glance:** download BV-BRC AMR data → k-mer features →
> out-of-core XGBoost → evaluation → top k-mers → BLAST annotation →
> biological report. Steps are numbered `00 → 09` and run in order.

---

## 0. Install on a fresh machine (laptop or HPC)

The whole stack — Python packages **and** the external tools (KMC, BLAST+,
Nextflow) — is reproducible with conda in one command. This is the recommended
path on an HPC where you pull from GitHub.

```bash
# 1. Get the code
git clone https://github.com/demirbase/ML_AMR_Prediction_v2.git
cd ML_AMR_Prediction_v2

# 2. Create + activate the environment (installs KMC / BLAST+ / Nextflow too)
#    On HPC: load conda first, e.g.  `module load anaconda3`  or  `module load miniconda`
conda env create -f environment.yml
conda activate amr-prediction

# 3. Sanity check (must all resolve)
python -c "import pandas, numpy, scipy, sklearn, xgboost, optuna, Bio, yaml, certifi; print('python deps OK')"
kmc       --help    >/dev/null 2>&1 && echo "kmc OK"
blastn    -version  | head -1
nextflow  -version  | head -2
```

**No conda?** Use pip for the Python side and install the tools yourself:

```bash
python -m venv .venv && source .venv/bin/activate   # Python 3.10–3.12
pip install -r requirements.txt
# then install KMC / BLAST+ / Nextflow (see "External tools" below)
```

### External tools (only if NOT using conda)

| Tool | Used by | Install |
|---|---|---|
| **KMC** ≥3.2 | steps 02, 03 | `conda install -c bioconda kmc` · Linux release binary from [KMC releases](https://github.com/refresh-bio/KMC/releases) · macOS `brew install brewsci/bio/kmc` |
| **BLAST+** ≥2.12 | step 08 | `conda install -c bioconda blast` · `apt-get install ncbi-blast+` · `brew install blast` |
| **Nextflow** ≥22.10 | step 08 | `curl -s https://get.nextflow.io \| bash` (needs Java) |

**Tool discovery is automatic and cross-platform.** Each script finds `kmc` /
`blastn` on your `PATH` (conda/module). The repo also ships a macOS binary under
`bin/bin/` used only as a fallback **on macOS**. To point at a custom build,
export an override, e.g. `export AMR_KMC_BIN=/path/to/kmc`.

> **HPC note:** there is no BV-BRC desktop app on a cluster — that's fine, data
> acquisition (step 00a) uses the HTTP API by default (`--backend api`). See §2.

---

## 1. Configure the target

Edit `config/config.yaml`:

```yaml
project:
  organism: "ecoli"              # registry slug (config/registry/organisms.yaml)
  target_antibiotic: "ampicillin"
```

Other tunables live here too: `preprocessing` (k_length=21, min_support,
`chunk_size=200`), `training` (n_trials, test/validation fractions), `analysis`
(top_n_features), `blast`, and `ncbi.entrez_email` (set this before step 09, or
`export AMR_ENTREZ_EMAIL=…`).

> **Adding an organism:** add a block to `config/registry/organisms.yaml`
> (`enabled: true`) — no code changes. Antibiotic name variants are normalised
> via `config/registry/antibiotics.yaml`.

---

## 2. Acquire the data (steps 00a → 00)

Downloads the BV-BRC AMR table for the organism, cleans it, and fetches the
matching genome assemblies as `{genome_id}.fna`.

```bash
# Quick dry run first (samples a few AMR genomes via the API — finishes in seconds)
python scripts/00a_download_bvbrc.py --organism ecoli --max-genomes 10

# Full acquisition (no cap). --backend api is the portable default (works on HPC).
python scripts/00a_download_bvbrc.py --organism ecoli --backend api

# Build the binary phenotype matrix from the genomes that actually downloaded
python scripts/00_prepare_metadata.py --organism ecoli
```

Backends: **`api`** (HTTP, portable, fast — default) · **`cli`** (BV-BRC `p3-*`
tools, macOS app only) · **`--raw-csv path.csv`** (skip the network; use the
website "DOWNLOAD" of the AMR Phenotypes tab, including the *Testing Standard*
column). Retry only failures with `--retry-failed`.

After this you'll have:
`data/raw/ecoli/genomes/*.fna` and
`data/external/ecoli/metadata/amr_phenotypes.csv`
(`Genome ID` column + one 0/1 column per antibiotic; blank = untested).
Review the cleaning counts in `logs/ecoli/cleaning_report.json`.

> **Already have genomes + a phenotype CSV?** Skip step 00 entirely — just place
> `.fna` files in `data/raw/{organism}/genomes/` and the CSV at
> `data/external/{organism}/metadata/amr_phenotypes.csv`.

---

## 3. Run the pipeline (steps 01 → 09)

Run in order; each reads `config.yaml` for the organism/antibiotic.

```bash
python scripts/01_data_validation.py      # metadata validation + EDA plots
python scripts/02_kmer_extraction.py      # KMC k-mer counting
python scripts/03_matrix_construction.py  # sparse binary .npz matrix chunks
python scripts/04_optimization.py         # Optuna HPO -> config/experiments/{org}/config_{ab}.yaml
python scripts/05_model_training.py       # out-of-core XGBoost -> model + manifest
python scripts/06_evaluation.py           # metrics, ROC/PR, calibration, bootstrap CIs
python scripts/07b_feature_stability.py   # 5-seed stability (AUC ±, Jaccard, stable k-mers)
python scripts/07_explainability.py       # top-N gain k-mers ∪ stable set -> CSV + FASTA
python scripts/08_blast_annotation.py     # BLAST vs CARD + NCBI (Nextflow; needs internet)
python scripts/09_biological_summary.py   # tiered report + recovery rate / composite / novel fraction
python scripts/10_kmer_background_frequency.py  # resistant-vs-susceptible prevalence + discriminativeness (Fisher)
python scripts/11_variant_snp_check.py    # CARD variant-model SNP allele check (optional; needs full CARD download)
```

You can also override the target per-invocation without editing config:
`--organism ecoli --antibiotic ampicillin` (or env vars `AMR_ORGANISM` /
`AMR_ANTIBIOTIC`).

### Running on HPC (SLURM sketch)

Steps 02–05 are the heavy ones (k-mer counting + training). Submit them in a job:

```bash
#!/bin/bash
#SBATCH --job-name=amr-ecoli
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
module load anaconda3            # or miniconda — cluster-specific
conda activate amr-prediction
cd $SLURM_SUBMIT_DIR
python scripts/02_kmer_extraction.py
python scripts/03_matrix_construction.py
python scripts/04_optimization.py
python scripts/05_model_training.py
```

Tune `preprocessing.kmc_mem` (GB) / `preprocessing.threads` and the XGBoost
`n_jobs` in `config.yaml` to match `--mem` / `--cpus-per-task`.
`chunk_size` keeps RAM bounded (one chunk in memory at a time) — leave at 200
unless you hit memory limits.

---

## 4. Where the results land

```
runs/{org}/{ab}/{run_id}/   run_metadata.json (git hash, versions, seed), metrics.json
models/{org}/{ab}/          trained model + manifest.json
results/{org}/{ab}/
  01_data_exploration/      class balance, missingness
  02_matrix_qc/             sparsity, prevalence
  03_model_optimization/    Optuna history & importance
  04_evaluation/            confusion matrix, ROC/PR, calibration, metrics CSV, bootstrap CIs
  05_explainability/        top-feature CSV/FASTA, CARD/NCBI BLAST TSVs, 05_final_biological_report.md
```

Every plot is saved next to the `.csv` it was drawn from. `results/`, `logs/`
and `runs/` are generated and **not** version-controlled.

---

## 5. Verify without a full run

```bash
pytest                  # smoke + unit tests (seconds)
pytest -m integration   # tiny synthetic 02->07b end-to-end (minutes; needs xgboost + KMC)
```

The integration test is the primary real-environment validator — run it first
on a new machine/HPC to confirm the tool chain works. See `tests/README.md`.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `KMC executable not found` | `conda install -c bioconda kmc`, or `export AMR_KMC_BIN=/full/path/kmc` |
| `CERTIFICATE_VERIFY_FAILED` (step 00a) | `certifi` is missing — `pip install certifi` (in `environment.yml`/`requirements.txt`) |
| step 00a returns 0 rows on a dry run | increase `--max-genomes`, or run full (no cap), or use `--raw-csv` |
| `base_score must be in (0,1)` | upgrade — fixed in step 04/05 (pinned `base_score=0.5`) |
| step 04 `division by zero` / 0 train chunks | too few genomes for the split — acquire more, or lower `chunk_size` for tiny test sets |
| step 09 NCBI rate-limit warning | set `config.yaml → ncbi.entrez_email` or `export AMR_ENTREZ_EMAIL=…` |
| step 08 NCBI remote BLAST `SIGXCPU` / empty output | the public NCBI server kills `blastn-short`/word7 over `nt`; the NCBI pass is decoupled to `blastn`/word11 + taxid `-entrez_query` (handled in `08_blast_annotation.py`) |
| step 08 Nextflow stalls under `nohup` (process state `T`) | ANSI console + no tty → SIGTTOU; 08 sets `NXF_ANSI_LOG=false` automatically |
