# METHODOLOGY.md — ML_AMR_Prediction_v2

> **A rigorous technical exposition of the biological, mathematical, and statistical foundations of the AMR prediction pipeline.**

> ⚠️ **Canonical pipeline (2026-07, schema 0.4.0) vs this document.** After the
> literature-review pivot, the *canonical* feature unit is
> the **unitig** (compacted de Bruijn graph, `bcalm2`/`unitig-caller`; step 03u),
> the *canonical* validation is **lineage-aware CV** (PopPUNK + StratifiedGroupKFold,
> step 07b), stability is **CPSS** (B=100, π≥0.6, PFER-bounded; step 13) with
> **TreeSHAP** importance, and biomarkers additionally carry pyseer-LMM significance
> (step 14), CARD SNP-allele checks (step 11), genome QC (CheckM2+QUAST; step 02d)
> and external AMRFinderPlus/ResFinder concordance (step 16, M13). **Sections 1–3
> below describe the raw-k-mer / Gain / single-split *baseline* (still a valid,
> runnable path) — read them as the foundation, not the current canonical method;**
> §4.4 is authoritative where they differ. A full rewrite
> to unitig-first is tracked (audit Issue 1).

---

## Table of Contents

1. [Biological Foundations](#1-biological-foundations)
2. [Feature Space Mathematics: The Curse of Dimensionality](#2-feature-space-mathematics-the-curse-of-dimensionality)
3. [Statistical & ML Architecture](#3-statistical--ml-architecture)
   - [3.1 Binary Histogram Quantization (`max_bin = 2`)](#31-binary-histogram-quantization-max_bin--2)
   - [3.2 Optuna HPO and the Square Root Heuristic](#32-optuna-hpo-and-the-square-root-heuristic)
   - [3.3 Stratified Linspace Chunk Selection](#33-stratified-linspace-chunk-selection)
   - [3.4 Full-Data Boosting over a Streaming `QuantileDMatrix`](#34-full-data-boosting-over-a-streaming-quantiledmatrix)
4. [Explainable AI and Biological Validation](#4-explainable-ai-and-biological-validation)
   - [4.1 Feature Importance Mapping](#41-feature-importance-mapping)
   - [4.2 Automated Nextflow BLAST Pipeline](#42-automated-nextflow-blast-pipeline)
   - [4.3 Automated Biological Reporting](#43-automated-biological-reporting)
   - [4.4 Statistical Signal ≠ Biological Causation (reification safeguard)](#44-statistical-signal--biological-causation-reification-safeguard)

---

## 1. Biological Foundations

### 1.1 Whole-Genome Sequencing (WGS)

Whole-Genome Sequencing (WGS) is the process of determining the complete nucleotide sequence of an organism's genome in a single laboratory run. For bacterial samples, modern short-read platforms (e.g., Illumina) produce millions of reads — short DNA fragments of 150–300 base pairs — that are assembled into a draft genome. The result is a FASTA file representing the full genetic blueprint of a bacterial isolate.

In the context of Antimicrobial Resistance (AMR), the key insight is:

> **Resistance is encoded in the genome.** Whether a bacterium survives exposure to an antibiotic is determined by specific mutations, insertions, deletions, or acquired horizontal gene transfer events — all of which are directly observable in the WGS data.

### 1.2 K-mers as Alignment-Free Genomic Features

#### Definition

A **k-mer** is any contiguous subsequence of length $k$ extracted from a DNA string. For a genome of length $L$, the number of k-mers is:

$$N_{\text{kmers}} = L - k + 1$$

For a typical bacterial genome of $L \approx 5 \times 10^6$ bp and $k = 21$:

$$N_{\text{kmers}} \approx 5{,}000{,}000 - 21 + 1 \approx 4.999 \times 10^6$$

#### Canonical K-mers

DNA is double-stranded. For any k-mer on the forward strand, its **reverse complement** appears on the reverse strand encoding the same biological sequence. To avoid redundancy, we use **canonical k-mers**: the lexicographically smaller of a k-mer and its reverse complement:

$$k_{\text{canonical}} = \min(k, \overline{k})$$

where $\overline{k}$ denotes the reverse complement of $k$.

#### Why $k = 21$?

The choice $k = 21$ is biologically and statistically motivated:

| Criterion | Value for $k=21$ |
|-----------|-----------------|
| Uniqueness probability (random genome) | $\approx 1 - e^{-L/4^k} \approx 1 - e^{-1.19 \times 10^{-6}} \approx 0$ (near-unique) |
| Sensitivity to single-nucleotide mutations | One SNP generates $k$ altered k-mers |
| Ability to span a resistance gene codon | 21 bp covers 7 codons — sufficient for most AMR-relevant mutations |

#### From K-mers to AMR Features — Without Alignment

Traditional AMR pipelines (e.g., ARIBA, ResFinder) require alignment of reads to a reference database of known resistance genes. This approach has critical limitations:
- It misses **novel resistance mutations** not catalogued in databases.
- It is sensitive to reference quality and database completeness.

Our approach is **alignment-free**: we treat every genome as a **bag of k-mers** and learn which k-mers co-occur with resistance directly from phenotypic labels. Any resistance-conferring SNP, insertion, or deletion creates a unique set of k-mers that do not appear in susceptible genomes — and the XGBoost model discovers this signal in the high-dimensional k-mer space.

---

## 2. Feature Space Mathematics: The Curse of Dimensionality

### 2.1 Theoretical Feature Space

The DNA alphabet is $\Sigma = \{A, C, G, T\}$, so $|\Sigma| = 4$. The total number of distinct k-mers over this alphabet is:

$$|\mathcal{F}| = 4^k$$

For $k = 21$:

$$|\mathcal{F}_{21}| = 4^{21} = 2^{42} \approx 4.398 \times 10^{12}$$

Accounting for canonical k-mers (which halve the space):

$$|\mathcal{F}_{21}^{\text{canonical}}| = \frac{4^{21} + 4^{\lceil 21/2 \rceil}}{2} \approx 2.2 \times 10^{12}$$

This is a ~2.2 trillion dimensional feature space.

### 2.2 The Observed Feature Space: Sparsity

In practice, only a small fraction of this theoretical space is observed across real bacterial genomes. Given a dataset of $n$ genomes (isolates), the observed feature matrix $X \in \{0, 1\}^{n \times p}$ has:

- **Rows** = individual bacterial genomes ($n$ samples)
- **Columns** = unique k-mers ever observed in any genome ($p$ features)
- **Values** = binary presence/absence indicator

For a typical AMR dataset:

$$n \sim 10^2\text{–}10^3 \quad \text{and} \quad p \sim 10^6\text{–}10^7$$

This establishes a severe **$p \gg n$ regime** (ultra-high dimensional, small-sample setting).

### 2.3 Sparsity Structure

The matrix $X$ is extremely sparse. Empirically, for $k = 21$:

$$\text{Sparsity} = 1 - \frac{\text{nnz}(X)}{n \cdot p} \approx 0.97\text{–}0.999$$

This is exploited by storing $X$ in **Compressed Sparse Row (CSR)** format (SciPy `csr_matrix`), which stores only the non-zero entries:

$$\text{Storage}_{\text{CSR}} = \mathcal{O}(\text{nnz}(X))$$

compared to a dense matrix requiring $\mathcal{O}(n \cdot p)$ bytes. For our problem, this represents a **100x–1000x memory reduction**.

### 2.4 Prevalence Filtering and Matrix Dimensionality Reduction

Before model training (step 03), uninformative k-mers are removed by filtering on
prevalence across genomes. A k-mer present in (nearly) all genomes carries no
discriminative signal, and one present in too few is rare noise / lineage-specific.
For feature $j$ over $n$ genomes we keep:

$$\text{keep}_j = \mathbf{1}\left[\; s_{\min} \le \sum_{i=1}^{n} X_{ij} \le n-1 \;\right]$$

The upper bound $n-1$ drops zero-variance core-genome k-mers. The lower bound
$s_{\min}$ (minimum support) is **data-adaptive** rather than a fixed count, so the
same configuration behaves sensibly across antibiotics/organisms of very different
sizes:

$$s_{\min} = \max\!\left(\, s_{\text{floor}},\; \lceil \rho \cdot n \rceil \,\right)$$

with an absolute noise floor $s_{\text{floor}} = 5$ (removes singleton/sequencing-error
k-mers regardless of $n$) and a prevalence fraction $\rho = 0.01$ (1%). Thus a large
dataset (e.g. $n=4373$ → $s_{\min}=44$) gets aggressive de-confounding of the rare,
often lineage-specific tail, while a small one ($n \le 500 → s_{\min}=5$) falls back to
the floor and keeps all real markers. The 1% floor is deliberately far below the
$\sim$10%-prevalence a k-mer needs to reach the discriminativeness criterion of the
background-frequency analysis ($|\Delta\text{prev}| \ge 0.10$), so no individually
informative marker is discarded — only the noise/confounder tail. An explicit
`preprocessing.min_support` integer overrides the formula when a fixed value is
desired. This reduces $p$ from tens of millions to a smaller, more informative set
while preserving discriminative features.

---

## 3. Statistical & ML Architecture

### 3.1 Binary Histogram Quantization (`max_bin = 2`)

XGBoost's histogram-based tree learning algorithm discretizes continuous features into bins before split finding. For a dense continuous feature, a typical default setting is `max_bin = 256`, creating 256 potential split points per feature and storing an 8-bit histogram.

Our k-mer features are **binary** ($X_{ij} \in \{0, 1\}$). A binary feature has only one meaningful split point: $X_{ij} < 0.5$ (i.e., absent) vs. $X_{ij} \geq 0.5$ (i.e., present). Therefore, we set:

$$\texttt{max\\_bin} = 2$$

This has a profound impact on memory:

**Memory per feature in XGBoost histogram:**

| `max_bin` | Bits per bin | Total bits/feature |
|-----------|-------------|-------------------|
| 256 | 8 bits | 8 × 256 = 2048 bits |
| 2 | 1 bit | 1 × 2 = 2 bits |

**Memory reduction factor:**

$$\frac{\text{Memory}_{256}}{\text{Memory}_{2}} = \frac{256}{2} = 128\times$$

For $p = 5 \times 10^6$ features, this reduces histogram memory from ~1.28 GB to ~10 MB per tree node — a critical enabler for training on an 8 GB machine.

### 3.2 Optuna HPO and the Square Root Heuristic

#### Hyperparameter Optimization Framework

Optuna performs **Bayesian optimization** using a Tree-structured Parzen Estimator (TPE). For each trial $t$ with parameters $\boldsymbol{\theta}_t$, Optuna fits a probabilistic model over the objective function $f(\boldsymbol{\theta})$ (validation AUC-ROC) and proposes the next trial by maximizing the **Expected Improvement (EI)**:

$$\boldsymbol{\theta}_{t+1} = \arg\max_{\boldsymbol{\theta}} \text{EI}(\boldsymbol{\theta}) = \mathbb{E}\left[\max(f(\boldsymbol{\theta}) - f^*, 0)\right]$$

where $f^*$ is the current best observed value.

#### The Square Root Heuristic for Feature Subsampling

In the $p \gg n$ setting, selecting all $p$ features per tree split is both computationally prohibitive and statistically harmful (overfitting). A well-established heuristic from random forests suggests using $m \approx \sqrt{p}$ features per split. In XGBoost's `colsample_bytree` parameter, this is expressed as a fraction:

$$\texttt{colsample\\_bytree} = \frac{m}{p} = \frac{\sqrt{p}}{p} = \frac{1}{\sqrt{p}} = p^{-1/2}$$

For $p = 5 \times 10^6$:

$$\texttt{colsample\\_bytree} = \frac{1}{\sqrt{5 \times 10^6}} \approx \frac{1}{2236} \approx 4.5 \times 10^{-4}$$

This means each tree sees only ~0.045% of features — a massive regularization effect that simultaneously reduces computation from $\mathcal{O}(p \cdot n)$ per split to $\mathcal{O}(\sqrt{p} \cdot n)$.

> **Implementation note.** The $1/\sqrt{p}$ value is used as the *anchor* of the
> Optuna search space rather than a single fixed value. `04_optimization.py`
> reads the actual feature count $p$ from `features.txt` and searches
> `colsample_bytree` over a **log-scale window bracketing $1/\sqrt{p}$**
> (`compute_colsample_range()`: roughly $[0.5/\sqrt{p},\, 20/\sqrt{p}]$). This
> keeps the search consistent with the square-root heuristic while letting the
> optimizer fine-tune around it. (Earlier versions hardcoded a fixed
> `[0.05, 0.30]` window — ~100× larger than $1/\sqrt{p}$ — which contradicted
> this derivation; that discrepancy has been removed.)

#### Early Stopping for `n_estimators`

Rather than letting Optuna randomly search over `n_estimators`, we **fix `num_boost_round = 1000`** and use XGBoost's built-in early stopping (patience = `early_stopping_rounds`). The optimal number of trees is determined empirically:

$$n_{\text{trees}}^* = \arg\min_{t \leq 1000} \mathcal{L}_{\text{val}}(t)$$

This is captured from `model.best_iteration` and stored as a trial user attribute, then merged into the final configuration. This prevents the Optuna anti-pattern of **random search conflicting with early stopping**, which otherwise leads to overfitting via the selection of unnecessarily large `n_estimators`.

### 3.3 Stratified Linspace Chunk Selection

#### Problem: Imbalanced Mini-Batches

When the full dataset is stored in $C$ chunks on disk and the model is trained on a subset of $k < C$ chunks per trial, naive random chunk selection risks drawing a subset with:
- All resistant samples (minority class dominates)
- Almost no resistant samples (majority class dominates)

This creates **biased gradient updates** that misrepresent the true class distribution.

#### Solution: Stratified Linspace Sampling

Each chunk $c$ has an associated **resistance ratio**:

$$r_c = \frac{|\{i \in c : y_i = 1\}|}{|c|}$$

To select $k$ chunks that collectively preserve the global resistance ratio $\bar{r}$, we sort chunks by $r_c$ and select indices using `numpy.linspace`:

$$\text{selected\\_indices} = \text{round}\left(\text{linspace}(0,\, C-1,\, k)\right)$$

applied to the **sorted** array of $(c, r_c)$ pairs. This ensures selected chunks are spread uniformly across the resistance distribution, providing a balanced sample regardless of which $k$ chunks are selected.

**Formal property:** Let $S = \{c_1, \ldots, c_k\}$ be the selected chunks with resistance ratios $\{r_{c_1}, \ldots, r_{c_k}\}$. The stratified selection minimizes:

$$\left| \frac{1}{k} \sum_{j=1}^{k} r_{c_j} - \bar{r} \right|$$

compared to random selection, by ensuring the selected ratios span the full observed range of $r_c$ values.

### 3.4 Full-Data Boosting over a Streaming `QuantileDMatrix`

#### Problem: the matrix does not fit in RAM, but per-chunk training is weak

The genome × k-mer matrix is far too large to hold densely in memory (e.g. ~109 GB decompressed for ~4.4k genomes × 50.8M k-mers, 21.8B non-zeros). An earlier version handled this by **incremental warm-started boosting** — one tree per chunk via repeated `xgb.train(num_boost_round=1, xgb_model=...)` over shuffled chunks. While this bounded memory, it had two drawbacks: **(i)** each tree was fit to the residuals of only a single ~200-genome chunk, so no tree ever saw the full training distribution (a weaker fit than standard boosting); and **(ii)** the work was dominated by serial chunk decompression with tiny per-tree compute, leaving HPC cores idle (a TRUBA low-efficiency warning at ~13% utilisation).

#### Solution: stream chunks into one quantised DMatrix, then boost normally

XGBoost's external-data iterator API lets us build a **single** in-core, quantised `QuantileDMatrix` by pulling one chunk at a time, without ever materialising the full sparse matrix. We implement `ChunkDMatrixIter` (`scripts/lib/xgb_data.py`), an `xgb.DataIter` whose `next()` loads chunk $c$, optionally applies a sample-level row mask (used by 07b's seed splits), and feeds $(X_c, y_c, w_c)$ to XGBoost. Because the data are binary, `max_bin = 2` makes the quantised histogram ~1 byte per non-zero, so the resulting DMatrix is compact (~22 GB here) and peak memory stays at roughly **one chunk + the histogram**.

Training is then ordinary gradient boosting on the whole training set:

$$\mathcal{L}^{(t)} = \sum_{i=1}^{N_{\text{train}}} \ell\!\left(y_i,\, \hat{y}_i^{(t-1)} + f_t(x_i)\right) + \Omega(f_t)$$

where the sum now runs over **all** $N_{\text{train}}$ training rows for every tree $f_t$, not a single chunk. The number of trees $T_{\text{total}}$ is the `n_estimators` budget found by early stopping during HPO (Section 3.3); we keep `num_boost_round = T_{\text{total}}` over the full DMatrix.

#### Class imbalance

Imbalance is corrected **once**, globally: positive rows receive instance weight $w^{+} = N^{-}_{\text{train}} / N^{+}_{\text{train}}$ (negatives weighted 1.0). HPO (Section 3.3) deliberately leaves `scale_pos_weight` untuned so the correction is applied a single time at training, never double-counted. The operating threshold is fixed at $0.5$ and is **not** tuned on the test set (leakage prevention; Section 4 / `06_evaluation.py` only reports Youden's J).

This same regime is reused by `07b_feature_stability.py`: each of the 5 seeds builds its own train-split `QuantileDMatrix` via a sample-level `row_mask`, so the stability analysis is methodologically identical to the final model. Both are organism/antibiotic-agnostic — the iterator simply streams whatever chunk files it is given.

### 3.5 Reproducibility & MLOps Best Practices

#### Timestamp Versioning and Artifact Provenance
To safeguard high-cost computational artifacts, the optimization (`04_optimization.py`) and training (`05_model_training.py`) scripts employ strict timestamp versioning. Each Optuna study database and XGBoost model binary is backed up with a precise timestamp upon creation. This prevents accidental overwriting during hyperparameter tuning iterations and creates a clear, reproducible lineage for every model deployed.

#### Publication-Ready Source Data Extraction
For maximum scientific transparency and reproducibility, all visualization modules (e.g., `06_evaluation.py`) are engineered to export the exact numerical arrays underlying any generated plot. Alongside every `.png` figure, a corresponding raw `.csv` file is exported, providing the source data required for researchers to independently redraw and modify figures in third-party software such as GraphPad Prism or R.

---

## 4. Explainable AI and Biological Validation

### 4.1 Feature Importance Mapping

Our methodology ensures that the machine learning models remain entirely interpretable. Once the XGBoost model is trained, we extract the top k-mer sequences using the **Gain** metric. In tree-based models, Gain calculates the fractional contribution of each feature to the model's overall predictive power, essentially quantifying how much a specific 21-mer improves the classification of resistance. High-gain k-mers represent critical biological signals. These top features are subsequently converted back into the `.fasta` format to facilitate downstream biological querying.

### 4.2 Automated Nextflow BLAST Pipeline

To translate mathematical importance into biological relevance, we employ a dual-pronged biological validation strategy via an automated Nextflow pipeline:

1. **CARD Local BLAST:** We query the top features against a rigorously curated local installation of the Comprehensive Antibiotic Resistance Database (CARD). This step rapidly identifies acquired resistance mechanisms such as horizontal gene transfer events, specific efflux pumps (e.g., *msbA*), and plasmid-mediated resistance determinants (e.g., *OXA* variants).
2. **NCBI Remote BLAST (`nt` database):** For zero-alignment, reference-free discovery of chromosomal mutations, features are queried against the massive NCBI `nt` database. This alignment-free approach enables the autonomous discovery of core-genome point mutations (SNPs), seamlessly identifying phenomena such as the Quinolone Resistance-Determining Region (QRDR) mutations within the *gyrA* and *parC* genes.

### 4.3 Automated Biological Reporting

To bridge the gap between raw alignment metrics (BLAST `outfmt 6` TSV format) and final biological discovery, an automated reporting mechanism (`09_biological_summary.py`) distills the pipeline's outputs into a synthesized summary. By enforcing strict FASTA header ID matching (e.g., `Rank_1|Score_154.4288|Feature_...`) and implementing regex-based text mining, the script filters low-quality alignments and extracts precise AMR determinants. It cleanly isolates specific resistance symbols (like `OXA-909` or `msbA`) from CARD and unambiguous species/strain identifiers from NCBI, ultimately generating a human-readable, publication-ready Markdown report.

### 4.4 Statistical Signal ≠ Biological Causation (reification safeguard)

Every claim the knowledge base makes about a unitig is **associational, not causal**, and the pipeline is worded and structured to keep that distinction explicit (cf. Takefuji 2025 on the over-interpretation of feature-importance). A high XGBoost gain, a stable CPSS selection frequency, or a large SHAP value means only that a feature *carries predictive signal for the resistance label in this population* — it does **not** establish that the sequence *causes*, *confers*, or *mechanistically determines* resistance. Reporting therefore uses associational verbs ("is associated with", "is predictive of", "co-occurs with") and never causal ones ("causes", "confers", "is responsible for") outside of determinants independently confirmed by an orthogonal, mechanism-aware line of evidence.

Three structural safeguards operationalise this:

1. **Layered, orthogonal evidence rather than a single score.** A biomarker is only elevated to a strong tier when independent lines agree: importance/stability (gain, CPSS π≥0.6, PFER bound), sequence identity to a curated determinant (CARD/ARO BLAST tier), discriminativeness (R-vs-S Δprevalence + Fisher/BH-FDR), allele-level confirmation for SNP mechanisms (step 11 variant model), population-structure-corrected association (pyseer LMM + Bonferroni), and external genotype–phenotype concordance (AMRFinderPlus/ResFinder, M13). Any one alone is treated as a hypothesis, not a fact.

2. **Confounding is measured, not assumed away.** Lineage-aware GroupKFold (PopPUNK) and pyseer's kinship random effect quantify how much of a signal is clonal/lineage co-carriage versus mechanism-driven. The cross-antibiotic H3 result is reported honestly as a *negative* finding (ampicillin=TEM vs cefotaxime=CTX-M/CMY share no gene family; the only β-lactamase overlap is cross-class CTX-M co-carriage) precisely because the method refuses to manufacture a mechanistic story the data does not support.

3. **Provenance over assertion.** Each `validation_evidence` row records the *evidence type, source (tool + version), score, and pipeline run* — so a KB consumer sees *why* a claim is tiered as it is and can re-derive it, rather than trusting an unqualified label. "Novel candidate" strictly means "no curated homolog found", an explicit knowledge gap, not a discovery claim.

This makes the KB's epistemic status auditable: it is a ranked, confidence-tiered, reproducible catalogue of *statistical AMR signals with biological context*, not a claim of causal mechanism.

---

## Summary of Design Decisions

| Component | Problem Solved | Technical Solution |
|-----------|---------------|-------------------|
| k-mer features, $k=21$ | Alignment-free genomic representation | Canonical 21-mer presence/absence vectors |
| CSR sparse matrix + `.npz` chunks | $p \gg n$ high-dimensional sparsity | SciPy CSR format; chunked disk storage |
| `max_bin = 2` | RAM exhaustion from large histograms | 1-bit histograms for binary features (128× reduction) |
| `colsample_bytree = 1/√p` | Overfitting + computational cost | Square Root Heuristic for column subsampling |
| Optuna TPE + early stopping | Conflicting HPO and early stopping | Fixed `num_boost_round`; `best_iteration` captured |
| Stratified linspace chunk selection | Biased mini-batch resistance ratios | Sorted-by-ratio linspace chunk indexing |
| Streaming `QuantileDMatrix` + full-data boosting | Matrix too large for RAM, yet per-chunk training is weak/inefficient | `ChunkDMatrixIter` streams chunks into one quantised DMatrix; standard boosting sees all train rows per tree |
| Nextflow Dual-BLAST (CARD + NCBI) | Black-box ML lack of biological interpretability | Asynchronous pipeline mapping mathematical Gain scores back to known physical AMR mechanisms (SNPs & Plasmids). |
| MLOps Artifact Versioning | Accidental loss of high-cost optimization and model binaries | Strict timestamp-based backup system protecting historical Optuna studies and models. |
| Source Data Extraction | Opaque, irreproducible numerical plots | Automated parallel export of plot arrays to `.csv` for transparent third-party rendering. |
| Automated Biological Reporting | Raw BLAST TSV outputs are unreadable and cluttered | Regex-based `09_biological_summary.py` script distills raw data into synthesized Markdown reports. |

---

*Document version: July 2026 (§4.4 + canonical-pipeline banner added; KB schema 0.4.0). Sections 1–3 document the raw-k-mer baseline; the canonical unitig / lineage-CV / CPSS pipeline (§4.4) is authoritative. Maintained alongside `scripts/` as the mathematical reference for the pipeline.*
