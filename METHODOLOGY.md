# METHODOLOGY.md — ML_AMR_Prediction_v2

> **A rigorous technical exposition of the biological, mathematical, and statistical foundations of the AMR prediction pipeline.**

> ⚠️ **Canonical pipeline (schema 0.7.1) vs this document.** After the
> literature-review pivot, the *canonical* feature unit is
> the **unitig** (compacted de Bruijn graph over k=21, `bcalm2`/`unitig-caller`;
> step 03u), the *canonical* validation is **lineage-aware CV** (PopPUNK +
> StratifiedGroupKFold, step 07b), stability is **CPSS** (B=100, π≥0.6,
> PFER-bounded; step 13) with **TreeSHAP** importance, and biomarkers additionally
> carry pyseer-LMM significance (step 14) and genome QC (CheckM2+QUAST; step 02d).
> **Sections 1–3 below describe the raw-k-mer / Gain / single-split *baseline*
> (still a valid, runnable path) — read them as the foundation, not the current
> canonical method;** §4 and §5 are authoritative where they differ. A full rewrite
> to unitig-first is tracked (audit Issue 1).
>
> **§5 records what the delivered 45-model run actually did**, including the two
> validation layers that produced nothing (step 11 SNP, step 12 MDA) and the two
> that never ran (step 16 M13, step 15 cross-antibiotic). Any Methods text drawn
> from this document must reconcile with §5, not with the design intent in §4.

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
   - [4.2 Dual BLAST Annotation (CARD local + NCBI remote)](#42-dual-blast-annotation-card-local--ncbi-remote)
   - [4.3 Automated Biological Reporting](#43-automated-biological-reporting)
   - [4.4 Statistical Signal ≠ Biological Causation (reification safeguard)](#44-statistical-signal--biological-causation-reification-safeguard)
5. [What the delivered run actually did](#5-what-the-delivered-run-actually-did)
   - [5.1 Panel and parameters as executed](#51-panel-and-parameters-as-executed)
   - [5.2 The evidence ladder as executed](#52-the-evidence-ladder-as-executed)
   - [5.3 Limitations that must be stated](#53-limitations-that-must-be-stated)

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

> **What the delivered run used.** The formula above governs the **raw-k-mer**
> path (`preprocessing.min_support: null` → auto-derive). The canonical **unitig**
> path has its own knob, `unitig.min_support`, and it is set to a **fixed 10** —
> so all 45 delivered models recorded `pipeline_runs.min_support = 10`, across
> dataset sizes from $n=518$ to $n=4495$. The adaptive derivation was therefore
> *not* in force for any delivered model. Two reasons the fixed floor is defensible
> here: unitigs are already ~25× less redundant than raw 21-mers (a compacted
> de Bruijn path replaces the dozens of overlapping k-mers a single SNP generates),
> and the organism-level unitig store is built once at `db_min_support: 2` and then
> subset per antibiotic, so the per-antibiotic floor must stay low enough not to
> pre-drop what a subset would keep. Methods text must quote **10**, not the formula.

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

Our methodology ensures that the machine learning models remain entirely interpretable. Once the XGBoost model is trained, we extract the top features using the **Gain** metric. In tree-based models, Gain calculates the fractional contribution of each feature to the model's overall predictive power, essentially quantifying how much a specific sequence improves the classification of resistance. High-gain features represent critical biological signals. These top features are subsequently converted back into the `.fasta` format to facilitate downstream biological querying.

In the canonical pipeline the queried unit is the **unitig** (variable length; 31–1424 bp in the delivered run), not the fixed 21-mer — this is what makes BLAST meaningful at all, since a 21 bp query yields no interpretable E-value. Step 07 writes `02_top_{N}_features_{ab}.fasta` as the **union of two sets**: the single model's Gain top-$N$ ($N=50$), *plus* the seed-reproducible stable features from 07b that did not make that top-$N$. Both downstream annotators (step 08 BLAST, step 11 SNP check) read this one FASTA. Note that this union is **not** the same candidate set that the knowledge base scores and tiers — the KB's universe comes from step 10 (background frequency) plus step 13 (CPSS); see §5.2 for the consequence.

### 4.2 Dual BLAST Annotation (CARD local + NCBI remote)

To translate mathematical importance into biological relevance, step 08 runs a dual-pronged annotation. Both passes are driven from plain Python via `subprocess` — an earlier Nextflow orchestration was removed in the M9 review, and the pipeline carries no JVM/Nextflow dependency. (The `.nf` file had been deleted while `main()` still shelled out to it, so step 08 could not run at all until this was fixed; the two facts are recorded together so the correction is not silently re-introduced.)

1. **CARD Local BLAST:** We query the FASTA from §4.1 against a local installation of the Comprehensive Antibiotic Resistance Database (CARD). This identifies acquired resistance determinants — horizontally transferred genes, efflux components, and plasmid-mediated β-lactamases such as the *OXA* and *CTX-M* families. The BLAST task is chosen from the **median query length** (`blastn-short` with `word_size 7` below 50 bp, otherwise `blastn`), because the config's general `word_size` truncated or missed full-length hits on 30–50 bp unitigs.
2. **NCBI Remote BLAST (`nt` database):** The same queries are searched against NCBI `nt`. This pass is **decoupled** from the local one: the public NCBI server kills `blastn-short` + `word_size 7` over `nt` with SIGXCPU, so it uses `blastn` + `word_size 11`, and it is **restricted to the study organism** via an `entrez_query` derived from the registry taxid, with `max_target_seqs = 50`.

**These two passes do not have equal standing in the knowledge base, and Methods must say so.** Only the CARD pass enters the KB: `populate_database.py` writes `blast_annotations.source_db = 'card'` and nothing else, so all 3611 annotation rows in the delivered KB are CARD. The NCBI results (45 files, ~231k alignments) live on disk and feed the human-readable step-09 report and the step-18 novel-biomarker context analysis.

This asymmetry is deliberate and should be defended, not apologised for: CARD is a *curated resistance catalogue* carrying the ARO ontology, so a CARD hit is evidence about resistance. Organism-restricted `nt` is a *sequence archive*, so an `nt` hit establishes where a sequence sits in that species — a locus, never a mechanism. Consequently the `blast` evidence layer means "known CARD determinant", and `strong_novel` means "no curated CARD determinant" — an explicit knowledge gap, with the `nt` placement supplied separately by step 18.

### 4.3 Automated Biological Reporting

To bridge the gap between raw alignment metrics (BLAST `outfmt 6` TSV format) and final biological discovery, an automated reporting mechanism (`09_biological_summary.py`) distills the pipeline's outputs into a synthesized summary. By enforcing strict FASTA header ID matching (e.g., `Rank_1|Score_154.4288|Feature_...`) and implementing regex-based text mining, the script filters low-quality alignments and extracts precise AMR determinants. It cleanly isolates specific resistance symbols (like `OXA-909` or `msbA`) from CARD and unambiguous species/strain identifiers from NCBI, ultimately generating a human-readable Markdown report.

This report is a **convenience artefact, not a KB input** — nothing downstream consumes it, and the KB is unaffected by its state. In the delivered run the Markdown is incomplete for 7 of the 45 models because the NCBI Entrez `efetch` calls dropped mid-run. The underlying BLAST TSVs for those models are intact.

### 4.4 Statistical Signal ≠ Biological Causation (reification safeguard)

Every claim the knowledge base makes about a unitig is **associational, not causal**, and the pipeline is worded and structured to keep that distinction explicit (cf. Takefuji 2025 on the over-interpretation of feature-importance). A high XGBoost gain, a stable CPSS selection frequency, or a large SHAP value means only that a feature *carries predictive signal for the resistance label in this population* — it does **not** establish that the sequence *causes*, *confers*, or *mechanistically determines* resistance. Reporting therefore uses associational verbs ("is associated with", "is predictive of", "co-occurs with") and never causal ones ("causes", "confers", "is responsible for") outside of determinants independently confirmed by an orthogonal, mechanism-aware line of evidence.

Three structural safeguards operationalise this:

1. **Layered, orthogonal evidence rather than a single score.** A biomarker is only elevated to a strong tier when independent lines agree: importance/stability (gain, CPSS π≥0.6, PFER bound), sequence identity to a curated determinant (CARD/ARO BLAST tier), discriminativeness (R-vs-S Δprevalence + Fisher/BH-FDR), allele-level confirmation for SNP mechanisms (step 11 variant model), population-structure-corrected association (pyseer LMM + Bonferroni), and external genotype–phenotype concordance (AMRFinderPlus/ResFinder, M13). Any one alone is treated as a hypothesis, not a fact. **This is the design; §5.2 records which of these layers actually contributed to a tier in the delivered run — two of them did not, and one never ran.**

2. **Confounding is measured, not assumed away.** Lineage-aware GroupKFold (PopPUNK) and pyseer's kinship random effect quantify how much of a signal is clonal/lineage co-carriage versus mechanism-driven. The starkest instance is *A. baumannii* ceftazidime: single-split ROC-AUC 0.985, random-CV 0.867, **lineage-CV 0.429** — held out across lineages the model is worse than chance, so what looked like a determinant was a memorised clone. Across the whole panel this inflation is systematic (45/45 models random > lineage; Wilcoxon $p = 5.7\times10^{-14}$; mean inflation $+0.088$) and scales with clonality.

   The cross-antibiotic H3 analysis is reported at the level the data supports: same-class drug pairs share **more ARO gene families** than cross-class pairs (mean overlap 0.84 vs 0.29; Mann–Whitney $p = 0.0015$). The claim rests on that **contrast**, not on individual pairs — 0 of 138 pairs survive Benjamini–Yekutieli correction, and only 5 within-class pairs exist because the panel was curated for class coverage. Components are ARO gene families, not raw unitigs, with the universe defined as the families reachable by that organism.

3. **Provenance over assertion.** Each `validation_evidence` row records the *evidence type, source (tool + version), score, and pipeline run* — so a KB consumer sees *why* a claim is tiered as it is and can re-derive it, rather than trusting an unqualified label. "Novel candidate" strictly means "no curated homolog found", an explicit knowledge gap, not a discovery claim.

This makes the KB's epistemic status auditable: it is a ranked, confidence-tiered, reproducible catalogue of *statistical AMR signals with biological context*, not a claim of causal mechanism.

---

## 5. What the delivered run actually did

Sections 1–4 describe design. This section records the **executed** 45-model run
behind KB schema 0.7.1 (Zenodo DOI `10.5281/zenodo.21789464`). Every figure below
was read back out of the shipped `amrk.db` or the run outputs. Where design and
execution differ, **this section wins**, and thesis Methods must follow it.

### 5.1 Panel and parameters as executed

| Quantity | Value |
|---|---|
| Models · organisms · antibiotic classes | 45 · 6 · 14 |
| Genome–phenotype pairs | 78,556 (per-model $n$ from 518 to 4495) |
| Cross-validation | `lineage_group_kfold_5fold` on **45/45** models — no fallback to random CV anywhere |
| Headline metric | `models.auc_mean_seeds` (lineage-CV, 5 seeds). `models.roc_auc` is a single split and is retained **only** as the contrast that exposes clonal inflation |
| Panel mean lineage-CV ROC-AUC | 0.842 (min 0.429 *A. baumannii* ceftazidime; max 0.975 *E. faecium* teicoplanin) |
| Per organism | *E. faecium* 0.933 (6) · *E. coli* 0.917 (8) · *S. aureus* 0.836 (8) · *K. pneumoniae* 0.834 (11) · *A. baumannii* 0.759 (8) · *P. aeruginosa* 0.755 (4) |
| Unitig support floor | `min_support = 10`, fixed, on all 45 runs (see §2.4) |
| de Bruijn $k$ | 21 for every unitig; unitig length 31–1424 bp |
| Random seed | 42 on all 45 runs |
| CARD snapshot | one version across all 45 runs |

### 5.2 The evidence ladder as executed

`classify_evidence_tier()` folds a CARD BLAST hit plus five statistical layers
(`prevalence`, `snp`, `mda`, `cpss`, `pyseer`) into one grade per (unitig, model).
Three different layer counts circulate in this project's documents and they must
not be conflated:

* **7** — the number of orthogonal analyses the pipeline *produces* (`blast`,
  `background_frequency`, `permutation_mda`, `stability_selection`, `snp`,
  `pyseer_lmm`, `label_permutation`), as listed in `docs/KB_KAVRAMLAR.md`.
* **6** — the number the tier function *counts* (`label_permutation` is
  model-level, not per-unitig, so it grades no biomarker).
* **5** — the number that *actually fires* in the delivered KB (`blast`,
  `prevalence`, `snp`, `cpss`, `pyseer`). The maximum observed
  `n_evidence_layers` is **4**, reached by 9 biomarkers.

Delivered tier distribution over 3571 (unitig, model) pairs: `weak` 1915 ·
`candidate` 947 · `confirmed` 349 · `none` 337 · **`strong_novel` 23**.

One designed layer contributes **nothing**:

**`snp` — was a load-time defect, repaired 2026-09-01.** Step 11 reports its hits
by the FASTA header it queried (`Rank_n|Score_x|Feature_f...`), not by sequence.
`populate_database.populate_snp` fell back to that header column when no `kmer`
column was present and passed it to `unitig_id()`, which registered the
identifier string itself as a unitig. The consequence was **not** that the two
candidate sets were disjoint — they overlap — but that every SNP row was attached
to a freshly minted pseudo-unitig, so the layer joined to nothing while appearing
to have run, and 335 non-DNA rows accumulated in `unitigs`.

The loader now resolves the header back to its k-mer through the same FASTA step
11 queried, and refuses any value that is not DNA. After the repair all 21
`resistant_allele` calls reach the graded universe, the layer fires on **18** of
3571 pairs, and 5 biomarkers move `weak` → `candidate`. `unitigs` drops from 3844
to **3509**. The 18 are the canonical target-site substitutions the homolog-model
BLAST cannot distinguish: *gyrA* S83L (*E. coli*), *gyrA* T83I (*P. aeruginosa*),
*gyrA* S84L and *parC* S80F (*S. aureus*), *parC* S84L (*A. baumannii*).
Step 11 itself was **not** re-run; its outputs were correct all along.
⚠️ The KB therefore differs from the archived v0.7.1 and needs a new version.

**`mda` — a genuine null result (not fixable by rewiring).** Here the sets match
perfectly (2409/2409 candidates overlap the scored set), the
`permutation_significant` column is present in all 45 output files, and its sum
across all 2409 rows is **0**: no unitig survives BH-FDR at $q<0.05$. This must be
reported as a negative finding *with its power caveat* — the test uses $R=100$
permutations against ~2400 candidates, so the attainable $q$ floor is coarse and
the procedure is underpowered. "No MDA-significant unitigs" therefore does **not**
license "no feature matters"; it means this particular test could not resolve one.

Both tables that once shipped empty have since been filled, without a schema
change. `external_concordance` holds the step 16 / M13 comparison against
AMRFinderPlus 4.2.7 and ResFinder 4.5.0: 83 rows over 44 of the 45 models, run
2026-09-01. `unitig_antibiotic_overlap` holds step 15's cross-antibiotic output,
29 (unitig, pair) records over 160 organism-internal pairs, so the API's
`/overlap` route now returns data.

For the annotation layer, all 3611 `blast_annotations` rows are CARD (§4.2), and
they are **not interchangeable**: 3007 sit at `tier='none'` (mean query coverage
0.38, E-values up to 9.3 — spurious homology from short queries), against 480
`confirmed` (99.9% mean identity, full coverage), 117 `weak` and 7 `candidate`.
All 2035 rows whose `gene_symbol` is the literal string `"nan"` fall in
`tier='none'`. **Any biological statement drawn from this table must filter
`tier IN ('confirmed','candidate')`**; without that filter the unfiltered join
returns, for example, staphylococcal *mecA* under an *A. baumannii* model.
Filtered, the panel reproduces textbook biology: *S. aureus* cefoxitin/oxacillin →
*mecA*+*mecR1*; *E. faecium* teicoplanin → all seven *vanA*-cluster genes;
*E. coli* cefotaxime → CTX-M; *K. pneumoniae* carbapenems → KPC/NDM;
*A. baumannii* carbapenems → OXA; and the models with no CARD hit are exactly the
mutation-driven/intrinsic ones.

Finally, the 23 `strong_novel` biomarkers were placed against organism-restricted
`nt` (step 18): **23/23 align at 100% identity over full query length**, so none is
an assembly artefact. By replicon, using an ≥80% majority over all retained
alignments: 10 chromosomal, 5 plasmid, 8 mixed. The mixed class is informative
rather than a failed call — a sequence occurring at full identity on both
replicon types within one species is a mobile-element signature.

### 5.3 Limitations that must be stated

1. **External validation is lineage hold-out only.** Temporal validation is
   *impossible* with this data: BV-BRC AMR phenotypes end in 2021 (≥2023 isolates:
   *E. coli* 28, *K. pneumoniae* 13, *A. baumannii* 11, *S. aureus* 0). Geographic
   validation was *not performed*, and the collections are country-dominated
   (*E. coli* 58% Norway, *A. baumannii* 63% USA). M13 concordance HAS now run
   (2026-09-01), but it scores the tools on the model's own held-out split, which
   is a chunk split rather than a lineage-aware one — a design that favours the
   model by the margin §5.3 item 4 measures.
2. **Labels are BV-BRC as published.** No MIC re-interpretation was attempted —
   raw MIC completeness is as low as 9% (*S. aureus*) and units are mixed.
3. **PFER exceeds 1 in 21 of the 45 models** (max 12.9): in those stable sets the
   expected number of false positives is above one. Both numbers are read from
   `kb_overview.csv` at draw time by figure 02 and by `limitations.csv`, so the
   claim and the artefacts cannot drift apart. (This item previously read "max ≈14",
   which no delivered artefact supported — figure 02 has always rendered 12.9.)
4. **Co-carriage is linkage, not causation.** `sul`/`qacEdelta1` co-occurrence is a
   class-1 integron; a novel *K. pneumoniae* gentamicin unitig maps to an
   MCR-1-carrying plasmid without being *mcr-1*. All KB claims are associational.
5. **`git_dirty = 1` on all 45 runs.** Each run records its commit, but the working
   tree was modified at execution time, so the commit hash alone does not
   reconstruct the exact code state. Seed, config hash, CARD version and tool
   versions are recorded and consistent; `bcalm` reports an honest NULL version
   (no version CLI).
6. **`nt` proportions are bounded samples.** The remote pass used
   `max_target_seqs = 50`, so replicon proportions in §5.2 are proportions of the
   alignments BLAST retained, not a census of `nt`.
7. **The `organisms` table holds 7 rows, the panel has 6.** *Enterobacter cloacae*
   is registered in the reference table but has no trained model.
8. **Genome QC enforces two of the four criteria it measures.** CheckM2/QUAST ran on
   all 17,742 assemblies and 98.7% (17,516) pass, but the enforced gate is
   completeness ≥95% and contamination ≤5% only. N50 ≥50 kb and contigs ≤500 are
   computed and reported, then deliberately not applied: an N50 gate removes **1,305
   of 2,078 *E. faecium* genomes (63%)**, which selects on assembly provenance rather
   than genome quality for a species whose BV-BRC entries are routinely short-contig
   drafts. Figure 12 plots all four and marks which two are enforced. The consequence
   is that assembly fragmentation is an uncontrolled covariate, most heavily in
   *E. faecium* — the organism with the highest panel AUC (0.933).
9. **pyseer p-values are strongly inflated relative to a uniform null** (genomic
   inflation λ = 0.7–79.3, median 2.6). Some of this is real: thousands of unitigs in
   tight LD tag the same locus. The QQ plot cannot separate that from stratification
   the kinship term did not absorb, so figure 29 is a diagnostic of the scan, not
   evidence that population structure was fully controlled.

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
| Dual BLAST (CARD local + NCBI remote), pure-Python `subprocess` | Black-box ML lack of biological interpretability | Maps Gain/CPSS-selected unitigs back to curated determinants (CARD/ARO → the KB) and to genomic locus (organism-restricted `nt` → context only). |
| MLOps Artifact Versioning | Accidental loss of high-cost optimization and model binaries | Strict timestamp-based backup system protecting historical Optuna studies and models. |
| Source Data Extraction | Opaque, irreproducible numerical plots | Automated parallel export of plot arrays to `.csv` for transparent third-party rendering. |
| Automated Biological Reporting | Raw BLAST TSV outputs are unreadable and cluttered | Regex-based `09_biological_summary.py` script distills raw data into synthesized Markdown reports. |

---

*Document version: August 2026 (KB schema 0.7.1, 45-model delivered run). Sections 1–3 document the raw-k-mer baseline; the canonical unitig / lineage-CV / CPSS pipeline (§4) is authoritative for design, and §5 is authoritative for what was executed. Maintained alongside `scripts/` as the mathematical reference for the pipeline.*
