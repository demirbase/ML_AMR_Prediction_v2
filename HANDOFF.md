# AMR k-mer Knowledge Base — Project Handoff Document

> **Repo:** `ML_AMR_Prediction_v2` · **Active branch now `main`** (HEAD `3ddc476`, pushed to `github.com/demirbase/ML_AMR_Prediction_v2`). (Earlier work was on `fix/amr-audit-remediation`, since merged to `main`.)
> **Local (Mac) path:** `~/Desktop/IU_master/projects/ML_project_kopyasi`
> **Last full local run:** 2026-06-15, *E. coli* / ampicillin, 1788 genomes — `00a → 11` end-to-end.
> **LIVE: real production run on TRUBA (ARF) HPC, 5470 genomes — in progress (see §0 below).**

---

# 0. TRUBA (ARF) deployment — LIVE STATE (resume here)

> A full real run is in progress on the **TRUBA ARF** cluster (user `edemirbas`). New session: continue from here.

**Where we are:** Data acquired (5470 *E. coli* genomes). Pipeline on TRUBA: `00a✓ 00✓ 01✓ 02 KMC✓ 02b QC✓ 03 matrix✓` — **FEATURES DONE**: `data/processed/ecoli/ampicillin/matrix/` has **22 `.npz` chunks + `features.txt` (1.27 GB)** = **4446 ampicillin genomes × 50.8M k-mers** (~90% sparse; full matrix ≈ **109 GB** decompressed / 21.8B nnz).

**Training regime refactored to full-data boosting (DONE & re-run on TRUBA 2026-06-21).** The old run (2026-06-18, incremental 1-tree/chunk) gave test ROC-AUC 0.903 / MCC 0.693 / acc 0.84. Refactored to **standard full-data gradient boosting** (`scripts/lib/xgb_data.py`); 04/05/07b updated and pushed to `main` (HEAD `a5b9ddc`). Key engineering hurdles solved on TRUBA: 04 HPO now runs **parallel Optuna trials** over a `QuantileDMatrix` subset (`training.optuna_threads_per_trial`, default 2) to use cores without OOM; 05/07b use **`ExtMemQuantileDMatrix`** (external memory, pages spilled to scratch) because an in-core full-train QuantileDMatrix peaked **>400 GB** and OOM-killed the 384 GB barbun node.

**New result (2026-06-21, 4373 genomes, 50 HPO trials, full-data boosting):** test **ROC-AUC 0.930 (CI 0.914–0.945), PR-AUC 0.965, MCC 0.739, acc 0.866, balanced-acc 0.885** (threshold 0.5). Improved over the old regime across the board (MCC 0.693→0.739). NB: 04+05+06 jobs ran to completion; the TRUBA **low-efficiency warning DID fire** during ExtMem training (Eff ~10%, I/O-bound) but **did NOT kill the job** — confirms the warning is a nag, not an auto-killer (every actual death this project was OOM or manual `scancel`).

**UNDERFITTING HYPOTHESIS TESTED → DISPROVEN (2026-06-21/22).** Added full-data early stopping to 05 (`training.max_boost_rounds`, `validation_fraction` for the ES split) and ran it: early stopping found **best tree count = 29** (≈ the HPO subset's 30), val AUC peaks there, more trees overfit. So the model is **well-fit, not underfit**; ROC-AUC ~0.93 / MCC ~0.73 is the **realistic signal ceiling** for ampicillin k-mers (consistent with literature), not a tuning bug. The compressed probability range is just low-lr × few-trees, not capacity-limited. **No ML tweak (early stopping, more trees, min_support) will raise this ceiling — it lives in the feature representation.** Note also: `min_support` removes rare k-mers but barely reduces nnz (dominated by common k-mers), so it does NOT shrink RAM (only fewer GENOMES via `training.max_train_chunks` does); its real value is lineage/noise de-confounding of candidate k-mers, a KB-quality lever, not a prediction lever.

**STRATEGIC PIVOT (agreed): stop optimising the prediction engine; move to the thesis contribution.** Per `docs/ROADMAP.md` the research question is the **queryable, biologically-validated AMR k-mer Knowledge Base**, NOT prediction accuracy — AUC 0.93 is already more than enough to rank k-mers for the KB. ROADMAP M2 accepts **5-seed repeated holdout OR 5-fold CV** (we use 5-seed = 07b); no must-have requires a prediction-accuracy bar. So the training regime is not a thesis lever; we keep **standard full-data boosting** (more conventional/defensible than the old incremental 1-tree/chunk, which a reviewer would question).

**GPU evaluated and REJECTED (2026-06-22).** Built `amr-gpu.sif` (CUDA xgboost, `amr-gpu.def`, USE_CUDA True) and smoke-tested on a Tesla **V100 16 GB** (akya-cuda/barbun-cuda, AllowAccounts=ALL). The 50.8M-feature ultra-wide sparse matrix is GPU-hostile: in-core GPU OOMs (3 chunks alone wanted 11 GB — GPU ELLPACK ≈4 bytes/nnz → full set ≈88 GB ≫ 16 GB) and ExtMem-GPU did not even finish building in 25 min. **CPU is the path.** (`amr-gpu.def` kept in repo for the record; don't pursue GPU on this hardware for this data.)

**ENGINE SETTLED — feature filter made DATA-ADAPTIVE (the real lever for speed + KB quality + generality).** Since the AUC ceiling is feature-bound and full-data training is slow only because of the 50.8M features, 03 now derives min_support adaptively: `min_support = max(min_support_floor=5, ceil(min_prevalence=0.01 * n_genomes))` (config `preprocessing`, with `min_support: null` = auto, or an int to force). This scales across antibiotics/organisms (ampicillin 4373 → 44; 1788 → 18; ≤500 → floor 5), so a small dataset is never over-filtered and a large one gets de-confounding + ~2× faster training. Biologically safe: 1% prevalence is far below step-10's ~10%-prevalence discriminativeness threshold, so no individually-strong marker is dropped (only the rare/lineage/error tail — exactly the confounders ROADMAP §⚠️/S3 flag).

**LITERATURE REVIEW DONE (2026-06-22) → MAJOR METHODOLOGICAL PIVOT (see `docs/ROADMAP.md` §0).** Two rounds of systematic literature research were completed and distilled into binding decisions now in ROADMAP **§0** (authoritative; supersedes older ROADMAP sections). The current raw-k-mer + adaptive-min_support engine WORKS and is a valid baseline, BUT the literature mandates a publication-grade overhaul. Key decisions:
- **Unitigs replace raw k-mers** (`bcalm2` + `unitig-caller`): ~10M→~730k features, ~212→~18 GB, ~7h→~50min, BLAST-mappable, GWAS-standard. Downstream XGBoost unchanged (binary matrix). **This dissolves the min_support/speed/GPU pain we fought all session** — so the raw-k-mer 03 rebuild + GPU work are now mostly moot.
- **Lineage-aware CV** (PopPUNK clusters → `GroupKFold`) replaces random/chunk split. Random CV inflates AUC 20-30%. **Biggest reviewer-blocker.** Final model still trained on all data.
- **Stability: CPSS (B=100, 50%, π≥0.6) + Chi²/MI staged prefilter + SHAP** replaces 5-seed 07b. Importance: **SHAP** not Gain.
- **MUST add:** external validation (temporal/geographic hold-out + AMRFinderPlus/ResFinder concordance: Kappa/McNemar/bACC); pyseer LMM+Bonferroni; CheckM2+QUAST QC; BH-FDR on step 10; ARO/CARD ontology mapping in the KB; reification-safe wording.
- **Confirmed (no change):** k=21, binary+max_bin=2, class-weight (no SMOTE), BV-BRC+EUCAST/CLSI, 4373 genomes adequate, AUC ~0.93 consistent with literature.
- **Novelty reframe:** NOT "first ML AMR DB" (BV-BRC exists) → "first PFER-bounded (CPSS), lineage-validated, k-mer/unitig-resolution, transparent+FAIR open AMR biomarker KB." Target: Database(Oxford)/Briefings in Bioinformatics. Showcase: ciprofloxacin (gyrA/parC SNP) + a β-lactam (acquired gene).

**TRUBA state at handoff:** raw-k-mer 03 matrix rebuild (job 5952577, adaptive min_support=44 → 21.4M features) was running/likely done — **but it is now superseded by the unitig pivot, so do NOT keep building on it.** KB schema (`scripts/lib/kb_schema.py`) + `scripts/populate_database.py` (M8 foundation) are built, tested, committed (HEAD `main`). CARD 4.0.1 full bundle downloaded (step 11 active). TRUBA scratch cleaned (GPU `.sif`, kmc_curve, smoke removed); pending: `rm -rf $AMR_WORK/data/interim/ecoli/kmc_outputs/tmp/*` after 03.

**IMMEDIATE NEXT (resume here):** ROADMAP §0 order: (1) ✅ **unitig pipeline**; (2) ✅ **PopPUNK lineage + GroupKFold**; (1.5) ✅ **unitig ML run 04→05→06→07b** (RESULTS in §0.1); (Block 1) ✅ **07→11 biology on unitigs DONE & reproducible** (NCBI remote BLAST fixed — see §0.1 "Block 1 + KB DONE"); (M8) ✅ **KB populated** (`results/ecoli/kb/amrk.db`, schema 0.2.0); (M16) ✅ **ARO ontology in KB**; (M6) ✅ **CARD 4.0.1 recorded**. **→ NOW START Block 2:** (a) **M9 MDA permutation** (next, fastest); (b) M14 pyseer LMM + Bonferroni (`amr-tools.sif`); (c) M4 CPSS+SHAP rewrite of 07b; (d) M13 external validation (temporal/geo + AMRFinderPlus/ResFinder); (e) M15 CheckM2/QUAST QC (`amr-checkm2.sif`); then M10 provenance (run_metadata/git hash — currently MISSING in KB), S8 FastAPI. S10 reification + BH-FDR(10) already coded.

## §0.1 Unitig + lineage-CV + biology — STATUS (2026-06-24)

### Block 1 + KB — DONE & REPRODUCIBLE (2026-06-24)
Block 1 biology ran end-to-end on unitigs and the **KB is populated**. Commits `1e09219`→`853143f` on `main`.
- **NCBI remote BLAST fixed (the session's main fight).** The public NCBI server **kills `blastn-short` + `word_size 7` over nt with SIGXCPU** (CPU-usage limit) — short 7-base seeds explode across nt **even when restricted to one species** (tested, still SIGXCPU). So the remote pass is now **DECOUPLED from CARD** (`08_blast_pipeline.nf` + `08_blast_annotation.py`): CARD pass keeps `blastn-short`/word7 (local, fine); **NCBI pass uses `blastn` + `word_size 11` + `-entrez_query txid<taxid>[Organism:exp]` + `-max_target_seqs 50`**. The entrez_query is **taxid-based from the registry** (not the scientific name — a space breaks the Nextflow CLI launcher → "Illegal option --"). Result: CARD 3605 + NCBI 4522 hits (E. coli-restricted). `09` got an `AMR_ENTREZ_EMAIL` env override (no config.yaml edit).
- **Nextflow-under-nohup gotcha:** Nextflow's ANSI console does terminal ioctls; a backgrounded JVM gets **SIGTTOU and STOPS (state `T`)** before submitting processes. Fix baked in: 08 sets **`NXF_ANSI_LOG=false`** (+ strips `NXF_OPTS` to kill the benign "Illegal option --"). HPC background runs now unblocked.
- **M8 KB populated** (`results/ecoli/kb/amrk.db`, **schema 0.2.0**): 65 kmers, 60 model-scores, 60 blast_annotations, 60 background-freq, 11 SNP, 131 validation_evidence. `populate_database.py` prefers `10_kmer_background_frequency` (superset of `07_kb_candidates` — has all candidate cols + prevalence/fisher/discriminative).
- **M16 ARO in KB** (schema bumped 0.1.0→0.2.0): `blast_annotations` gained `aro_accession/aro_gene_family/aro_drug_class/aro_resistance_mechanism`; **13/60 ARO-mapped** (= the 13 confirmed CARD hits: CMY-198 `ARO:3008132`, OXA-1042, CTX-M-260/278, TEM-258 `ARO:3009077`, sul1).
- **M6 CARD version recorded** via `AMR_CARD_VERSION` env override → `kb_metadata.card_version = 4.0.1`.
- **Validation metrics (M7/H2):** 13 confirmed / 47 none; **recovery 32% → H2 FALSE** (<40%); novel fraction 68% (19 stable-novel unitigs → H4). H2 failing is OK per ROADMAP §0.4 reframe (KB value = PFER-bound + lineage-validated + transparent + novel, not recovery %); revisit under M4 CPSS+SHAP.
- **Env overrides added this session:** `AMR_ENTREZ_EMAIL`, `AMR_CARD_VERSION` (both bypass config.yaml). To repro Block 1: `export AMR_FEATURE_REPR=unitig AMR_ENTREZ_EMAIL=… NXF_ANSI_LOG=false` then `08_blast_annotation.py` (UI node, internet) → `09` → `populate_database.py` (`AMR_CARD_VERSION=4.0.1`).
- **Known KB gaps (→ M10):** `run_metadata.json` MISSING (04 HPO was cut short) → `pipeline_runs` has no git_commit/timestamp (run_id `…__unknown`); `blast_annotations.coverage`/`description` NULL (not emitted to candidate CSV); `delta_prevalence` not in CSV.

### RESULTS (the thesis headline)
- **Unitig matrix:** 4,938,938 unitigs × 4373 genomes (median len 34 bp, min 31, p90 54, max 10030), 22 chunks, ~8 GB (vs 21.4M k-mers / 49 GB).
- **06 chunk-split test (unitig):** ROC-AUC **0.9534** (CI 0.936–0.969), MCC **0.8185**, bal-acc 0.913 — **beats k-mer baseline (0.930 / 0.739)**.
- **07b lineage-aware 5-fold GroupKFold (HONEST headline):** ROC-AUC **0.9505 ± 0.0102** (folds 0.934/0.943/0.961/0.957/0.957), **28 stable unitigs** (freq≥0.6), Jaccard 0.236. Lineage-CV ≈ chunk-split → **near-zero lineage leakage → signal is mechanism-driven, generalises across lineages** (answers the Yu-2024 reviewer-blocker).
- **05 early stopping:** best tree count = **8** (lr 0.024, depth 10). `n_estimators: 8` is in the experiment config.
- **Biology (08 fixed):** candidate unitigs map full-length (cov=1, 100% id) to β-lactamases **TEM-258/257, CTX-M-260/278, OXA-1042, CMY-198** (+ catA1, sul1, aadA24) — the real ampicillin mechanisms. rank-1 (32 bp, best 14 bp) is likely novel (H4).

### Containers on TRUBA (`$AMR_WORK/containers/`) — all built
- `amr.sif` — core + `unitig-caller 1.3.2 (Bifrost)` + `bcalm`. (NOTE: rebuild pulled **Nextflow 26.04.4**, strict parser — see 08 fix below.)
- `amr-pp.sif` — core + `poppunk 2.7.8` (env pins `setuptools<81`; PopPUNK needs legacy `pkg_resources`).
- `amr-tools.sif` — `pyseer 1.4.1 + quast 5.3.0 + ncbi-amrfinderplus 4.2.7 + resfinder` (M13/M14). DBs NOT baked: `amrfinder -u`, ResFinder DB clone — download to scratch when needed.
- `amr-checkm2.sif` — `checkm2 1.1.0` (separate: it pins python<3.9). DB: `checkm2 database --download --path $AMR_WORK/data/external/checkm2_db`.
- Build on a **debug node**: `unset APPTAINER_BINDPATH` (build sandbox can't bind `/arf`); login node hits a CPU-time ulimit on `mksquashfs`. `def`s defined in `amr.def`/`amr-pp.def`/`amr-tools.def`/`amr-checkm2.def`.

### Scripts / wiring (all on `main`, HEAD ~`400657f`)
- **`03u_unitig_matrix.py`** — `unitig-caller --call --rtab` → 03's exact chunk contract → `matrix_unitig/`. `--build-db` = organism-level store (`processed/{org}/unitig_all/`); per-antibiotic then SUBSETS it (no re-run). Config `unitig:` (`out_subdir`, `min_support: 10`, `db_min_support: 2`, `threads`).
- **`02c_lineage_poppunk.py`** + **`lib/lineage.py`** — PopPUNK **dbscan** (bgmm degenerate→refine NaN) → `processed/{org}/lineage/poppunk_clusters.csv`. `group_kfold_masks` → StratifiedGroupKFold sample masks. Config `lineage:` (`model: dbscan`, `refine: false`, `n_splits: 5`).
- **`07b`** — `build_cv_splits()`: lineage GroupKFold if `poppunk_clusters.csv` exists, else 5-seed fallback. n_estimators (=8) read from experiment config best_params.
- **`08`/`08_blast_pipeline.nf`** — FIXED: (a) `def OUTFMT`→`params.outfmt` (strict Nextflow 26 rejected top-level `def`); (b) BLAST task by **median query length** (`blastn-short` if median<50, with **word_size 7** — word_size 11 / `blastn` truncated short-unitig hits to ~14 bp noise). `blast.task` overrides.
- **`09`** — ARO mapping (M16): `aro_from_sseqid` + `load_aro_index` (`data/external/card/aro_index.tsv`) → KB cols `aro_accession/aro_gene_family/aro_drug_class/aro_resistance_mechanism`; coverage now = aln/`qlen` (08 emits qlen); reification note (S10).
- **`10`** — BH-FDR (§0.2): `fisher_q` + `discriminative_fdr`.
- **ENV OVERRIDES (no config edit on TRUBA):** `AMR_FEATURE_REPR=unitig` (matrix_dir→matrix_unitig), `AMR_EXTERNAL_MEMORY=false` (in-core 05/07b — faster but ~211 GB RAM for 4.9M feats), `AMR_OPTUNA_PATIENCE=15` (04 early stop), `AMR_POPPUNK_BIN`, `AMR_UNITIG_CALLER_BIN`.

### TRUBA artefacts present (verified)
matrix_unitig (features.txt + 22 X chunks + y + genomes + 70 GB unitigs.rtab) · lineage/poppunk_clusters.csv (324) · models/ecoli/ampicillin/xgboost_ampicillin_final_v2.json (8-tree unitig model) · config/experiments/ecoli/config_ampicillin.yaml (MANUAL — built from Optuna trial 20 after HPO was cut short; `n_estimators: 8`) · full CARD at data/external/card/ (card.json, aro_index.tsv, variant_model fasta) · CARD homolog DB at data/external/blast_db/card_nt/.

### Manual-config caveat (why)
The unitig 04 HPO (50 trials, ~5 h on 4.9M feats, slow) was **cancelled at trial 35**; best = **trial 20** (subset AUC 0.9903). `config_ampicillin.yaml` was written BY HAND (snippet) with trial-20 params + a manual_linspace chunk split (test = parts 0,7,14,21). 04 itself was NOT completed for unitigs. To redo properly later: re-run 04 with `AMR_OPTUNA_PATIENCE=15` (now coded) so it stops ~trial 35 cleanly.

### TRUBA rules that bit us
- **Compute nodes have NO outbound internet** → `blastn -remote` (NCBI) + Entrez FAIL there (errorStrategy 'ignore' keeps CARD alive). Run **08 + 09 on the UI node (arf-ui1)** for NCBI/Entrez; 10/11 on compute (offline OK). Local nt (~200 GB, N6) is the alternative.
- **barbun = whole-node feel:** `-c40 --mem 300G` waits for a fully-free node + fairshare (heavy daily usage lowers priority); `-c20 --mem 100G` (half node) schedules instantly. In-core 05 used **211 GB** (4.9M feats) — request ≥256 G if `AMR_EXTERNAL_MEMORY=false`.
- Pull code with **targeted `git checkout origin/main -- <file>`**; NEVER `config.yaml` (manual HPC tuning) or `02p`/`02b`/`03` (manual parallel patches). The TRUBA `config.yaml` has hand-added `unitig:`/`lineage:` sections, but env overrides make config edits unnecessary.

## §0.2 SLURM templates (copy-paste; submit from `$AMR_WORK`)

Common header rules: submit from `/arf/scratch` (`cd $AMR_WORK && sbatch ...`); `export APPTAINER_BINDPATH=/arf`; env overrides set the behaviour (no config edit); `2>&1` so Optuna/Nextflow stderr lands in the `.out`. `H=/arf/home/edemirbas/ML_AMR_Prediction_v2`, `SIF=$AMR_WORK/containers/amr.sif`. **Half node (`-c20 --mem 100G`) schedules instantly; full node (`-c40 --mem 300G`) waits.** Internet-needing steps (08-NCBI, 09-Entrez) run on the **UI node** (not SLURM).

**A) Unitig ML chain 04→05→06→07b** (`-c40 --mem 300G` in-core, or `-c20 --mem 120G` + `AMR_EXTERNAL_MEMORY=true` ExtMem):
```bash
#SBATCH -J amr-ml -p barbun -N1 -c20 --mem=120G --time=1-00:00:00 -o amr-ml-%j.out -e amr-ml-%j.err
set -euo pipefail; export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig AMR_EXTERNAL_MEMORY=true AMR_OPTUNA_PATIENCE=15
cd $H; for s in 04_optimization 05_model_training 06_evaluation 07b_feature_stability; do apptainer exec $SIF python -u scripts/$s.py 2>&1; done
```
**B) 07b only (in-core, fast)** — `-c40 --mem 300G`, `AMR_EXTERNAL_MEMORY=false`; needs the experiment config from 04.
**C) Biology 10+11 (compute, offline)** — `-c20 --mem 100G`, `AMR_FEATURE_REPR=unitig`; run `10_kmer_background_frequency` + `11_variant_snp_check`.
**D) Biology 07+08+09 on the UI node (interactive, internet for NCBI/Entrez):** `export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig; cd $H; apptainer exec $SIF python -u scripts/08_blast_annotation.py` then `09_biological_summary.py`.
**E) Organism unitig store (once/org, long)** — `-c40 --mem 300G`: `apptainer exec $SIF python -u scripts/03u_unitig_matrix.py --build-db --db-min-support 10`.
**F) PopPUNK lineage (once/org)** — `-c20 --mem 120G`, `SIF=amr-pp.sif`: `python -u scripts/02c_lineage_poppunk.py` (dbscan; `--reuse-db` to skip re-sketch).
**G) Container build (debug node)**: `unset APPTAINER_BINDPATH; export APPTAINER_TMPDIR=/tmp/apptmp APPTAINER_CACHEDIR=/tmp/apcache; cd $H; apptainer build --fakeroot $AMR_WORK/containers/<name>.sif <name>.def`.
**H) Block-2 tool runs** (`SIF=amr-tools.sif`/`amr-checkm2.sif`): pyseer LMM (M14), CheckM2+QUAST QC (M15), AMRFinderPlus/ResFinder concordance (M13) — scripts not written yet.

**Connection / layout (all on TRUBA):**
- Login: OpenVPN → `ssh edemirbas@172.16.6.11` (UI = `arf-ui1`; transfer hosts `arf-ui4/5` = .14/.15).
- `$AMR_HOME=/arf/home/edemirbas/ML_AMR_Prediction_v2` (git clone of `main` = **code**).
- `$AMR_WORK=/arf/scratch/edemirbas/amr` (**data/outputs/container**). `$SIF=$AMR_WORK/containers/amr.sif`.
- All three + `APPTAINER_BINDPATH=/arf` are in `~/.bashrc`; `~/.bash_profile` does `source ~/.bashrc`.
- In the repo, `data/ results/ logs/ runs/ models/` are **symlinks → `$AMR_WORK/…`**. CARD homolog DB copied to `$AMR_WORK/data/external/blast_db/card_nt/`.

**Environment = Apptainer container** (TRUBA forbids conda/pip on the shared FS):
- `apptainer` is at `/usr/bin/apptainer` (v1.3.6, **no `module load` needed**).
- `amr.sif` built from `$AMR_HOME/amr.def` (Bootstrap docker `condaforge/miniforge3` + `environment.yml`) on an **interactive debug node** with `apptainer build --fakeroot` (set `APPTAINER_TMPDIR/CACHEDIR=/tmp/...`). Contains python+xgboost+sklearn+biopython+certifi, KMC 3.2.4, BLAST+ 2.17, Nextflow.
- Run everything as `apptainer exec $SIF python scripts/...`.

**Hard TRUBA rules learned (critical):**
1. **Submit jobs from `/arf/scratch`** (`cd $AMR_WORK` before `sbatch`/`srun`) — else `srun: error: Lutfen islerinizi /arf/scratch/ ...`.
2. **`APPTAINER_BINDPATH=/arf`** required, else the container can't see scratch (symlinks → `mkdir`/FileNotFound errors).
3. **`ftp.bv-brc.org` is FIREWALL-BLOCKED** on TRUBA. Genome FASTAs are fetched from the **BV-BRC Data API** (`www.bv-brc.org/api/genome_sequence`, dna+fasta) — already the repo default in `00a` (commit `5d0c9a3`).
4. **Queues:** `barbun` **min 20 cores/node** (hamsi 28, orfoz 56). **MS-student limit = 40 cores.** → use `barbun -c 20` (or `-c 40` for the parallel ML step). `debug` ≤4h for tests.
5. **Low CPU efficiency → TRUBA warns (`Eff:%…`) and may auto-cancel + cut your core quota.** barbun's 20-core minimum means any single-threaded step looks ~5%. Mitigations applied & on `main`: **`scripts/02p_kmer_parallel.py`** (parallel KMC, 5470 genomes in ~2.5 min), **parallel 02b** spectra extraction (`1cd0119`), **parallel 03** per-genome KMC dump (`3ddc476`). **Caveat:** 03's per-genome *parse* (matching 5M k-mers/genome against the ~8 GB `kmer_to_index` dict) is **GIL-bound**, so thread parallelism only partly helps → 03 still ran ~5% and got warned. **03 is resume-safe** (skips existing `*.npz` + `features.txt`), so a kill is harmless — just resubmit. **True 03 fix (future):** 2-bit integer-encode k-mers + numpy `searchsorted` (drop the Python dict) or multiprocessing. **ML steps (04/05) use cores via `n_jobs=40`** → expected high efficiency.
6. Quotas (banner): `/arf/home` 100 GB / 100K inode; `/arf/scratch` 1 TB / 200K inode; **no backup**; `/arf` is NVMe Lustre (fast — no separate `/tmp` needed for KMC).

**`config.yaml` on TRUBA (NOT committed — TRUBA-specific tuning):** `kmc_mem:128`, `threads:20`, `n_jobs:40`, `chunk_size:200`, `n_trials:30`, `target_antibiotic:ampicillin`.

**SLURM scripts in `$AMR_HOME/slurm/`:** `00_test.slurm` (debug sanity ✓), `run_features.slurm` (02p→02b→03, barbun `-c 20 --mem 120G`, **done**), `run_matrix.slurm` (03-only spare), `run_ml.slurm` (04→05→06, `-c 40 --mem 300G --time 3-00:00:00`, **current**). All use `export APPTAINER_BINDPATH=/arf`, `set -euo pipefail`, mail to `eren0demirbas@gmail.com`, submit from `$AMR_WORK`.

**Data state:** 5470 genomes downloaded (395 had no API sequence, dropped). `amr_phenotypes.csv` = 5470×72. ampicillin 4446 tested (recommended target). KMC dbs for all 5470 in `$AMR_WORK/data/interim/ecoli/kmc_outputs/`.

**Commits made during deployment (all on `main`):** `bvbrc` NaN-safe filters (`17b1301`); `00a` FTP→API download + new `02p` parallel KMC (`5d0c9a3`); `02b` str-path KMC check (`c352c48`); **parallel 02b spectra (`1cd0119`)**; **parallel 03 dump (`3ddc476`)**. **Do NOT `git pull` on TRUBA** — its working copy is manually patched (02p, parallel 02b/03) + TRUBA-specific `config.yaml`; pull would conflict. `main` already contains all the code fixes.

**Immediate next steps (new session):**
1. Check `run_ml.slurm` (`amr-ml`): `squeue -u $USER`; `tail $AMR_WORK/amr-ml-*.out`. 04 HPO (30 Optuna trials × out-of-core over 22 chunks × 50.8M features) is the long phase. On success: `config/experiments/ecoli/config_ampicillin.yaml` + model + `06_evaluation` metrics appear.
2. If `run_ml` is killed for low efficiency: 04 is **not** resume-safe (restarts HPO). Re-submit; if it keeps getting flagged, note that 04/05's per-chunk DMatrix *load* is serial between trees — acceptable for a one-off, or reduce scope.
3. After ML → **biology job** (`barbun -c 20`): `07b → 07 → 09 → 10`. `08` CARD-local BLAST works on compute nodes; `08` NCBI-remote + `09` Entrez need internet → run those on the **UI** (compute nodes may lack outbound internet) or skip.
4. Then download results: `tar` `results/ models/ runs/` from scratch → home or `rsync` to laptop (scratch auto-purges in 30 days).
5. Token: pushes use a fine-grained PAT pasted in chat (expires fast; ask for a fresh one with Contents:write). Full step-by-step + real-run corrections in `docs/TRUBA_Kurulum_ve_Calistirma_Rehberi.md` (Appendix).

---

# 1. Project Overview — what & why

- **Goal.** Predict antimicrobial resistance (AMR) in *E. coli* from whole-genome assemblies using **alignment-free k-mer features + out-of-core XGBoost**, then *reverse-translate* the most important k-mers back into biology (BLAST vs CARD/NCBI). The thesis-level goal is a **queryable, confidence-tiered, cross-antibiotic AMR k-mer Knowledge Base (AMRK-DB)**.
- **The real research gap (why this is novel).** Many papers predict AMR from k-mers; almost none turn the ML feature-importance output into a *reproducible, stability-filtered, biologically validated, discriminativeness-checked* knowledge base. The contribution chain is: **Gain → seed stability → BLAST confidence tier → discriminativeness (R vs S) → SNP-allele check → cross-antibiotic overlap → KB**. See `docs/ROADMAP.md` for the full thesis framing (hypotheses H1–H4, must-haves M1–M11).
- **Current stage.** Pipeline is implemented end-to-end (`00a → 11` + `07b`), multi-organism/antibiotic, audit-clean, **cross-environment (macOS/Linux/HPC)**, and **verified on the real 1788-genome dataset**. The KB persistence layer (SQLite/Postgres + API, M8/M10/M11) and the cross-antibiotic hypergeometric test (S1) are **not yet built** — those are the next big items.

---

# 2. Repository Structure

```
config/
  config.yaml                       # global config: organism, target_antibiotic, params, tiers, tool/data paths
  registry/organisms.yaml           # organism -> taxid, data paths, antibiotic set (single source of truth)
  registry/antibiotics.yaml         # antibiotic classes + ALIASES (name normalisation single source)
  experiments/{organism}/config_{antibiotic}.yaml   # AUTO-generated by step 04 (data split + best HPO params)
scripts/
  00a_download_bvbrc.py             # BV-BRC download (API/CLI/--raw-csv) + clean + parallel {id}.fna fetch
  00_prepare_metadata.py            # cleaned long table -> wide binary amr_phenotypes.csv (∩ present .fna)
  01_data_validation.py / 01b_*     # phenotype validation + EDA plots; ML-target recommendation
  02_kmer_extraction.py             # KMC k-mer counting (k=21) -> per-genome .kmc_pre/.kmc_suf
  02b_global_qc_analysis.py         # global QC: complexity outliers (IQR) + min_support elbow advisory
  03_matrix_construction.py         # global vocab (KMC) -> sparse binary CSR .npz chunks (+ y, genomes, features.txt)
  03u_unitig_matrix.py              # ROADMAP §0 M12: unitig-caller rtab -> SAME chunked matrix (-> matrix_unitig/); replaces raw k-mers downstream
  03b_matrix_validation_qc.py       # matrix QC (sparsity/prevalence)
  04_optimization.py                # Optuna HPO -> experiment config + run_metadata.json
  05_model_training.py              # full-data boosting over streaming QuantileDMatrix -> model + manifest.json + threshold
  lib/xgb_data.py                   # ChunkDMatrixIter + build_quantile_dmatrix (streaming DMatrix; used by 05/07b)
  06_evaluation.py                  # metrics, ROC/PR, calibration, bootstrap 95% CIs, error analysis
  07b_feature_stability.py          # 5-seed repeated holdout: AUC mean±std, selection freq, Jaccard, stable set
  07_explainability.py              # Gain top-N ∪ 07b stable set -> CSV + FASTA (flagged)
  08_blast_annotation.py / .nf      # Nextflow: CARD local + NCBI remote BLAST (blastn-short)
  09_biological_summary.py          # tiered report + KB candidates + recovery/composite/novel metrics
  10_kmer_background_frequency.py   # R-vs-S prevalence + Fisher exact -> discriminativeness
  11_variant_snp_check.py           # CARD variant-model SNP allele check (resistant vs wildtype)
  lib/                              # shared package (config, registry, chunking, io_utils, run_metadata, bvbrc)
  constants.py, utils.py            # thin backward-compat shims -> lib/ (kept for old imports/tests)
  run_pipeline.py                   # orchestrator: runs the numbered steps in order (subprocess + logging)
  lib/logging_utils.py              # standard timestamped logger factory (orchestrator + new code)
  migrate_to_organism_layout.py     # reversible data-layout migration (already applied)
tests/                              # pytest: smoke / unit / integration + README
docs/
  TECHNICAL_REVIEW.md               # consolidated audit findings + resolution status
  SCALE_MLOPS_PLAN.md               # multi-organism + KB + MLOps plan
  ROADMAP.md                        # thesis roadmap (hypotheses, M1-M11, 6-month plan)
data/  models/  results/  logs/  runs/   # generated — only the CARD homolog BLAST DB is version-controlled
# Research-software-engineering scaffolding:
LICENSE (MIT) · CITATION.cff · pyproject.toml (PEP 621 + ruff/mypy/pytest config)
.github/workflows/ci.yml (ruff + unit/smoke on py3.10-3.12) · .pre-commit-config.yaml
Makefile · CONTRIBUTING.md · CHANGELOG.md
requirements.txt, environment.yml, pytest.ini, README.md, QUICKSTART.md, METHODOLOGY.md
```

---

# 3. Pipeline — step by step (what it does, why, how)

Run order (config-driven; each reads `config.yaml` for organism/antibiotic):

```
00a  download + clean BV-BRC AMR  ─►  00  binary phenotype matrix
01   validate / pick target           02  KMC k-mers     02b  global QC
03   sparse matrix chunks             03b  matrix QC
04   Optuna HPO  ─►  05  train  ─►  06  evaluate
07b  5-seed stability  ─►  07  candidate k-mers (gain ∪ stable)
08   BLAST (CARD + NCBI)  ─►  09  tiered report + KB candidates + metrics
10   discriminativeness (R vs S)       11  variant-model SNP allele check
```

- **00a `download_bvbrc`** — *why:* get a reproducible, cleaned AMR label set + the matching assemblies. *How:* fetches the BV-BRC `genome_amr` table (`--backend api` default, `cli` via `p3-*`, or `--raw-csv` from the website), cleans it via `lib/bvbrc.py` (EUCAST/CLSI only, Lab-Method evidence, R/S→1/0, antibiotic-name normalisation, duplicate-conflict resolution), then downloads each surviving genome as `{genome_id}.fna` in parallel (retry/resume/`--retry-failed`/`--max-genomes`). Writes `amr_cleaned_long.csv`, `download_manifest.json`, logs + reports.
- **00 `prepare_metadata`** — *why:* the numbered pipeline needs a wide binary label matrix. *How:* pivots the cleaned long table to `amr_phenotypes.csv` (`Genome ID` + one 0/1 column per antibiotic, blank = untested), intersected with the `.fna` files actually present. Genomes with AMR labels but no downloadable assembly are dropped (a few % is normal, e.g. `562.1`).
- **01 / 01b `data_validation`** — class balance, missingness, EDA plots, and a scientific ML-target recommendation (minority count/ratio per antibiotic). Uses antibiotic classes from the registry.
- **02 `kmer_extraction`** — KMC counts canonical 21-mers per genome (`min_count=1`); outputs binary KMC DBs. Re-runnable (skips genomes already counted).
- **02b `global_qc_analysis`** — scans all KMC DBs for genome-complexity outliers (IQR on unique-k-mer count) and computes a `min_support` "elbow" advisory (does **not** change config).
- **03 `matrix_construction`** — builds the **global k-mer vocabulary** with one KMC pass over all genomes (`-ci min_support` rare filter, `-cx max_support` drops core-genome k-mers), dumps it to `features.txt`, then writes the genome×k-mer presence/absence matrix as **CSR `.npz` chunks** of `chunk_size=200` genomes + `y_{ab}.csv` + `genomes_{ab}.csv`.
- **04 `optimization`** — Optuna HPO (`n_trials=25`, `eval_metric=auc`). Stratified chunk split into train/test/optuna subsets. `colsample_bytree` searched on a **√p-anchored log range**; `n_estimators = best_iteration+1`; `base_score=0.5` pinned. Writes `config/experiments/{organism}/config_{antibiotic}.yaml` (the data split + best params) and `run_metadata.json`.
- **05 `model_training`** — **standard full-data gradient boosting** over a single, streaming **`ExtMemQuantileDMatrix`** (external memory): `lib/xgb_data.ChunkDMatrixIter` feeds the chunks to XGBoost one at a time and the quantised pages are **spilled to fast scratch** (`cache_prefix`), so the matrix never has to fit in RAM. (An in-core `QuantileDMatrix` of the full train set peaked **>400 GB** and OOM-killed the 384 GB node — external memory keeps RAM bounded to ~one page + histograms.) Trains the Optuna-tuned `n_estimators` trees on the whole training set (every tree sees all training rows), with a single **global `neg/pos` instance weight** for class imbalance, `base_score=0.5`. Saves the model + `manifest.json`, and an **operating threshold fixed at 0.5** (global weighting; no test-set tuning → no leakage). The disk cache (`models/.../_xgb_cache_train`) is removed after training. *Replaces the previous 1-tree-per-chunk incremental regime (weaker fit + very low HPC CPU efficiency).* **07b uses the same external-memory regime per seed.**
- **06 `evaluation`** — single stratified split; ROC/PR curves, calibration, confusion matrix, **bootstrap 95% CIs**, MCC/κ, error analysis. **Does not overwrite** the config threshold (leakage fix).
- **07b `feature_stability`** — runs **before** 07. 5 seeds `[42,123,777,1024,2025]`, stratified 80/20, **fixed HPO across seeds** (no per-seed retune → leakage-safe; "repeated holdout", Mahé 2018 resampling, not true k-fold because of the out-of-core constraint). Each seed trains with the **same full-data boosting regime as 05** — a streaming `QuantileDMatrix` built from the seed's train rows (sample-level `row_mask` over the chunks), so it stays out-of-core (one chunk in RAM at a time) while every tree sees the whole train split. Reports AUC mean±std, per-k-mer **selection_frequency** (`stable` if ≥ `stability_threshold=0.6`), and mean pairwise **Jaccard** of the top-N sets.
- **07 `explainability`** — extracts the single model's Gain top-N k-mers, **then merges in the 07b stable set** (k-mers reproducible across seeds but not in the gain top-N), flagging each row `in_gain_topN` / `stable` / `selection_frequency`. Emits the candidate CSV + FASTA. So BLAST/biology (08–11) covers **both** the gain and stability views.
- **08 `blast_annotation`** — Nextflow runs two parallel BLASTs of the candidate FASTA: **CARD local** and **NCBI nt remote**, both with `-task blastn-short -dust no` (correct for 21-mer queries). Records the CARD DB version (`blastdbcmd -info`).
- **09 `biological_summary`** — grades every hit into **confirmed / candidate / weak / none** using **identity + coverage** (alignment length / k) as the primary, DB-size-independent criteria, with E-value a loose secondary gate (E-value is *not* comparable between CARD and NCBI). Joins 07b stability, resolves NCBI gene names via Entrez, and writes: the Markdown report (with quantitative summary + reification-fallacy + gyrA/SNP caveats), `07_kb_candidates_{ab}.csv` (per-k-mer KB record incl. **composite_score = stability × log10(1/E) × identity**), and `08_validation_metrics_{ab}.json` (**M7 known-mechanism recovery rate**, **H2 pass/fail ≥40%**, **H4 novel-candidate fraction**, tier counts).
- **10 `kmer_background_frequency`** — *why:* BLAST says *which gene*, not *whether the k-mer discriminates*. Streams the matrix once and computes each candidate's prevalence in **resistant vs susceptible** genomes + **Fisher's exact** + a `discriminative` flag (|Δprev| ≥ 0.10 AND p < 0.05). Flags BLAST-confirmed-but-ubiquitous k-mers (likely wildtype/lineage). Output `10_kmer_background_frequency_{ab}.csv`.
- **11 `variant_snp_check`** — *why:* a homolog hit to gyrA only proves "gyrA region present", not "resistance SNP present". BLASTs candidates against CARD's **protein-variant-model** sequences, parses each model's resistance SNPs from `card.json` (protein position, e.g. `S83L`), maps protein pos → CDS codon, reads the k-mer's strand-aware codon, translates, and classifies **resistant_allele / wildtype / other_variant / ambiguous**. Needs the full CARD download; **skips cleanly** with instructions if absent. Output `11_variant_snp_check_{ab}.csv`.

---

# 4. Last real-data run (2026-06-15, ampicillin, 1788 genomes)

Proof the whole chain works on real data, not just synthetic smoke:

| Step | Result |
|---|---|
| 00a / 00 | 1788 genomes; ampicillin 758 R / 828 S / 202 untested |
| 02 KMC | 1788/1788, 0 failures |
| 03 matrix | 1552 genomes × **30,082,953** k-mers, 8 chunks, `features.txt` ≈ 757 MB |
| 04 HPO | best ROC-AUC **0.935** (10-trial test run; config now back to 25) |
| 06 eval | **test ROC-AUC 0.862** (CI 0.81–0.91), PR-AUC 0.90, MCC 0.64, acc 0.82 |
| 07b | AUC **0.750 ± 0.051**, Jaccard 0.23, **5 stable** k-mers |
| 07 | 10 gain ∪ 3 added stable = 13 candidates |
| 08 BLAST | CARD 440 hits (after the blastn-short fix; was 23), NCBI 1.4 MB |
| 09 | tiers: 2 confirmed / 4 weak / 7 none; recovery 20%, novel 20% |
| 10 | **9/13 discriminative**; APH(6)-Id real marker (prev_R 0.53 vs prev_S 0.11, p=4e-76); **OXA-1238 confirmed-but-ubiquitous** (0.63 vs 0.53) |
| 11 | 0 resistant-allele SNPs (expected — ampicillin = β-lactamase **acquisition**, a homolog mechanism, not a point mutation; the SNP machinery's payoff is on ciprofloxacin/gyrA) |

---

# 5. Cross-environment design (macOS / Linux / HPC)

The project must run on the user's Mac **and** a remote HPC pulled from GitHub. Key mechanisms:

- **PATH-aware tool resolution** — `lib.config.resolve_tool(config_key, command, …)` finds KMC/BLAST in this order: **(1)** env override `AMR_<TOOL>_BIN`, **(2)** `shutil.which` (conda/module on PATH), **(3)** the bundled macOS binary under `bin/bin/` **only on Darwin** (it is a Mach-O arm64 build that would mis-fire on Linux). Used by **02, 02b, 03** (`kmc`/`kmc_tools`) and **11** (`blastn`/`makeblastdb`); **08** prepends the conda bin to PATH so Nextflow's `blastn` resolves. *Why:* the old hardcoded `bin/bin/kmc` path broke every Linux/HPC run.
- **HTTPS / SSL** — `00a` builds an SSL context from `certifi` for all HTTPS (BV-BRC API + FTP); conda Python otherwise fails cert verification.
- **No CWD assumptions** — every script anchors paths to `PROJECT_ROOT = Path(__file__).resolve().parent.parent` and imports `lib` because Python puts the script's dir on `sys.path[0]`; running `python scripts/XX.py` works from any directory.
- **`.gitignore` + data hygiene** — all generated data is ignored (`data/raw`, `data/interim`, `data/processed`, `data/external/**/metadata`, the full CARD bundle `data/external/card`, `*.npz`, `features.txt`, models, results, logs, runs, nextflow `work/`). The **only** committed data is the CARD homolog BLAST DB (`data/external/blast_db/card_nt/card.*`, ~8.5 MB) so step 08 works out-of-the-box. A fresh clone reproduces everything by running `00a → …`.
- **Dependencies** — `environment.yml` (conda, installs KMC/BLAST/Nextflow via bioconda) or `requirements.txt` (pip Python deps; tools installed separately). `certifi` is listed. `QUICKSTART.md` has the full fresh-machine + HPC (SLURM) setup.

---

# 6. Key config knobs (`config/config.yaml`)

- `project.organism` = `ecoli`, `project.target_antibiotic` = **`ampicillin`** (was gentamicin; changed because gentamicin is all-susceptible in the current sample).
- `preprocessing`: `k_length=21`, `min_support=5`, `chunk_size=200`, `kmc_mem`, `threads`.
- `training.n_trials=25` (was temporarily 10 for a fast test run; restored).
- `analysis`: `top_n_features=50`, `stability_threshold=0.6`, and **`confidence_tiers`** (identity+coverage+evalue per tier) + `report_max_evalue`.
- `blast`: `card_db_dir/name`, `evalue`, `word_size`, `threads`, and **`card_variant_fasta` / `card_json`** (step 11; full CARD download).
- `ncbi.entrez_email` — empty in config; set it before step 09 to avoid NCBI rate-limit warnings, **or** `export AMR_ENTREZ_EMAIL=…` (env override).
- Tool overrides via env: `AMR_KMC_BIN`, `AMR_KMC_TOOLS_BIN`, `AMR_BLASTN_BIN`, etc.

---

# 7. Important decisions (rationale)

1. **Conflict resolution** for duplicate (genome, antibiotic): majority vote → tie: newest `testing_standard_year` → still tied: drop (NaN) + log.
2. **Antibiotic normalisation** via registry aliases (`co-trimoxazole`→`trimethoprim/sulfamethoxazole`; canonical = registry spelling).
3. **Cleaning filters:** EUCAST/CLSI standards, Lab-Method evidence, R/S phenotypes only.
4. **HPO once, fixed across seeds** in 07b (leakage-safe repeated holdout, not k-fold — out-of-core constraint).
5. **Confidence tiers grade on identity + coverage**, not raw E-value, because E-value depends on DB size and is not comparable across CARD vs NCBI. **Weak hits are kept and flagged**, never silently dropped (so e.g. a gyrA-type partial hit stays visible).
6. **07 carries both** the gain top-N and the 07b stable set forward so biology covers both.
7. **Discriminativeness (step 10) is separate from BLAST** — a k-mer can hit a known ARG yet be non-discriminative (ubiquitous); both facts are recorded.
8. **Organism-scoping** of every path incl. the auto-generated experiment config.
9. **Generated outputs are not version-controlled**; only the CARD homolog DB is.
10. **FASTA filenames must be `{genome_id}.fna`** (pipeline globs `*.fna`, uses the stem as Genome ID).

---

# 8. Audit findings — all resolved (see `docs/TECHNICAL_REVIEW.md`)

- **P-01** data leakage (Youden's J fit on test set) — FIXED (threshold on train/val only; 06 doesn't overwrite config).
- **P-02** shell injection (`shell=True`) — FIXED (`lib.io_utils.run_command`, shlex, never shell=True).
- **P-03** `eval_metric` hardcode — FIXED (from config).
- **P-04** colsample √p mismatch — FIXED (dynamic √p-anchored range).
- R/S double-count (01), hardcoded `top_50` filenames (08/09), O(N²) Gram SVD (03b), Entrez email/api_key, PR-AUC inconsistency, KMC resume, NCBI errorStrategy, CARD-version record, duplicated code → `lib/` — ALL FIXED.
- Integration-test-caught runtime bugs: `base_score must be in (0,1)`, `n_estimators=0`, empty feature-importance crash, empty-stability table — ALL FIXED.

---

# 9. Recent work log (newest first)

- **2026-06-24 Block 1 biology made reproducible + KB populated (M8/M16/M6) — see §0.1 "Block 1 + KB DONE".** Commits `1e09219`→`853143f` on `main`. Fixed NCBI remote BLAST (SIGXCPU): decoupled the remote pass from CARD → `blastn`/word11 + taxid `-entrez_query` + `-max_target_seqs`; taxid (not scientific name — space breaks Nextflow launcher); `NXF_ANSI_LOG=false` so backgrounded Nextflow doesn't stall (SIGTTOU); `AMR_ENTREZ_EMAIL`/`AMR_CARD_VERSION` env overrides. Re-ran 08 (CARD 3605 + NCBI 4522, E. coli) → 09 → `populate_database.py` → `amrk.db` schema **0.2.0** (M8). Added ARO ontology cols to `blast_annotations` (**M16**, 13/60 mapped) + recorded CARD 4.0.1 (**M6**). H2 still FALSE (recovery 32%) — fine per §0.4. **Next: Block 2, starting M9 MDA permutation.**
- **2026-06-23 FULL unitig pivot executed end-to-end (ROADMAP §0) — see §0.1/§0.2.** Unitig matrix (4.94M unitigs) + PopPUNK dbscan lineages (324) + lineage-aware 07b → **honest headline ROC-AUC 0.9505 ± 0.0102** (5-fold GroupKFold; ≈ chunk-split 0.9534 → near-zero lineage leakage; beats k-mer 0.930). 04 HPO cut short at trial 35 → experiment config written by hand from trial-20 params (`n_estimators: 8`). Biology fixes: Nextflow 26 strict-parser (`def`→`params.outfmt`), BLAST task by median length + **word_size 7** → candidate unitigs hit **TEM/CTX-M/OXA/CMY** β-lactamases at cov=1. Added 09 ARO mapping (M16), 10 BH-FDR, S10 reification, env overrides (`AMR_FEATURE_REPR`/`AMR_EXTERNAL_MEMORY`/`AMR_OPTUNA_PATIENCE`), 4 containers (amr/amr-pp/amr-tools/amr-checkm2). Commits `0909589`→`400657f` on `main`. **IN PROGRESS:** Block 1 biology (08/09 on UI for NCBI, 10/11 on compute).
- **2026-06-22 unitig pipeline — step 1 started (ROADMAP §0 / IMMEDIATE NEXT #1):** added `bcalm`+`unitig-caller` to `environment.yml`; wrote `03u_unitig_matrix.py`; added the config `unitig:` section. (Superseded by the 2026-06-23 entry above.)
- **2026-06-22 literature review (2 rounds) → ROADMAP §0 methodological pivot:** systematic review (Sections A–F + an implementation "how" round) concluded that the raw-k-mer + random-CV approach, while a working baseline, is not publication-grade. Binding decisions written to `docs/ROADMAP.md` **§0** (and must-have table M2/M4 revised + M12–M16 added): switch to **unitigs** (bcalm2/unitig-caller), **lineage-aware CV** (PopPUNK+GroupKFold), **CPSS stability selection + SHAP**, add **external validation + AMRFinderPlus/ResFinder concordance**, **pyseer LMM**, **CheckM2/QUAST QC**, **BH-FDR (step 10)**, **ARO/CARD ontology mapping**, reification-safe wording; novelty reframed away from "first ML AMR DB". Confirmed unchanged: k=21, binary+max_bin=2, class-weight/no-SMOTE, BV-BRC, AUC~0.93. **This pivot makes the day's raw-k-mer/min_support/GPU work largely moot (unitigs dissolve the speed/memory problem).** Next session implements §0.
- **2026-06-22 KB layer started (M8)** (`d49377f`, `54178a1`): `scripts/lib/kb_schema.py` (SQLite DDL, 11 tables per ROADMAP §1.1, stdlib `sqlite3` — no new dep / no container rebuild) + `scripts/populate_database.py` (loads run_metadata, manifest, 06 metrics, 07b holdout, 09/10 candidate+background, 11 SNP → `results/{org}/kb/amrk.db`; idempotent, multi-antibiotic, graceful on missing inputs; writes `validation_evidence` per result (M11) + `kb_metadata` FAIR row incl. CARD 4.0.1). Functionally tested on synthetic inputs. **Next:** run the pipeline to produce real outputs, then `populate_database.py` for real + validate column mapping; then FastAPI endpoints (S8) + cross-antibiotic overlap (S1) + permutation test (M9). **Also done 2026-06-22:** data-adaptive min_support (03), GPU evaluated+rejected (V100 16 GB), TRUBA scratch cleaned (~2.4 GB), full CARD 4.0.1 downloaded (step 11 now active for all antibiotics).
- **2026-06-22 prediction-engine consolidation/audit:** after many incremental patches, audited 04/05/07b/lib for consistency, reproducibility, no-hardcoding. Fixes: seeded the Optuna TPE sampler (`random_seed`) for reproducible HPO; moved the 05 early-stopping val fraction to config (`training.validation_fraction`, was hardcoded 0.15); 07b now respects `training.external_memory` (in-core vs ExtMem) like 05; manifest records the actual `n_trees` + corrected stale `threshold_type` (global neg/pos weight, not per-chunk "Dynamic Instance Weighting"); fixed stale "incremental/epoch" docstrings in 05/07b. No hardcoded abs paths remain in the engine. 57 unit/smoke + 1 integration pass. Config knobs now: `optuna_threads_per_trial`, `max_boost_rounds`, `external_memory`, `max_train_chunks` (all documented in `config.yaml`).
- **2026-06-18→22 full-data boosting refactor** (`a5b9ddc`→`7e5c22b`): replaced incremental 1-tree/chunk with standard boosting over a streaming `(Ext)QuantileDMatrix` (`lib/xgb_data.py`); parallel Optuna trials (04); ExtMem to survive the >400 GB in-core peak; early stopping in 05. Result on 4373 genomes: ROC-AUC 0.930 / MCC 0.739. Early stopping disproved the underfit hypothesis (best ≈ 29–30 trees) → metrics are at the k-mer signal ceiling; pivot to the KB contribution (see §0).
- **2026-06-15 RSE / open-science scaffolding:** added `LICENSE` (MIT), `CITATION.cff`, `pyproject.toml` (packaging metadata + ruff/mypy/pytest config), GitHub Actions CI (`ruff` + unit/smoke on py3.10–3.12), `.pre-commit-config.yaml`, `CONTRIBUTING.md`, `CHANGELOG.md`, a `Makefile`, a `run_pipeline.py` orchestrator, and `lib/logging_utils.py`. Added type hints to `lib/config.py`. Untracked the auto-generated experiment configs (now gitignored). Smoke test now covers steps 10/11 → **56 tests pass**.
- **2026-06-15 cross-env cleanup** (`647c0d9`): `resolve_tool` adopted in 03 + 02b (was hardcoded `bin/bin/kmc`); `.gitignore` hardened + generated artifacts untracked (only CARD homolog DB kept); README pipeline list updated to `00a→11`. Verified: all scripts py_compile (3.10), 54 unit/smoke + 1 integration test pass, no hardcoded/CWD-dependent paths.
- **2026-06-15 step 11** (`306a18f`): `11_variant_snp_check.py` — CARD variant-model SNP allele check (pure helpers unit-tested incl. ± strand; fixed a `card.json` parse bug where `param_value` entries are bare strings).
- **2026-06-15 step 10** (`541f90d`): `10_kmer_background_frequency.py` — R-vs-S prevalence + Fisher + discriminativeness.
- **2026-06-15 BLAST/tier fix** (`e7e2e62`): CARD `blastn-short -dust no` (CARD hits 23→440); tiers → identity+coverage; hits sorted by E-value; gyrA/SNP + cross-DB caveats.
- **2026-06-14 biology pass** (`84dc94d`): M7 recovery rate, composite score, H4 novel fraction in 09; tiers moved to config; 07 merges 07b stable set; order 07b→07.
- **2026-06-13 portability + BV-BRC fixes** (`64b0fab`, `ab68d02`): `resolve_tool`, certifi SSL, API URL-encode (HTTP 400 fix), dry-run AMR sampling, batched CLI fetch, `base_score=0.5` in 04, rewritten QUICKSTART.

---

# 10. Known issues / technical debt

- **`config ncbi.entrez_email` is empty** → step 09 Entrez warns / may rate-limit (does NOT crash — falls back to stitle parsing). Set it via `config.yaml` **or** `export AMR_ENTREZ_EMAIL=…` (env override, added 2026-06-24; no config edit needed on HPC). Same for `AMR_ENTREZ_API_KEY`.
- **Step 11 needs the full CARD download** (`data/external/card/card.json` + `nucleotide_fasta_protein_variant_model.fasta`); not shipped (ignored). Skips cleanly if absent.
- **BV-BRC API deep-pagination cap** on the full 243k-row table — `--backend api` works and is fast, but if truncated fall back to `--backend cli` (batched, slow) or website `--raw-csv`.
- **No pipeline orchestrator** (`run_pipeline.py`) — steps run manually.
- ~~**04 vs 05 training regimes differ** (single early-stopped fit vs 1-tree/chunk incremental)~~ — **RESOLVED 2026-06-18**: 05 (and 07b) now use the same standard full-data boosting as 04's HPO (streaming `QuantileDMatrix`, `lib/xgb_data.py`). 04 still tunes `n_estimators` on a representative chunk subset; 05 trains that many trees on the full set. Minor remaining nuance: the budget is tuned on a subset, not the full train set (acceptable, keeps HPO fast).
- **BV-BRC env footgun:** `source /Applications/BV-BRC.app/user-env.sh` shadows `python` with BV-BRC's Python 2.7 → run scripts with `python3` or the explicit conda path while it's sourced (so `p3-*` stay on PATH).
- **Only ampicillin** has been run on real data so far; cefotaxime/ciprofloxacin/gentamicin matrices not regenerated at full scale.

---

# 11. Next priority tasks

**High priority (thesis novelty)**
- `10_cross_antibiotic_analysis.py` + **hypergeometric / Fisher test** (S1 / H3) — uses 07b stable-k-mer sets across antibiotics (β-lactam-internal vs cross-class overlap).
- **Permutation test** (M9): null Gain distribution from shuffled labels vs the real model.
- Run the full chain on **ciprofloxacin** to exercise step 11 on a genuine SNP mechanism (gyrA/parC) and validate `resistant_allele` calls.

**KB layer (M8/M10/M11)**
- `populate_database.py` → SQLite (then Postgres) with tables incl. `kmer_background_frequency`, `validation_evidence`, `pipeline_runs`; `kb_schema_version` + Zenodo DOI; minimal FastAPI endpoints (`/kmers`, `/overlap`, `/metadata`).

**Infra / optional**
- `run_pipeline.py` orchestrator + Makefile; temporal validation (S2); MLST/phylogenetic bias flagging (S3).

---

# 12. Working conventions & notes for the next session

- **Branch** `fix/amr-audit-remediation` (NOT main); **Conventional Commits**. The user has asked that pushes be attributed to them only — **do not add a `Co-Authored-By: Claude` trailer** to commits that will be pushed.
- **Push** only when asked; uses a GitHub PAT pasted in chat (fine-grained, `Contents: Read and write`). Remind the user to revoke/rotate exposed tokens.
- **Verification discipline:** the assistant sandbox is Python 3.13 but the user runs **3.10** — always `py_compile` with `/opt/anaconda3/envs/bitirme_vol2/bin/python`. The user runs xgboost/KMC/BLAST/BV-BRC locally; **`pytest` (54 tests) + `pytest -m integration` (synthetic 02→07b) are the primary validators** and have caught several real runtime bugs.
- **Single sources of truth:** `lib/` (helpers), `config/registry/` (organisms + antibiotics + aliases), `config/config.yaml` (params/tiers/paths).
- **Naming:** organism slug `ecoli` ↔ taxid 562; antibiotic ids lowercase with `/` preserved; genome FASTA `{genome_id}.fna`; experiment config `config/experiments/{organism}/config_{antibiotic}.yaml`; run_id `{org}__{ab}__{UTC}__{git7}`.
- **Do not rewrite working code.** Keep changes incremental, organism-scoped, cross-environment, and tested. Read `docs/ROADMAP.md`, `docs/TECHNICAL_REVIEW.md`, `config/registry/*`, `lib/bvbrc.py`, and `lib/config.py` (esp. `resolve_path` / `resolve_tool`) before starting.
