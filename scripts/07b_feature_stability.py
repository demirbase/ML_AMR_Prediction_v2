#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Step 07b — Feature stability via 5-seed repeated holdout.

PURPOSE (to be implemented):
    Quantify how reproducible the top k-mer features are across resampling, and
    provide an out-of-core-friendly reliability estimate for the reported
    ROC-AUC. This is the methodological backbone of the Knowledge Base claim
    (ROADMAP §1.5 / §4; presentation slides "Stability").

DESIGN (agreed — content to be written later):
    SEEDS = [42, 123, 777, 1024, 2025]                      # ROADMAP §4 / sunum.js
    For each seed:
        - stratified 80/20 split at the GENOME (sample) level (StratifiedShuffleSplit)
        - train XGBoost with the FIXED hyperparameters from
          config/experiments/{organism}/config_{antibiotic}.yaml
          (HPO is done ONCE in step 04 and held fixed across seeds — it is NOT
           re-tuned per seed and never sees the seed's test split: this keeps the
           repeated-holdout estimate honest / leakage-free)
        - evaluate ROC-AUC on that seed's held-out 20%
        - extract the top-N (config: analysis.top_n_features) Gain k-mers

    Aggregate and report:
        - ROC-AUC mean ± std across the 5 seeds
        - per-k-mer selection_frequency = (#seeds in top-N) / len(SEEDS);
          "stable" k-mers = selection_frequency >= 0.6  (H1 acceptance criterion)
        - mean pairwise Jaccard similarity of the 5 top-N k-mer sets

    Outputs (organism/antibiotic + seed scoped):
        models/{organism}/{antibiotic}/seed{S}/xgboost_{antibiotic}_seed{S}.json
        results/{organism}/{antibiotic}/04_evaluation/
            10_repeated_holdout_summary_{antibiotic}.csv   # per-seed AUC + mean/std + mean Jaccard
        results/{organism}/{antibiotic}/05_explainability/
            06_feature_stability_{antibiotic}.csv          # k-mer, selection_frequency, mean_gain

OPEN IMPLEMENTATION DECISION (to confirm before writing content):
    Memory model for the sample-level split. The current matrix is chunked
    (.npz parts) and the existing split (step 04) is at the chunk level. A fresh
    80/20 split at the sample level requires either:
      (A) vstack all chunks once into RAM, then slice rows per seed
          (simplest; same approach 04/06 already use on subsets; fine at the
           E. coli scale of ~5k genomes), or
      (B) keep it out-of-core by mapping global sample indices back to chunk
          offsets and streaming (lower peak RAM, more code).
    Default recommendation: (A).

NOTE: This is a SKELETON. Implementation pending confirmation of the memory
model above.
"""

import sys


def main() -> None:
    """TODO: implement 5-seed repeated holdout + stability + Jaccard."""
    print("07b_feature_stability.py is a placeholder — implementation pending.")
    print("Design: 5-seed repeated holdout (SEEDS=[42,123,777,1024,2025]),")
    print("ROC-AUC mean±std, selection_frequency>=0.6, mean pairwise Jaccard.")
    sys.exit(0)


if __name__ == "__main__":
    main()
