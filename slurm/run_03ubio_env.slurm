#!/bin/bash
#SBATCH -p barbun -N1 -c20 --mem=180G --time=1-00:00:00 --exclude=barbun45,barbun46
#SBATCH -o /arf/scratch/edemirbas/amr/%x-%j.out -e /arf/scratch/edemirbas/amr/%x-%j.err
set -euo pipefail
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig
: "${AMR_ORGANISM:?}" "${AMR_ANTIBIOTIC:?}"
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2; SIF=$W/containers/amr.sif
o=$AMR_ORGANISM; a=$AMR_ANTIBIOTIC; cd $H
echo "== 03u(rtab) $o/$a =="; apptainer exec --no-home $SIF python -u scripts/03u_unitig_matrix.py --threads 20 2>&1 | tail -3
if [ -s "$H/results/$o/$a/05_explainability/13_stability_selection_$a.csv" ]; then echo "== bio zaten var, atla =="; else
  echo "== bio 10-13b =="
  apptainer exec --no-home $SIF python -u scripts/10_kmer_background_frequency.py 2>&1 | tail -1
  apptainer exec --no-home $SIF python -u scripts/11_variant_snp_check.py 2>&1 | tail -1
  apptainer exec --no-home $SIF python -u scripts/12_permutation_test.py --n-permutations 100 2>&1 | tail -1
  apptainer exec --no-home $SIF python -u scripts/12b_label_permutation_test.py --n-permutations 50 2>&1 | tail -1
  apptainer exec --no-home $SIF python -u scripts/13_stability_selection.py --n-candidates 5000 --B 100 --pi 0.6 --base-trees 10 2>&1 | tail -1
  apptainer exec --no-home $SIF python -u scripts/13b_stable_annotation.py 2>&1 | tail -1
fi
echo "== 03ubio ($o/$a) BITTI =="
