#!/bin/bash
#SBATCH -p barbun -N1 -c40 --mem=300G --time=1-00:00:00
#SBATCH --exclude=barbun45,barbun46 -o %x-%j.out -e %x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -euo pipefail
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig
: "${AMR_ORGANISM:?}" "${AMR_ANTIBIOTIC:?}"
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2; SIF=$W/containers/amr.sif
cd $H
apptainer exec --no-home $SIF python -u scripts/10_kmer_background_frequency.py 2>&1
apptainer exec --no-home $SIF python -u scripts/11_variant_snp_check.py 2>&1
apptainer exec --no-home $SIF python -u scripts/12_permutation_test.py --n-permutations 100 2>&1
apptainer exec --no-home $SIF python -u scripts/12b_label_permutation_test.py --n-permutations 50 2>&1
apptainer exec --no-home $SIF python -u scripts/13_stability_selection.py --n-candidates 5000 --B 100 --pi 0.6 --base-trees 10 2>&1
apptainer exec --no-home $SIF python -u scripts/13b_stable_annotation.py 2>&1
