#!/bin/bash
#SBATCH -p barbun -N1 -c20 --mem=180G --time=1-00:00:00 --exclude=barbun45,barbun46
#SBATCH -o /arf/scratch/edemirbas/amr/%x-%j.out -e /arf/scratch/edemirbas/amr/%x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -euo pipefail
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig AMR_EXTERNAL_MEMORY=false AMR_OPTUNA_PATIENCE=15
: "${AMR_ORGANISM:?}" "${AMR_ANTIBIOTIC:?}"
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2; SIF=$W/containers/amr.sif
cd $H
echo "===== 03u ($AMR_ORGANISM/$AMR_ANTIBIOTIC) ====="
apptainer exec --no-home $SIF python -u scripts/03u_unitig_matrix.py --threads 20 2>&1
for s in 04_optimization 05_model_training 06_evaluation 07b_feature_stability; do
  echo "===== $s ($AMR_ORGANISM/$AMR_ANTIBIOTIC) ====="
  apptainer exec --no-home $SIF python -u scripts/$s.py 2>&1
done
echo "===== 03u+ML ($AMR_ORGANISM/$AMR_ANTIBIOTIC) BITTI ====="
