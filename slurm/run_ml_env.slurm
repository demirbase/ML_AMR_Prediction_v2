#!/bin/bash
#SBATCH -p barbun -N1 -c40 --mem=300G --time=1-00:00:00
#SBATCH --exclude=barbun45,barbun46 -o %x-%j.out -e %x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -euo pipefail
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig AMR_EXTERNAL_MEMORY=false AMR_OPTUNA_PATIENCE=15
: "${AMR_ORGANISM:?}" "${AMR_ANTIBIOTIC:?}"
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2; SIF=$W/containers/amr.sif
cd $H
for s in 04_optimization 05_model_training 06_evaluation 07b_feature_stability; do
  apptainer exec --no-home $SIF python -u scripts/$s.py 2>&1
done
