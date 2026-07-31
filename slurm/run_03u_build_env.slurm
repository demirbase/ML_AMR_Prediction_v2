#!/bin/bash
#SBATCH -p barbun -N1 -c40 --mem=300G --time=1-00:00:00 --exclude=barbun45,barbun46
#SBATCH -o /arf/scratch/edemirbas/amr/%x-%j.out -e /arf/scratch/edemirbas/amr/%x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -euo pipefail
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig
: "${AMR_ORGANISM:?}"
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2
cd $H
apptainer exec --no-home $W/containers/amr.sif python -u scripts/03u_unitig_matrix.py --build-db --organism "$AMR_ORGANISM" --threads 40 2>&1
