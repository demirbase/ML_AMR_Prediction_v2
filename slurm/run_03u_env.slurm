#!/bin/bash
#SBATCH -p barbun -N1 -c40 --mem=300G --time=12:00:00
#SBATCH --exclude=barbun45,barbun46 -o %x-%j.out -e %x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -euo pipefail
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig
: "${AMR_ORGANISM:?}" "${AMR_ANTIBIOTIC:?}"
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2
cd $H
apptainer exec --no-home $W/containers/amr.sif python -u scripts/03u_unitig_matrix.py --threads 40 2>&1
