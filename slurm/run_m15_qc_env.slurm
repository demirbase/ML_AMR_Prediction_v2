#!/bin/bash
#SBATCH -p barbun -N1 -c40 --mem=300G --time=1-00:00:00 --exclude=barbun45,barbun46
#SBATCH -o /arf/scratch/edemirbas/amr/%x-%j.out -e /arf/scratch/edemirbas/amr/%x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -uo pipefail
export APPTAINER_BINDPATH=/arf
: "${AMR_ORGANISM:?}"
H=/arf/home/edemirbas/ML_AMR_Prediction_v2
W=/arf/scratch/edemirbas/amr
AMR=$W/containers/amr.sif; CHECKM2=$W/containers/amr-checkm2.sif; TOOLS=$W/containers/amr-tools.sif
DMND=$W/data/external/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd
export TMPDIR=$W/tmp_m15_$AMR_ORGANISM; mkdir -p "$TMPDIR"
cd $H
echo "== prep ($AMR_ORGANISM) =="
apptainer exec --no-home $AMR python -u scripts/02d_genome_qc.py --mode prep --organism "$AMR_ORGANISM" 2>&1
source $(find $H/results/$AMR_ORGANISM -name "02d_qc_paths_$AMR_ORGANISM.sh" | head -1)
echo "GENOMES=$GENOMES_DIR | CHECKM2_OUT=$CHECKM2_OUT | QUAST_OUT=$QUAST_OUT | N=$N_GENOMES"
echo "== CheckM2 (long) =="
apptainer exec --no-home $CHECKM2 checkm2 predict --input "$GENOMES_DIR" -x fna \
  --output-directory "$CHECKM2_OUT" --threads 40 --database_path "$DMND" --tmpdir "$TMPDIR" --force 2>&1 | tail -25
echo "== QUAST =="
apptainer exec --no-home $TOOLS bash -c "cd '$GENOMES_DIR' && quast.py *.fna -o '$QUAST_OUT' --threads 40 --no-plots --no-html --no-icarus" 2>&1 | tail -15
echo "== post ($AMR_ORGANISM) =="
apptainer exec --no-home $AMR python -u scripts/02d_genome_qc.py --mode post --organism "$AMR_ORGANISM" 2>&1
echo "===== M15 ($AMR_ORGANISM) BITTI ====="
