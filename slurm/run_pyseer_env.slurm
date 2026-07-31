#!/bin/bash
#SBATCH -p barbun -N1 -c20 --mem=80G --time=06:00:00
#SBATCH --exclude=barbun45,barbun46 -o %x-%j.out -e %x-%j.err
#SBATCH --mail-type=END,FAIL --mail-user=eren0demirbas@gmail.com
set -euo pipefail
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig
: "${AMR_ORGANISM:?}" "${AMR_ANTIBIOTIC:?}"
W=/arf/scratch/edemirbas/amr; H=/arf/home/edemirbas/ML_AMR_Prediction_v2
AMR=$W/containers/amr.sif; TOOLS=$W/containers/amr-tools.sif
cd $H
echo "== prep =="; apptainer exec $AMR python -u scripts/14_pyseer_lmm.py --mode prep
P=$H/results/$AMR_ORGANISM/$AMR_ANTIBIOTIC/05_explainability/14_pyseer_paths_$AMR_ANTIBIOTIC.sh
source "$P"; D=$(dirname "$ASSOC")
echo "== cand + sub rtab =="
apptainer exec $AMR python -c "import pandas as pd;d=pd.read_csv('$D/13_stability_selection_${AMR_ANTIBIOTIC}.csv');open('$D/14_cand_seqs.txt','w').write('\n'.join(d.kmer.astype(str))+'\n')"
awk -F'\t' -v C="$D/14_cand_rtab.rtab" -v S="$D/14_sub_rtab.rtab" \
  'NR==FNR{seq[$1]=1;next} FNR==1{print>C;print>S;next} ($1 in seq){print>C} (FNR%100==1){print>S}' \
  "$D/14_cand_seqs.txt" "$PRES"
echo "== similarity =="; [ -s "$SIM" ] || apptainer exec $TOOLS similarity_pyseer "$SAMPLES" --pres "$D/14_sub_rtab.rtab" > "$SIM"
echo "== pyseer lmm =="; [ -s "$ASSOC" ] || apptainer exec $TOOLS pyseer --lmm --phenotypes "$PHENO" --pres "$D/14_cand_rtab.rtab" --similarity "$SIM" --output-patterns "$PATTERNS" --cpu 20 > "$ASSOC"
echo "== post =="
PYSEER_VER=$(apptainer exec $TOOLS pyseer --version 2>&1 | head -1)
AMR_PYSEER_VERSION="$PYSEER_VER" apptainer exec $AMR python -u scripts/14_pyseer_lmm.py --mode post
echo "===== 14 ($AMR_ORGANISM/$AMR_ANTIBIOTIC) BITTI ====="
