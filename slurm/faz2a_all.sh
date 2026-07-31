#!/bin/bash
export AMR_HOME=/arf/home/edemirbas/ML_AMR_Prediction_v2 AMR_WORK=/arf/scratch/edemirbas/amr
export APPTAINER_BINDPATH=/arf AMR_FEATURE_REPR=unitig
cd $AMR_HOME; SIF=$AMR_WORK/containers/amr.sif; LOG=$AMR_WORK/faz2a_all.log
declare -A PANEL=(
  [ecoli]="ampicillin amoxicillin_clavulanic_acid cefotaxime ciprofloxacin gentamicin trimethoprim_sulfamethoxazole tetracycline chloramphenicol"
  [kpneumoniae]="gentamicin meropenem imipenem cefoxitin ceftazidime ciprofloxacin trimethoprim_sulfamethoxazole tetracycline tigecycline aztreonam colistin"
  [staphylococcus_aureus]="cefoxitin oxacillin ciprofloxacin erythromycin clindamycin gentamicin tetracycline trimethoprim_sulfamethoxazole"
  [acinetobacter_baumannii]="imipenem meropenem amikacin ampicillin_sulbactam ceftazidime ciprofloxacin tetracycline trimethoprim_sulfamethoxazole"
  [pseudomonas_aeruginosa]="meropenem ceftazidime ciprofloxacin amikacin"
  [enterococcus_faecium]="vancomycin teicoplanin gentamicin ampicillin clindamycin tetracycline"
)
echo "FAZ2A START $(date)" > $LOG
for o in "${!PANEL[@]}"; do for a in ${PANEL[$o]}; do
  echo "===== $o/$a $(date) =====" >> $LOG
  if AMR_ORGANISM=$o AMR_ANTIBIOTIC=$a apptainer exec --no-home $SIF python -u scripts/07_explainability.py  >>$LOG 2>&1 \
  && AMR_ORGANISM=$o AMR_ANTIBIOTIC=$a apptainer exec --no-home $SIF python -u scripts/08_blast_annotation.py >>$LOG 2>&1 \
  && AMR_ORGANISM=$o AMR_ANTIBIOTIC=$a apptainer exec --no-home $SIF python -u scripts/09_biological_summary.py >>$LOG 2>&1
  then echo "OK   $o/$a" >> $LOG; else echo "FAIL $o/$a (rc=$?)" >> $LOG; fi
done; done
echo "FAZ2A DONE $(date)" >> $LOG
