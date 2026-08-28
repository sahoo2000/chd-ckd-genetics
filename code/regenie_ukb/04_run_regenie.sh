#!/usr/bin/env bash
# 04_run_regenie.sh — step 1 (whole-genome model) then step 2 (burden), quantitative + binary.
#
# The primary analysis is the POOLED gene set against eGFR. Everything else is secondary.
set -euo pipefail

# ---------------------------------------------------------------- configure
GENO_ARRAY=${GENO_ARRAY:-ukb_cal_allChrs}          # array genotypes for step 1 (bed/bim/fam prefix)
EXOME_BGEN=${EXOME_BGEN:-ukb_exome.bgen}           # exome dosages for step 2
EXOME_SAMPLE=${EXOME_SAMPLE:-ukb_exome.sample}
PHENO_Q=${PHENO_Q:-../phenotypes/pheno_quant.tsv}
PHENO_B=${PHENO_B:-../phenotypes/pheno_binary.tsv}
COVAR=${COVAR:-../phenotypes/covariates.tsv}
KEEP=${KEEP:-../phenotypes/keep.tsv}
SETS=${SETS:-sets_tier1}                            # output dir from 03_build_regenie_sets.py
THREADS=${THREADS:-16}
OUT=${OUT:-../results}

mkdir -p "$OUT"
CATCOV="centre,batch"
QCOV="age,age2,sex,age_sex,$(seq -s, -f 'PC%g' 1 20)"

# ---------------------------------------------------------------- step 1
# Ridge prediction on array genotypes. --lowmem is essential at UKB scale.
step1 () {
  local trait_flag=$1 pheno=$2 tag=$3
  regenie --step 1 \
    --bed "$GENO_ARRAY" \
    --keep "$KEEP" \
    --phenoFile "$pheno" \
    --covarFile "$COVAR" \
    --covarColList "$QCOV" \
    --catCovarList "$CATCOV" \
    $trait_flag \
    --bsize 1000 --lowmem --lowmem-prefix "$OUT/tmp_$tag" \
    --threads "$THREADS" \
    --out "$OUT/step1_$tag"
}

# ---------------------------------------------------------------- step 2
step2 () {
  local trait_flag=$1 pheno=$2 tag=$3 setlist=$4 label=$5 extra=${6:-}
  regenie --step 2 \
    --bgen "$EXOME_BGEN" --sample "$EXOME_SAMPLE" \
    --keep "$KEEP" \
    --phenoFile "$pheno" \
    --covarFile "$COVAR" \
    --covarColList "$QCOV" \
    --catCovarList "$CATCOV" \
    $trait_flag \
    --pred "$OUT/step1_${tag}_pred.list" \
    --anno-file  "$SETS/annotations.txt" \
    --set-list   "$setlist" \
    --mask-def   "$SETS/masks.txt" \
    --aaf-bins 0.0001,0.001,0.01 \
    --vc-tests skato,acatv,acato \
    --build-mask max --write-mask-snplist \
    --minMAC 3 --bsize 200 --threads "$THREADS" \
    $extra \
    --out "$OUT/burden_${label}"
}

echo "=============================================================="
echo " PRIMARY  — quantitative eGFR, pooled gene set (1 test)"
echo "=============================================================="
step1 "--qt --apply-rint" "$PHENO_Q" qt
step2 "--qt --apply-rint" "$PHENO_Q" qt "$SETS/setlist_pooled.txt"  egfr_pooled

echo "=============================================================="
echo " SECONDARY — quantitative eGFR, per gene"
echo "=============================================================="
step2 "--qt --apply-rint" "$PHENO_Q" qt "$SETS/setlist_pergene.txt" egfr_pergene

echo "=============================================================="
echo " SECONDARY — binary CHD, pooled gene set"
echo "  --firth --approx is REQUIRED: with ~3,400 cases and rare masks,"
echo "  unpenalised logistic regression produces separation artefacts"
echo "  (the OR~1e5 singleton estimates seen in the earlier run)."
echo "=============================================================="
step1 "--bt" "$PHENO_B" bt
step2 "--bt --firth --approx --pThresh 0.05" "$PHENO_B" bt \
      "$SETS/setlist_pooled.txt" chd_pooled

echo "done. results in $OUT"
