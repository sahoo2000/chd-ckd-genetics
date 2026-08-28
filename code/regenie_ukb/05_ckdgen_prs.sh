#!/usr/bin/env bash
# 05_ckdgen_prs.sh — CKD polygenic score from CKDGen, tested in UK Biobank CHD cases.
#
# This is the design that directly tests the comorbidity hypothesis at the COMMON-variant
# level, and it needs no exome data. With ~3,385 CHD cases there is 80% power to detect
# OR/SD = 1.05, which is inside the range where real cross-trait PRS effects live.
#
# CRITICAL: use a CKDGen stratum that EXCLUDES UK Biobank, or the test is circular.
# Stanzick 2021 reports CKDGen and UKB separately — take the CKDGen-only file.
set -euo pipefail

WORK=${WORK:-ckdgen}
TARGET=${TARGET:-ukb_imputed}          # PLINK2 pgen/pvar/psam prefix for UKB
KEEP=${KEEP:-../phenotypes/keep.tsv}
PHENO_B=${PHENO_B:-../phenotypes/pheno_binary.tsv}
COVAR=${COVAR:-../phenotypes/covariates.tsv}
THREADS=${THREADS:-16}
mkdir -p "$WORK"

# ---------------------------------------------------------------- 1. summary stats
# Browse https://ckdgen.imbi.uni-freiburg.de/datasets and take the CKDGen-only
# (non-UKB) eGFRcrea meta-analysis. Direct download, no application needed.
if [ ! -f "$WORK/egfr.txt.gz" ]; then
  echo "Place the CKDGen eGFRcrea summary statistics at $WORK/egfr.txt.gz, then re-run."
  echo "  https://ckdgen.imbi.uni-freiburg.de/datasets"
  exit 1
fi

# ---------------------------------------------------------------- 2. munge
# Expected CKDGen columns: RSID Chr Pos_b37 Allele1 Allele2 Freq1 Effect StdErr P-value n
# Harmonise to build 38 first if your target is 38 (CrossMap / LiftOver).
zcat "$WORK/egfr.txt.gz" \
  | awk 'NR==1{print "SNP A1 A2 BETA P"; next} {print $1, toupper($4), toupper($5), $7, $9}' \
  > "$WORK/egfr.munged"

# ---------------------------------------------------------------- 3. clump
plink2 --pfile "$TARGET" --keep "$KEEP" \
  --clump "$WORK/egfr.munged" --clump-p1 1 --clump-r2 0.1 --clump-kb 250 \
  --clump-field P --clump-snp-field SNP \
  --threads "$THREADS" --out "$WORK/egfr"

awk 'NR>1 && $3!=""{print $3}' "$WORK/egfr.clumps" > "$WORK/keep_snps.txt"
echo "[prs] independent variants retained: $(wc -l < "$WORK/keep_snps.txt")"

# ---------------------------------------------------------------- 4. score
# Several p-value thresholds; pick the best-performing in a held-out set, or use
# all-variant scores from PRS-CS / LDpred2 instead if you have the LD reference.
for T in 5e-8 1e-5 1e-3 0.05 1.0; do
  awk -v t="$T" 'NR==1 || $5+0 <= t' "$WORK/egfr.munged" > "$WORK/sc_$T.txt"
  plink2 --pfile "$TARGET" --keep "$KEEP" --extract "$WORK/keep_snps.txt" \
    --score "$WORK/sc_$T.txt" 1 2 4 header cols=+scoresums \
    --threads "$THREADS" --out "$WORK/prs_$T"
done

# ---------------------------------------------------------------- 5. test
python3 06_test_prs.py \
  --prs-glob "$WORK/prs_*.sscore" \
  --pheno "$PHENO_B" --covar "$COVAR" \
  --outcomes CHD CAKUT CHD_CAKUT \
  --out ../results/prs_results.tsv
