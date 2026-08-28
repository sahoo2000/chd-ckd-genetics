# 08_run_magma.sh
#
# MAGMA turns a GWAS (which gives a p-value for every SNP) into a p-value
# for every gene, and then tests whether our candidate genes have smaller
# p-values than genes in general.
#
# To run:  bash 08_run_magma.sh
#
# WARNING: step 2 is slow. It took about 4 hours on a laptop.

RAW=../rawdata
OUT=../data/results

# ---------------------------------------------------------------
# Step 0: prepare the two input files MAGMA needs from the GWAS
# ---------------------------------------------------------------
# MAGMA wants one file with SNP positions and one with SNP p-values.
# The CKDGen file columns are:
#   4=n  6=effect  7=se  8=pvalue  11=chr  12=pos  15=rsid

echo "preparing input files..."
gzip -cd $RAW/ckdgen_egfr.gz | awk -F'\t' 'NR>1 && $15 ~ /^rs/ {
    print $15"\t"$11"\t"$12 > "egfr.snploc"
    print $15"\t"$8"\t"$4  > "egfr.pval.body"
}'
echo -e "SNP\tP\tN" > egfr.pval
cat egfr.pval.body >> egfr.pval
rm egfr.pval.body

# ---------------------------------------------------------------
# Step 1: decide which SNPs belong to which gene
# ---------------------------------------------------------------
echo "step 1: annotating SNPs to genes"
$RAW/magma --annotate \
    --snp-loc egfr.snploc \
    --gene-loc $RAW/NCBI37.3.gene.loc \
    --out egfr

# ---------------------------------------------------------------
# Step 2: work out a p-value for every gene (SLOW)
# ---------------------------------------------------------------
echo "step 2: gene analysis - this takes a few hours"
$RAW/magma \
    --bfile $RAW/g1000_eur/g1000_eur \
    --pval egfr.pval use=SNP,P ncol=N \
    --gene-annot egfr.genes.annot \
    --out egfr

# ---------------------------------------------------------------
# Step 3: test our candidate gene sets
# ---------------------------------------------------------------
echo "step 3: gene set analysis"
$RAW/magma \
    --gene-results egfr.genes.raw \
    --set-annot ../data/genes/magma_genesets.txt \
    --out egfr_genesets

mv egfr.genes.out egfr_genesets.gsa.out $OUT/ 2>/dev/null
echo "done - results are in $OUT"
