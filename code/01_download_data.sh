# 01_download_data.sh
#
# Downloads every public dataset this project uses.
# Everything here is free and needs no application or login.
#
# To run:   bash 01_download_data.sh
#
# It makes a folder called "rawdata" next to the code folder and puts
# everything in there. If a file is already downloaded it is skipped,
# so it is safe to run this script again.

mkdir -p ../rawdata
cd ../rawdata

echo "=========================================="
echo "1. Human Phenotype Ontology"
echo "=========================================="

# The ontology itself: all the phenotype terms and how they relate to each other
if [ ! -f hp.obo ]; then
    curl -L -o hp.obo "https://purl.obolibrary.org/obo/hp.obo"
else
    echo "hp.obo already downloaded, skipping"
fi

# Which gene causes which phenotype in which disease
if [ ! -f genes_to_phenotype.txt ]; then
    curl -L -o genes_to_phenotype.txt \
        "https://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt"
else
    echo "genes_to_phenotype.txt already downloaded, skipping"
fi

echo ""
echo "=========================================="
echo "2. Gene information"
echo "=========================================="

# Official gene names and symbols
if [ ! -f hgnc_complete_set.tsv ]; then
    curl -L -o hgnc_complete_set.tsv \
        "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.tsv"
else
    echo "hgnc_complete_set.tsv already downloaded, skipping"
fi

# Where each gene sits on the chromosome (genome build GRCh38)
if [ ! -f MANE.summary.txt.gz ]; then
    curl -L -o MANE.summary.txt.gz \
        "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.5/MANE.GRCh38.v1.5.summary.txt.gz"
else
    echo "MANE.summary.txt.gz already downloaded, skipping"
fi

# How intolerant each gene is to being broken (pLI and LOEUF scores)
if [ ! -f gnomad_constraint.tsv ]; then
    curl -L -o gnomad_constraint.tsv \
        "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv"
else
    echo "gnomad_constraint.tsv already downloaded, skipping"
fi

echo ""
echo "=========================================="
echo "3. GWAS summary statistics"
echo "=========================================="

# Kidney function (eGFR) from the CKDGen consortium, about 1 million people
if [ ! -f ckdgen_egfr.gz ]; then
    curl -L -o ckdgen_egfr.gz \
        "https://ckdgen.imbi.uni-freiburg.de/files/Stanzick2021/metal_eGFR_meta_ea1.TBL.map.annot.gc.gz"
else
    echo "ckdgen_egfr.gz already downloaded, skipping"
fi

# FinnGen release 10. These files are about 800 MB each, so this part is slow.
# CONGEN_HEART_ARTER      = congenital malformations of heart and great arteries
# N14_CHRONKIDNEYDIS      = chronic kidney disease
# Q17_CYSTIC_KIDNEY_DISEA = cystic kidney disease
for PHENOTYPE in CONGEN_HEART_ARTER N14_CHRONKIDNEYDIS Q17_CYSTIC_KIDNEY_DISEA
do
    if [ ! -f finngen_R10_$PHENOTYPE.gz ]; then
        echo "downloading $PHENOTYPE ..."
        curl -L -O \
          "https://storage.googleapis.com/finngen-public-data-r10/summary_stats/finngen_R10_$PHENOTYPE.gz"
    else
        echo "$PHENOTYPE already downloaded, skipping"
    fi
done

echo ""
echo "=========================================="
echo "4. Reference files for MAGMA and LD score regression"
echo "=========================================="

# MAGMA program and its helper files
if [ ! -f magma ]; then
    curl -L -o magma_mac.zip "https://vu.data.surfsara.nl/index.php/s/1M1d9vHtVidEwvU/download"
    unzip -o magma_mac.zip
    chmod +x magma
fi

if [ ! -f NCBI37.3.gene.loc ]; then
    curl -L -o geneloc37.zip "https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download"
    unzip -o geneloc37.zip
fi

# 1000 Genomes European reference panel (about 3 GB, this takes a while)
if [ ! -d g1000_eur ]; then
    curl -L -o g1000_eur.zip "https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download"
    unzip -o g1000_eur.zip -d g1000_eur
fi

# Pre-computed LD scores for LD score regression
if [ ! -d eur_w_ld_chr ]; then
    curl -L -o eur_w_ld_chr.tar.gz \
        "https://zenodo.org/api/records/18749273/files/eur_w_ld_chr.tar.gz/content"
    tar -xzf eur_w_ld_chr.tar.gz
fi

echo ""
echo "All downloads finished. Files are in the rawdata folder:"
ls -lh
