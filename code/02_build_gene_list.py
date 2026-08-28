# 02_build_gene_list.py
#
# Builds the candidate gene list used by every later script.
#
# The idea: take every gene that causes BOTH a heart malformation and a
# kidney malformation, then clean the list up so that only genes we can
# actually test are left.
#
# Three cleaning steps matter:
#
#  1. Some "genes" in these databases are not testable - they are
#     microRNAs, pseudogenes, or mitochondrial genes with no nuclear
#     protein-coding sequence.
#
#  2. Some genes only cause diseases that kill babies. Those people can
#     never appear in a biobank that recruits adults aged 40 to 69, so
#     including them just adds noise.
#
#  3. Some genes get in through autoimmune disease (HLA, PTPN22 and so on).
#     That is a different biological mechanism, and the HLA region has
#     such complicated genetics that burden tests there are unreliable.
#
# To run:   python3 02_build_gene_list.py

import os
import pandas as pd
import numpy as np

RAW = "../rawdata"
OUT_GENES = "../data/genes"

# The two functions below are the same ones used in script 03. They are
# copied here rather than imported so that this script works on its own.


def read_ontology(filename):
    """Read hp.obo and return term names, parent-child links, and dead terms."""
    term_name = {}
    children_of = {}
    obsolete_terms = set()
    current_term = None
    for line in open(filename):
        line = line.rstrip("\n")
        if line == "[Term]":
            current_term = None
        elif line.startswith("id: HP:"):
            current_term = line[4:]
        elif line.startswith("name: ") and current_term is not None:
            term_name[current_term] = line[6:]
        elif line.startswith("is_a: ") and current_term is not None:
            parent = line[6:].split(" ! ")[0].strip()
            if parent not in children_of:
                children_of[parent] = []
            children_of[parent].append(current_term)
        elif line.startswith("is_obsolete: true") and current_term is not None:
            obsolete_terms.add(current_term)
    return term_name, children_of, obsolete_terms


def get_all_terms_below(start_term, children_of, obsolete_terms):
    """Collect a term and every more specific term underneath it."""
    found = set([start_term])
    to_visit = [start_term]
    while len(to_visit) > 0:
        term = to_visit.pop()
        for child in children_of.get(term, []):
            if child not in found:
                found.add(child)
                to_visit.append(child)
    return found - obsolete_terms


# ======================================================================
print("Reading the ontology and annotations...")
term_name, children_of, obsolete_terms = read_ontology(os.path.join(RAW, "hp.obo"))
annotations = pd.read_csv(os.path.join(RAW, "genes_to_phenotype.txt"), sep="\t", dtype=str)

heart_terms = get_all_terms_below("HP:0001627", children_of, obsolete_terms)
kidney_terms = get_all_terms_below("HP:0012210", children_of, obsolete_terms)
ckd_terms = get_all_terms_below("HP:0012622", children_of, obsolete_terms)
ckd_terms = ckd_terms | get_all_terms_below("HP:0000083", children_of, obsolete_terms)


def count_diseases_per_gene(term_set):
    """For each gene, how many different diseases link it to these terms?"""
    rows = annotations[annotations["hpo_id"].isin(term_set)]
    return rows.groupby("gene_symbol")["disease_id"].nunique()


heart_counts = count_diseases_per_gene(heart_terms)
kidney_counts = count_diseases_per_gene(kidney_terms)
ckd_counts = count_diseases_per_gene(ckd_terms)

# genes that appear in both lists
overlap_genes = sorted(set(heart_counts.index) & set(kidney_counts.index))
print("genes causing both a heart and a kidney malformation:", len(overlap_genes))

genes = pd.DataFrame({"hpo_symbol": overlap_genes})
genes["n_dis_card"] = heart_counts.reindex(genes["hpo_symbol"]).fillna(0).astype(int).values
genes["n_dis_kid"] = kidney_counts.reindex(genes["hpo_symbol"]).fillna(0).astype(int).values
genes["n_dis_ckd"] = ckd_counts.reindex(genes["hpo_symbol"]).fillna(0).astype(int).values


# ----------------------------------------------------------------------
# CLEANING STEP 1: fix gene names and keep only protein-coding autosomes
# ----------------------------------------------------------------------

print("")
print("Cleaning step 1: checking gene names against HGNC...")

hgnc = pd.read_csv(os.path.join(RAW, "hgnc_complete_set.tsv"), sep="\t",
                   dtype=str, low_memory=False)
hgnc = hgnc[hgnc["status"] == "Approved"]

# build a lookup from old or alternative names to the current official name
rename_map = {}
for i in hgnc.index:
    official = hgnc.loc[i, "symbol"]
    for column in ["prev_symbol", "alias_symbol"]:
        value = hgnc.loc[i, column]
        if isinstance(value, str):
            for old_name in value.split("|"):
                if old_name not in rename_map:
                    rename_map[old_name] = official

official_symbols = set(hgnc["symbol"])
current_names = []
for name in genes["hpo_symbol"]:
    if name in official_symbols:
        current_names.append(name)
    else:
        current_names.append(rename_map.get(name, name))
genes["gene_symbol"] = current_names

hgnc_info = hgnc.set_index("symbol")[["locus_type", "location", "ensembl_id"
                                      if "ensembl_id" in hgnc.columns
                                      else "ensembl_gene_id"]]
hgnc_info.columns = ["locus_type", "location", "ensembl_id"]
genes = genes.join(hgnc_info, on="gene_symbol")

is_protein_coding = genes["locus_type"] == "gene with protein product"
print("  protein coding:", int(is_protein_coding.sum()), "of", len(genes))
genes = genes[is_protein_coding]

# keep autosomes only (chromosomes 1 to 22)
autosome_names = [str(number) for number in range(1, 23)]
chromosome = genes["location"].fillna("").str.extract(r"^(\d+|X|Y|mitochondria)")[0]
genes["chrom"] = chromosome
genes = genes[genes["chrom"].isin(autosome_names)]
print("  autosomal and protein coding:", len(genes))


# ----------------------------------------------------------------------
# CLEANING STEP 2: drop genes whose diseases kill in infancy
# ----------------------------------------------------------------------

print("")
print("Cleaning step 2: removing genes only linked to lethal infant disease...")

lethal_terms = ["HP:0001522",   # death in infancy
                "HP:0003811",   # neonatal death
                "HP:0003826",   # stillbirth
                "HP:0034241",   # prenatal death
                "HP:0003819"]   # death in childhood

adult_terms = ["HP:0003581", "HP:0003584", "HP:0003596", "HP:0011462",
               "HP:0003621", "HP:0011463", "HP:0025708", "HP:0025709",
               "HP:0025710", "HP:0033763", "HP:0033764", "HP:0033765",
               "HP:0100613"]

lethal_counts = annotations[annotations["hpo_id"].isin(lethal_terms)]
lethal_counts = lethal_counts.groupby("gene_symbol")["disease_id"].nunique()

adult_counts = annotations[annotations["hpo_id"].isin(adult_terms)]
adult_counts = adult_counts.groupby("gene_symbol")["disease_id"].nunique()

genes["n_lethal"] = lethal_counts.reindex(genes["hpo_symbol"]).fillna(0).astype(int).values
genes["n_adult"] = adult_counts.reindex(genes["hpo_symbol"]).fillna(0).astype(int).values

# keep a gene if it has any adult-onset disease, or if it has no lethal record
keep_gene = (genes["n_adult"] > 0) | (genes["n_lethal"] == 0)
print("  removing", int((~keep_gene).sum()), "genes as implausible in an adult cohort")
genes = genes[keep_gene]


# ----------------------------------------------------------------------
# CLEANING STEP 3: add gene positions and constraint scores
# ----------------------------------------------------------------------

print("")
print("Cleaning step 3: adding coordinates and constraint scores...")

mane = pd.read_csv(os.path.join(RAW, "MANE.summary.txt.gz"), sep="\t", dtype=str)
mane_positions = mane.drop_duplicates("symbol").set_index("symbol")[["chr_start", "chr_end"]]
genes = genes.join(mane_positions, on="gene_symbol")
genes["gene_start"] = pd.to_numeric(genes["chr_start"], errors="coerce")
genes["gene_end"] = pd.to_numeric(genes["chr_end"], errors="coerce")

constraint = pd.read_csv(os.path.join(RAW, "gnomad_constraint.tsv"), sep="\t",
                         low_memory=False)
if "mane_select" in constraint.columns:
    constraint = constraint[constraint["mane_select"] == True]
constraint = constraint.groupby("gene").agg(pLI=("lof.pLI", "max"),
                                            LOEUF=("lof.oe_ci.upper", "min"))
genes = genes.join(constraint, on="gene_symbol")
genes = genes[genes["gene_start"].notna()]


# ----------------------------------------------------------------------
# Put the genes into tiers, and pull out the ciliary genes
# ----------------------------------------------------------------------

# Tier 1 = heart + kidney + a documented chronic kidney disease link
genes["tier"] = np.where(genes["n_dis_ckd"] > 0, "1", "2")

# Immune genes come in through autoimmune disease, which is a different story
immune_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DRB5", "HLA-DQB1",
                "HLA-DQA1", "HLA-DQA2", "HLA-DQB2", "PTPN22", "STAT4", "MEFV",
                "TNXB", "IL6", "IL10", "TLR4", "SAA1", "HFE"]
genes["is_immune"] = genes["gene_symbol"].isin(immune_genes)

# The ciliary genes. These are the ones the analysis eventually points to.
ciliary_genes = ["AHI1", "ALG9", "ANKS6", "BBIP1", "BBS1", "BBS10", "BBS12", "BBS2",
                 "BBS4", "BBS5", "BBS7", "BBS9", "CC2D2A", "CEP290", "DNAJB11",
                 "DYNC2I1", "GANAB", "IFT122", "IFT140", "IFT172", "INVS", "LZTFL1",
                 "MKKS", "MKS1", "NEK8", "NPHP1", "NPHP3", "PKD1", "PKD2",
                 "RPGRIP1L", "SDCCAG8", "TMEM216", "TMEM231", "TMEM237", "TMEM67",
                 "TTC8", "ZNF423"]

# These are already known to cause cystic kidney disease. Keeping them in a
# pooled test guarantees a positive result that tells us nothing new, so we
# also write out a version with them taken out.
known_cystogenic = ["PKD1", "PKD2", "IFT140", "ALG9", "GANAB", "DNAJB11"]

tier1 = genes[(genes["tier"] == "1") & (~genes["is_immune"])]
ciliary = genes[genes["gene_symbol"].isin(ciliary_genes)]
ciliary_novel = ciliary[~ciliary["gene_symbol"].isin(known_cystogenic)]
ciliary_control = ciliary[ciliary["gene_symbol"].isin(known_cystogenic)]

if not os.path.exists(OUT_GENES):
    os.makedirs(OUT_GENES)

columns_to_save = ["gene_symbol", "ensembl_id", "chrom", "gene_start", "gene_end",
                   "n_dis_card", "n_dis_kid", "n_dis_ckd", "pLI", "LOEUF"]

print("")
print("Final gene sets:")
for set_name, gene_set in [("all_candidates", genes),
                           ("tier1", tier1),
                           ("ciliary", ciliary),
                           ("ciliary_novel", ciliary_novel),
                           ("ciliary_cystogenic_control", ciliary_control)]:
    table = gene_set[columns_to_save].sort_values(["chrom", "gene_start"])
    table.to_csv(os.path.join(OUT_GENES, set_name + ".tsv"), sep="\t", index=False)

    symbol_file = open(os.path.join(OUT_GENES, set_name + "_symbols.txt"), "w")
    for symbol in sorted(gene_set["gene_symbol"]):
        symbol_file.write(symbol + "\n")
    symbol_file.close()

    print("  %-28s %4d genes" % (set_name, len(gene_set)))

print("")
print("Done. Gene lists are in " + OUT_GENES)
