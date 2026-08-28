# 05_genebass_analysis.py
#
# Looks at the results that script 04 downloaded from Genebass and works
# out whether there is any real signal in them.
#
# The problem we have is that each gene on its own is underpowered. One
# gene tested against one phenotype in 394,841 people usually cannot
# reach significance unless the effect is very large. So instead of
# asking "is this one gene significant", we ask three questions about
# the whole set of results:
#
#   1. Are there MORE small p-values than we would expect by chance?
#   2. Do the results that ARE significant point the same way?
#   3. Which group of genes do they belong to?
#
# Question 2 is the clever one. If all our results were just noise, then
# among the ones with p < 0.05 about half would say "kidney gets worse"
# and half would say "kidney gets better". If nearly all of them say
# "worse", that is very unlikely to be luck.
#
# To run:   python3 05_genebass_analysis.py


import pandas as pd
import numpy as np
from scipy.stats import binomtest, mannwhitneyu


# ======================================================================
# SETTINGS - change these instead of using command line options
# ======================================================================

# which results file to look at
RESULTS_FILE = "ciliary_renal.tsv"

# the usual significance cut-off for gene-based tests
EXOME_WIDE_CUTOFF = 0.0000025      # this is 2.5e-6

# Genes we already know cause cystic kidney disease. If we leave these in,
# any test will come out positive just because of them, which tells us
# nothing new. So we remove them and see if anything is left.
ALREADY_KNOWN_GENES = ["PKD1", "PKD2", "IFT140", "ALG9", "GANAB", "DNAJB11"]

# Which part of the cilium each gene belongs to. The cilium is built and
# maintained by some of these groups, while others only sort cargo that
# passes through it, so it is interesting to see which group matters.
GENE_GROUPS = {}

for gene in ["BBS1", "BBS2", "BBS4", "BBS5", "BBS7", "BBS9", "BBIP1", "TTC8"]:
    GENE_GROUPS[gene] = "BBSome"

for gene in ["MKKS", "BBS10", "BBS12"]:
    GENE_GROUPS[gene] = "Chaperonin (BBS6/10/12)"

for gene in ["IFT122", "IFT140", "WDR19", "WDR35"]:
    GENE_GROUPS[gene] = "IFT-A"

for gene in ["IFT172", "IFT27", "IFT74"]:
    GENE_GROUPS[gene] = "IFT-B"

for gene in ["MKS1", "TMEM67", "TMEM216", "TMEM231", "TMEM237",
             "CC2D2A", "RPGRIP1L", "AHI1"]:
    GENE_GROUPS[gene] = "Transition zone / MKS"

for gene in ["NPHP1", "NPHP3", "NPHP4", "INVS", "NEK8", "ANKS6", "SDCCAG8"]:
    GENE_GROUPS[gene] = "NPHP module"

for gene in ["DYNC2H1", "DYNC2LI1", "DYNC2I1"]:
    GENE_GROUPS[gene] = "Dynein-2"

for gene in ["PKD1", "PKD2", "GANAB", "DNAJB11", "ALG9"]:
    GENE_GROUPS[gene] = "ADPKD / cystogenic"


# ======================================================================
# STEP 1: read the results file
# ======================================================================

print("Reading", RESULTS_FILE)
results = pd.read_csv(RESULTS_FILE, sep="\t")

# The p-value column is called P_SKATO. Make a shorter name for it, and
# also make a short name for the effect size column.
results["p_value"] = pd.to_numeric(results["P_SKATO"], errors="coerce")
results["effect"] = pd.to_numeric(results["BETA"], errors="coerce")

# Throw away any rows where the p-value did not read properly.
results = results.dropna(subset=["p_value"])

number_of_tests = len(results)
number_of_genes = results["gene"].nunique()

# If we do many tests we must make the cut-off stricter. The simplest way
# is to divide 0.05 by the number of tests. This is called the Bonferroni
# correction.
bonferroni_cutoff = 0.05 / number_of_tests

print("")
print("tests done      :", number_of_tests)
print("genes tested    :", number_of_genes)
print("Bonferroni cutoff: %.3g" % bonferroni_cutoff)
print("exome-wide cutoff: %.3g" % EXOME_WIDE_CUTOFF)


# ======================================================================
# STEP 2: which results are properly significant?
# ======================================================================

print("")
print("=" * 62)
print("STEP 2 - results that pass the exome-wide cutoff")
print("=" * 62)

significant = results[results["p_value"] < EXOME_WIDE_CUTOFF]
significant = significant.sort_values("p_value")

print("found", len(significant), "results in",
      significant["gene"].nunique(), "genes")
print("")

# print them one at a time
for row_number in significant.index:
    gene = results.loc[row_number, "gene"]
    phenotype = results.loc[row_number, "phenotype"]
    p_value = results.loc[row_number, "p_value"]
    effect = results.loc[row_number, "effect"]
    print("  %-9s %-42s p = %-10.3g beta = %.4g"
          % (gene, phenotype[:42], p_value, effect))


# ======================================================================
# STEP 3: is anything left after removing the genes we already knew?
# ======================================================================

print("")
print("=" * 62)
print("STEP 3 - is there signal left after removing known genes?")
print("=" * 62)

# Build a version of the table without the already-known genes.
remaining = results[~results["gene"].isin(ALREADY_KNOWN_GENES)]

print("removed", len(ALREADY_KNOWN_GENES), "known genes")
print("that leaves", len(remaining), "tests in",
      remaining["gene"].nunique(), "genes")
print("")
print("  %-12s %10s %10s %8s" % ("cutoff", "observed", "expected", "ratio"))

# For each cut-off, count how many results beat it, and compare that with
# how many we would expect just by chance. If p-values were pure noise
# then 5% of them would be below 0.05, 1% below 0.01, and so on.
for cutoff in [0.05, 0.01, 0.001, 0.0001]:
    observed = 0
    for p_value in remaining["p_value"]:
        if p_value < cutoff:
            observed = observed + 1

    expected = cutoff * len(remaining)
    ratio = observed / expected

    print("  p < %-8g %10d %10.1f %7.1fx" % (cutoff, observed, expected, ratio))


# ======================================================================
# STEP 4: do the significant results point the same way?
# ======================================================================

print("")
print("=" * 62)
print("STEP 4 - do the results point in the same direction?")
print("=" * 62)
print("")
print("A positive effect means rare damage to the gene makes kidney")
print("function WORSE, which is what we would expect if the gene matters.")
print("If our results were noise, half would be positive and half negative.")
print("")

# We check three different sets of results
sets_to_check = [("all tests", results),
                 ("nominally significant (p < 0.05)", results[results["p_value"] < 0.05]),
                 ("nominally significant, known genes removed",
                  remaining[remaining["p_value"] < 0.05])]

for set_name, table in sets_to_check:
    # count how many effects are positive and how many are not
    positive_count = 0
    negative_count = 0
    for effect in table["effect"]:
        if pd.notna(effect):
            if effect > 0:
                positive_count = positive_count + 1
            else:
                negative_count = negative_count + 1

    total_count = positive_count + negative_count
    if total_count == 0:
        continue

    percent_positive = 100.0 * positive_count / total_count

    # A binomial test asks: if the true chance were 50/50, how surprising
    # is this many positives?
    test_result = binomtest(positive_count, total_count, 0.5, alternative="greater")

    print("  %-44s %3d of %3d positive (%3.0f%%)  p = %.3g"
          % (set_name, positive_count, total_count, percent_positive,
             test_result.pvalue))


# ======================================================================
# STEP 5: which part of the cilium carries the signal?
# ======================================================================

print("")
print("=" * 62)
print("STEP 5 - which group of ciliary genes carries the signal?")
print("=" * 62)

# Label every row with the group its gene belongs to
group_labels = []
for gene in results["gene"]:
    if gene in GENE_GROUPS:
        group_labels.append(GENE_GROUPS[gene])
    else:
        group_labels.append("Other")
results["group"] = group_labels

print("")
print("  %-26s %6s %7s %8s %10s" % ("group", "genes", "tests", "nominal", "enrichment"))

group_rows = []
for group_name in sorted(set(results["group"])):
    group_table = results[results["group"] == group_name]

    # count nominally significant results in this group
    nominal_count = 0
    for p_value in group_table["p_value"]:
        if p_value < 0.05:
            nominal_count = nominal_count + 1

    expected_count = 0.05 * len(group_table)
    enrichment = nominal_count / expected_count

    group_rows.append({"group": group_name,
                       "genes": group_table["gene"].nunique(),
                       "tests": len(group_table),
                       "nominal": nominal_count,
                       "enrichment": enrichment})

# sort the groups so the most enriched comes first
group_rows.sort(key=lambda row: row["enrichment"], reverse=True)

for row in group_rows:
    print("  %-26s %6d %7d %8d %9.2fx"
          % (row["group"], row["genes"], row["tests"],
             row["nominal"], row["enrichment"]))


# ======================================================================
# STEP 6: is the signal just because some genes have more variants?
# ======================================================================

print("")
print("=" * 62)
print("STEP 6 - is the signal just about how well powered each gene is?")
print("=" * 62)

# A gene with more rare variants in it is easier to detect. We should
# check that our result is not simply picking out the biggest genes.
if "n_variants" in results.columns:

    # For every gene, find its variant count and how many nominal hits it got
    variants_per_gene = results.groupby("gene")["n_variants"].max()

    hits_per_gene = {}
    for gene in results["gene"].unique():
        gene_table = results[results["gene"] == gene]
        count = 0
        for p_value in gene_table["p_value"]:
            if p_value < 0.05:
                count = count + 1
        hits_per_gene[gene] = count

    # split the genes into two groups
    genes_with_hits = []
    genes_without_hits = []
    for gene in hits_per_gene:
        if hits_per_gene[gene] >= 2:
            genes_with_hits.append(variants_per_gene[gene])
        else:
            genes_without_hits.append(variants_per_gene[gene])

    print("  genes with 2 or more nominal hits: %3d, median %5.0f variants"
          % (len(genes_with_hits), np.median(genes_with_hits)))
    print("  genes with fewer than 2 hits     : %3d, median %5.0f variants"
          % (len(genes_without_hits), np.median(genes_without_hits)))

    test_result = mannwhitneyu(genes_with_hits, genes_without_hits,
                               alternative="greater")
    print("  Mann-Whitney p = %.3g" % test_result.pvalue)

    if test_result.pvalue > 0.05:
        print("  -> the signal is NOT explained by variant counts")
    else:
        print("  -> better powered genes do show more hits, so be careful here")

print("")
print("Finished.")
