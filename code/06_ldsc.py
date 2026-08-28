# 06_ldsc.py
#
# LD score regression. This measures two things:
#
#   1. Heritability: how much of the variation in a trait is explained
#      by common genetic variants.
#
#   2. Genetic correlation: whether two traits are influenced by the
#      same genetic variants. This is the number we actually want,
#      because it tells us directly whether congenital heart disease
#      and kidney disease share a genetic basis.
#
# HOW IT WORKS
#
# Some parts of the genome are inherited in big blocks, so one variant
# carries information about its neighbours. The "LD score" of a variant
# counts how many neighbours it speaks for. A variant with a high LD
# score is standing in for a lot of the genome.
#
# If a trait is genuinely heritable, variants with high LD scores will
# show stronger associations, because they are tagging more real signal.
# So we plot association strength against LD score and look at the slope:
#
#     chi-squared  =  1  +  (N * heritability / M) * LD_score
#
# The slope gives us heritability. For two traits at once we multiply
# their z-scores together and do the same thing, and the slope then gives
# the genetic covariance.
#
# The intercept is useful too. It picks up confounding, and in the
# two-trait version it absorbs any overlap between the two samples.
#
# NOTE: the original LDSC software is written in Python 2, which no
# longer installs. This script does the same calculation. It was checked
# by confirming that (a) the heritability of kidney function comes out
# close to the published value and (b) the positive control pair
# behaves as expected. It is not the original program.
#
# To run:   python3 06_ldsc.py


import pandas as pd
import numpy as np
from scipy.stats import norm


# ======================================================================
# SETTINGS
# ======================================================================

LDSCORE_FILE = "ldscores.tsv.gz"

# M is the number of variants used to build the LD scores. It comes from
# the .M_5_50 files that ship with the LD score reference data.
M = 1173569

# The traits to analyse.
#   name        : what to call it in the output
#   file        : the filtered summary statistics
#   n_cases     : leave as 0 for a continuous trait like eGFR
#   n_controls  : leave as 0 for a continuous trait
TRAITS = [
    {"name": "eGFR_CKDGen",   "file": "egfr.hm3.tsv",
     "n_cases": 0,     "n_controls": 0},
    {"name": "CHD_FinnGen",   "file": "fg_CONGEN_HEART_ARTER.hm3.tsv",
     "n_cases": 4270,  "n_controls": 407911},
    {"name": "CKD_FinnGen",   "file": "fg_N14_CHRONKIDNEYDIS.hm3.tsv",
     "n_cases": 10039, "n_controls": 396706},
    {"name": "Cystic_FinnGen", "file": "fg_Q17_CYSTIC_KIDNEY_DISEA.hm3.tsv",
     "n_cases": 947,   "n_controls": 410449},
]

OUTPUT_FILE = "ldsc_results.tsv"

# We split the genome into this many blocks to work out error bars.
NUMBER_OF_BLOCKS = 200

# Variants where the two alleles are A/T or C/G are ambiguous, because we
# cannot tell which strand of the DNA the study reported. We drop them.
AMBIGUOUS_PAIRS = [("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")]


# ======================================================================
# READING AND CLEANING A SET OF SUMMARY STATISTICS
# ======================================================================

def read_and_clean(filename, trait_name, n_cases, n_controls):
    """
    Read one summary statistics file and tidy it up.

    The file has no header. Its columns are:
        SNP, effect allele, other allele, beta, standard error,
        then either (N, frequency) for a continuous trait
        or (frequency) for a binary trait.
    """
    column_names = ["SNP", "A1", "A2", "BETA", "SE", "COL6", "COL7"]
    table = pd.read_csv(filename, sep="\t", header=None, names=column_names)

    starting_rows = len(table)

    # ---- work out the sample size ----
    if n_cases > 0:
        # For a yes/no trait, the useful sample size is smaller than the
        # total, because a study with 100 cases and 100,000 controls is
        # not as powerful as one with 50,000 of each. This formula gives
        # the "effective" sample size.
        effective_n = 4.0 / (1.0 / n_cases + 1.0 / n_controls)
        table["N"] = effective_n
        table["FREQ"] = pd.to_numeric(table["COL6"], errors="coerce")
    else:
        # For a continuous trait the file already has a sample size column
        table["N"] = pd.to_numeric(table["COL6"], errors="coerce")
        table["FREQ"] = pd.to_numeric(table["COL7"], errors="coerce")

    # ---- turn the effect into a z-score ----
    table["BETA"] = pd.to_numeric(table["BETA"], errors="coerce")
    table["SE"] = pd.to_numeric(table["SE"], errors="coerce")
    table = table.dropna(subset=["BETA", "SE", "N"])
    table = table[table["SE"] > 0]
    table["Z"] = table["BETA"] / table["SE"]

    # ---- remove variants we cannot use ----

    # any SNP listed more than once
    table = table[~table["SNP"].duplicated(keep=False)]

    # ambiguous strand variants
    keep_row = []
    for allele1, allele2 in zip(table["A1"], table["A2"]):
        if (allele1, allele2) in AMBIGUOUS_PAIRS:
            keep_row.append(False)
        else:
            keep_row.append(True)
    table = table[keep_row]

    # anything that is not a plain single letter
    table = table[table["A1"].isin(["A", "C", "G", "T"])]
    table = table[table["A2"].isin(["A", "C", "G", "T"])]

    # very rare variants, which are unreliable
    table = table[table["FREQ"] > 0.01]
    table = table[table["FREQ"] < 0.99]

    # A handful of variants have enormous effects. They are usually real
    # but they would dominate the regression, so LDSC removes them.
    largest_allowed = max(80.0, 0.001 * table["N"].max())
    table = table[table["Z"] * table["Z"] <= largest_allowed]

    print("  %-16s %8d variants kept out of %8d   (mean sample size %s)"
          % (trait_name, len(table), starting_rows,
             "{:,.0f}".format(table["N"].mean())))

    return table[["SNP", "A1", "A2", "Z", "N"]]


# ======================================================================
# THE REGRESSION ITSELF
# ======================================================================

def weighted_straight_line(x_values, y_values, weights):
    """
    Fit y = intercept + slope * x, where some points count more than
    others. Returns the intercept and the slope.
    """
    # Build the design matrix: a column of ones, then the x values
    ones_column = np.ones(len(x_values))
    design = np.column_stack([ones_column, x_values])

    # Applying weights is the same as multiplying both sides by the
    # square root of the weight before doing ordinary least squares.
    root_weights = np.sqrt(weights)
    weighted_design = design * root_weights[:, None]
    weighted_y = y_values * root_weights

    answer = np.linalg.lstsq(weighted_design, weighted_y, rcond=None)[0]
    intercept = answer[0]
    slope = answer[1]
    return intercept, slope


def slope_with_error_bar(x_values, y_values, weights, want="slope"):
    """
    Fit the line, then work out an error bar by leaving out one block of
    the genome at a time and seeing how much the answer moves. This is
    called a jackknife.
    """
    def fit(x, y, w):
        intercept, slope = weighted_straight_line(x, y, w)
        if want == "slope":
            return slope
        else:
            return intercept

    full_answer = fit(x_values, y_values, weights)

    # split the data into blocks
    all_positions = np.arange(len(x_values))
    blocks = np.array_split(all_positions, min(NUMBER_OF_BLOCKS, len(x_values)))

    answers_without_one_block = []
    for block in blocks:
        keep = np.ones(len(x_values), dtype=bool)
        keep[block] = False
        answers_without_one_block.append(fit(x_values[keep],
                                             y_values[keep],
                                             weights[keep]))

    answers_without_one_block = np.array(answers_without_one_block)
    number_of_blocks = len(answers_without_one_block)

    # standard jackknife formula
    pseudo_values = (number_of_blocks * full_answer
                     - (number_of_blocks - 1) * answers_without_one_block)
    standard_error = pseudo_values.std(ddof=1) / np.sqrt(number_of_blocks)

    return full_answer, standard_error


def estimate_heritability(trait_table, ld_scores):
    """How much of this trait is explained by common variants?"""
    merged = trait_table.merge(ld_scores, on="SNP")

    ld_score = merged["L2"].values
    chi_squared = merged["Z"].values ** 2
    sample_size = merged["N"].values

    # The model is:  chi_squared = 1 + (N * h2 / M) * LD_score
    # so if we regress chi_squared on (LD_score * N / M) the slope is h2.
    x_values = ld_score * sample_size / M

    # Variants that tag a lot of the genome are correlated with each
    # other, so they should not each count fully. Weighting by 1/LD score
    # corrects for that.
    ld_weights = 1.0 / np.maximum(ld_score, 1.0)

    # First pass, to get a rough heritability
    rough_intercept, rough_h2 = weighted_straight_line(x_values, chi_squared,
                                                       ld_weights)
    if rough_h2 < 0:
        rough_h2 = 0.0
    if rough_h2 > 1:
        rough_h2 = 1.0

    # Second pass. Points with bigger expected values are noisier, so we
    # down-weight them using the rough answer from the first pass.
    spread = (1.0 + rough_h2 * x_values) ** 2
    final_weights = 1.0 / (np.maximum(ld_score, 1.0) * spread)

    heritability, h2_error = slope_with_error_bar(x_values, chi_squared,
                                                  final_weights, "slope")
    intercept, intercept_error = slope_with_error_bar(x_values, chi_squared,
                                                      final_weights, "intercept")

    return {"h2": heritability, "h2_se": h2_error,
            "intercept": intercept, "intercept_se": intercept_error,
            "mean_chi2": chi_squared.mean(), "n_snp": len(merged)}


def estimate_correlation(table_one, table_two, ld_scores, h2_one, h2_two):
    """Do two traits share the same genetic variants?"""
    merged = table_one.merge(table_two, on="SNP", suffixes=("_1", "_2"))
    merged = merged.merge(ld_scores, on="SNP")

    # The two studies may have written the alleles the other way round.
    # If they match, keep the z-score as it is. If they are swapped, flip
    # the sign so both refer to the same allele.
    flipped_z = []
    keep_row = []
    for i in merged.index:
        allele1_first = merged.loc[i, "A1_1"]
        allele2_first = merged.loc[i, "A2_1"]
        allele1_second = merged.loc[i, "A1_2"]
        allele2_second = merged.loc[i, "A2_2"]
        z_second = merged.loc[i, "Z_2"]

        if allele1_first == allele1_second and allele2_first == allele2_second:
            keep_row.append(True)
            flipped_z.append(z_second)
        elif allele1_first == allele2_second and allele2_first == allele1_second:
            keep_row.append(True)
            flipped_z.append(-z_second)
        else:
            # alleles do not match at all, so drop this variant
            keep_row.append(False)
            flipped_z.append(np.nan)

    merged["Z_2_aligned"] = flipped_z
    merged = merged[keep_row]

    ld_score = merged["L2"].values
    z_one = merged["Z_1"].values
    z_two = merged["Z_2_aligned"].values
    n_one = merged["N_1"].values
    n_two = merged["N_2"].values

    average_n = np.sqrt(n_one * n_two)
    x_values = ld_score * average_n / M
    y_values = z_one * z_two

    ld_weights = 1.0 / np.maximum(ld_score, 1.0)
    rough_intercept, rough_covariance = weighted_straight_line(x_values, y_values,
                                                               ld_weights)

    # second pass weights, same idea as for heritability
    spread_one = h2_one * ld_score * n_one / M + 1.0
    spread_two = h2_two * ld_score * n_two / M + 1.0
    final_weights = 1.0 / (np.maximum(ld_score, 1.0)
                           * (spread_one * spread_two
                              + (rough_covariance * x_values) ** 2))

    covariance, covariance_error = slope_with_error_bar(x_values, y_values,
                                                        final_weights, "slope")
    intercept, intercept_error = slope_with_error_bar(x_values, y_values,
                                                      final_weights, "intercept")

    # Correlation is covariance divided by the two heritabilities.
    bottom = np.sqrt(max(h2_one, 0.000001) * max(h2_two, 0.000001))
    correlation = covariance / bottom
    correlation_error = covariance_error / bottom

    return {"rg": correlation, "rg_se": correlation_error,
            "covariance": covariance, "intercept": intercept,
            "n_snp": len(merged)}


# ======================================================================
# MAIN
# ======================================================================

print("Reading the LD scores...")
ld_scores = pd.read_csv(LDSCORE_FILE, sep="\t")
ld_scores = ld_scores[["SNP", "L2"]]
print("  %s variants with an LD score" % "{:,}".format(len(ld_scores)))

print("")
print("Reading and cleaning the summary statistics...")
cleaned_traits = {}
for trait in TRAITS:
    cleaned_traits[trait["name"]] = read_and_clean(trait["file"],
                                                   trait["name"],
                                                   trait["n_cases"],
                                                   trait["n_controls"])

# ---------------------------------------------------------------- h2
print("")
print("=" * 70)
print("HERITABILITY")
print("=" * 70)

heritabilities = {}
rows_to_save = []

for trait in TRAITS:
    name = trait["name"]
    answer = estimate_heritability(cleaned_traits[name], ld_scores)
    heritabilities[name] = answer["h2"]

    # The z-score tells us how solid the estimate is. Below 4 is shaky.
    z_score = answer["h2"] / answer["h2_se"]

    print("  %-16s h2 = %.4f (SE %.4f)   z = %5.2f   intercept = %.3f"
          % (name, answer["h2"], answer["h2_se"], z_score, answer["intercept"]))

    if z_score < 4:
        print("      WARNING: z below 4, so correlations involving this")
        print("      trait are underpowered and should be read carefully.")

    rows_to_save.append({"analysis": "h2", "trait1": name, "trait2": "",
                         "estimate": answer["h2"], "se": answer["h2_se"],
                         "z": z_score, "intercept": answer["intercept"],
                         "n_snp": answer["n_snp"]})

# ---------------------------------------------------------------- rg
print("")
print("=" * 70)
print("GENETIC CORRELATION")
print("=" * 70)

for i in range(len(TRAITS)):
    for j in range(i + 1, len(TRAITS)):
        name_one = TRAITS[i]["name"]
        name_two = TRAITS[j]["name"]

        answer = estimate_correlation(cleaned_traits[name_one],
                                      cleaned_traits[name_two],
                                      ld_scores,
                                      heritabilities[name_one],
                                      heritabilities[name_two])

        z_score = answer["rg"] / answer["rg_se"]
        p_value = 2 * norm.sf(abs(z_score))
        lower = answer["rg"] - 1.96 * answer["rg_se"]
        upper = answer["rg"] + 1.96 * answer["rg_se"]

        print("")
        print("  %s  vs  %s" % (name_one, name_two))
        print("      rg = %.4f (SE %.4f)   p = %.3g   95%% CI %.3f to %.3f"
              % (answer["rg"], answer["rg_se"], p_value, lower, upper))

        rows_to_save.append({"analysis": "rg", "trait1": name_one,
                             "trait2": name_two, "estimate": answer["rg"],
                             "se": answer["rg_se"], "z": z_score,
                             "intercept": answer["intercept"],
                             "n_snp": answer["n_snp"]})

pd.DataFrame(rows_to_save).to_csv(OUTPUT_FILE, sep="\t", index=False,
                                  float_format="%.6g")
print("")
print("Saved results to", OUTPUT_FILE)
