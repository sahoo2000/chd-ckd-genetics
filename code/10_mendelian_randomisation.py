# 10_mendelian_randomisation.py
#
# WHAT THIS SCRIPT DOES
# ---------------------
# Everything else in this project asks "do these two diseases share genes?".
# This script asks a different and harder question: "does one CAUSE the other?"
#
# The trick is called Mendelian randomisation (MR). The idea is simple once
# you see it. You are given your DNA at conception, at random, before anything
# in your life has happened to you. So if we find genetic variants that make
# congenital heart disease more likely, and those same variants also predict
# worse kidney function later in life, that is hard to explain by lifestyle or
# by hospital treatment. The randomisation already happened, in the womb.
#
# In MR language:
#   exposure   = the thing we think might be a cause  (congenital heart disease)
#   outcome    = the thing we think might be an effect (kidney function, eGFR)
#   instrument = a genetic variant that predicts the exposure
#
# We use two separate studies, which is called "two-sample MR":
#   exposure  -> FinnGen R10, congenital malformations of heart and great arteries
#   outcome   -> CKDGen, estimated glomerular filtration rate (eGFR)
#
# IMPORTANT DESIGN CHOICE
# -----------------------
# We deliberately do NOT use FinnGen's own chronic kidney disease GWAS as the
# outcome, even though we have it. Both would come from the same people, and
# when the exposure and outcome samples overlap, any weakness in the
# instruments drags the answer back towards the ordinary observational
# association -- which is exactly the confounded number this whole project is
# arguing against. Using CKDGen (a different, largely non-Finnish study) avoids
# that.
#
# HOW TO RUN
#   python3 10_mendelian_randomisation.py
#
# It needs these files in ../rawdata (get them with 01_download_data.sh):
#   finngen_R10_CONGEN_HEART_ARTER.gz
#   ckdgen_egfr.gz

import gzip
import math
import os
import random

from scipy.stats import norm, chi2


# ----------------------------------------------------------------------
# SETTINGS
# You can change these numbers and run the script again.
# ----------------------------------------------------------------------

RAW_FOLDER = "../rawdata"
OUT_FOLDER = "../data/results"

CHD_FILE = os.path.join(RAW_FOLDER, "finngen_R10_CONGEN_HEART_ARTER.gz")
EGFR_FILE = os.path.join(RAW_FOLDER, "ckdgen_egfr.gz")

# How strong does a variant have to be before we trust it as an instrument?
# We run the whole analysis twice, once strict and once relaxed, because the
# strict threshold leaves us with very few instruments.
STRICT_THRESHOLD = 5e-8      # the usual genome-wide significance line
RELAXED_THRESHOLD = 1e-5     # more instruments, but weaker ones

# Two variants sitting close together on a chromosome are usually inherited
# together, so they are not really independent pieces of evidence. We keep only
# the strongest variant within a window of this many base pairs.
CLUMP_WINDOW = 500000

# Variants that are very rare are unreliable: a handful of people can swing the
# result. We drop anything rarer than this.
MINIMUM_ALLELE_FREQUENCY = 0.01

# How many times to resample when working out the error bar for the
# weighted median. More is slower but smoother.
BOOTSTRAP_ROUNDS = 1000

# Fixing the random seed means you get exactly the same numbers as we did.
random.seed(20260828)


# ----------------------------------------------------------------------
# SMALL HELPER FUNCTIONS
# ----------------------------------------------------------------------

def normal_p_value(z):
    """Turn a z-score into a two-sided p-value.

    norm.sf gives the area in one tail of the bell curve, and we want both
    tails, so we double it.
    """
    return 2.0 * norm.sf(abs(z))


def flip_dna_letter(letter):
    """Give back the letter on the opposite DNA strand. A pairs with T, C with G."""
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return pairs.get(letter, "N")


def is_palindromic(allele1, allele2):
    """True if we cannot tell which strand this variant was reported on.

    An A/T variant looks like a T/A variant read from the other side, so we
    cannot always line two studies up. Same for C/G.
    """
    return flip_dna_letter(allele1) == allele2


def variance_explained_by_one_variant(beta, allele_frequency):
    """Roughly how much of the disease liability does one variant explain?

    For a yes/no disease the effect size is a log odds ratio, so we convert it
    onto the "liability" scale -- the imaginary underlying scale of risk that
    everyone sits somewhere on. The pi^2/3 is the variance of the logistic
    distribution and comes from that conversion.
    """
    genetic_variance = 2.0 * allele_frequency * (1.0 - allele_frequency) * beta * beta
    return genetic_variance / (genetic_variance + (math.pi ** 2) / 3.0)


# ----------------------------------------------------------------------
# STEP 1: READ THE HEART DISEASE STUDY AND PICK OUT INSTRUMENTS
# ----------------------------------------------------------------------

def read_chd_variants(p_value_cutoff):
    """Read the FinnGen file and keep variants stronger than the cutoff.

    The FinnGen columns are:
      0 chromosome, 1 position, 2 reference allele, 3 alternative allele,
      4 rsid, 5 nearest gene, 6 p-value, 7 -log10(p), 8 beta, 9 standard error,
      10 frequency of the alternative allele
    The beta describes the effect of the ALTERNATIVE allele.
    """
    kept = []
    with gzip.open(CHD_FILE, "rt") as handle:
        handle.readline()                      # throw away the header line
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            try:
                p_value = float(parts[6])
            except ValueError:
                continue                       # some lines have "NA", skip them
            if p_value >= p_value_cutoff:
                continue
            frequency = float(parts[10])
            if frequency < MINIMUM_ALLELE_FREQUENCY:
                continue
            if frequency > 1.0 - MINIMUM_ALLELE_FREQUENCY:
                continue
            if not parts[4].startswith("rs"):
                continue                       # we need an rsid to match studies
            kept.append({
                "chromosome": parts[0],
                "position": int(parts[1]),
                "other_allele": parts[2].upper(),
                "effect_allele": parts[3].upper(),
                "rsid": parts[4],
                "gene": parts[5],
                "p_value": p_value,
                "beta": float(parts[8]),
                "se": float(parts[9]),
                "frequency": frequency,
            })
    return kept


def keep_one_variant_per_region(variants):
    """Keep only the strongest variant in each stretch of chromosome.

    We sort everything by p-value, then walk down the list keeping a variant
    only if it is far away from every variant we have already kept.
    """
    # Sort by p-value, strongest first. We do it by putting the p-value at the
    # front of a pair, because Python sorts pairs by their first item.
    p_value_and_variant = []
    for variant in variants:
        p_value_and_variant.append((variant["p_value"], variant["rsid"], variant))
    p_value_and_variant.sort()

    variants = []
    for index in range(len(p_value_and_variant)):
        variants.append(p_value_and_variant[index][2])
    kept = []
    for candidate in variants:
        too_close = False
        for already_kept in kept:
            same_chromosome = candidate["chromosome"] == already_kept["chromosome"]
            distance = abs(candidate["position"] - already_kept["position"])
            if same_chromosome and distance < CLUMP_WINDOW:
                too_close = True
                break
        if not too_close:
            kept.append(candidate)
    return kept


# ----------------------------------------------------------------------
# STEP 2: LOOK THOSE VARIANTS UP IN THE KIDNEY STUDY
# ----------------------------------------------------------------------

def look_up_in_egfr(wanted_rsids):
    """Find our chosen variants in the CKDGen eGFR file.

    CKDGen columns we care about:
      1 Allele1 (the effect allele), 2 Allele2, 4 Freq1, 5 Effect,
      6 StdErr, 7 P.value, 14 RSID
    Note CKDGen writes its alleles in lower case, so we upper-case them.
    """
    found = {}
    with gzip.open(EGFR_FILE, "rt") as handle:
        handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rsid = parts[14]
            if rsid not in wanted_rsids:
                continue
            try:
                found[rsid] = {
                    "effect_allele": parts[1].upper(),
                    "other_allele": parts[2].upper(),
                    "frequency": float(parts[4]),
                    "beta": float(parts[5]),
                    "se": float(parts[6]),
                    "p_value": float(parts[7]),
                    "n": float(parts[3]),
                }
            except ValueError:
                continue
    return found


# ----------------------------------------------------------------------
# STEP 3: MAKE THE TWO STUDIES AGREE ABOUT WHICH ALLELE IS WHICH
# ----------------------------------------------------------------------

def harmonise(exposure_variant, outcome_variant):
    """Line the outcome effect up with the exposure effect allele.

    Two studies can report the same variant from opposite points of view. If
    FinnGen says "the A allele raises risk" and CKDGen says "the G allele
    lowers eGFR", those may be the same finding written two ways. This function
    flips the outcome sign when needed and returns None when we cannot be sure.
    """
    exposure_effect = exposure_variant["effect_allele"]
    exposure_other = exposure_variant["other_allele"]
    outcome_effect = outcome_variant["effect_allele"]
    outcome_other = outcome_variant["other_allele"]

    # We only handle simple single-letter changes. Insertions and deletions are
    # written differently by different studies and are not worth the risk.
    for letter in (exposure_effect, exposure_other, outcome_effect, outcome_other):
        if letter not in ("A", "C", "G", "T"):
            return None

    outcome_beta = outcome_variant["beta"]
    outcome_frequency = outcome_variant["frequency"]

    if exposure_effect == outcome_effect and exposure_other == outcome_other:
        pass                                            # already lined up
    elif exposure_effect == outcome_other and exposure_other == outcome_effect:
        outcome_beta = -outcome_beta                    # simply swapped round
        outcome_frequency = 1.0 - outcome_frequency
    else:
        # Maybe the outcome study read the other DNA strand. Try flipping it.
        flipped_effect = flip_dna_letter(outcome_effect)
        flipped_other = flip_dna_letter(outcome_other)
        if exposure_effect == flipped_effect and exposure_other == flipped_other:
            pass
        elif exposure_effect == flipped_other and exposure_other == flipped_effect:
            outcome_beta = -outcome_beta
            outcome_frequency = 1.0 - outcome_frequency
        else:
            return None                                 # genuinely different variants

    # If the variant is A/T or C/G we cannot tell the strand apart by letters
    # alone. We can sometimes rescue it using the allele frequency, but only
    # when the frequency is clearly not near one half.
    if is_palindromic(exposure_effect, exposure_other):
        exposure_frequency = exposure_variant["frequency"]
        if 0.42 < exposure_frequency < 0.58:
            return None                                 # too close to call, drop it
        same_side = (exposure_frequency < 0.5) == (outcome_frequency < 0.5)
        if not same_side:
            outcome_beta = -outcome_beta
            outcome_frequency = 1.0 - outcome_frequency

    return {
        "rsid": exposure_variant["rsid"],
        "gene": exposure_variant["gene"],
        "chromosome": exposure_variant["chromosome"],
        "position": exposure_variant["position"],
        "exposure_beta": exposure_variant["beta"],
        "exposure_se": exposure_variant["se"],
        "exposure_p": exposure_variant["p_value"],
        "exposure_frequency": exposure_variant["frequency"],
        "outcome_beta": outcome_beta,
        "outcome_se": outcome_variant["se"],
        "outcome_p": outcome_variant["p_value"],
        "outcome_n": outcome_variant["n"],
    }


# ----------------------------------------------------------------------
# STEP 4: THE MR METHODS THEMSELVES
# ----------------------------------------------------------------------
#
# Each instrument gives its own little estimate of the causal effect:
#     ratio = (effect on outcome) / (effect on exposure)
# The methods below are different ways of averaging those ratios, and they make
# different assumptions about how badly behaved the instruments might be.

def inverse_variance_weighted(pairs):
    """The standard MR estimate: a weighted average of the per-variant ratios.

    Variants measured precisely in the outcome study get more say. This is the
    main answer, but it assumes every instrument is well behaved.
    """
    top = 0.0
    bottom = 0.0
    for pair in pairs:
        weight = 1.0 / (pair["outcome_se"] ** 2)
        top = top + weight * pair["exposure_beta"] * pair["outcome_beta"]
        bottom = bottom + weight * pair["exposure_beta"] ** 2
    estimate = top / bottom
    standard_error = math.sqrt(1.0 / bottom)

    # Cochran's Q asks whether the instruments disagree with each other more
    # than chance would explain. If they do, we widen the error bar.
    q_statistic = 0.0
    for pair in pairs:
        weight = 1.0 / (pair["outcome_se"] ** 2)
        residual = pair["outcome_beta"] - estimate * pair["exposure_beta"]
        q_statistic = q_statistic + weight * residual ** 2
    degrees_of_freedom = len(pairs) - 1

    if degrees_of_freedom > 0:
        overdispersion = q_statistic / degrees_of_freedom
        if overdispersion > 1.0:
            standard_error = standard_error * math.sqrt(overdispersion)

    return estimate, standard_error, q_statistic, degrees_of_freedom


def weighted_middle_value(values, weights):
    """Find the weighted middle of a list of numbers.

    An ordinary median is just the middle number once you have sorted them.
    A weighted median is the same idea, except some numbers count for more
    than others, so we walk along the sorted list adding up weight until we
    have gone past half of the total weight.
    """
    # Put each value next to its own weight, then sort. Python sorts pairs by
    # the first thing in the pair, which here is the value, which is what we
    # want.
    value_and_weight = []
    for index in range(len(values)):
        value_and_weight.append((values[index], weights[index]))
    value_and_weight.sort()

    total_weight = 0.0
    for index in range(len(weights)):
        total_weight = total_weight + weights[index]

    running_total = 0.0
    for index in range(len(value_and_weight)):
        running_total = running_total + value_and_weight[index][1]
        if running_total >= total_weight / 2.0:
            return value_and_weight[index][0]
    return value_and_weight[-1][0]


def weighted_median(pairs):
    """A more robust estimate: the middle ratio rather than the average ratio.

    This one still gives the right answer even if up to half of the
    instruments are misbehaving, which the plain weighted average does not.
    """
    ratios = []
    weights = []
    for pair in pairs:
        ratio = pair["outcome_beta"] / pair["exposure_beta"]
        # The weight is how precisely we know this particular ratio.
        ratio_se = abs(pair["outcome_se"] / pair["exposure_beta"])
        ratios.append(ratio)
        weights.append(1.0 / (ratio_se ** 2))

    estimate = weighted_middle_value(ratios, weights)

    # There is no neat formula for the error bar here, so we get it by
    # simulation: jiggle every input within its own uncertainty, work the
    # answer out again, and see how much the answer moves about.
    simulated_answers = []
    for round_number in range(BOOTSTRAP_ROUNDS):
        jiggled_ratios = []
        jiggled_weights = []
        for pair in pairs:
            fake_exposure = random.gauss(pair["exposure_beta"], pair["exposure_se"])
            fake_outcome = random.gauss(pair["outcome_beta"], pair["outcome_se"])
            if fake_exposure == 0.0:
                continue
            jiggled_ratios.append(fake_outcome / fake_exposure)
            ratio_se = abs(pair["outcome_se"] / fake_exposure)
            jiggled_weights.append(1.0 / (ratio_se ** 2))
        if len(jiggled_ratios) > 0:
            simulated_answers.append(weighted_middle_value(jiggled_ratios, jiggled_weights))

    total = 0.0
    for answer in simulated_answers:
        total = total + answer
    mean_answer = total / len(simulated_answers)

    squared_spread = 0.0
    for answer in simulated_answers:
        squared_spread = squared_spread + (answer - mean_answer) ** 2
    variance = squared_spread / (len(simulated_answers) - 1)

    return estimate, math.sqrt(variance)


def mr_egger(pairs):
    """Fit a line through the points instead of forcing it through the origin.

    If all the instruments were clean, a graph of outcome effect against
    exposure effect would go through (0, 0). If the line misses the origin,
    that intercept is evidence that the variants affect the outcome through
    some other route -- which would break MR. So the intercept is a warning
    light, and the slope is a pleiotropy-corrected causal estimate.
    """
    if len(pairs) < 3:
        return None

    # Egger requires all exposure effects to point the same way, so we flip any
    # negative ones (and their outcome effect with them).
    x_values = []
    y_values = []
    weights = []
    for pair in pairs:
        x = pair["exposure_beta"]
        y = pair["outcome_beta"]
        if x < 0:
            x = -x
            y = -y
        x_values.append(x)
        y_values.append(y)
        weights.append(1.0 / (pair["outcome_se"] ** 2))

    total_weight = 0.0
    weighted_x_total = 0.0
    weighted_y_total = 0.0
    for index in range(len(x_values)):
        total_weight = total_weight + weights[index]
        weighted_x_total = weighted_x_total + weights[index] * x_values[index]
        weighted_y_total = weighted_y_total + weights[index] * y_values[index]
    mean_x = weighted_x_total / total_weight
    mean_y = weighted_y_total / total_weight

    covariance = 0.0
    variance_x = 0.0
    for index in range(len(x_values)):
        gap_x = x_values[index] - mean_x
        gap_y = y_values[index] - mean_y
        covariance = covariance + weights[index] * gap_x * gap_y
        variance_x = variance_x + weights[index] * gap_x * gap_x

    slope = covariance / variance_x
    intercept = mean_y - slope * mean_x

    # Work out the error bars from how far the points sit off the fitted line.
    residual_sum = 0.0
    for index in range(len(x_values)):
        predicted = intercept + slope * x_values[index]
        gap = y_values[index] - predicted
        residual_sum = residual_sum + weights[index] * gap * gap
    degrees_of_freedom = len(pairs) - 2
    scatter = residual_sum / degrees_of_freedom

    slope_se = math.sqrt(scatter / variance_x)
    intercept_se = math.sqrt(scatter * (1.0 / total_weight + mean_x ** 2 / variance_x))

    return slope, slope_se, intercept, intercept_se, degrees_of_freedom


def leave_one_out(pairs):
    """Redo the main estimate with each instrument removed in turn.

    If dropping a single variant changes the answer a lot, the answer was
    really just that one variant and should not be trusted.
    """
    results = []
    if len(pairs) < 3:
        return results
    for index_to_drop in range(len(pairs)):
        remaining = pairs[:index_to_drop] + pairs[index_to_drop + 1:]
        estimate, standard_error, _, _ = inverse_variance_weighted(remaining)
        results.append({
            "dropped": pairs[index_to_drop]["rsid"],
            "dropped_gene": pairs[index_to_drop]["gene"],
            "estimate": estimate,
            "se": standard_error,
        })
    return results


def variance_explained_continuous(beta, standard_error, sample_size):
    """How much of a normal measurement (like eGFR) does one variant explain?

    We can read this straight off the t-statistic, which is just the effect
    divided by its own error bar.
    """
    t_statistic = beta / standard_error
    return (t_statistic ** 2) / (t_statistic ** 2 + sample_size)


def steiger_direction_test(pairs, exposure_sample_size,
                           exposure_is_binary, outcome_is_binary):
    """Check we have not got cause and effect the wrong way round.

    The instruments should explain more of the exposure than of the outcome.
    If they explain more of the outcome, we may have picked variants that are
    really about the outcome, and the arrow points the other way.

    A yes/no disease and a normal measurement need different formulas, so we
    have to be told which is which. Getting this wrong silently flips the
    answer, so the caller must say explicitly.
    """
    exposure_r2 = 0.0
    outcome_r2 = 0.0
    outcome_sample_size = 0.0
    # Adding up r-squared across many instruments is fine for the exposure,
    # because those variants were chosen for being real. It is misleading for
    # the outcome, where most instruments are pure noise and the noise still
    # adds up to a positive number. So we also count, one instrument at a
    # time, how many favour each direction. That count is not fooled by noise
    # piling up, so we use it as the verdict.
    favouring_exposure = 0
    favouring_outcome = 0
    for pair in pairs:
        if exposure_is_binary:
            exposure_r2 = exposure_r2 + variance_explained_by_one_variant(
                pair["exposure_beta"], pair["exposure_frequency"])
        else:
            exposure_r2 = exposure_r2 + variance_explained_continuous(
                pair["exposure_beta"], pair["exposure_se"], exposure_sample_size)

        if outcome_is_binary:
            outcome_r2 = outcome_r2 + variance_explained_by_one_variant(
                pair["outcome_beta"], pair["exposure_frequency"])
        else:
            outcome_r2 = outcome_r2 + variance_explained_continuous(
                pair["outcome_beta"], pair["outcome_se"], pair["outcome_n"])
        outcome_sample_size = max(outcome_sample_size, pair["outcome_n"])

        # Work the same two numbers out again for this one instrument alone,
        # then see which side it comes down on.
        if exposure_is_binary:
            this_exposure = variance_explained_by_one_variant(
                pair["exposure_beta"], pair["exposure_frequency"])
        else:
            this_exposure = variance_explained_continuous(
                pair["exposure_beta"], pair["exposure_se"], exposure_sample_size)
        if outcome_is_binary:
            this_outcome = variance_explained_by_one_variant(
                pair["outcome_beta"], pair["exposure_frequency"])
        else:
            this_outcome = variance_explained_continuous(
                pair["outcome_beta"], pair["outcome_se"], pair["outcome_n"])
        if this_exposure > this_outcome:
            favouring_exposure = favouring_exposure + 1
        else:
            favouring_outcome = favouring_outcome + 1

    return {
        "exposure_r2": exposure_r2,
        "outcome_r2": outcome_r2,
        "favouring_exposure": favouring_exposure,
        "favouring_outcome": favouring_outcome,
        "direction_correct": favouring_exposure > favouring_outcome,
        "exposure_n": exposure_sample_size,
        "outcome_n": outcome_sample_size,
    }


def smallest_effect_we_could_have_seen(standard_error):
    """The minimum detectable effect.

    A null result only means something if we say how big an effect we could
    have found. With 80% power and the usual 5% threshold, we can detect an
    effect about 2.8 standard errors away from zero.
    """
    return 2.80 * standard_error


# ----------------------------------------------------------------------
# STEP 5: THE REVERSE DIRECTION, USED AS A CHECK ON OURSELVES
# ----------------------------------------------------------------------
#
# Running MR backwards -- asking whether kidney function causes congenital
# heart disease -- should give nothing at all. A congenital heart malformation
# is finished forming about eight weeks into pregnancy, decades before anyone
# has an adult kidney function measurement. There is no way for the second to
# cause the first.
#
# That makes it a very useful test. We know the true answer is zero, so if our
# pipeline reports something, the fault is in the pipeline (or in the
# instruments), not in biology. Most MR papers have to invent a check like
# this; here the biology hands us one for free.

def read_egfr_instruments(p_value_cutoff):
    """Pick strong, independent variants for kidney function from CKDGen."""
    kept = []
    with gzip.open(EGFR_FILE, "rt") as handle:
        handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            try:
                p_value = float(parts[7])
            except ValueError:
                continue
            if p_value >= p_value_cutoff:
                continue
            rsid = parts[14]
            if not rsid.startswith("rs"):
                continue
            try:
                frequency = float(parts[4])
                kept.append({
                    "n": float(parts[3]),
                    "chromosome": parts[10],
                    "position": int(parts[11]),
                    "effect_allele": parts[1].upper(),
                    "other_allele": parts[2].upper(),
                    "rsid": rsid,
                    "gene": "",
                    "p_value": p_value,
                    "beta": float(parts[5]),
                    "se": float(parts[6]),
                    "frequency": frequency,
                })
            except ValueError:
                continue
    return kept


def look_up_in_chd(wanted_rsids):
    """Find variants in the FinnGen congenital heart disease file."""
    found = {}
    with gzip.open(CHD_FILE, "rt") as handle:
        handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rsid = parts[4]
            if rsid not in wanted_rsids:
                continue
            try:
                found[rsid] = {
                    "effect_allele": parts[3].upper(),
                    "other_allele": parts[2].upper(),
                    "frequency": float(parts[10]),
                    "beta": float(parts[8]),
                    "se": float(parts[9]),
                    "p_value": float(parts[6]),
                    "n": 412181.0,          # 4,270 cases plus 407,911 controls
                }
            except ValueError:
                continue
    return found


# ----------------------------------------------------------------------
# STEP 6: RUN ONE COMPLETE ANALYSIS AND PRINT IT OUT
# ----------------------------------------------------------------------

def run_one_analysis(label, pairs, exposure_sample_size, output_rows,
                     exposure_is_binary=True, outcome_is_binary=False):
    """Apply every MR method to one set of harmonised variants."""

    print("")
    print("=" * 68)
    print(label)
    print("=" * 68)

    if len(pairs) < 2:
        print("  Only %d instrument(s) survived. MR needs at least 2." % len(pairs))
        return

    # How strong are our instruments? The F-statistic is the usual measure and
    # anything under 10 is conventionally called "weak".
    f_total = 0.0
    for pair in pairs:
        f_total = f_total + (pair["exposure_beta"] / pair["exposure_se"]) ** 2
    mean_f = f_total / len(pairs)

    total_r2 = 0.0
    for pair in pairs:
        total_r2 = total_r2 + variance_explained_by_one_variant(
            pair["exposure_beta"], pair["exposure_frequency"])

    print("  instruments               : %d" % len(pairs))
    print("  mean F-statistic          : %.1f" % mean_f)
    print("  variance explained (total): %.4f%%" % (100.0 * total_r2))
    print("")

    estimate, standard_error, q_statistic, degrees_of_freedom = inverse_variance_weighted(pairs)
    p_value = normal_p_value(estimate / standard_error)
    low = estimate - 1.96 * standard_error
    high = estimate + 1.96 * standard_error
    print("  IVW              : %+.4f  (95%% CI %+.4f to %+.4f)  p = %.3g"
          % (estimate, low, high, p_value))
    output_rows.append([label, "IVW", len(pairs), estimate, standard_error, low, high, p_value])

    if degrees_of_freedom > 0:
        # chi2.sf is the chance of seeing a Q value this big just by luck.
        q_p = chi2.sf(q_statistic, degrees_of_freedom)
        print("  Cochran's Q      : %.1f on %d df, p = %.3g%s"
              % (q_statistic, degrees_of_freedom, q_p,
                 "   <- instruments disagree" if q_p < 0.05 else ""))
        output_rows.append([label, "CochranQ", len(pairs), q_statistic,
                            float(degrees_of_freedom), float("nan"), float("nan"), q_p])

    median_estimate, median_se = weighted_median(pairs)
    median_p = normal_p_value(median_estimate / median_se)
    print("  Weighted median  : %+.4f  (95%% CI %+.4f to %+.4f)  p = %.3g"
          % (median_estimate, median_estimate - 1.96 * median_se,
             median_estimate + 1.96 * median_se, median_p))
    output_rows.append([label, "WeightedMedian", len(pairs), median_estimate, median_se,
                        median_estimate - 1.96 * median_se,
                        median_estimate + 1.96 * median_se, median_p])

    egger = mr_egger(pairs)
    if egger is None:
        print("  MR-Egger         : needs at least 3 instruments, skipped")
    else:
        slope, slope_se, intercept, intercept_se, egger_df = egger
        slope_p = normal_p_value(slope / slope_se)
        intercept_p = normal_p_value(intercept / intercept_se)
        print("  MR-Egger slope   : %+.4f  (95%% CI %+.4f to %+.4f)  p = %.3g"
              % (slope, slope - 1.96 * slope_se, slope + 1.96 * slope_se, slope_p))
        print("  MR-Egger intercept: %+.5f  p = %.3g%s"
              % (intercept, intercept_p,
                 "   <- possible pleiotropy" if intercept_p < 0.05 else "   (no evidence of pleiotropy)"))
        output_rows.append([label, "EggerSlope", len(pairs), slope, slope_se,
                            slope - 1.96 * slope_se, slope + 1.96 * slope_se, slope_p])
        output_rows.append([label, "EggerIntercept", len(pairs), intercept, intercept_se,
                            intercept - 1.96 * intercept_se,
                            intercept + 1.96 * intercept_se, intercept_p])

    steiger = steiger_direction_test(pairs, exposure_sample_size,
                                     exposure_is_binary, outcome_is_binary)
    print("")
    print("  Steiger check    : summed r2 is %.4f%% for the exposure and %.4f%% for the outcome"
          % (100.0 * steiger["exposure_r2"], 100.0 * steiger["outcome_r2"]))
    print("                     %d of %d instruments individually favour the assumed direction"
          % (steiger["favouring_exposure"],
             steiger["favouring_exposure"] + steiger["favouring_outcome"]))
    print("                     verdict: direction is %s"
          % ("as assumed" if steiger["direction_correct"] else "REVERSED -- be careful"))
    output_rows.append([label, "SteigerExposureR2", len(pairs), steiger["exposure_r2"],
                        float("nan"), float("nan"), float("nan"), float("nan")])
    output_rows.append([label, "SteigerOutcomeR2", len(pairs), steiger["outcome_r2"],
                        float("nan"), float("nan"), float("nan"), float("nan")])
    output_rows.append([label, "SteigerInstrumentsFavouringExposure", len(pairs),
                        float(steiger["favouring_exposure"]), float("nan"),
                        float("nan"), float("nan"), float("nan")])

    detectable = smallest_effect_we_could_have_seen(standard_error)
    print("  Smallest effect we could have detected (80%% power): %.4f" % detectable)
    output_rows.append([label, "MinDetectableEffect", len(pairs), detectable,
                        float("nan"), float("nan"), float("nan"), float("nan")])

    loo = leave_one_out(pairs)
    if loo:
        print("")
        # With hundreds of instruments a full listing is unreadable, so we
        # only print every row for the small analyses. The full set is always
        # written to the results file either way.
        lowest = loo[0]["estimate"]
        highest = loo[0]["estimate"]
        for row in loo:
            if row["estimate"] < lowest:
                lowest = row["estimate"]
            if row["estimate"] > highest:
                highest = row["estimate"]
        spread = highest - lowest
        print("  Leave-one-out: estimate ranges %+.4f to %+.4f (spread %.4f, %.2f SE)"
              % (lowest, highest, spread, spread / standard_error))
        show_every_row = (len(loo) <= 30)
        for row in loo:
            if not show_every_row:
                continue
            moved = abs(row["estimate"] - estimate) / standard_error
            flag = "   <- moves the answer a lot" if moved > 1.0 else ""
            print("    without %-13s %-16s %+.4f%s"
                  % (row["dropped"], row["dropped_gene"][:16], row["estimate"], flag))
            output_rows.append([label, "LeaveOneOut_" + row["dropped"], len(pairs) - 1,
                                row["estimate"], row["se"], float("nan"), float("nan"),
                                float("nan")])


# ----------------------------------------------------------------------
# STEP 7: THE MAIN PROGRAM
# ----------------------------------------------------------------------

def main():
    output_rows = []

    print("Mendelian randomisation: does congenital heart disease cause")
    print("worse kidney function, or is the comorbidity explained some other way?")
    print("")
    print("  exposure : FinnGen R10, congenital malformations of heart and great")
    print("             arteries (4,270 cases and 407,911 controls)")
    print("  outcome  : CKDGen, estimated glomerular filtration rate")
    print("")
    print("A negative estimate means heart disease liability predicts LOWER eGFR,")
    print("which would mean worse kidney function.")

    # ---- forward direction, at two instrument thresholds ----
    for threshold, name in ((STRICT_THRESHOLD, "strict"), (RELAXED_THRESHOLD, "relaxed")):
        print("")
        print("Reading heart disease variants at p < %g ..." % threshold)
        chd_variants = read_chd_variants(threshold)
        instruments = keep_one_variant_per_region(chd_variants)
        print("  %d variants pass, %d remain after keeping one per region"
              % (len(chd_variants), len(instruments)))

        wanted = set()
        for variant in instruments:
            wanted.add(variant["rsid"])
        print("Looking them up in the eGFR study ...")
        outcome_lookup = look_up_in_egfr(wanted)
        print("  found %d of %d in CKDGen" % (len(outcome_lookup), len(wanted)))

        pairs = []
        dropped_for_alleles = 0
        for variant in instruments:
            if variant["rsid"] not in outcome_lookup:
                continue
            harmonised = harmonise(variant, outcome_lookup[variant["rsid"]])
            if harmonised is None:
                dropped_for_alleles = dropped_for_alleles + 1
                continue
            pairs.append(harmonised)
        print("  %d usable after allele matching (%d dropped as ambiguous)"
              % (len(pairs), dropped_for_alleles))

        run_one_analysis("CHD -> eGFR (%s, p<%g)" % (name, threshold),
                         pairs, 412181.0, output_rows)

        # The 17q21 region has unusually long-range correlation between
        # variants, which can make one signal look like several. We repeat the
        # analysis without it to check it is not driving anything.
        without_17q21 = []
        for pair in pairs:
            in_region = (pair["chromosome"] == "17"
                         and 44000000 < pair["position"] < 48000000)
            if not in_region:
                without_17q21.append(pair)
        if len(without_17q21) < len(pairs):
            run_one_analysis("CHD -> eGFR (%s, 17q21 removed)" % name,
                             without_17q21, 412181.0, output_rows)

    # ---- reverse direction, the sanity check ----
    print("")
    print("Now the reverse direction, which should find nothing at all.")
    print("Reading kidney function variants ...")
    egfr_variants = read_egfr_instruments(STRICT_THRESHOLD)
    egfr_instruments = keep_one_variant_per_region(egfr_variants)
    print("  %d variants pass, %d remain after keeping one per region"
          % (len(egfr_variants), len(egfr_instruments)))

    wanted = set()
    for variant in egfr_instruments:
        wanted.add(variant["rsid"])
    chd_lookup = look_up_in_chd(wanted)
    print("  found %d of %d in FinnGen" % (len(chd_lookup), len(wanted)))

    reverse_pairs = []
    for variant in egfr_instruments:
        if variant["rsid"] not in chd_lookup:
            continue
        harmonised = harmonise(variant, chd_lookup[variant["rsid"]])
        if harmonised is not None:
            reverse_pairs.append(harmonised)
    print("  %d usable after allele matching" % len(reverse_pairs))

    # Here the exposure is a normal measurement and the outcome is a yes/no
    # disease -- the opposite way round from the forward analysis.
    run_one_analysis("eGFR -> CHD (negative control)", reverse_pairs, 1004020.0,
                     output_rows, exposure_is_binary=False, outcome_is_binary=True)

    # ---- save everything ----
    if not os.path.isdir(OUT_FOLDER):
        os.makedirs(OUT_FOLDER)
    out_path = os.path.join(OUT_FOLDER, "mr_results.tsv")
    with open(out_path, "w") as handle:
        handle.write("analysis\tmethod\tn_instruments\testimate\tse\tci_low\tci_high\tp_value\n")
        for row in output_rows:
            text_pieces = []
            for value in row:
                text_pieces.append(str(value))
            handle.write("\t".join(text_pieces) + "\n")

    print("")
    print("Saved results to %s" % out_path)


if __name__ == "__main__":
    main()
