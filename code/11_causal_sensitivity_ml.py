# 11_causal_sensitivity_ml.py
#
# WHAT THIS SCRIPT DOES
# ---------------------
# In 03_hpo_analysis.py we adjusted for "syndromic pleiotropy" using ordinary
# logistic regression with two extra columns: how many organ systems a disease
# affects, and how many phenotype terms it has.
#
# That adjustment carries a lot of weight in this project, because it is the
# step that turns a huge unadjusted odds ratio (5.81) into a modest adjusted
# one. So it is worth asking two awkward questions about it.
#
#   Question 1. Logistic regression assumes the effect of "number of organ
#               systems" is a straight line on the log-odds scale. What if the
#               real relationship bends?
#
#   Question 2. "How many organ systems" is a crude summary. Two diseases can
#               both affect 8 systems and yet be completely different. What if
#               we adjust for WHICH systems, not just how many?
#
# This script answers both. It uses a machine learning model (gradient
# boosting) that can bend and interact however it likes, and it feeds that
# model one column per organ system instead of a single count.
#
# IMPORTANT: the machine learning here is NOT used to predict anything for its
# own sake. Prediction accuracy is not the goal and we never report it. The
# model is only a flexible way to describe the confounders, and it sits inside
# a formula (AIPW, explained below) that still targets the same causal
# quantity as before. Swapping in a fancier model does not change the
# question; it only removes an assumption about the shape of the answer.
#
# TWO THINGS WE MUST NOT ADJUST FOR
# ---------------------------------
# We deliberately leave out the cardiovascular system column and the
# genitourinary system column.
#   - Cardiovascular is what the exposure IS. Adjusting for it would adjust
#     away the very thing we are studying.
#   - Genitourinary is where the outcome lives. Adjusting for it would block
#     the effect we are trying to measure.
# Putting either one in would guarantee a null result for the wrong reason.
#
# To run:   python3 11_causal_sensitivity_ml.py

import os

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import norm
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold

RAW = "../rawdata"
OUT = "../data/results"

# How many pieces to cut the data into for cross-fitting (explained below).
NUMBER_OF_FOLDS = 5

# Propensity scores very close to 0 or 1 make the AIPW formula explode,
# because we end up dividing by a tiny number. We pull them back inside
# these limits, which is standard practice.
SMALLEST_ALLOWED_PROPENSITY = 0.01
LARGEST_ALLOWED_PROPENSITY = 0.99

# Fixing the seed means the folds are cut the same way every time, so you get
# the same numbers we did.
RANDOM_SEED = 20260828

# The two organ systems we must not adjust for, as explained above.
CARDIOVASCULAR_SYSTEM = "HP:0001626"
GENITOURINARY_SYSTEM = "HP:0000119"


# ======================================================================
# PART 1: read the ontology
# ======================================================================
#
# These two functions are the same ones used in 03_hpo_analysis.py. They are
# repeated here so this script can be read on its own.

def read_ontology(filename):
    """Read the HPO file and work out which term sits under which.

    Returns three things: the name of every term, a dictionary saying which
    terms are directly below each term, and the set of retired terms.
    """
    names = {}
    children_of = {}
    obsolete_terms = set()

    current_id = None
    current_name = None
    current_is_obsolete = False
    current_parents = []

    handle = open(filename)
    for raw_line in handle:
        line = raw_line.strip()

        if line == "[Term]":
            # We have reached a new term, so file away the one we just read.
            if current_id is not None:
                names[current_id] = current_name
                if current_is_obsolete:
                    obsolete_terms.add(current_id)
                for parent in current_parents:
                    if parent not in children_of:
                        children_of[parent] = set()
                    children_of[parent].add(current_id)
            current_id = None
            current_name = None
            current_is_obsolete = False
            current_parents = []

        elif line.startswith("id: HP:"):
            current_id = line[4:]
        elif line.startswith("name: "):
            current_name = line[6:]
        elif line.startswith("is_obsolete: true"):
            current_is_obsolete = True
        elif line.startswith("is_a: "):
            # A line looks like "is_a: HP:0001671 ! Abnormal cardiac septum".
            # We only want the code, so we cut at the exclamation mark.
            parent_part = line[6:]
            parent_id = parent_part.split("!")[0].strip()
            current_parents.append(parent_id)
    handle.close()

    # The very last term in the file still needs filing away.
    if current_id is not None:
        names[current_id] = current_name
        if current_is_obsolete:
            obsolete_terms.add(current_id)
        for parent in current_parents:
            if parent not in children_of:
                children_of[parent] = set()
            children_of[parent].add(current_id)

    return names, children_of, obsolete_terms


def get_all_terms_below(start_term, children_of, obsolete_terms):
    """Collect a term and everything underneath it, however deep.

    "Abnormality of the kidney" has children, and those have children, and so
    on. We want all of them. We keep a list of terms still to look at, and
    keep going until the list is empty.
    """
    found = set()
    still_to_check = [start_term]
    while len(still_to_check) > 0:
        term = still_to_check.pop()
        if term in found:
            continue
        if term in obsolete_terms:
            continue
        found.add(term)
        if term in children_of:
            for child in children_of[term]:
                still_to_check.append(child)
    return found


# ======================================================================
# PART 2: the AIPW estimator
# ======================================================================
#
# AIPW stands for "augmented inverse probability weighting". The name is
# horrible but the idea is not.
#
# We want to know: if a disease affects the heart, how much more likely is it
# to affect the kidney, once we have accounted for how syndromic it is?
#
# We build two helpers:
#   - the OUTCOME model guesses, for each disease, whether it would affect the
#     kidney if it did affect the heart, and if it did not.
#   - the PROPENSITY model guesses how likely each disease was to affect the
#     heart in the first place.
#
# AIPW combines them. If either one of the two is right, the answer is right.
# That is why it is called "doubly robust" -- you get two chances.
#
# The last ingredient is CROSS-FITTING. If we train a flexible model on the
# same rows we then use it on, it will have partly memorised them, and that
# memorisation leaks into the final number. So we cut the data into 5 pieces,
# and for each piece we train on the OTHER 4 and only make guesses about the
# piece the model has never seen.

def cross_fitted_predictions(features, treatment, outcome):
    """Get honest out-of-sample guesses from both helper models."""
    number_of_rows = len(treatment)
    propensity = np.zeros(number_of_rows)
    outcome_if_treated = np.zeros(number_of_rows)
    outcome_if_untreated = np.zeros(number_of_rows)

    splitter = StratifiedKFold(n_splits=NUMBER_OF_FOLDS, shuffle=True,
                               random_state=RANDOM_SEED)

    # We stratify on treatment and outcome together so every fold contains
    # some of the rare "heart AND kidney" diseases.
    stratify_on = treatment * 2 + outcome

    for train_rows, test_rows in splitter.split(features, stratify_on):
        train_features = features[train_rows]
        test_features = features[test_rows]

        # --- helper 1: how likely was this disease to affect the heart? ---
        propensity_model = HistGradientBoostingClassifier(
            max_iter=200, learning_rate=0.1, max_depth=4,
            random_state=RANDOM_SEED)
        propensity_model.fit(train_features, treatment[train_rows])
        propensity[test_rows] = propensity_model.predict_proba(test_features)[:, 1]

        # --- helper 2: how likely is a kidney problem, given the heart? ---
        # We add the treatment as one more column, then ask the fitted model
        # twice: once pretending every disease affects the heart, once
        # pretending none of them do.
        train_with_treatment = np.column_stack([train_features, treatment[train_rows]])
        outcome_model = HistGradientBoostingClassifier(
            max_iter=200, learning_rate=0.1, max_depth=4,
            random_state=RANDOM_SEED)
        outcome_model.fit(train_with_treatment, outcome[train_rows])

        all_treated = np.ones(len(test_rows))
        all_untreated = np.zeros(len(test_rows))
        test_as_treated = np.column_stack([test_features, all_treated])
        test_as_untreated = np.column_stack([test_features, all_untreated])

        outcome_if_treated[test_rows] = outcome_model.predict_proba(test_as_treated)[:, 1]
        outcome_if_untreated[test_rows] = outcome_model.predict_proba(test_as_untreated)[:, 1]

    return propensity, outcome_if_treated, outcome_if_untreated


def aipw_estimate(treatment, outcome, propensity,
                  outcome_if_treated, outcome_if_untreated):
    """Combine the two helper models into one causal estimate.

    The number we get back is a RISK DIFFERENCE: how many percentage points
    more likely a kidney problem is when the heart is affected. We also work
    out an error bar from how much the per-disease contributions vary.
    """
    # Keep the propensity away from 0 and 1 so we never divide by almost zero.
    safe_propensity = np.clip(propensity, SMALLEST_ALLOWED_PROPENSITY,
                              LARGEST_ALLOWED_PROPENSITY)

    # Each disease contributes one number. The first half of each line is the
    # model's guess; the second half is a correction based on how wrong the
    # guess turned out to be for that disease.
    treated_part = (outcome_if_treated
                    + treatment * (outcome - outcome_if_treated) / safe_propensity)
    untreated_part = (outcome_if_untreated
                      + (1 - treatment) * (outcome - outcome_if_untreated)
                      / (1 - safe_propensity))

    per_disease_contribution = treated_part - untreated_part

    estimate = float(np.mean(per_disease_contribution))
    standard_error = float(np.std(per_disease_contribution, ddof=1)
                           / np.sqrt(len(per_disease_contribution)))
    return estimate, standard_error, safe_propensity


# ======================================================================
# PART 3: the main program
# ======================================================================

def main():
    print("Reading the ontology ...")
    names, children_of, obsolete_terms = read_ontology(os.path.join(RAW, "hp.obo"))

    # The same phenotype definitions used in 03_hpo_analysis.py.
    heart_terms = get_all_terms_below("HP:0001627", children_of, obsolete_terms)
    kidney_terms = get_all_terms_below("HP:0012210", children_of, obsolete_terms)
    # Chronic kidney disease is spread over two branches of the ontology, so
    # we take everything under both and join the two sets together with "|".
    ckd_terms = get_all_terms_below("HP:0012622", children_of, obsolete_terms)
    ckd_terms = ckd_terms | get_all_terms_below("HP:0000083", children_of, obsolete_terms)

    print("Reading the gene to phenotype annotations ...")
    annotations = pd.read_csv(os.path.join(RAW, "genes_to_phenotype.txt"), sep="\t")

    terms_per_disease = annotations.groupby("disease_id")["hpo_id"].apply(set)
    print("  %d diseases" % len(terms_per_disease))

    # The big organ systems are the direct children of "Phenotypic abnormality".
    organ_system_ids = children_of["HP:0000118"]
    organ_system_terms = {}
    for system_id in organ_system_ids:
        organ_system_terms[system_id] = get_all_terms_below(system_id, children_of,
                                                            obsolete_terms)

    # ---- build one row per disease ----
    print("Building the table, one disease at a time ...")

    disease_ids = []
    has_chd_list = []
    has_kidney_list = []
    has_ckd_list = []
    n_terms_list = []
    n_systems_list = []

    # We keep the organ systems we are allowed to adjust for in a fixed order,
    # so every disease gets its columns in the same places.
    usable_system_ids = []
    for system_id in organ_system_ids:
        if system_id == CARDIOVASCULAR_SYSTEM:
            continue
        if system_id == GENITOURINARY_SYSTEM:
            continue
        usable_system_ids.append(system_id)
    usable_system_ids.sort()

    system_columns = []

    for disease in terms_per_disease.index:
        terms = terms_per_disease[disease]

        disease_ids.append(disease)
        has_chd_list.append(int(len(terms & heart_terms) > 0))
        has_kidney_list.append(int(len(terms & kidney_terms) > 0))
        has_ckd_list.append(int(len(terms & ckd_terms) > 0))
        n_terms_list.append(len(terms))

        # Count all systems (for the count-only comparison) but record the
        # individual yes/no flags only for the ones we are allowed to use.
        count = 0
        for system_id in organ_system_ids:
            if len(terms & organ_system_terms[system_id]) > 0:
                count = count + 1
        n_systems_list.append(count)

        this_disease_flags = []
        for system_id in usable_system_ids:
            if len(terms & organ_system_terms[system_id]) > 0:
                this_disease_flags.append(1)
            else:
                this_disease_flags.append(0)
        system_columns.append(this_disease_flags)

    treatment = np.array(has_chd_list)
    n_terms = np.array(n_terms_list, dtype=float)
    n_systems = np.array(n_systems_list, dtype=float)
    system_flags = np.array(system_columns, dtype=float)

    print("  %d diseases, %d organ-system columns kept"
          % (len(treatment), len(usable_system_ids)))
    print("  (cardiovascular and genitourinary deliberately left out)")
    print("")

    results = []

    for outcome_name, outcome_values in [("kidney anomaly", np.array(has_kidney_list)),
                                         ("chronic kidney disease", np.array(has_ckd_list))]:

        print("=" * 68)
        print("OUTCOME: %s" % outcome_name)
        print("=" * 68)

        # ---- the plain unadjusted difference, for comparison ----
        risk_if_chd = float(np.mean(outcome_values[treatment == 1]))
        risk_if_not = float(np.mean(outcome_values[treatment == 0]))
        print("  unadjusted risk with a heart malformation   : %.3f" % risk_if_chd)
        print("  unadjusted risk without one                 : %.3f" % risk_if_not)
        print("  unadjusted risk difference                  : %+.4f"
              % (risk_if_chd - risk_if_not))
        print("")

        # ---- the logistic adjustment used in the main analysis ----
        logistic_predictors = pd.DataFrame({
            "has_chd": treatment.astype(float),
            "n_systems": n_systems,
            "n_terms": n_terms,
        })
        logistic_predictors = sm.add_constant(logistic_predictors)
        logistic_model = sm.Logit(outcome_values.astype(float),
                                  logistic_predictors).fit(disp=0)
        logistic_odds_ratio = float(np.exp(logistic_model.params["has_chd"]))
        print("  logistic adjustment (counts only, as in 03):")
        print("    odds ratio for has_chd                    : %.3f  (p = %.3g)"
              % (logistic_odds_ratio, float(logistic_model.pvalues["has_chd"])))
        print("")

        # ---- the machine learning adjustment ----
        # First with the same two count columns, so the ONLY thing that has
        # changed is the shape of the model. Then with one column per organ
        # system, so the adjustment itself is richer.
        count_only_features = np.column_stack([n_systems, n_terms])
        full_features = np.column_stack([n_systems, n_terms, system_flags])

        for feature_name, features in [("counts only", count_only_features),
                                       ("counts plus each organ system", full_features)]:
            propensity, guess_if_treated, guess_if_untreated = cross_fitted_predictions(
                features, treatment, outcome_values)

            estimate, standard_error, safe_propensity = aipw_estimate(
                treatment, outcome_values, propensity,
                guess_if_treated, guess_if_untreated)

            low = estimate - 1.96 * standard_error
            high = estimate + 1.96 * standard_error
            # Turn the estimate into a p-value. norm.sf is the area in one
            # tail of the bell curve, and we want both tails, so we double it.
            z_score = estimate / standard_error
            p_value = float(2.0 * norm.sf(abs(z_score)))

            print("  AIPW with gradient boosting (%s):" % feature_name)
            print("    risk difference                           : %+.4f  (95%% CI %+.4f to %+.4f)"
                  % (estimate, low, high))
            print("    p-value                                   : %.3g" % p_value)

            # Positivity check. If some diseases had essentially no chance of
            # affecting the heart, the formula is being asked to imagine
            # something that never happens, and we should say so.
            fraction_clipped = float(np.mean((propensity < SMALLEST_ALLOWED_PROPENSITY)
                                             | (propensity > LARGEST_ALLOWED_PROPENSITY)))
            print("    propensity range                          : %.4f to %.4f"
                  % (float(np.min(propensity)), float(np.max(propensity))))
            print("    fraction needing clipping                 : %.4f" % fraction_clipped)
            print("")

            results.append({
                "outcome": outcome_name,
                "method": "AIPW_" + feature_name.replace(" ", "_"),
                "estimate_risk_difference": estimate,
                "se": standard_error,
                "ci_low": low,
                "ci_high": high,
                "p_value": p_value,
                "fraction_clipped": fraction_clipped,
            })

        results.append({
            "outcome": outcome_name,
            "method": "logistic_counts_only_OR",
            "estimate_risk_difference": logistic_odds_ratio,
            "se": float(logistic_model.bse["has_chd"]),
            "ci_low": float("nan"),
            "ci_high": float("nan"),
            "p_value": float(logistic_model.pvalues["has_chd"]),
            "fraction_clipped": float("nan"),
        })
        results.append({
            "outcome": outcome_name,
            "method": "unadjusted_risk_difference",
            "estimate_risk_difference": risk_if_chd - risk_if_not,
            "se": float("nan"), "ci_low": float("nan"), "ci_high": float("nan"),
            "p_value": float("nan"), "fraction_clipped": float("nan"),
        })

    if not os.path.isdir(OUT):
        os.makedirs(OUT)
    output_table = pd.DataFrame(results)
    output_path = os.path.join(OUT, "causal_ml_sensitivity.tsv")
    output_table.to_csv(output_path, sep="\t", index=False)
    print("Saved results to %s" % output_path)


if __name__ == "__main__":
    main()
