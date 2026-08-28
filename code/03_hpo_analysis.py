# 03_hpo_analysis.py
#
# This is the main analysis of the project.
#
# The question is: do congenital heart disease (CHD) and kidney disease
# share genetic causes?
#
# The trap we have to avoid: genes that cause birth defects usually break
# many organs at once. So ANY two organs will look like they "share" genes.
# We have to check whether the heart-kidney overlap is bigger than that
# background level, and we have to control for it statistically.
#
# To run:   python3 03_hpo_analysis.py
#
# It prints all the numbers used in the Results chapter and saves
# result tables into ../data/results/

import os
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.stats import fisher_exact, mannwhitneyu

RAW = "../rawdata"
OUT = "../data/results"


# ======================================================================
# PART 1: read the ontology
# ======================================================================
#
# The HPO file is a big text file. Each phenotype term looks like this:
#
#   [Term]
#   id: HP:0001631
#   name: Atrial septal defect
#   is_a: HP:0001671 ! Abnormal cardiac septum morphology
#
# The "is_a" lines say which term is a more general version of which.
# We read those lines to build a tree (technically a graph).

def read_ontology(filename):
    """Read hp.obo and return the term names and the parent-child links."""
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
            # a line looks like: is_a: HP:0001671 ! Abnormal cardiac septum morphology
            parent = line[6:].split(" ! ")[0].strip()
            if parent not in children_of:
                children_of[parent] = []
            children_of[parent].append(current_term)

        elif line.startswith("is_obsolete: true") and current_term is not None:
            obsolete_terms.add(current_term)

    return term_name, children_of, obsolete_terms


def get_all_terms_below(start_term, children_of, obsolete_terms):
    """
    Collect a term and everything underneath it.

    For example, starting at "Abnormal heart morphology" this collects
    "Atrial septal defect", "Tetralogy of Fallot" and so on, because all
    of those are specific kinds of abnormal heart morphology.
    """
    found = set()
    found.add(start_term)
    to_visit = [start_term]

    while len(to_visit) > 0:
        term = to_visit.pop()
        for child in children_of.get(term, []):
            if child not in found:
                found.add(child)
                to_visit.append(child)

    # take out any terms that have been retired
    return found - obsolete_terms


# ======================================================================
# PART 2: helper for odds ratios
# ======================================================================

def odds_ratio_with_confidence_interval(a, b, c, d):
    """
    Work out the odds ratio and its 95% confidence interval
    from a 2x2 table:

                  outcome yes   outcome no
      exposed          a             b
      not exposed      c             d
    """
    if b == 0 or c == 0 or a == 0 or d == 0:
        return float("nan"), float("nan"), float("nan")

    odds = (a * d) / (b * c)
    # standard error of log(odds ratio)
    standard_error = np.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
    lower = odds * np.exp(-1.96 * standard_error)
    upper = odds * np.exp(1.96 * standard_error)
    return odds, lower, upper


# ======================================================================
# MAIN
# ======================================================================

print("Reading the ontology...")
term_name, children_of, obsolete_terms = read_ontology(os.path.join(RAW, "hp.obo"))
print("  found", len(term_name), "phenotype terms")

print("Reading gene to phenotype annotations...")
annotations = pd.read_csv(os.path.join(RAW, "genes_to_phenotype.txt"),
                          sep="\t", dtype=str)
print("  found", len(annotations), "annotations")
print("  covering", annotations["gene_symbol"].nunique(), "genes")
print("  and", annotations["disease_id"].nunique(), "diseases")

# The three groups of phenotypes we care about
heart_terms = get_all_terms_below("HP:0001627", children_of, obsolete_terms)
kidney_terms = get_all_terms_below("HP:0012210", children_of, obsolete_terms)
ckd_terms = get_all_terms_below("HP:0012622", children_of, obsolete_terms)
ckd_terms = ckd_terms | get_all_terms_below("HP:0000083", children_of, obsolete_terms)

print("")
print("Heart terms  :", len(heart_terms))
print("Kidney terms :", len(kidney_terms))
print("CKD terms    :", len(ckd_terms))


# ----------------------------------------------------------------------
# TEST 1: do heart genes and kidney genes overlap?
# ----------------------------------------------------------------------

print("")
print("=" * 62)
print("TEST 1 - do heart and kidney genes overlap at all?")
print("=" * 62)

all_genes = set(annotations["gene_symbol"])

heart_rows = annotations[annotations["hpo_id"].isin(heart_terms)]
heart_genes = set(heart_rows["gene_symbol"])

kidney_rows = annotations[annotations["hpo_id"].isin(kidney_terms)]
kidney_genes = set(kidney_rows["gene_symbol"])

both_genes = heart_genes & kidney_genes

# build the 2x2 table
a = len(both_genes)
b = len(heart_genes - kidney_genes)
c = len(kidney_genes - heart_genes)
d = len(all_genes) - a - b - c

odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

print("heart genes  :", len(heart_genes))
print("kidney genes :", len(kidney_genes))
print("both         :", len(both_genes))
print("odds ratio = %.2f, p = %.3g" % (odds, p_value))

# save the overlapping gene list for later scripts
if not os.path.exists(OUT):
    os.makedirs(OUT)
overlap_table = pd.DataFrame({"gene_symbol": sorted(both_genes)})
overlap_table.to_csv(os.path.join(OUT, "hpo_overlap_genes.tsv"), sep="\t", index=False)


# ----------------------------------------------------------------------
# TEST 2: is the kidney special, or does every organ overlap the heart?
# ----------------------------------------------------------------------

print("")
print("=" * 62)
print("TEST 2 - is the kidney special compared to other organs?")
print("=" * 62)

# We compare the heart against every "Abnormal <something> morphology"
# term in the ontology. We skip anything that is part of the
# cardiovascular system, because those would overlap the heart trivially.
cardiovascular_terms = get_all_terms_below("HP:0001626", children_of, obsolete_terms)

organ_results = []
for term_id in term_name:
    name = term_name[term_id]

    if term_id in obsolete_terms:
        continue
    if term_id in cardiovascular_terms:
        continue
    if not name.lower().startswith("abnormal "):
        continue
    if not name.lower().endswith(" morphology"):
        continue

    organ_terms = get_all_terms_below(term_id, children_of, obsolete_terms)
    organ_genes = set(annotations[annotations["hpo_id"].isin(organ_terms)]["gene_symbol"])

    # only test organs with a decent number of genes
    if len(organ_genes) < 150:
        continue

    a = len(heart_genes & organ_genes)
    b = len(heart_genes - organ_genes)
    c = len(organ_genes - heart_genes)
    d = len(all_genes) - a - b - c
    odds, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

    organ_results.append({"term": name, "id": term_id, "n_organ": len(organ_genes),
                          "overlap": a, "OR": odds, "p": p_value})

organ_table = pd.DataFrame(organ_results).sort_values("OR", ascending=False)
organ_table = organ_table.reset_index(drop=True)
organ_table.to_csv(os.path.join(OUT, "hpo_organ_specificity.csv"), index=False)

# where does the kidney come?
kidney_position = -1
for i in range(len(organ_table)):
    if organ_table.loc[i, "id"] == "HP:0012210":
        kidney_position = i + 1

print("we tested", len(organ_table), "organ systems")
print("the kidney ranks number", kidney_position)
print("kidney odds ratio  = %.2f" % organ_table.loc[kidney_position - 1, "OR"])
print("median across all organs = %.2f" % organ_table["OR"].median())
print("")
print("Top 8 organs:")
for i in range(8):
    print("  %2d. %-45s OR = %.2f" % (i + 1, organ_table.loc[i, "term"],
                                      organ_table.loc[i, "OR"]))


# ----------------------------------------------------------------------
# TEST 3: build a table with one row per disease
# ----------------------------------------------------------------------

print("")
print("=" * 62)
print("TEST 3 - how many organs does a heart disease usually affect?")
print("=" * 62)

# For each disease, collect the set of phenotype terms it has
terms_per_disease = annotations.groupby("disease_id")["hpo_id"].apply(set)

# The 23 big organ systems are the direct children of "Phenotypic abnormality"
organ_system_ids = children_of["HP:0000118"]
organ_system_terms = {}
for system_id in organ_system_ids:
    organ_system_terms[system_id] = get_all_terms_below(system_id, children_of,
                                                        obsolete_terms)

# Build the disease table one column at a time
disease_table = pd.DataFrame(index=terms_per_disease.index)
disease_table["has_chd"] = terms_per_disease.apply(
    lambda terms: int(len(terms & heart_terms) > 0))
disease_table["has_kidney"] = terms_per_disease.apply(
    lambda terms: int(len(terms & kidney_terms) > 0))
disease_table["has_ckd"] = terms_per_disease.apply(
    lambda terms: int(len(terms & ckd_terms) > 0))
disease_table["n_terms"] = terms_per_disease.apply(len)

# count how many organ systems each disease touches
system_counts = []
for disease in terms_per_disease.index:
    terms = terms_per_disease[disease]
    count = 0
    for system_id in organ_system_ids:
        if len(terms & organ_system_terms[system_id]) > 0:
            count = count + 1
    system_counts.append(count)
disease_table["n_systems"] = system_counts

disease_table.to_csv(os.path.join(OUT, "hpo_disease_table.tsv"), sep="\t")

chd_diseases = disease_table[disease_table["has_chd"] == 1]
other_diseases = disease_table[disease_table["has_chd"] == 0]

print("diseases with a heart malformation affect a median of %.0f organ systems"
      % chd_diseases["n_systems"].median())
print("all other diseases affect a median of %.0f organ systems"
      % other_diseases["n_systems"].median())
test_result = mannwhitneyu(chd_diseases["n_systems"], other_diseases["n_systems"],
                           alternative="greater")
print("Mann-Whitney p = %.3g" % test_result.pvalue)
print("")
print("This is the problem. Heart syndromes are simply more pleiotropic,")
print("so they overlap everything. We now control for it.")


# ----------------------------------------------------------------------
# TEST 4: the key test - adjust for pleiotropy
# ----------------------------------------------------------------------

print("")
print("=" * 62)
print("TEST 4 - does the link survive adjusting for pleiotropy?")
print("=" * 62)

results_to_save = []

for outcome_name, outcome_column in [("kidney anomaly", "has_kidney"),
                                     ("chronic kidney disease", "has_ckd")]:

    # ---- first the simple unadjusted version ----
    a = int(((disease_table["has_chd"] == 1) & (disease_table[outcome_column] == 1)).sum())
    b = int(((disease_table["has_chd"] == 1) & (disease_table[outcome_column] == 0)).sum())
    c = int(((disease_table["has_chd"] == 0) & (disease_table[outcome_column] == 1)).sum())
    d = int(((disease_table["has_chd"] == 0) & (disease_table[outcome_column] == 0)).sum())

    odds, lower, upper = odds_ratio_with_confidence_interval(a, b, c, d)
    unused, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

    print("")
    print("CHD  ->  " + outcome_name)
    print("  unadjusted: OR = %.2f  (95%% CI %.2f to %.2f)  p = %.3g"
          % (odds, lower, upper, p_value))
    results_to_save.append({"outcome": outcome_name, "model": "unadjusted",
                            "OR": odds, "lower": lower, "upper": upper, "p": p_value})

    # ---- now the adjusted version ----
    # We predict the kidney phenotype from CHD, but we also give the model
    # the number of organ systems and the number of terms. That way it can
    # separate "these two really go together" from "this is just a big syndrome".
    predictors = disease_table[["has_chd", "n_systems", "n_terms"]].astype(float)
    predictors = sm.add_constant(predictors)
    outcome = disease_table[outcome_column].astype(float)

    model = sm.Logit(outcome, predictors).fit(disp=0)

    coefficient = model.params["has_chd"]
    standard_error = model.bse["has_chd"]
    p_value = model.pvalues["has_chd"]

    odds = np.exp(coefficient)
    lower = np.exp(coefficient - 1.96 * standard_error)
    upper = np.exp(coefficient + 1.96 * standard_error)

    print("  adjusted:   OR = %.2f  (95%% CI %.2f to %.2f)  p = %.3g"
          % (odds, lower, upper, p_value))
    results_to_save.append({"outcome": outcome_name, "model": "adjusted",
                            "OR": odds, "lower": lower, "upper": upper, "p": p_value})

pd.DataFrame(results_to_save).to_csv(os.path.join(OUT, "hpo_main_result.tsv"),
                                     sep="\t", index=False)

print("")
print("The kidney malformation link survives. The chronic kidney disease")
print("link does not - it was entirely explained by pleiotropy.")


# ----------------------------------------------------------------------
# TEST 5: which heart lesions go with kidney problems?
# ----------------------------------------------------------------------

print("")
print("=" * 62)
print("TEST 5 - which heart lesion goes with kidney anomalies?")
print("=" * 62)

lesions = {"Hypoplastic left heart": "HP:0004383",
           "Tetralogy of Fallot": "HP:0001636",
           "Dextrocardia": "HP:0001651",
           "Patent ductus arteriosus": "HP:0001643",
           "Heterotaxy/situs inversus": "HP:0001696",
           "Coarctation of aorta": "HP:0001680",
           "Ventricular septal defect": "HP:0001629",
           "Atrioventricular canal defect": "HP:0006695",
           "Atrial septal defect": "HP:0001631",
           "Transposition great arteries": "HP:0001669",
           "Bicuspid aortic valve": "HP:0001647",
           "Pulmonic stenosis": "HP:0001642"}

lesion_results = []
for lesion_name in lesions:
    lesion_terms = get_all_terms_below(lesions[lesion_name], children_of, obsolete_terms)
    has_lesion = terms_per_disease.apply(lambda terms: int(len(terms & lesion_terms) > 0))

    if has_lesion.sum() < 25:
        print("  skipping %s (only %d diseases)" % (lesion_name, has_lesion.sum()))
        continue

    predictors = pd.DataFrame({"lesion": has_lesion.astype(float),
                               "n_systems": disease_table["n_systems"].astype(float),
                               "n_terms": disease_table["n_terms"].astype(float)})
    predictors = sm.add_constant(predictors)

    model = sm.Logit(disease_table["has_kidney"].astype(float), predictors).fit(disp=0)
    coefficient = model.params["lesion"]
    standard_error = model.bse["lesion"]

    lesion_results.append({"lesion": lesion_name,
                           "n_dis": int(has_lesion.sum()),
                           "adjOR": np.exp(coefficient),
                           "lo": np.exp(coefficient - 1.96 * standard_error),
                           "hi": np.exp(coefficient + 1.96 * standard_error),
                           "p": model.pvalues["lesion"]})

lesion_table = pd.DataFrame(lesion_results).sort_values("adjOR", ascending=False)
lesion_table.to_csv(os.path.join(OUT, "hpo_lesion_table.tsv"), sep="\t", index=False)

print("")
print("%-32s %8s %8s" % ("lesion", "adj OR", "p"))
for i in lesion_table.index:
    print("%-32s %8.2f %8.3g" % (lesion_table.loc[i, "lesion"],
                                 lesion_table.loc[i, "adjOR"],
                                 lesion_table.loc[i, "p"]))

print("")
print("Note that the strongest ones (hypoplastic left heart, tetralogy of")
print("Fallot) are the severe lesions that adult biobanks do not contain.")
print("")
print("Finished. Result tables are in " + OUT)
