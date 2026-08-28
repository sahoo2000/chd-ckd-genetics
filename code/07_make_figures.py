# 07_make_figures.py
#
# This script draws all the figures for the thesis.
# Run it after the analysis scripts have made their result files.
#
# To run:   python3 07_make_figures.py
#
# Every figure is saved into the "figures" folder as a PNG file at 300 dpi,
# which is good enough for printing in a thesis.

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import chi2, binomtest

import palette
from palette import INK, MUTED, HAIRLINE

# ----------------------------------------------------------------------
# Settings that are the same for every figure
# ----------------------------------------------------------------------

DATA = "../data/results"
FIGS = "../figures"

# Pastel scheme shared with the schematic figures (see palette.py).
# The darker "line" colour is used for points and bar outlines so that
# the marks stay readable; the pale "fill" colour is used for bar bodies.
RED = palette.ROSE[1]          # heart / negative result
RED_FILL = palette.ROSE[0]
TEAL = palette.TEAL[1]         # kidney
TEAL_FILL = palette.TEAL[0]
GREEN = palette.SAGE[1]        # confirmed / positive
GREEN_FILL = palette.SAGE[0]
GREY = palette.STONE[1]        # neutral
LIGHT = palette.STONE[0]

palette.apply_style()


def save(fig, filename):
    """Save a figure into the figures folder and print a message."""
    if not os.path.exists(FIGS):
        os.makedirs(FIGS)
    path = os.path.join(FIGS, filename)
    fig.savefig(path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print("saved " + path)


# ----------------------------------------------------------------------
# Figure 1 - odds ratios before and after adjusting for pleiotropy
# ----------------------------------------------------------------------

def figure1_forest():
    # These four numbers come from 02_hpo_analysis.py
    labels = ["CHD to kidney anomaly\n(unadjusted)",
              "CHD to kidney anomaly\n(adjusted)",
              "CHD to chronic kidney disease\n(unadjusted)",
              "CHD to chronic kidney disease\n(adjusted)"]
    odds = [4.21, 1.77, 1.57, 0.95]
    low = [3.69, 1.52, 1.25, 0.73]
    high = [4.81, 2.07, 1.97, 1.23]
    colours = [RED, RED, TEAL, TEAL]

    fig, ax = plt.subplots(figsize=(7, 3.2))
    y_positions = [3, 2, 1, 0]

    for i in range(len(odds)):
        y = y_positions[i]
        # the horizontal line is the 95% confidence interval
        ax.plot([low[i], high[i]], [y, y], color=colours[i], linewidth=2)
        # small vertical caps at each end
        ax.plot([low[i], low[i]], [y - 0.12, y + 0.12], color=colours[i], linewidth=2)
        ax.plot([high[i], high[i]], [y - 0.12, y + 0.12], color=colours[i], linewidth=2)
        # the dot is the point estimate; hollow if it crosses 1 (not significant)
        if low[i] < 1 < high[i]:
            ax.plot(odds[i], y, "o", color="white", markeredgecolor=colours[i],
                    markersize=8, markeredgewidth=2)
        else:
            ax.plot(odds[i], y, "o", color=colours[i], markersize=8)
        ax.text(6.0, y, str(odds[i]), va="center", fontsize=9)

    ax.axvline(1, color="black", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_xticks([0.6, 1, 2, 3, 4, 5])
    ax.set_xticklabels(["0.6", "1", "2", "3", "4", "5"])
    ax.set_xlim(0.55, 6.6)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.set_ylim(-0.6, 3.6)
    ax.set_xlabel("Odds ratio (95% confidence interval)")
    ax.set_title("Adjusting for how many organs a syndrome affects\n"
                 "removes the CHD to CKD association completely", loc="left")
    save(fig, "fig01_forest_adjusted_odds_ratios.png")


# ----------------------------------------------------------------------
# Figure 2 - how special is the kidney compared with other organs?
# ----------------------------------------------------------------------

def figure2_organ_specificity():
    table = pd.read_csv(os.path.join(DATA, "hpo_organ_specificity.csv"))
    table = table.sort_values("OR", ascending=False).reset_index(drop=True)
    top = table.head(25).copy()

    # find where the kidney sits
    kidney_rank = -1
    for i in range(len(table)):
        if table.loc[i, "id"] == "HP:0012210":
            kidney_rank = i
            break
    median_or = table["OR"].median()

    fig, ax = plt.subplots(figsize=(7, 6))
    colours = []
    for i in range(len(top)):
        if top.loc[i, "id"] == "HP:0012210":
            colours.append(TEAL_FILL)
        else:
            colours.append(LIGHT)

    positions = np.arange(len(top))[::-1]
    ax.barh(positions, top["OR"], height=0.68, color=colours,
            edgecolor=GREY, linewidth=0.7)
    ax.grid(axis="y", visible=False)
    names = []
    for i in range(len(top)):
        names.append(top.loc[i, "term"].replace("Abnormal ", "").replace(" morphology", ""))
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=8)
    ax.axvline(median_or, color=RED, linestyle="--", linewidth=1)
    ax.text(median_or + 0.15, 1, "median across\nall 194 organs\n(OR = %.2f)" % median_or,
            color=RED, fontsize=8, va="bottom")
    ax.set_xlabel("Odds ratio for overlap with congenital heart disease genes")
    ax.set_title("The kidney ranks %d of %d organ systems\n"
                 "Many unrelated organs overlap the heart just as much"
                 % (kidney_rank + 1, len(table)), loc="left")
    save(fig, "fig02_organ_specificity.png")


# ----------------------------------------------------------------------
# Figure 3 - which heart lesions go with kidney anomalies?
# ----------------------------------------------------------------------

def figure3_lesions():
    table = pd.read_csv(os.path.join(DATA, "hpo_lesion_table.tsv"), sep="\t")
    table = table.sort_values("adjOR", ascending=True).reset_index(drop=True)

    # Which lesions actually appear in an adult cohort like UK Biobank?
    # Severe lesions are missing because those patients rarely reach age 40-69.
    rare_in_adults = ["Hypoplastic left heart", "Tetralogy of Fallot",
                      "Atrioventricular canal defect", "Transposition great arteries"]

    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    for i in range(len(table)):
        if table.loc[i, "lesion"] in rare_in_adults:
            colour = RED
        else:
            colour = TEAL
        ax.plot([table.loc[i, "lo"], table.loc[i, "hi"]], [i, i], color=colour, linewidth=2)
        ax.plot(table.loc[i, "adjOR"], i, "o", color=colour, markersize=7)

    ax.axvline(1, color="black", linestyle="--", linewidth=0.8)
    ax.set_xscale("log")
    ax.set_yticks(range(len(table)))
    ax.set_yticklabels(table["lesion"], fontsize=8.5)
    ax.set_xlabel("Adjusted odds ratio for a structural kidney anomaly")
    ax.set_title("The lesions most linked to kidney anomalies are the ones\n"
                 "an adult cohort does not contain", loc="left")
    # simple legend
    ax.plot([], [], "o-", color=RED, label="rare in adults (poor survival)")
    ax.plot([], [], "o-", color=TEAL, label="present in adult cohorts")
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    save(fig, "fig03_lesion_odds_ratios.png")


# ----------------------------------------------------------------------
# Figure 4 - Genebass results for the ciliary genes (volcano-style)
# ----------------------------------------------------------------------

def figure4_genebass_volcano():
    table = pd.read_csv(os.path.join(DATA, "ciliary_renal.tsv"), sep="\t")
    table["P"] = pd.to_numeric(table["P_SKATO"], errors="coerce")
    table["BETA"] = pd.to_numeric(table["BETA"], errors="coerce")
    table = table.dropna(subset=["P", "BETA"])
    table["logp"] = -np.log10(table["P"])

    cystogenic = ["PKD1", "PKD2", "IFT140", "ALG9", "GANAB", "DNAJB11"]

    fig, ax = plt.subplots(figsize=(7, 4.6))
    for i in table.index:
        if table.loc[i, "gene"] in cystogenic:
            colour = RED
            size = 34
        else:
            colour = TEAL
            size = 18
        ax.scatter(table.loc[i, "BETA"], table.loc[i, "logp"],
                   color=colour, s=size, alpha=0.75, edgecolors="none")

    ax.axhline(-np.log10(2.5e-6), color="black", linestyle="--", linewidth=0.9)
    ax.text(0.42, -np.log10(2.5e-6) + 2, "exome-wide significance (p = 2.5e-6)", fontsize=8)
    ax.axvline(0, color=GREY, linewidth=0.6)

    # label the strongest few points
    strongest = table.nlargest(5, "logp")
    for i in strongest.index:
        ax.annotate(table.loc[i, "gene"],
                    (table.loc[i, "BETA"], table.loc[i, "logp"]),
                    textcoords="offset points", xytext=(6, 3), fontsize=8)

    ax.set_xlabel("Effect of rare loss-of-function burden (positive = worse kidney measure)")
    ax.set_ylabel("-log10 p-value")
    ax.set_title("Only the dominant cystogenic genes reach significance\n"
                 "in 394,841 UK Biobank exomes", loc="left")
    ax.plot([], [], "o", color=RED, label="known cystogenic genes")
    ax.plot([], [], "o", color=TEAL, label="other ciliary genes")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    save(fig, "fig04_genebass_volcano.png")


# ----------------------------------------------------------------------
# Figure 5 - QQ plot comparing ciliary and non-ciliary genes
# ----------------------------------------------------------------------

def qq_points(pvalues):
    """Return the expected and observed -log10 p-values for a QQ plot."""
    pvalues = np.sort(np.asarray(pvalues))
    n = len(pvalues)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(pvalues)
    return expected, observed


def figure5_qq():
    tier1 = pd.read_csv(os.path.join(DATA, "tier1_renal_pLoF.tsv"), sep="\t")
    tier1["P"] = pd.to_numeric(tier1["P_SKATO"], errors="coerce")
    tier1 = tier1.dropna(subset=["P"])

    ciliary_genes = []
    for line in open("../data/genes/ciliary_symbols.txt"):
        if line.strip() != "":
            ciliary_genes.append(line.strip())

    # The known cystogenic genes give p-values around 1e-95, which squashes
    # everything else. So we draw two panels: one with them, one without.
    cystogenic = ["PKD1", "PKD2", "IFT140", "ALG9", "GANAB", "DNAJB11"]

    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.6))

    for panel in [0, 1]:
        ax = axes[panel]
        if panel == 0:
            data = tier1
            title = "All genes"
        else:
            data = tier1[~tier1["gene"].isin(cystogenic)]
            title = "Known cystogenic genes removed"

        is_ciliary = data["gene"].isin(ciliary_genes)
        ciliary_p = data.loc[is_ciliary, "P"].values
        other_p = data.loc[~is_ciliary, "P"].values

        n_cil = data.loc[is_ciliary, "gene"].nunique()
        n_oth = data.loc[~is_ciliary, "gene"].nunique()

        for pvals, colour, label in [
                (other_p, GREY, "non-ciliary (n=%d genes)" % n_oth),
                (ciliary_p, TEAL, "ciliary (n=%d genes)" % n_cil)]:
            expected, observed = qq_points(pvals)
            ax.scatter(expected, observed, s=10, color=colour, alpha=0.7,
                       edgecolors="none", label=label)

        biggest = max(qq_points(ciliary_p)[0].max(), qq_points(other_p)[0].max())
        ax.plot([0, biggest], [0, biggest], color="black", linestyle="--", linewidth=0.8)
        ax.set_xlabel("Expected -log10 p")
        if panel == 0:
            ax.set_ylabel("Observed -log10 p")
        ax.set_title(title, loc="left", fontsize=9.5)
        ax.legend(fontsize=8, frameon=False, loc="upper left")

    fig.suptitle("Ciliary genes depart from the null even after the known "
                 "cystogenic genes are removed", x=0.02, ha="left", fontsize=10.5)
    fig.tight_layout()
    save(fig, "fig05_qq_ciliary_vs_other.png")


# ----------------------------------------------------------------------
# Figure 6 - enrichment by ciliary subcomplex
# ----------------------------------------------------------------------

def figure6_subcomplex():
    # These come from 05_genebass_analysis.py
    families = ["ADPKD / cystogenic", "IFT-A", "IFT-B", "Dynein-2",
                "NPHP module", "Transition zone / MKS",
                "Chaperonin (BBS6/10/12)", "BBSome"]
    enrichment = [6.67, 5.5, 5.0, 5.0, 2.5, 1.88, 1.67, 1.0]

    fig, ax = plt.subplots(figsize=(7.2, 3.6))

    # Darker teal for the strongly enriched groups, pale stone for the
    # ones sitting at or near no enrichment.
    colours = []
    edges = []
    for value in enrichment:
        if value >= 4:
            colours.append(TEAL_FILL)
            edges.append(TEAL)
        elif value >= 2:
            colours.append(GREEN_FILL)
            edges.append(GREEN)
        else:
            colours.append(LIGHT)
            edges.append(GREY)

    positions = np.arange(len(families))[::-1]
    bars = ax.barh(positions, enrichment, height=0.62,
                   color=colours, linewidth=0.9)
    for i in range(len(bars)):
        bars[i].set_edgecolor(edges[i])

    # write the value at the end of each bar
    for i in range(len(enrichment)):
        ax.text(enrichment[i] + 0.12, positions[i], "%.2fx" % enrichment[i],
                va="center", fontsize=8, color=INK)

    ax.axvline(1, color=INK, linestyle="--", linewidth=0.9)
    # put the label above the plot area so it cannot sit on top of a bar
    ax.text(1.15, len(families) - 0.35, "no enrichment", ha="left",
            fontsize=8, color=MUTED)

    ax.set_yticks(positions)
    ax.set_yticklabels(families, fontsize=8.5)
    ax.set_xlim(0, 7.6)
    ax.set_ylim(-0.7, len(families) - 0.1)
    ax.set_xlabel("Nominal hits divided by hits expected by chance")
    ax.set_title("The signal follows the machinery that builds the cilium.\n"
                 "The BBSome, which only sorts cargo, shows nothing.", loc="left")
    ax.grid(axis="y", visible=False)
    save(fig, "fig06_subcomplex_enrichment.png")


# ----------------------------------------------------------------------
# Figure 7 - direction of effect, loss-of-function versus missense
# ----------------------------------------------------------------------

def figure7_direction():
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.4))

    for ax, filename, title in [
            (axes[0], "ciliary_renal.tsv", "Loss-of-function variants"),
            (axes[1], "ciliary_renal_missense.tsv", "Missense variants (control)")]:

        table = pd.read_csv(os.path.join(DATA, filename), sep="\t")
        table["P"] = pd.to_numeric(table["P_SKATO"], errors="coerce")
        table["BETA"] = pd.to_numeric(table["BETA"], errors="coerce")
        table = table.dropna(subset=["P", "BETA"])
        hits = table[table["P"] < 0.05]

        positive = int((hits["BETA"] > 0).sum())
        negative = int((hits["BETA"] <= 0).sum())
        pvalue = binomtest(positive, positive + negative, 0.5, alternative="greater").pvalue

        ax.bar([0, 1], [positive, negative], color=[TEAL_FILL, LIGHT],
               edgecolor=GREY, linewidth=0.8, width=0.55)
        ax.grid(axis="x", visible=False)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["worse kidney\nfunction", "better kidney\nfunction"], fontsize=8.5)
        ax.set_ylabel("Number of nominally significant results")
        ax.set_title("%s\n%d of %d in expected direction (p = %.2g)"
                     % (title, positive, positive + negative, pvalue),
                     loc="left", fontsize=9)
        ax.text(0, positive + 1.5, str(positive), ha="center", fontsize=9)
        ax.text(1, negative + 1.5, str(negative), ha="center", fontsize=9)
        ax.set_ylim(0, max(positive, negative) * 1.25 + 4)

    fig.suptitle("Loss-of-function results point one way; unfiltered missense results do not",
                 x=0.02, ha="left", fontsize=10)
    fig.tight_layout()
    save(fig, "fig07_direction_of_effect.png")


# ----------------------------------------------------------------------
# Figure 8 - genetic correlations from LD score regression
# ----------------------------------------------------------------------

def figure8_ldsc():
    table = pd.read_csv(os.path.join(DATA, "ldsc_results.tsv"), sep="\t")
    correlations = table[table["analysis"] == "rg"].copy()
    correlations["label"] = correlations["trait1"] + "\nvs " + correlations["trait2"]
    correlations = correlations.sort_values("estimate").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(7, 3.8))
    for i in range(len(correlations)):
        estimate = correlations.loc[i, "estimate"]
        se = correlations.loc[i, "se"]
        low = estimate - 1.96 * se
        high = estimate + 1.96 * se
        # colour by whether the interval excludes zero
        if low > 0 or high < 0:
            colour = RED
        else:
            colour = GREY
        ax.plot([low, high], [i, i], color=colour, linewidth=2)
        ax.plot(estimate, i, "o", color=colour, markersize=7)

    ax.axvline(0, color="black", linestyle="--", linewidth=0.8)
    ax.set_yticks(range(len(correlations)))
    ax.set_yticklabels(correlations["label"], fontsize=7.5)
    ax.set_xlabel("Genetic correlation (rg) with 95% confidence interval")
    ax.set_title("Kidney function correlates with chronic kidney disease,\n"
                 "but not with congenital heart disease", loc="left")
    save(fig, "fig08_ldsc_genetic_correlation.png")


# ----------------------------------------------------------------------
# Figure 9 - why the comorbidity phenotype cannot be studied in UK Biobank
# ----------------------------------------------------------------------

def figure9_sample_sizes():
    names = ["Kidney biomarkers\n(creatinine, cystatin C)",
             "Chronic kidney disease\n(N18)",
             "Congenital heart disease\n(Q20-Q26)",
             "Kidney malformation\n(Q60-Q64)",
             "CHD and kidney\nmalformation together"]
    counts = [376624, 14086, 3385, 1386, 18]

    fig, ax = plt.subplots(figsize=(7, 3.6))
    colours = [TEAL_FILL, TEAL_FILL, RED_FILL, RED_FILL, palette.ROSE[1]]
    positions = np.arange(len(names))[::-1]
    ax.barh(positions, counts, height=0.62, color=colours,
            edgecolor=GREY, linewidth=0.7)
    ax.grid(axis="y", visible=False)
    ax.set_xscale("log")
    ax.set_yticks(positions)
    ax.set_yticklabels(names, fontsize=8.5)
    ax.set_xlabel("Number of people in UK Biobank (log scale)")
    for i in range(len(counts)):
        ax.text(counts[i] * 1.25, positions[i], "{:,}".format(counts[i]),
                va="center", fontsize=8.5)
    ax.set_xlim(5, 2500000)
    ax.set_title("The comorbidity phenotype has about 18 people.\n"
                 "A continuous kidney measurement has 376,624.", loc="left")
    save(fig, "fig09_sample_sizes.png")


# ----------------------------------------------------------------------
# Figure 10 - statistical power of the different study designs
# ----------------------------------------------------------------------

def power_for_quantitative_trait(n_carriers, effect_size, alpha, total_n=450000):
    """Power to detect an effect of a given size in a continuous trait."""
    from scipy.stats import ncx2
    non_centrality = effect_size ** 2 * (n_carriers * (total_n - n_carriers) / total_n)
    threshold = chi2.isf(alpha, 1)
    return 1 - ncx2.cdf(threshold, 1, non_centrality)


def figure10_power():
    effect_sizes = np.linspace(0.01, 0.5, 120)

    fig, ax = plt.subplots(figsize=(6.6, 4))
    settings = [(3000, 0.05, "31 genes pooled into one test", RED),
                (3000, 0.05 / 31, "31 genes tested separately", TEAL),
                (450, 0.05 / 31, "one gene on its own", GREY)]

    for n_carriers, alpha, label, colour in settings:
        powers = []
        for effect in effect_sizes:
            powers.append(power_for_quantitative_trait(n_carriers, effect, alpha))
        ax.plot(effect_sizes, powers, color=colour, linewidth=2, label=label)

    ax.axhline(0.8, color="black", linestyle="--", linewidth=0.8)
    ax.text(0.42, 0.82, "80% power", fontsize=8)
    ax.set_xlabel("True effect size (standard deviations of eGFR)")
    ax.set_ylabel("Power")
    ax.set_ylim(0, 1.02)
    ax.set_title("Pooling genes into a single test is what makes the study possible",
                 loc="left")
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    save(fig, "fig10_power_curves.png")


# ----------------------------------------------------------------------
# Run everything
# ----------------------------------------------------------------------

if __name__ == "__main__":
    print("Making figures...")
    figure1_forest()
    figure2_organ_specificity()
    figure3_lesions()
    figure4_genebass_volcano()
    figure5_qq()
    figure6_subcomplex()
    figure7_direction()
    figure8_ldsc()
    figure9_sample_sizes()
    figure10_power()
    print("Done.")


# ----------------------------------------------------------------------
# Figure 11 - rare variants versus common variants
# ----------------------------------------------------------------------

def figure11_rare_vs_common():
    """
    Compare what we learn from rare variants (Genebass) with what we learn
    from common variants (MAGMA on the CKDGen GWAS), for the same genes.
    """
    burden = pd.read_csv(os.path.join(DATA, "ciliary_renal.tsv"), sep="\t")
    burden["P"] = pd.to_numeric(burden["P_SKATO"], errors="coerce")
    # best rare-variant p-value per gene
    best_rare = burden.groupby("gene")["P"].min()

    common = pd.read_csv(os.path.join(DATA, "magma_genes_egfr.tsv"), sep="\t")
    common = common.dropna(subset=["SYMBOL"])
    best_common = common.set_index("SYMBOL")["P"]

    shared_genes = []
    for gene in best_rare.index:
        if gene in best_common.index:
            shared_genes.append(gene)

    rare_values = -np.log10(best_rare[shared_genes].values)
    common_values = -np.log10(best_common[shared_genes].values)

    cystogenic = ["PKD1", "PKD2", "IFT140", "ALG9", "GANAB", "DNAJB11"]

    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    for i in range(len(shared_genes)):
        gene = shared_genes[i]
        if gene in cystogenic:
            colour = RED
        else:
            colour = TEAL
        ax.scatter(common_values[i], rare_values[i], color=colour, s=40,
                   alpha=0.8, edgecolors="none")
        # label the interesting ones
        if rare_values[i] > 5 or common_values[i] > 7:
            ax.annotate(gene, (common_values[i], rare_values[i]),
                        textcoords="offset points", xytext=(6, 3), fontsize=8)

    genome_wide = -np.log10(0.05 / 18030)
    ax.axvline(genome_wide, color=GREY, linestyle="--", linewidth=0.9)
    ax.axhline(-np.log10(2.5e-6), color=GREY, linestyle="--", linewidth=0.9)
    ax.text(genome_wide + 0.3, 1, "genome-wide\n(common)", fontsize=7.5, color=GREY)
    ax.text(0.4, -np.log10(2.5e-6) + 2, "exome-wide (rare)", fontsize=7.5, color=GREY)

    ax.set_xlabel("Common variant evidence: -log10 p (MAGMA, CKDGen eGFR)")
    ax.set_ylabel("Rare variant evidence: -log10 p (Genebass burden)")
    ax.set_title("Rare and common variant evidence pick out\ndifferent ciliary genes",
                 loc="left")
    ax.plot([], [], "o", color=RED, label="known cystogenic")
    ax.plot([], [], "o", color=TEAL, label="other ciliary")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    save(fig, "fig11_rare_vs_common.png")
