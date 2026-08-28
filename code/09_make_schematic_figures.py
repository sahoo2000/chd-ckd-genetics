# 09_make_schematic_figures.py
#
# Draws the explanatory diagrams: the graphical abstract, the two
# introduction figures and the two methods figures.
#
# These are drawings rather than plots of data, so they are built out of
# rectangles, arrows and text rather than from a results table.
#
# To run:   python3 09_make_schematic_figures.py

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Ellipse

FIGS = "../figures"

RED = "#A8203A"
TEAL = "#12707F"
GREY = "#8A929B"
LIGHT = "#E4E9ED"
DARK = "#141C25"
AMBER = "#B07A20"

plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["figure.dpi"] = 300


def save(fig, filename):
    if not os.path.exists(FIGS):
        os.makedirs(FIGS)
    path = os.path.join(FIGS, filename)
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print("saved " + path)


def box(ax, x, y, width, height, text, facecolour, edgecolour,
        textcolour="white", fontsize=9, bold=False):
    """Draw a rounded box with text in the middle."""
    patch = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0.012,rounding_size=0.02",
                           facecolor=facecolour, edgecolor=edgecolour, linewidth=1.2)
    ax.add_patch(patch)
    weight = "bold" if bold else "normal"
    ax.text(x + width / 2.0, y + height / 2.0, text, ha="center", va="center",
            fontsize=fontsize, color=textcolour, weight=weight, linespacing=1.45)


def arrow(ax, x1, y1, x2, y2, colour=GREY, width=1.6, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                                 mutation_scale=13, color=colour,
                                 linewidth=width, shrinkA=2, shrinkB=2))


def blank_axes(ax):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


# ======================================================================
# GRAPHICAL ABSTRACT
# ======================================================================

def graphical_abstract():
    fig, ax = plt.subplots(figsize=(10.0, 5.6))
    blank_axes(ax)

    ax.text(0.5, 0.975, "Towards a Mechanistic Explanation of CHD-CKD",
            ha="center", fontsize=15, weight="bold", color=DARK)

    # ---------------- column 1: the question ----------------
    ax.text(0.15, 0.895, "THE ASSUMPTION", ha="center", fontsize=8.5,
            weight="bold", color=GREY)
    box(ax, 0.035, 0.775, 0.105, 0.085, "Congenital\nheart disease", RED, RED, "white", 7.5)
    box(ax, 0.160, 0.775, 0.105, 0.085, "Chronic\nkidney disease", TEAL, TEAL, "white", 7.5)
    arrow(ax, 0.143, 0.8175, 0.157, 0.8175, DARK, 1.4, "<|-|>")
    ax.text(0.15, 0.735, "co-occur clinically", ha="center", fontsize=8, color=DARK)
    ax.text(0.15, 0.700, "assumed to share genes", ha="center", fontsize=8,
            color=DARK, style="italic")

    # ---------------- column 2: what we did ----------------
    ax.text(0.50, 0.895, "THE TEST", ha="center", fontsize=8.5, weight="bold", color=GREY)
    box(ax, 0.325, 0.800, 0.35, 0.058,
        "9,155 curated Mendelian diseases (HPO)", LIGHT, GREY, DARK, 8)
    box(ax, 0.325, 0.712, 0.35, 0.058,
        "adjust for how many organs a syndrome breaks", "#F3E7D2", AMBER, DARK, 8)
    box(ax, 0.325, 0.624, 0.35, 0.058,
        "replicate in 394,841 exomes + 1.4M GWAS", LIGHT, GREY, DARK, 8)
    arrow(ax, 0.50, 0.798, 0.50, 0.772, GREY, 1.3)
    arrow(ax, 0.50, 0.710, 0.50, 0.684, GREY, 1.3)

    # ---------------- column 3: the answer ----------------
    ax.text(0.855, 0.895, "THE ANSWER", ha="center", fontsize=8.5,
            weight="bold", color=GREY)
    box(ax, 0.725, 0.785, 0.26, 0.072,
        "CHD  vs  CKD\nOR 0.95   (p = 0.68)", "#F4DEDB", RED, RED, 8.5, True)
    ax.text(0.855, 0.757, "no shared genetics", ha="center", fontsize=7.5, color=RED)

    box(ax, 0.725, 0.640, 0.26, 0.072,
        "CHD  vs  kidney malformation\nOR 1.77   (p = 3e-13)",
        "#DCEBE1", "#2E6B45", "#2E6B45", 8.0, True)
    ax.text(0.855, 0.598, "genuine shared axis", ha="center", fontsize=7.5, color="#2E6B45")

    arrow(ax, 0.272, 0.818, 0.320, 0.828, GREY, 1.8)
    arrow(ax, 0.680, 0.745, 0.720, 0.800, GREY, 1.8)
    arrow(ax, 0.680, 0.660, 0.720, 0.672, GREY, 1.8)

    # ---------------- bottom band: the mechanism ----------------
    ax.plot([0.035, 0.985], [0.545, 0.545], color=LIGHT, linewidth=1.5)
    ax.text(0.5, 0.500, "THE MECHANISM", ha="center", fontsize=8.5,
            weight="bold", color=GREY)

    # cilium drawing, kept well clear of its caption
    ax.plot([0.095, 0.095], [0.300, 0.420], color=TEAL, linewidth=5,
            solid_capstyle="round")
    ax.add_patch(Ellipse((0.095, 0.283), 0.105, 0.075, facecolor=LIGHT,
                         edgecolor=GREY, linewidth=1.2))
    ax.add_patch(Circle((0.095, 0.278), 0.013, facecolor=GREY, edgecolor="none"))
    ax.text(0.095, 0.205, "primary cilium", ha="center", fontsize=8.5,
            weight="bold", color=DARK)
    ax.text(0.095, 0.170, "built by IFT + dynein-2", ha="center", fontsize=7.5, color=TEAL)

    arrow(ax, 0.155, 0.300, 0.205, 0.335, RED, 1.5)
    arrow(ax, 0.155, 0.270, 0.205, 0.215, TEAL, 1.5)

    box(ax, 0.210, 0.310, 0.185, 0.070, "left-right patterning\nof the heart tube",
        "#F4DEDB", RED, RED, 8)
    box(ax, 0.210, 0.175, 0.185, 0.070, "nephron induction\nin the kidney",
        "#D6EAED", TEAL, TEAL, 8)

    arrow(ax, 0.400, 0.345, 0.455, 0.305, GREY, 1.5)
    arrow(ax, 0.400, 0.210, 0.455, 0.255, GREY, 1.5)

    box(ax, 0.460, 0.235, 0.185, 0.090,
        "One ciliary lesion,\ntwo organs affected", LIGHT, GREY, DARK, 8.5, True)

    arrow(ax, 0.650, 0.280, 0.710, 0.280, GREY, 1.8)

    box(ax, 0.715, 0.310, 0.185, 0.062,
        "IFT / dynein-2\n5.0-5.5x enriched", "#DCEBE1", "#2E6B45", "#2E6B45", 7.5)
    box(ax, 0.715, 0.205, 0.185, 0.062,
        "BBSome\n1.0x (no signal)", LIGHT, GREY, DARK, 7.5)
    ax.text(0.8075, 0.135, "assembly matters,\ncargo sorting does not",
            ha="center", fontsize=7.5, color=DARK, style="italic", linespacing=1.4)

    ax.text(0.5, 0.045,
            "Rare loss-of-function burden, 394,841 UK Biobank exomes   |   "
            "9.7x excess at p < 0.001   |   59/60 in the predicted direction",
            ha="center", fontsize=7.5, color=GREY)

    save(fig, "fig00_graphical_abstract.png")


# ======================================================================
# INTRODUCTION FIGURE 1 - the two competing explanations
# ======================================================================

def intro_figure_hypotheses():
    fig, ax = plt.subplots(figsize=(8.6, 4.0))
    blank_axes(ax)

    ax.text(0.5, 0.95, "Two explanations for the same clinical observation",
            ha="center", fontsize=11.5, weight="bold", color=DARK)

    # ---- left: shared genetic cause ----
    ax.text(0.25, 0.845, "A.  Shared genetic cause", ha="center",
            fontsize=9.5, weight="bold", color=DARK)
    box(ax, 0.175, 0.66, 0.15, 0.085, "one gene\n(e.g. a cilium gene)", LIGHT, GREY, DARK, 8)
    box(ax, 0.06, 0.40, 0.15, 0.085, "heart\nmalformation", RED, RED, "white", 8)
    box(ax, 0.29, 0.40, 0.15, 0.085, "kidney\ndisease", TEAL, TEAL, "white", 8)
    arrow(ax, 0.22, 0.655, 0.14, 0.49, GREY, 1.6)
    arrow(ax, 0.28, 0.655, 0.36, 0.49, GREY, 1.6)
    ax.text(0.25, 0.32, "both organs affected\nfrom conception", ha="center",
            fontsize=8, color=DARK, style="italic")

    # dividing line
    ax.plot([0.5, 0.5], [0.12, 0.86], color=LIGHT, linewidth=1.5)

    # ---- right: acquired ----
    ax.text(0.75, 0.845, "B.  Kidney damage acquired later", ha="center",
            fontsize=9.5, weight="bold", color=DARK)
    box(ax, 0.575, 0.66, 0.15, 0.085, "heart\nmalformation", RED, RED, "white", 8)
    box(ax, 0.775, 0.66, 0.16, 0.085,
        "cyanosis, surgery,\ncontrast, drugs", "#F3E7D2", AMBER, DARK, 7.5)
    box(ax, 0.675, 0.40, 0.15, 0.085, "kidney\ndisease", TEAL, TEAL, "white", 8)
    arrow(ax, 0.725, 0.70, 0.773, 0.70, GREY, 1.6)
    arrow(ax, 0.85, 0.655, 0.78, 0.49, GREY, 1.6)
    ax.text(0.75, 0.32, "kidney healthy at birth,\ndamaged over decades",
            ha="center", fontsize=8, color=DARK, style="italic")

    ax.text(0.5, 0.135,
            "These make different genetic predictions. Under A the two conditions share "
            "genetic architecture;\nunder B they do not. This project tests which holds.",
            ha="center", fontsize=8.5, color=DARK)

    save(fig, "figI1_two_hypotheses.png")


# ======================================================================
# INTRODUCTION FIGURE 2 - why naive gene overlap is misleading
# ======================================================================

def intro_figure_pleiotropy():
    fig, ax = plt.subplots(figsize=(8.6, 3.9))
    blank_axes(ax)

    ax.text(0.5, 0.95, "Why counting shared genes is not enough",
            ha="center", fontsize=11.5, weight="bold", color=DARK)

    box(ax, 0.03, 0.55, 0.16, 0.10, "a developmental\ngene breaks",
        LIGHT, GREY, DARK, 8.5)

    organs = ["heart", "kidney", "thumb", "palate", "brain", "ear"]
    colours = [RED, TEAL, GREY, GREY, GREY, GREY]
    y_top = 0.80
    for i in range(len(organs)):
        y = y_top - i * 0.115
        box(ax, 0.30, y, 0.115, 0.075, organs[i], colours[i], colours[i], "white", 8)
        arrow(ax, 0.195, 0.60, 0.297, y + 0.037, GREY, 1.0)

    ax.text(0.3575, 0.055, "one gene, many organs", ha="center",
            fontsize=8, color=DARK, style="italic")

    ax.plot([0.48, 0.48], [0.08, 0.88], color=LIGHT, linewidth=1.5)

    ax.text(0.745, 0.86, "Consequence", ha="center", fontsize=9.5,
            weight="bold", color=DARK)
    ax.text(0.745, 0.70,
            "Every pair of organs now shares this gene.\n"
            "Heart and kidney share it. So do heart and thumb.\n"
            "An overlap test cannot tell these apart.",
            ha="center", fontsize=8.5, color=DARK, linespacing=1.6)

    box(ax, 0.565, 0.40, 0.36, 0.13,
        "Observed: CHD syndromes affect a median of\n"
        "11 organ systems; other diseases affect 5",
        "#F3E7D2", AMBER, DARK, 8.5)

    box(ax, 0.565, 0.16, 0.36, 0.155,
        "Solution used here:\n"
        "model the kidney phenotype while adjusting for\n"
        "the number of organ systems the disease affects",
        "#DCEBE1", "#2E6B45", DARK, 8.5)
    arrow(ax, 0.745, 0.395, 0.745, 0.322, GREY, 1.6)

    save(fig, "figI2_pleiotropy_problem.png")


# ======================================================================
# METHODS FIGURE 1 - the analysis workflow
# ======================================================================

def methods_figure_workflow():
    fig, ax = plt.subplots(figsize=(9.0, 6.4))
    blank_axes(ax)

    ax.text(0.5, 0.975, "Analysis workflow", ha="center", fontsize=12,
            weight="bold", color=DARK)

    # ---- stage 1: data sources ----
    ax.text(0.5, 0.925, "STAGE 1  -  data sources (all public)", ha="center",
            fontsize=8.5, weight="bold", color=GREY)
    sources = [("Human Phenotype\nOntology 2026-06-23", 0.03),
               ("HGNC + MANE v1.5\n+ gnomAD v4.1", 0.265),
               ("Genebass\n394,841 exomes", 0.50),
               ("CKDGen + FinnGen R10\nGWAS summary stats", 0.735)]
    for label, x in sources:
        box(ax, x, 0.815, 0.225, 0.075, label, LIGHT, GREY, DARK, 7.5)

    # ---- stage 2: gene list ----
    ax.text(0.5, 0.765, "STAGE 2  -  build the candidate gene list", ha="center",
            fontsize=8.5, weight="bold", color=GREY)
    box(ax, 0.19, 0.665, 0.62, 0.068,
        "genes causing BOTH a heart and a kidney malformation   (n = 750)",
        "#D6EAED", TEAL, DARK, 8.5)
    arrow(ax, 0.14, 0.812, 0.30, 0.736, GREY, 1.3)
    arrow(ax, 0.378, 0.812, 0.42, 0.736, GREY, 1.3)

    filters = [("remove non-coding,\nnon-autosomal", "676"),
               ("remove infant-lethal\n(cannot reach a biobank)", "611"),
               ("keep Tier 1, drop\nimmune genes", "169"),
               ("ciliary subset", "37")]
    x = 0.045
    for label, remaining in filters:
        box(ax, x, 0.545, 0.20, 0.075, label, "white", GREY, DARK, 7.5)
        ax.text(x + 0.10, 0.522, "n = " + remaining, ha="center",
                fontsize=8, weight="bold", color=TEAL)
        x = x + 0.235
    for i in range(3):
        arrow(ax, 0.245 + i * 0.235, 0.583, 0.278 + i * 0.235, 0.583, GREY, 1.3)
    arrow(ax, 0.5, 0.662, 0.5, 0.625, GREY, 1.6)

    # ---- stage 3: the four analyses ----
    ax.text(0.5, 0.472, "STAGE 3  -  four independent analyses", ha="center",
            fontsize=8.5, weight="bold", color=GREY)
    analyses = [("Disease-level\nlogistic regression\n(adjusts for pleiotropy)",
                 0.03, "#F3E7D2", AMBER),
                ("Rare-variant burden\nlookup in Genebass\n(pLoF and missense)",
                 0.265, "#D6EAED", TEAL),
                ("Cross-trait LD score\nregression\n(CKDGen, FinnGen)",
                 0.50, "#D6EAED", TEAL),
                ("MAGMA common-variant\ngene-set analysis\n(CKDGen eGFR)",
                 0.735, "#D6EAED", TEAL)]
    for label, x, fill, edge in analyses:
        box(ax, x, 0.335, 0.225, 0.105, label, fill, edge, DARK, 7.5)
        arrow(ax, x + 0.1125, 0.51, x + 0.1125, 0.445, GREY, 1.2)

    # ---- stage 4: results ----
    ax.text(0.5, 0.285, "STAGE 4  -  what each analysis concluded", ha="center",
            fontsize=8.5, weight="bold", color=GREY)
    outcomes = [("CKD link vanishes\nOR 0.95, p = 0.68", 0.03, "#F4DEDB", RED),
                ("ciliary genes enriched\n9.7x at p < 0.001", 0.265, "#DCEBE1", "#2E6B45"),
                ("no correlation\nrg = -0.03, p = 0.67", 0.50, "#F4DEDB", RED),
                ("no common-variant\nenrichment, p = 0.48", 0.735, "#F4DEDB", RED)]
    for label, x, fill, edge in outcomes:
        box(ax, x, 0.155, 0.225, 0.085, label, fill, edge, DARK, 7.5)
        arrow(ax, x + 0.1125, 0.332, x + 0.1125, 0.245, GREY, 1.2)

    box(ax, 0.14, 0.025, 0.72, 0.075,
        "CONCLUSION:  CHD and CKD do not share Mendelian genetic architecture.\n"
        "CHD and congenital kidney malformation do, through ciliary biology.",
        DARK, DARK, "white", 9, True)

    save(fig, "figM1_workflow.png")


# ======================================================================
# METHODS FIGURE 2 - the gene filtering funnel, with reasons
# ======================================================================

def methods_figure_funnel():
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    blank_axes(ax)

    ax.text(0.5, 0.96, "Candidate gene list: what was removed and why",
            ha="center", fontsize=11.5, weight="bold", color=DARK)

    steps = [("Genes causing a heart malformation", 1660, ""),
             ("...and also a kidney malformation", 750, "910 removed"),
             ("...protein-coding and autosomal", 676, "74 non-coding / X / MT"),
             ("...survivable to adulthood", 611, "65 infant-lethal only"),
             ("...Tier 1, non-immune", 169, "442 no CKD link or immune"),
             ("...ciliary mechanism", 37, "132 other mechanisms"),
             ("...excluding known cystogenic genes", 31, "6 already established")]

    widest = 0.74
    y = 0.845
    for i in range(len(steps)):
        label, count, removed = steps[i]
        # bar width proportional to log of the count, so small sets stay visible
        width = widest * (np.log10(count) / np.log10(1660))
        x = 0.5 - width / 2.0

        if i < 2:
            colour = LIGHT
            textcolour = DARK
        elif i < 5:
            colour = "#9FC4CC"
            textcolour = DARK
        else:
            colour = TEAL
            textcolour = "white"

        box(ax, x, y, width, 0.072, "", colour, GREY)
        ax.text(0.5, y + 0.036, "%s   (n = %s)" % (label, "{:,}".format(count)),
                ha="center", va="center", fontsize=8.5, color=textcolour)
        if removed != "":
            ax.text(0.955, y + 0.036, removed, ha="right", va="center",
                    fontsize=7.5, color=RED)
        if i < len(steps) - 1:
            arrow(ax, 0.5, y, 0.5, y - 0.043, GREY, 1.3)
        y = y - 0.115

    ax.text(0.5, 0.045,
            "The 31-gene set at the bottom is the pre-specified primary analysis.\n"
            "The 6 removed cystogenic genes are kept separately as a positive control.",
            ha="center", fontsize=8, color=DARK, style="italic")

    save(fig, "figM2_gene_filtering.png")


# ======================================================================

if __name__ == "__main__":
    print("Drawing schematic figures...")
    graphical_abstract()
    intro_figure_hypotheses()
    intro_figure_pleiotropy()
    methods_figure_workflow()
    methods_figure_funnel()
    print("Done.")
