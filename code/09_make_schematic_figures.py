# 09_make_schematic_figures.py
#
# Draws the diagrams: the graphical abstract, two introduction figures
# and two methods figures.
#
# These are drawings rather than plots of data, so they are made of
# boxes, arrows and text instead of coming from a results table.
#
# Layout rule used throughout: boxes sit on a grid worked out by
# palette.row_positions(), so there is always a real gap between them.
# Each row of boxes gets its own horizontal band, and the arrows run in
# the empty space between bands, so labels never collide with arrows.
#
# To run:   python3 09_make_schematic_figures.py

import os
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse

import palette
from palette import (box, arrow, blank_axes, row_positions, fitted_text,
                     INK, MUTED, HAIRLINE,
                     ROSE, TEAL, SAGE, SAND, LILAC, STONE,
                     ROSE_LINE, TEAL_LINE, SAGE_LINE)

palette.apply_style()
FIGS = "../figures"


def save(fig, filename):
    if not os.path.exists(FIGS):
        os.makedirs(FIGS)
    path = os.path.join(FIGS, filename)
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print("saved " + path)


def band_label(ax, x, y, text):
    """A small left-aligned label sitting above a row of boxes."""
    ax.text(x, y, text, ha="left", va="center", fontsize=7.5,
            weight="bold", color=MUTED)


# ======================================================================
# GRAPHICAL ABSTRACT
# ======================================================================

def graphical_abstract():
    fig, ax = plt.subplots(figsize=(10.0, 5.3))
    blank_axes(ax)

    ax.text(0.5, 0.962, "Towards a Mechanistic Explanation of CHD-CKD",
            ha="center", fontsize=15, weight="bold", color=INK)

    columns = row_positions(3, left=0.03, right=0.97, gap=0.050)

    headings = ["THE ASSUMPTION", "WHAT WE DID", "WHAT WE FOUND"]
    for i in range(3):
        x, w = columns[i]
        ax.text(x + w / 2.0, 0.882, headings[i], ha="center", fontsize=8,
                weight="bold", color=MUTED)

    # ---- column 1 ----
    x, w = columns[0]
    box(ax, x, 0.772, w * 0.46, 0.072, "Congenital\nheart disease", ROSE,
        7.5, textcolour=ROSE_LINE)
    box(ax, x + w * 0.54, 0.772, w * 0.46, 0.072, "Chronic\nkidney disease", TEAL,
        7.5, textcolour=TEAL_LINE)
    ax.text(x + w / 2.0, 0.733, "co-occur clinically", ha="center",
            fontsize=8, color=INK)
    ax.text(x + w / 2.0, 0.701, "assumed to share genes", ha="center",
            fontsize=8, color=MUTED, style="italic")

    # ---- column 2 ----
    x, w = columns[1]
    steps = ["9,155 curated Mendelian diseases",
             "adjust for how many organs\na syndrome affects",
             "replicate in 394,841 exomes\nand 1.4M GWAS samples"]
    fills = [STONE, SAND, STONE]
    y = 0.786
    for i in range(3):
        box(ax, x, y, w, 0.056, steps[i], fills[i], 7.5)
        if i < 2:
            arrow(ax, x + w / 2.0, y - 0.003, x + w / 2.0, y - 0.025)
        y = y - 0.082

    # ---- column 3 ----
    x, w = columns[2]
    box(ax, x, 0.780, w, 0.066,
        "CHD and CKD\nOR 0.95   (p = 0.68)", ROSE, 8.5, bold=True,
        textcolour=ROSE_LINE)
    ax.text(x + w / 2.0, 0.757, "no shared genetics", ha="center",
            fontsize=7.5, color=ROSE_LINE)
    box(ax, x, 0.652, w, 0.066,
        "CHD and kidney malformation\nOR 1.77   (p = 3e-13)", SAGE, 8.0, bold=True,
        textcolour=SAGE_LINE)
    ax.text(x + w / 2.0, 0.629, "a genuine shared axis", ha="center",
            fontsize=7.5, color=SAGE_LINE)

    arrow(ax, columns[0][0] + columns[0][1] + 0.008, 0.808,
          columns[1][0] - 0.008, 0.814, width=1.6)
    arrow(ax, columns[1][0] + columns[1][1] + 0.008, 0.700,
          columns[2][0] - 0.008, 0.790, width=1.6)
    arrow(ax, columns[1][0] + columns[1][1] + 0.008, 0.678,
          columns[2][0] - 0.008, 0.690, width=1.6)

    # ---- mechanism band ----
    ax.plot([0.03, 0.97], [0.560, 0.560], color=HAIRLINE, linewidth=1.2)
    ax.text(0.5, 0.516, "THE MECHANISM", ha="center", fontsize=8,
            weight="bold", color=MUTED)

    ax.plot([0.082, 0.082], [0.320, 0.432], color=TEAL_LINE, linewidth=5,
            solid_capstyle="round")
    ax.add_patch(Ellipse((0.082, 0.302), 0.092, 0.060,
                         facecolor=TEAL[0], edgecolor=TEAL_LINE, linewidth=1.1))
    ax.add_patch(Circle((0.082, 0.298), 0.010, facecolor=TEAL_LINE, edgecolor="none"))
    ax.text(0.082, 0.232, "primary cilium", ha="center", fontsize=8.5,
            weight="bold", color=INK)
    ax.text(0.082, 0.200, "built by IFT and dynein-2", ha="center",
            fontsize=7.5, color=TEAL_LINE)

    arrow(ax, 0.136, 0.318, 0.190, 0.354, ROSE_LINE, 1.4)
    arrow(ax, 0.136, 0.286, 0.190, 0.246, TEAL_LINE, 1.4)

    box(ax, 0.194, 0.328, 0.178, 0.062, "left-right patterning\nof the heart tube",
        ROSE, 8, textcolour=ROSE_LINE)
    box(ax, 0.194, 0.210, 0.178, 0.062, "nephron induction\nin the kidney",
        TEAL, 8, textcolour=TEAL_LINE)

    arrow(ax, 0.378, 0.352, 0.432, 0.320)
    arrow(ax, 0.378, 0.244, 0.432, 0.276)

    box(ax, 0.436, 0.258, 0.178, 0.080,
        "one ciliary lesion,\ntwo organs affected", LILAC, 8.5, bold=True)

    arrow(ax, 0.620, 0.298, 0.672, 0.298, width=1.6)

    box(ax, 0.676, 0.322, 0.190, 0.056,
        "IFT and dynein-2\n5.0 to 5.5x enriched", SAGE, 7.5, textcolour=SAGE_LINE)
    box(ax, 0.676, 0.236, 0.190, 0.056,
        "BBSome\n1.0x, no signal", STONE, 7.5)
    ax.text(0.771, 0.170, "assembly matters,\ncargo sorting does not",
            ha="center", va="center", fontsize=7.5, color=MUTED,
            style="italic", linespacing=1.5)

    ax.text(0.5, 0.048,
            "Rare loss-of-function burden in 394,841 UK Biobank exomes   ·   "
            "9.7x excess at p < 0.001   ·   59 of 60 in the predicted direction",
            ha="center", fontsize=7.5, color=MUTED)

    save(fig, "fig00_graphical_abstract.png")


# ======================================================================
# INTRODUCTION FIGURE 1 - the two competing explanations
# ======================================================================

def intro_figure_hypotheses():
    fig, ax = plt.subplots(figsize=(8.4, 4.4))
    blank_axes(ax)

    ax.text(0.5, 0.950, "Two explanations for the same clinical observation",
            ha="center", fontsize=11.5, weight="bold", color=INK)

    ax.plot([0.5, 0.5], [0.20, 0.90], color=HAIRLINE, linewidth=1.2)

    # ---- panel A ----
    ax.text(0.25, 0.865, "A.  Shared genetic cause", ha="center",
            fontsize=9.5, weight="bold", color=INK)
    box(ax, 0.150, 0.700, 0.20, 0.072, "one gene, for example\na cilium gene",
        LILAC, 7.5)
    box(ax, 0.040, 0.482, 0.18, 0.072, "heart\nmalformation", ROSE, 8,
        textcolour=ROSE_LINE)
    box(ax, 0.270, 0.482, 0.18, 0.072, "kidney\nmalformation", TEAL, 8,
        textcolour=TEAL_LINE)
    arrow(ax, 0.212, 0.694, 0.148, 0.560)
    arrow(ax, 0.288, 0.694, 0.352, 0.560)
    ax.text(0.245, 0.412, "both organs affected from conception",
            ha="center", fontsize=8, color=MUTED, style="italic")
    box(ax, 0.070, 0.246, 0.35, 0.078,
        "prediction: the two conditions\nshare genetic architecture",
        SAGE, 8, textcolour=SAGE_LINE)

    # ---- panel B ----
    ax.text(0.75, 0.865, "B.  Kidney damage acquired later", ha="center",
            fontsize=9.5, weight="bold", color=INK)
    box(ax, 0.540, 0.700, 0.18, 0.072, "heart\nmalformation", ROSE, 8,
        textcolour=ROSE_LINE)
    box(ax, 0.775, 0.700, 0.185, 0.072, "cyanosis, surgery,\ncontrast, drugs",
        SAND, 7.5)
    box(ax, 0.655, 0.482, 0.19, 0.072, "chronic\nkidney disease", TEAL, 8,
        textcolour=TEAL_LINE)
    arrow(ax, 0.724, 0.736, 0.771, 0.736, width=1.4)
    arrow(ax, 0.855, 0.694, 0.775, 0.560)
    ax.text(0.75, 0.412, "kidney healthy at birth, damaged over decades",
            ha="center", fontsize=8, color=MUTED, style="italic")
    box(ax, 0.575, 0.246, 0.35, 0.078,
        "prediction: no shared\ngenetic architecture",
        ROSE, 8, textcolour=ROSE_LINE)

    ax.text(0.5, 0.085,
            "The two explanations make opposite genetic predictions,\n"
            "which is what makes the question testable.",
            ha="center", fontsize=8.5, color=INK, linespacing=1.6)

    save(fig, "figI1_two_hypotheses.png")


# ======================================================================
# INTRODUCTION FIGURE 2 - why naive gene overlap is misleading
# ======================================================================

def intro_figure_pleiotropy():
    fig, ax = plt.subplots(figsize=(9.0, 4.4))
    blank_axes(ax)

    ax.text(0.5, 0.955, "Why counting shared genes is not enough",
            ha="center", fontsize=11.5, weight="bold", color=INK)

    # ---- left half: one gene, many organs ----
    box(ax, 0.015, 0.470, 0.150, 0.100, "a developmental\ngene is disrupted",
        LILAC, 8)

    organs = ["heart", "kidney", "thumb", "palate", "brain", "ear"]
    fills = [ROSE, TEAL, STONE, STONE, STONE, STONE]
    text_colours = [ROSE_LINE, TEAL_LINE, INK, INK, INK, INK]

    top = 0.800
    step = 0.104
    for i in range(len(organs)):
        y = top - i * step
        box(ax, 0.245, y, 0.100, 0.066, organs[i], fills[i], 8,
            textcolour=text_colours[i])
        arrow(ax, 0.172, 0.520, 0.239, y + 0.033, width=0.9)

    fitted_text(ax, 0.295, 0.115, "one gene, many organs", 0.30,
                fontsize=8, colour=MUTED, italic=True)

    ax.plot([0.400, 0.400], [0.10, 0.87], color=HAIRLINE, linewidth=1.2)

    # ---- right half: what follows from it ----
    RIGHT_CENTRE = 0.700
    RIGHT_WIDTH = 0.545

    ax.text(RIGHT_CENTRE, 0.860, "The consequence", ha="center",
            fontsize=9.5, weight="bold", color=INK)

    fitted_text(ax, RIGHT_CENTRE, 0.745,
                "Every pair of organs now shares this gene.\n"
                "Heart and kidney share it. So do heart and thumb.\n"
                "A simple overlap test cannot tell these apart.",
                RIGHT_WIDTH, fontsize=8.5, linespacing=1.8)

    box(ax, 0.428, 0.470, RIGHT_WIDTH, 0.120,
        "What we observe:  diseases with a heart malformation\n"
        "affect a median of 11 organ systems;\n"
        "other diseases affect 5",
        SAND, 8.5)

    arrow(ax, RIGHT_CENTRE, 0.462, RIGHT_CENTRE, 0.408, width=1.4)

    box(ax, 0.428, 0.235, RIGHT_WIDTH, 0.165,
        "The solution used here:  model the kidney phenotype\n"
        "while adjusting for the number of organ systems\n"
        "affected and the number of annotated terms\n"
        "for each disease",
        SAGE, 8.5, textcolour=SAGE_LINE)

    save(fig, "figI2_pleiotropy_problem.png")


# ======================================================================
# METHODS FIGURE 1 - the analysis workflow
# ======================================================================

def methods_figure_workflow():
    fig, ax = plt.subplots(figsize=(9.4, 7.4))
    blank_axes(ax)

    ax.text(0.5, 0.980, "Analysis workflow", ha="center", fontsize=12.5,
            weight="bold", color=INK)

    LEFT = 0.040
    RIGHT = 0.985
    columns = row_positions(4, left=LEFT, right=RIGHT, gap=0.022)

    # ---------------- stage 1 ----------------
    band_label(ax, LEFT, 0.936, "STAGE 1      DATA SOURCES, ALL PUBLIC")
    sources = ["Human Phenotype\nOntology 2026-06-23",
               "HGNC, MANE v1.5\nand gnomAD v4.1",
               "Genebass\n394,841 exomes",
               "CKDGen and FinnGen R10\nsummary statistics"]
    for i in range(4):
        x, w = columns[i]
        box(ax, x, 0.856, w, 0.060, sources[i], STONE, 7.5)
        arrow(ax, x + w / 2.0, 0.852, 0.5, 0.836, width=0.8)

    # ---------------- stage 2 ----------------
    band_label(ax, LEFT, 0.812, "STAGE 2      BUILD THE CANDIDATE GENE LIST")
    arrow(ax, 0.5, 0.828, 0.5, 0.790, width=1.2)
    box(ax, 0.140, 0.712, 0.72, 0.058,
        "genes causing BOTH a heart and a kidney malformation        n = 750",
        TEAL, 9, bold=True, textcolour=TEAL_LINE)

    arrow(ax, 0.5, 0.712, 0.5, 0.676, width=1.4)

    filters = ["remove non-coding\nand non-autosomal",
               "remove infant-lethal\n(cannot reach a biobank)",
               "keep Tier 1,\ndrop immune genes",
               "take the ciliary subset"]
    counts = ["n = 676", "n = 611", "n = 169", "n = 37"]
    for i in range(4):
        x, w = columns[i]
        box(ax, x, 0.604, w, 0.066, filters[i], LILAC, 7.5)
        ax.text(x + w / 2.0, 0.581, counts[i], ha="center", fontsize=8,
                weight="bold", color=TEAL_LINE)
        if i < 3:
            arrow(ax, x + w + 0.003, 0.637, x + w + 0.019, 0.637, width=1.0)

    # ---------------- stage 3 ----------------
    band_label(ax, LEFT, 0.534, "STAGE 3      FOUR INDEPENDENT ANALYSES")
    analyses = ["Disease-level logistic\nregression, adjusting\nfor pleiotropy",
                "Rare-variant burden\nlookup in Genebass\n(pLoF and missense)",
                "Cross-trait LD score\nregression\n(CKDGen, FinnGen)",
                "MAGMA common-variant\ngene-set analysis\n(CKDGen eGFR)"]
    fills = [SAND, TEAL, TEAL, TEAL]
    for i in range(4):
        x, w = columns[i]
        box(ax, x, 0.398, w, 0.096, analyses[i], fills[i], 7.5)
        arrow(ax, x + w / 2.0, 0.518, x + w / 2.0, 0.500, width=1.0)

    # ---------------- stage 4 ----------------
    band_label(ax, LEFT, 0.356, "STAGE 4      WHAT EACH ANALYSIS CONCLUDED")
    outcomes = ["the CKD link vanishes\nOR 0.95,  p = 0.68",
                "ciliary genes enriched\n9.7x at p < 0.001",
                "no correlation\nrg = -0.03,  p = 0.67",
                "no common-variant\nenrichment,  p = 0.48"]
    result_fills = [ROSE, SAGE, ROSE, ROSE]
    result_text = [ROSE_LINE, SAGE_LINE, ROSE_LINE, ROSE_LINE]
    for i in range(4):
        x, w = columns[i]
        box(ax, x, 0.242, w, 0.076, outcomes[i], result_fills[i], 7.5,
            textcolour=result_text[i])
        arrow(ax, x + w / 2.0, 0.340, x + w / 2.0, 0.322, width=1.0)

    # ---------------- conclusion ----------------
    arrow(ax, 0.5, 0.238, 0.5, 0.208, width=1.6)
    box(ax, 0.070, 0.040, 0.86, 0.162,
        "CONCLUSION\n"
        "CHD and CKD do not share Mendelian genetic architecture.\n"
        "CHD and congenital kidney malformation do, through ciliary biology.",
        SAGE, 9.5, bold=True, textcolour=SAGE_LINE)

    save(fig, "figM1_workflow.png")


# ======================================================================
# METHODS FIGURE 2 - the gene filtering funnel
# ======================================================================

def methods_figure_funnel():
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    blank_axes(ax)

    ax.text(0.5, 0.960, "Building the candidate gene list",
            ha="center", fontsize=11.5, weight="bold", color=INK)
    ax.text(0.5, 0.915, "what was removed at each step, and why",
            ha="center", fontsize=8.5, color=MUTED)

    steps = [("genes causing a heart malformation", "1,660", "", STONE),
             ("also causing a kidney malformation", "750", "910 heart only", STONE),
             ("protein-coding and autosomal", "676", "74 non-coding, X or MT", LILAC),
             ("survivable to adulthood", "611", "65 infant-lethal only", LILAC),
             ("Tier 1 and non-immune", "169", "442 no CKD link or immune", LILAC),
             ("ciliary mechanism", "37", "132 other mechanisms", TEAL),
             ("excluding known cystogenic genes", "31", "6 already established", TEAL)]

    # The funnel narrows, but not so far that a label stops fitting.
    widest = 0.46
    narrowest = 0.32
    centre = 0.40
    y = 0.812
    row_height = 0.060
    spacing = 0.106

    for i in range(len(steps)):
        label, count, removed, colours = steps[i]
        fraction = i / float(len(steps) - 1)
        width = widest - (widest - narrowest) * fraction
        x = centre - width / 2.0

        box(ax, x, y, width, row_height, "", colours)
        ax.text(centre, y + row_height / 2.0, label,
                ha="center", va="center", fontsize=8, color=INK)

        ax.text(0.015, y + row_height / 2.0, "n = " + count,
                ha="left", va="center", fontsize=8.5, weight="bold",
                color=TEAL_LINE)
        if removed != "":
            ax.text(0.655, y + row_height / 2.0, "removed:  " + removed,
                    ha="left", va="center", fontsize=7.5, color=ROSE_LINE)

        if i < len(steps) - 1:
            arrow(ax, centre, y - 0.003, centre, y - 0.040, width=1.1)
        y = y - spacing

    ax.text(0.5, 0.038,
            "The 31-gene set is the pre-specified primary analysis. The six removed\n"
            "cystogenic genes are retained separately as a positive control.",
            ha="center", fontsize=8, color=MUTED, style="italic", linespacing=1.6)

    save(fig, "figM2_gene_filtering.png")



def blocking_mark(ax, x1, y1, x2, y2, position, colour):
    """Draw a short bar straight across an arrow, meaning "this path is shut".

    We work out the direction the arrow points, then draw a line at right
    angles to it. Swapping the two numbers over and flipping one sign is the
    standard trick for turning a direction into the direction at right angles
    to it.
    """
    along_x = x2 - x1
    along_y = y2 - y1
    length = (along_x ** 2 + along_y ** 2) ** 0.5
    across_x = -along_y / length
    across_y = along_x / length

    centre_x = x1 + along_x * position
    centre_y = y1 + along_y * position
    half = 0.045
    ax.plot([centre_x - across_x * half, centre_x + across_x * half],
            [centre_y - across_y * half, centre_y + across_y * half],
            color=colour, linewidth=2.4, solid_capstyle="round", zorder=8)


def causal_dag_figure():
    """Draws the causal diagram behind the pleiotropy adjustment.

    A causal diagram (a DAG) is just boxes with arrows saying what causes
    what. It is worth drawing because it makes clear WHY we adjust for the
    number of organ systems a disease affects, rather than that step looking
    like an arbitrary statistical tweak.

    Panel (a) is the problem. Panel (b) is what adjusting does about it.

    Layout note: the proxy box sits beside the confounder rather than below
    everything, so that nothing has to cross the arrow in the middle.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 3.9))

    for panel_index in range(len(axes)):
        ax = axes[panel_index]
        blank_axes(ax)
        ax.set_xlim(0, 1)
        # The drawing only uses the top part of the square, so we crop the
        # empty band underneath instead of leaving a big gap.
        ax.set_ylim(0.16, 1.0)

        adjusted = (panel_index == 1)

        if adjusted:
            title = "(b)  Adjusted: the shared path is blocked"
        else:
            title = "(a)  Unadjusted: both paths are open"
        ax.text(0.5, 0.99, title, ha="center", va="top",
                fontsize=9.5, color=INK, weight="bold")

        # --- the confounder, top left ---
        confounder_colours = STONE if adjusted else SAND
        box(ax, 0.02, 0.70, 0.45, 0.16,
            "Developmental gene disrupted\n(a syndromic, pleiotropic gene)",
            confounder_colours, fontsize=8.0)

        # --- what we actually measure, top right ---
        proxy_colours = SAGE if adjusted else STONE
        box(ax, 0.60, 0.70, 0.38, 0.16,
            "Measured proxy:\norgan systems affected,\nnumber of HPO terms",
            proxy_colours, fontsize=7.8)

        # Short dashed link between the proxy and the thing it stands for.
        ax.plot([0.595, 0.475], [0.78, 0.78], color=proxy_colours[1],
                linewidth=1.1, linestyle=(0, (2, 2.5)), zorder=1)
        if adjusted:
            link_words = "conditioned on"
        else:
            link_words = "stands in for"
        # The gap between the two boxes is narrower than the words, so the
        # label goes underneath them rather than between them.
        ax.text(0.535, 0.675, link_words, ha="center", va="top",
                fontsize=6.8, color=proxy_colours[1], style="italic")

        # --- the two organ outcomes, along the bottom ---
        box(ax, 0.02, 0.30, 0.40, 0.15, "Cardiac\nmalformation",
            ROSE, fontsize=8.5, bold=True)
        box(ax, 0.58, 0.30, 0.40, 0.15, "Kidney\nphenotype",
            TEAL, fontsize=8.5, bold=True)

        # --- the confounding path: one arrow down to each organ ---
        if adjusted:
            backdoor_colour = HAIRLINE
            backdoor_width = 1.0
        else:
            backdoor_colour = SAND[1]
            backdoor_width = 2.0

        arrow(ax, 0.16, 0.695, 0.20, 0.458, colour=backdoor_colour, width=backdoor_width)
        arrow(ax, 0.40, 0.695, 0.74, 0.458, colour=backdoor_colour, width=backdoor_width)

        if adjusted:
            blocking_mark(ax, 0.16, 0.695, 0.20, 0.458, 0.5, SAGE[1])
            blocking_mark(ax, 0.40, 0.695, 0.74, 0.458, 0.5, SAGE[1])
        else:
            ax.text(0.325, 0.60, "confounding\npath", ha="center", va="center",
                    fontsize=7.2, color=SAND[1], style="italic", linespacing=1.4)

        # --- the effect we actually want to estimate ---
        if adjusted:
            effect_colour = SAGE[1]
        else:
            effect_colour = HAIRLINE
        arrow(ax, 0.432, 0.375, 0.568, 0.375, colour=effect_colour, width=1.8)
        ax.text(0.50, 0.255, "the effect\nwe want", ha="center", va="top",
                fontsize=7.2, color=effect_colour, style="italic", linespacing=1.4)

    fig.text(0.5, 0.055,
             "The adjusted odds ratio is a causal quantity only if the measured proxy fully captures the shared developmental\n"
             "path. It is a stand-in, not the confounder itself, so some residual confounding should be expected.",
             ha="center", fontsize=7.8, color=MUTED, style="italic", linespacing=1.7)

    fig.subplots_adjust(bottom=0.20, top=0.96, wspace=0.09)
    save(fig, "figM3_causal_dag.png")


# ======================================================================

if __name__ == "__main__":
    print("Drawing schematic figures...")
    graphical_abstract()
    intro_figure_hypotheses()
    intro_figure_pleiotropy()
    methods_figure_workflow()
    methods_figure_funnel()
    causal_dag_figure()
    print("Done.")
