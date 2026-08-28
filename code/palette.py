# palette.py
#
# One place to keep the colours and the drawing helpers, so every figure
# in the thesis looks the same.
#
# The colours are pastels: soft fills with a slightly darker version of
# the same hue for the outline and the text. That keeps things readable
# when printed in black and white too, because the fills stay light.

import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# ---------------------------------------------------------------- colours

INK = "#33404A"          # main text, softer than pure black
MUTED = "#7C8894"        # captions and axis labels
HAIRLINE = "#C9D2D9"     # grid lines and dividers

# each entry is (light fill, darker line/text)
ROSE = ("#F6DCE0", "#C4737F")     # the heart
TEAL = ("#D9EBEF", "#4A93A3")     # the kidney
SAGE = ("#DFEDE1", "#5F9469")     # a positive or confirmed result
SAND = ("#F7EBD8", "#C0954E")     # a caution or an adjustment step
LILAC = ("#E8E4F0", "#8579A8")    # methods and processing steps
STONE = ("#EDF1F3", "#93A0AA")    # neutral background boxes

# short names for when only the line colour is wanted
ROSE_LINE = ROSE[1]
TEAL_LINE = TEAL[1]
SAGE_LINE = SAGE[1]
SAND_LINE = SAND[1]
STONE_LINE = STONE[1]


def apply_style():
    """Set the matplotlib defaults used by every figure."""
    plt.rcParams["font.family"] = "DejaVu Sans"
    plt.rcParams["font.size"] = 9
    plt.rcParams["text.color"] = INK
    plt.rcParams["axes.labelcolor"] = INK
    plt.rcParams["axes.edgecolor"] = HAIRLINE
    plt.rcParams["xtick.color"] = MUTED
    plt.rcParams["ytick.color"] = MUTED
    plt.rcParams["axes.spines.top"] = False
    plt.rcParams["axes.spines.right"] = False
    plt.rcParams["axes.grid"] = True
    plt.rcParams["grid.color"] = HAIRLINE
    plt.rcParams["grid.alpha"] = 0.55
    plt.rcParams["grid.linewidth"] = 0.5
    plt.rcParams["figure.dpi"] = 300


def blank_axes(ax):
    """Turn an axes into a plain drawing surface running 0 to 1 both ways."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def wrap(text, characters_per_line):
    """Break a line of text so it fits inside a box."""
    return "\n".join(textwrap.wrap(text, characters_per_line))


def box(ax, x, y, width, height, text, colours,
        fontsize=8, bold=False, textcolour=None, wrap_at=None):
    """
    Draw a soft rounded box with text inside it.

    x and y are the bottom-left corner. The text is centred. If wrap_at is
    given the text is wrapped to that many characters per line first, which
    is the easiest way to stop text spilling out of the box.
    """
    fill, line = colours
    patch = FancyBboxPatch((x, y), width, height,
                           boxstyle="round,pad=0,rounding_size=0.015",
                           facecolor=fill, edgecolor=line, linewidth=1.1)
    ax.add_patch(patch)

    if wrap_at is not None:
        text = wrap(text, wrap_at)
    if textcolour is None:
        textcolour = INK

    ax.text(x + width / 2.0, y + height / 2.0, text,
            ha="center", va="center", fontsize=fontsize, color=textcolour,
            weight="bold" if bold else "normal", linespacing=1.5)


def arrow(ax, x1, y1, x2, y2, colour=None, width=1.3):
    """Draw a thin arrow from one point to another."""
    if colour is None:
        colour = STONE_LINE
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=11, color=colour,
                                 linewidth=width, shrinkA=1, shrinkB=1))


def row_positions(n_boxes, left=0.06, right=0.98, gap=0.022):
    """
    Work out the x position and width for n boxes spread across the figure
    with an even gap between them. Returns a list of (x, width).
    """
    total_width = right - left
    box_width = (total_width - gap * (n_boxes - 1)) / n_boxes
    positions = []
    for i in range(n_boxes):
        positions.append((left + i * (box_width + gap), box_width))
    return positions


def stage_label(ax, y, text):
    """A small centred band label used between rows of a flow chart."""
    ax.text(0.5, y, text, ha="center", va="center", fontsize=8,
            weight="bold", color=MUTED)
