"""Shared figure style: Okabe-Ito colourblind-safe palette + matplotlib rcParams.

Okabe-Ito is the canonical CVD-safe categorical set (validated by construction),
so the categorical hues used across all figures pass the colour checks without
per-figure tuning. Sector hues and the alpha/control accents are drawn from it in
a fixed order and reused identically in every figure.
"""
import matplotlib as mpl

# Okabe-Ito
OKABE = {
    "orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00",
    "purple": "#CC79A7", "black": "#000000",
}

# fixed semantic assignments (reused in every figure)
SECTOR = {"occipital": OKABE["blue"], "parietal": OKABE["green"], "frontal": OKABE["orange"]}
ALPHA_ACCENT = OKABE["vermillion"]     # alpha phenotypes / target highlight
CONTROL_GREY = "#8a8f98"               # non-alpha control bands
CORTEX_GREY = "#d9dce1"                # non-target cortex
BOUNDARY = "#3a3d42"                   # parcel boundary lines
INK = "#222428"

# diverging pair for signed maps (CVD-friendly blue<->orange with light midpoint)
DIVERGING = [OKABE["blue"], "#f2f2f2", OKABE["orange"]]


def apply_rc():
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300,
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
        "axes.edgecolor": INK, "axes.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": "#e6e8eb", "grid.linewidth": 0.7,
        "xtick.color": INK, "ytick.color": INK, "text.color": INK, "axes.labelcolor": INK,
        "legend.frameon": False, "figure.facecolor": "white", "savefig.facecolor": "white",
    })
