"""Shared figure style for Phase 2 (colour-blind-safe palette; stability also encoded by line style)."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA, VIOLET, MAGENTA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#6250d6", "#e87ba4", "#8a8985"
INK, INK2 = "#0b0b0b", "#52514e"

plt.rcParams.update({
    "figure.dpi": 100, "savefig.dpi": 200, "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.grid": True, "grid.color": "#e4e3df",
    "grid.linewidth": 0.5, "legend.frameon": False, "legend.fontsize": 8, "lines.linewidth": 1.8,
})

K_LABEL = r"immune killing $k$ [1/(month$\cdot$E-unit)]"
DISCLAIMER = "Hypothetical illustrative parameters; mathematical results, not clinical predictions."


def baseline_note(p, extra=""):
    return (f"r={p.r}, q={p.q}, a={p.a}, $\\delta_D$={p.delta_D}, s={p.s}, $\\mu$={p.mu}, "
            f"$\\rho$={p.rho}, h={p.h}, $\\eta$={p.eta}, K={p.K} (rates 1/month)" + extra)


def footer(fig, text):
    fig.text(0.01, 0.005, text, fontsize=6.5, color=INK2, ha="left", va="bottom")
