#!/usr/bin/env python3
"""Figure 8 — what follows for experiment design.

a  the confound, drawn: what a flight-vs-ground contrast can and cannot see
b  the three-arm design that separates hardware from gravity
c  the archive gap this analysis ran into
"""
from __future__ import annotations
import os, sys
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
from paths import FIGURES  # noqa: E402


def box(ax, x, y, w, h, label, fc, ec, fontsize=6.9):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.045",
                                facecolor=fc, edgecolor=ec, linewidth=0.9))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=fontsize)


def panel_a(ax):
    style.panel(ax, "a", "What the standard contrast cannot see")
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off"); ax.grid(False)
    box(ax, 0.4, 4.3, 3.9, 1.3, "FLIGHT\nsealed, lit, microgravity", "#FBD5BC", "#8C5A3C")
    box(ax, 0.4, 2.1, 3.9, 1.3, "GROUND CONTROL\nsealed, lit, 1 g", "#E4E7EB", "#5A6169")
    ax.add_patch(FancyArrowPatch((4.5, 4.95), (5.6, 3.9), arrowstyle="-|>",
                                 mutation_scale=12, color=style.INK, linewidth=1.3))
    ax.add_patch(FancyArrowPatch((4.5, 2.75), (5.6, 3.6), arrowstyle="-|>",
                                 mutation_scale=12, color=style.INK, linewidth=1.3))
    box(ax, 5.7, 3.1, 4.0, 1.3, "the published\nflight-vs-ground contrast", "white",
        style.INK)
    ax.text(2.35, 6.05, "CO$_2$ drawn down in BOTH arms", ha="center", fontsize=7.2,
            color=style.VERMILION, fontweight="bold")
    ax.text(7.7, 2.75, "hardware effect\nsubtracts out —\n25$\\times$ the gravity\n"
                       "effect, invisible", ha="center", va="top", fontsize=6.6,
            color=style.VERMILION)
    ax.text(0.4, 1.55, "What survives is a persistent ~5–9 %\n"
                       "assimilation penalty, and downstream\n"
                       "of it a carbon-starvation signature —\n"
                       "which is what the omics show.",
            fontsize=6.7, va="top", color=style.INK)


def panel_b(ax):
    style.panel(ax, "b", "The design that separates them")
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off"); ax.grid(False)
    arms = [("ARM 1  flight, sealed", "#FBD5BC", 4.6),
            ("ARM 2  ground, sealed", "#E4E7EB", 2.9),
            ("ARM 3  ground, vented", "#D6ECF7", 1.2)]
    for label, fc, y in arms:
        box(ax, 0.4, y, 5.0, 1.25, label, fc, "#5A6169", fontsize=7.0)
    ax.add_patch(FancyArrowPatch((5.7, 5.22), (7.0, 4.6), arrowstyle="-|>",
                                 mutation_scale=11, color=style.BLUE, linewidth=1.2))
    ax.add_patch(FancyArrowPatch((5.7, 3.52), (7.0, 4.3), arrowstyle="-|>",
                                 mutation_scale=11, color=style.BLUE, linewidth=1.2))
    ax.text(7.15, 4.45, "1 $-$ 2 = gravity", fontsize=7.0, color=style.BLUE,
            va="center", fontweight="bold")
    ax.add_patch(FancyArrowPatch((5.7, 3.22), (7.0, 2.4), arrowstyle="-|>",
                                 mutation_scale=11, color=style.VERMILION, linewidth=1.2))
    ax.add_patch(FancyArrowPatch((5.7, 1.82), (7.0, 2.1), arrowstyle="-|>",
                                 mutation_scale=11, color=style.VERMILION, linewidth=1.2))
    ax.text(7.15, 2.25, "2 $-$ 3 = hardware", fontsize=7.0, color=style.VERMILION,
            va="center", fontweight="bold")
    ax.text(0.4, 0.95, "Two arms confound the two effects.\nThree separate them — and the third\n"
                       "arm needs no flight.",
            fontsize=6.7, va="top", color=style.INK)


def panel_c(ax):
    style.panel(ax, "c", "The gap this ran into")
    ax.set_xlim(0, 10); ax.set_ylim(0, 6.4); ax.axis("off"); ax.grid(False)
    items = [
        ("0 of 66", "OSDR plant studies with a metabolome", style.VERMILION),
        ("6 of 567", "OSDR studies that are metabolite profiling,\n"
                     "all mouse, human, rat or microbial", style.INK),
        ("0", "in-canister CO$_2$ logs, any BRIC-class flight", style.VERMILION),
    ]
    y = 6.0
    for big, small, colour in items:
        ax.text(0.3, y, big, fontsize=13.5, fontweight="bold", color=colour, va="top")
        ax.text(0.35, y - 0.72, small, fontsize=6.7, va="top", color=style.INK)
        y -= 2.05
    ax.text(0.3, y + 0.62, "A ~20-compound targeted panel on archived material\n"
                           "would be the first plant spaceflight metabolome.",
            fontsize=6.7, va="top", color=style.INK, fontweight="bold")


def main():
    style.apply()
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.1))
    panel_a(axes[0]); panel_b(axes[1]); panel_c(axes[2])
    fig.subplots_adjust(top=0.82, bottom=0.03, left=0.01, right=0.99, wspace=0.16)
    fig.suptitle("Designing the experiment that could tell hardware from gravity",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.985)
    style.save(fig, "fig08_implications", FIGURES)


if __name__ == "__main__":
    main()
