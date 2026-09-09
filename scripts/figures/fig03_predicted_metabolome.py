#!/usr/bin/env python3
"""Figure 3 — the metabolite layer OSDR does not contain, and where it comes from.

a  the 21 compounds, grouped by the modelled flux that drives them
b  the same compounds under the hardware contrast, which is 25x larger
c  what the model does NOT predict, and why

Every value here is model output. Nothing on this figure was measured.
"""
from __future__ import annotations
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
from paths import FIGURES, TABLES  # noqa: E402

DRIVER_LABEL = {"Vo": "C2 photorespiratory cycle\n(pool $\\propto$ oxygenation $V_o$)",
                "A": "Net photosynthate\n(pool $\\propto$ assimilation $A$)",
                "RuBP": "RuBP\n(rises as Rubisco is CO$_2$-limited)",
                "starvation": "Carbon starvation\n(pool $\\propto$ $-$assimilation)"}
DRIVER_COLOUR = {"Vo": style.GREY, "A": style.BLUE,
                 "RuBP": style.GREEN, "starvation": style.VERMILION}
ORDER = ["Vo", "A", "RuBP", "starvation"]

EXCLUDED = [
    ("TCA intermediates", "the model does not resolve respiratory flux, and carbon\n"
                          "limitation can deplete OR raise these pools"),
    ("Fermentation products", "fermentation needs hypoxia; a lit canister is releasing O$_2$"),
    ("Glycine / serine ratio", "both pools track the same flux, so the ratio cannot move —\n"
                               "reporting it would be circular"),
]


def bars(ax, df, col, title, letter):
    style.panel(ax, letter, title)
    df = df.copy()
    df["order"] = df["driver"].map({d: i for i, d in enumerate(ORDER)})
    df = df.sort_values(["order", "compound"], ascending=[True, False])
    y = np.arange(len(df))
    ax.barh(y, df[col], color=[DRIVER_COLOUR[d] for d in df["driver"]], height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels(df["compound"], fontsize=6.3)
    ax.axvline(0, color=style.INK, linewidth=0.9)
    ax.set_xlabel("predicted log$_2$ fold change")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    return df


def main():
    style.apply()
    df = pd.read_csv(os.path.join(TABLES, "T06_predicted_compounds.tsv"),
                     sep="\t", comment="#")
    fig = plt.figure(figsize=(7.5, 5.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 0.74], wspace=0.62)

    ax1 = fig.add_subplot(gs[0, 0])
    bars(ax1, df, "log2FC_FLT_vs_GC", "Flight vs ground (sealed)", "a")
    ax1.set_xlim(-0.16, 0.16)

    ax2 = fig.add_subplot(gs[0, 1])
    d2 = bars(ax2, df, "log2FC_BRIC_vs_VENTED", "Sealed vs vented hardware", "b")
    ax2.set_yticklabels([])
    ax2.set_xlim(-3.4, 3.4)
    ax2.text(0.5, -0.115, "note the axis: 25$\\times$ wider than panel a",
             transform=ax2.transAxes, ha="center", va="top", fontsize=6.6,
             color=style.VERMILION)

    handles = [plt.Rectangle((0, 0), 1, 1, color=DRIVER_COLOUR[d]) for d in ORDER]
    ax1.legend(handles, [DRIVER_LABEL[d] for d in ORDER], loc="upper left",
               bbox_to_anchor=(0.0, -0.14), fontsize=6.2, ncol=1, handlelength=1.1,
               labelspacing=0.85, borderaxespad=0)

    ax3 = fig.add_subplot(gs[0, 2])
    style.panel(ax3, "c", "Deliberately not predicted")
    ax3.axis("off"); ax3.grid(False)
    y = 0.97
    for what, why in EXCLUDED:
        ax3.text(0, y, what, fontsize=7.2, fontweight="bold", va="top", color=style.INK)
        ax3.text(0, y - 0.065, why, fontsize=6.3, va="top", color=style.GREY)
        y -= 0.235
    ax3.text(0, y + 0.04,
             "The prediction has only four distinct\nvalues, because every pool in a block\n"
             "inherits the same flux ratio. That is\nwhat the model knows: it resolves\n"
             "fluxes, not individual pool kinetics.",
             fontsize=6.3, va="top", color=style.INK)

    fig.suptitle("The metabolite layer had to be predicted — OSDR holds no plant metabolome",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.985)
    fig.subplots_adjust(top=0.88, bottom=0.30, left=0.16, right=0.985)
    style.predicted_note(fig, "EVERY VALUE ON THIS FIGURE IS MODEL OUTPUT. Nothing here was "
                              "measured. 0 of 66 OSDR plant studies has a metabolome.")
    style.save(fig, "fig03_predicted_metabolome", FIGURES)


if __name__ == "__main__":
    main()
