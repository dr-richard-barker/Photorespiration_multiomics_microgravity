#!/usr/bin/env python3
"""Figure 1 — the growth hardware, not the gravity, sets the atmosphere the plant sees.

a  the conductance chain, and where each term is set
b  enclosure CO2 over time, per hardware (LunarLeaf-CFD T7, all runs in microgravity)
c  12 h carbon gain as a percentage of Earth's (T11)
d  sealed-enclosure timescales (T8)
"""

from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
from paths import DATA, FIGURES  # noqa: E402

LL = os.path.join(DATA, "lunarleaf")
HW_COLOUR = {"BRIC light": style.VERMILION, "BRIC dark": "#8C2E0F",
             "CARA tape": style.ORANGE, "VEGGIE vented": style.BLUE,
             "open (ref)": style.GREY}


def panel_a(ax) -> None:
    style.panel(ax, "a", "The conductance chain")
    ax.set_xlim(0, 10); ax.set_ylim(0.42, 2.52); ax.axis("off"); ax.grid(False)

    boxes = [(0.35, "Bulk air\n$C_a$", style.GREY_LIGHT),
             (2.55, "Leaf surface\n$C_s$", "#FBE3D2"),
             (4.75, "Intercellular\n$C_i$", "#FBD5BC"),
             (6.95, "Chloroplast\n$C_c$, O$_2$", "#F7BE9B")]
    for x, label, fc in boxes:
        ax.add_patch(FancyBboxPatch((x, 1.35), 1.75, 0.85, boxstyle="round,pad=0.06",
                                    facecolor=fc, edgecolor="#8C5A3C", linewidth=0.8))
        ax.text(x + 0.875, 1.775, label, ha="center", va="center", fontsize=7.6)

    conductances = [(2.10, "$g_{bl}$", "boundary layer\nSET BY HARDWARE\n+ gravity",
                     style.VERMILION),
                    (4.30, "$g_s$", "stomata", style.GREY),
                    (6.50, "$g_m$", "mesophyll", style.GREY)]
    for x, sym, note, col in conductances:
        ax.add_patch(FancyArrowPatch((x - 0.30, 1.775), (x + 0.30, 1.775),
                                     arrowstyle="-|>", mutation_scale=11,
                                     color=col, linewidth=1.5))
        ax.text(x, 2.32, sym, ha="center", fontsize=8.6, color=col, fontweight="bold")
        ax.text(x, 1.16, note, ha="center", va="top", fontsize=6.1, color=col)

    ax.add_patch(FancyArrowPatch((7.83, 1.30), (7.83, 0.72), arrowstyle="-|>",
                                 mutation_scale=11, color="#8C5A3C", linewidth=1.4))
    ax.text(8.05, 1.00, "Rubisco:  $V_o/V_c = 2\\Gamma^*/C_c$", fontsize=7.4, va="center")
    ax.text(8.05, 0.55, "oxygenation rises as $C_c$ falls", fontsize=6.6,
            va="center", color=style.GREY)

    ax.text(0.35, 0.62, "A sealed, lit enclosure lowers $C_a$ itself —\n"
                        "upstream of every conductance in the chain.",
            fontsize=6.9, va="top", color=style.VERMILION)


def panel_b(ax) -> None:
    style.panel(ax, "b", "Enclosure CO$_2$ ($\\mu$g)")
    df = pd.read_csv(os.path.join(LL, "T7_hardware_timeseries.csv"))
    # VEGGIE and the open reference both sit flat at zero and would hide one another;
    # the vented trace is drawn slightly thicker and dashed over the reference.
    for case, g in df.groupby("case"):
        vented = case in ("VEGGIE vented", "open (ref)")
        ax.plot(g["step"] / 1000.0, g["dishmean_co2"],
                label=case, color=HW_COLOUR.get(case, style.GREY),
                linestyle="--" if vented else "-",
                linewidth=2.2 if case == "VEGGIE vented" else 1.6,
                alpha=0.95 if case == "VEGGIE vented" else 1.0,
                zorder=3 if case == "VEGGIE vented" else 2)
    ax.set_xlabel("model step ($\\times10^3$)")
    ax.set_ylabel("enclosure-mean CO$_2$ excess\n(model units)")
    ax.legend(loc="lower left", ncol=1, handlelength=1.6)
    ax.annotate("VEGGIE and the open\nreference both sit on zero",
                xy=(0.80, 0.0), xycoords=("axes fraction", "data"),
                xytext=(0.42, 0.80), textcoords="axes fraction",
                fontsize=6.2, color=style.GREY, ha="left", va="top",
                arrowprops=dict(arrowstyle="-", color=style.GREY, linewidth=0.6,
                                shrinkB=2))


def panel_c(ax) -> None:
    style.panel(ax, "c", "Carbon gained over 12 h")
    df = pd.read_csv(os.path.join(LL, "T11_photosynthesis_feedback.csv"))
    df = df.set_index("enclosure").loc[["Earth", "VEGGIE", "CARA", "BRIC"]].reset_index()
    cols = {"Earth": style.GREY, "VEGGIE": style.BLUE,
            "CARA": style.ORANGE, "BRIC": style.VERMILION}
    vals = df["12h carbon (% Earth)"]
    ax.bar(df["enclosure"], vals, color=[cols[e] for e in df["enclosure"]], width=0.62)
    for i, v in enumerate(vals):
        ax.text(i, v + 2.5, f"{v:.0f}%", ha="center", fontsize=8, fontweight="bold")
    ax.set_ylabel("% of Earth control")
    ax.set_ylim(0, 126)
    ax.tick_params(axis="x", labelrotation=0)


def panel_d(ax) -> None:
    style.panel(ax, "d", "How fast a sealed canister changes")
    ax.axis("off"); ax.grid(False)
    df = pd.read_csv(os.path.join(LL, "T8_enclosure_timescales.csv"))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    # Long descriptions wrap badly beside a right-aligned number, so the label column is
    # shortened here and the value column is given its own fixed x.
    short = {
        "BRIC light: CO2 400->~0 ppm (photosynthesis self-limits)":
            ("BRIC lit: CO$_2$ 400 $\\rightarrow$ ~0 ppm", True),
        "BRIC dark: CO2 400 ppm -> 1% (stress)":
            ("BRIC dark: CO$_2$ $\\rightarrow$ 1%", False),
        "BRIC dark: O2 21% -> 5% (hypoxia)":
            ("BRIC dark: O$_2$ 21% $\\rightarrow$ 5%", False),
        "CARA / VEGGIE: enclosure vents to cabin":
            ("CARA / VEGGIE: vents to cabin", False),
    }
    y = 0.96
    for _, r in df.iterrows():
        label, lit = short.get(str(r.iloc[0]), (str(r.iloc[0])[:38], False))
        value, unit = str(r.iloc[1]), str(r.iloc[2])
        colour = style.VERMILION if lit else style.INK
        ax.text(0.0, y, label, fontsize=6.8, va="top", color=colour)
        ax.text(0.03, y - 0.085, value if value == "bounded" else f"{value} {unit}",
                fontsize=8.6, va="top", ha="left", fontweight="bold", color=colour)
        y -= 0.235
    ax.text(0.0, y - 0.02,
            "The 7-minute drawdown is a mass balance, so it happens in the ground\n"
            "control too — which is why it cancels out of flight-vs-ground.",
            fontsize=6.6, va="top", color=style.GREY)


def main() -> None:
    style.apply()
    fig = plt.figure(figsize=(7.5, 4.7))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.62, 1.0],
                          hspace=0.46, wspace=0.44)
    panel_a(fig.add_subplot(gs[0, :]))
    panel_b(fig.add_subplot(gs[1, 0]))
    panel_c(fig.add_subplot(gs[1, 1]))
    panel_d(fig.add_subplot(gs[1, 2]))

    fig.suptitle("Growth hardware, not gravity, sets the atmosphere at the leaf",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.995)
    style.predicted_note(fig, "Gas transport from LunarLeaf-CFD (validated D2Q9 LBM solver); "
                              "all hardware runs are in microgravity.")
    style.save(fig, "fig01_hardware_atmosphere", FIGURES)


if __name__ == "__main__":
    main()
