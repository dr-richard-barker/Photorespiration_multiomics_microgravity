#!/usr/bin/env python3
"""Figure 7 — six spaceflight studies, and the axis that actually separates them.

a  the studies, laid out by enclosure and illumination
b  THE RESULT: starvation-over-photosynthesis separation, lit versus dark
c  the same six studies against the enclosure gradient the model predicted — which fails
d  gene-set responses per study, the whole matrix

The model made two predictions. Illumination separates the studies perfectly. The enclosure
gradient does not, and panel c shows the failure rather than hiding it.
"""

from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
from paths import FIGURES, TABLES  # noqa: E402

KEY = "starvation vs photosynthesis"
SET_ORDER = ["carbon_starvation_DIN", "photorespiration_core", "ath00630",
             "fermentation", "hypoxia_responsive", "ath00500", "ath00710",
             "photosynthesis_apparatus", "rubisco"]


def short(acc: str) -> str:
    return acc.replace("OSD-", "").replace("-light", " lit").replace("-dark", " dark")


def panel_a(ax, t12) -> None:
    style.panel(ax, "a", "The six studies")
    ax.set_xlim(-0.4, 2.6); ax.set_ylim(-0.6, 1.7)
    ax.set_xticks(range(3)); ax.set_xticklabels(["sealed\nBRIC", "tape\nCARA", "vented\nVEGGIE"])
    ax.set_yticks([0, 1]); ax.set_yticklabels(["dark", "lit"])
    ax.grid(False)
    xpos = {"sealed": 0, "tape": 1, "vented": 2}
    # Two studies share the sealed/dark cell, so cells are laid out first and their
    # members spread symmetrically about the cell centre.
    cells: dict[tuple[int, int], list] = {}
    for _, r in t12.iterrows():
        cells.setdefault((xpos[r["enclosure_class"]], 1 if r["light"] == "light" else 0),
                         []).append(r)
    for (x, y), members in cells.items():
        n = len(members)
        for i, r in enumerate(members):
            dx = 0.0 if n == 1 else (i - (n - 1) / 2) * 0.40
            ax.scatter(x + dx, y, s=330, color=style.LIGHT[r["light"]], alpha=0.9,
                       edgecolor="white", linewidth=1.2, zorder=3)
            ax.text(x + dx, y, short(r["accession"]).split()[0], ha="center",
                    va="center", fontsize=6.1, color="white", fontweight="bold",
                    zorder=4)
    ax.text(2, -0.42, "no dark vented study exists\nin OSDR", fontsize=6.2,
            color=style.GREY, ha="center", va="center")
    ax.text(-0.32, 1.55, "carbon retained by the enclosure $\\rightarrow$",
            fontsize=6.4, color=style.GREY)


def panel_b(ax, t12) -> None:
    style.panel(ax, "b", "Lit and dark separate completely")
    t = t12.sort_values(["light", "separation"])
    colours = [style.LIGHT[l] for l in t["light"]]
    ypos = np.arange(len(t))
    ax.barh(ypos, t["separation"], color=colours, height=0.62)
    ax.axvline(0, color=style.INK, linewidth=0.9)
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{short(a)}" for a in t["accession"]], fontsize=7)
    ax.set_xlabel("starvation $-$ photosynthesis\n(median log$_2$FC difference)")
    # Widen to the right so the statistics box has empty space of its own to sit in;
    # every dark bar runs left, so the lower-right quadrant stays clear.
    lo, hi = ax.get_xlim()
    ax.set_xlim(lo, hi + 3.4)

    lit = t[t["light"] == "light"]["separation"]
    dark = t[t["light"] == "dark"]["separation"]
    _u, p = stats.mannwhitneyu(lit, dark, alternative="greater")
    # Short enough to sit in the empty positive-x space beside the dark bars; the full
    # statement of the test belongs in the caption, not stacked on top of the data.
    ax.text(0.99, 0.03, f"every lit $>$ every dark\nMann–Whitney $p$ = {p:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.4,
            color=style.INK, linespacing=1.35)

    for spine in ("left",):
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)


def panel_c(ax, t12) -> None:
    style.panel(ax, "c", "The enclosure gradient does not hold")
    lit = t12[t12["light"] == "light"].dropna(subset=["separation"])
    for _, r in lit.iterrows():
        ax.scatter(r["carbon_pct_earth"], r["separation"], s=95,
                   color=style.ENCLOSURE[r["enclosure_class"]], zorder=3,
                   edgecolor="white", linewidth=1.0)
        ax.annotate(short(r["accession"]), (r["carbon_pct_earth"], r["separation"]),
                    textcoords="offset points", xytext=(0, 11), ha="center", fontsize=6.4)
    rho, prho = stats.spearmanr(lit["carbon_pct_earth"], lit["separation"])
    ax.set_xlabel("carbon the enclosure lets the plant keep\n(% of Earth, CFD prediction)")
    ax.set_ylabel("starvation $-$ photosynthesis")
    ax.set_xlim(-14, 128)
    ax.set_ylim(-0.35, 3.7)
    xs = np.linspace(0, 100, 10)
    ax.plot(xs, 2.6 - 0.022 * xs, color=style.GREY, linestyle="--", linewidth=1.0,
            label="what the model predicted")
    ax.legend(loc="upper right")
    ax.text(0.03, 0.03, f"Spearman $\\rho$ = {rho:+.2f}, $p$ = {prho:.2f}\n"
                        "three studies cannot rank three classes",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.4,
            color=style.VERMILION)


def panel_d(ax, t11) -> None:
    style.panel(ax, "d", "Every gene set, every study")
    pivot = t11.pivot_table(index="set", columns="accession", values="shift")
    order = [s for s in SET_ORDER if s in pivot.index]
    lit_first = [c for c in pivot.columns if "light" in c or c in ("OSD-522", "OSD-427")]
    rest = [c for c in pivot.columns if c not in lit_first]
    pivot = pivot.loc[order, lit_first + rest]

    vmax = float(np.nanmax(np.abs(pivot.values)))
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels([short(c) for c in pivot.columns], rotation=40, ha="right",
                       fontsize=6.6)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels([style.PRETTY.get(s, s).replace("\n", " ") for s in pivot.index],
                       fontsize=6.4)
    ax.grid(False)
    ax.axvline(len(lit_first) - 0.5, color=style.INK, linewidth=1.8)
    # The study labels carry the illumination themselves. Badges above the heatmap ran
    # into the panel title and badges below ran into the figure footer; colouring the
    # ticks says the same thing and cannot collide with anything.
    for tick, col in zip(ax.get_xticklabels(), lit_first + rest):
        tick.set_color(style.ORANGE if col in lit_first else style.BLUE)
        tick.set_fontweight("bold")
    ax.text(0.995, 1.035, "lit", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.2, fontweight="bold", color=style.ORANGE)
    ax.text(1.0, 1.035, "  |  dark", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=7.2, fontweight="bold", color=style.BLUE)
    cb = plt.colorbar(im, ax=ax, fraction=0.030, pad=0.02)
    cb.set_label("shift vs all other genes (log$_2$FC)", fontsize=6.4)
    cb.ax.tick_params(labelsize=6)


def main() -> None:
    style.apply()
    t11 = pd.read_csv(os.path.join(TABLES, "T11_hardware_ladder.tsv"), sep="\t")
    t12 = pd.read_csv(os.path.join(TABLES, "T12_ladder_contrasts.tsv"), sep="\t")
    t12 = t12[t12["contrast"] == KEY].copy()

    fig = plt.figure(figsize=(7.5, 6.4))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.85, 1.05], hspace=0.62, wspace=0.52)
    panel_a(fig.add_subplot(gs[0, 0]), t12)
    panel_b(fig.add_subplot(gs[0, 1]), t12)
    panel_c(fig.add_subplot(gs[0, 2]), t12)
    panel_d(fig.add_subplot(gs[1, :]), t11)

    fig.suptitle("Illumination, not sealing, is what orders the spaceflight response",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.988)
    style.predicted_note(fig, "Six matched flight-vs-ground contrasts from NASA OSDR; "
                              "gene sets from KEGG. Study selection in data/study_registry.tsv.")
    style.save(fig, "fig07_hardware_ladder", FIGURES)


if __name__ == "__main__":
    main()
