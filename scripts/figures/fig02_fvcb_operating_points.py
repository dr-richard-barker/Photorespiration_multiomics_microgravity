#!/usr/bin/env python3
"""Figure 2 — what a drawn-down enclosure does to Rubisco, and what gravity does not do.

a  oxygenation fraction against canister CO2, with the four operating points marked
b  THE REFUTATION: the gravity gap in oxygenation SHRINKS as CO2 falls, it does not amplify
c  the assimilation penalty that does survive, across three canopy scales
d  hardware effect against gravity effect, on the same axis

Panel b is the panel we did not expect to draw. The working hypothesis was that CO2
starvation would amplify the microgravity boundary-layer penalty. It does the opposite:
as assimilation falls toward the compensation point, the flux through the boundary layer
falls with it, so the CO2 drop across that layer (A/g_bl) tends to zero and the two
gravities converge.
"""

from __future__ import annotations

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
from fvcb import LeafParams, load_cfd_sweep, solve_operating_point  # noqa: E402
from paths import FIGURES, TABLES  # noqa: E402

CA_GRID = np.linspace(60, 420, 145)
SCALES = ["leaf", "rosette", "canopy"]
SCALE_COLOUR = {"leaf": style.SKY, "rosette": style.BLUE, "canopy": "#003D63"}


def gbl_by_scale() -> dict[str, tuple[float, float]]:
    """scale -> (Earth g_bl, microgravity g_bl) from the CFD export."""
    out: dict[str, dict[str, float]] = {}
    for r in load_cfd_sweep():
        key = "earth" if r["gravity_g"] > 5.0 else ("ug" if r["gravity_g"] < 0.5 else None)
        if key:
            out.setdefault(r["scale"], {}).setdefault(key, r["g_bl"])
    return {s: (v["earth"], v["ug"]) for s, v in out.items() if {"earth", "ug"} <= set(v)}


def sweep(g_bl: float, p: LeafParams) -> pd.DataFrame:
    rows = [solve_operating_point(g_bl, p, Ca=float(ca)) for ca in CA_GRID]
    return pd.DataFrame(rows, index=CA_GRID)


def panel_a(ax, gbl, p) -> None:
    style.panel(ax, "a", "Oxygenation rises as the canister empties")
    earth, ug = gbl["rosette"]
    for g, colour, label in ((earth, style.GREY, "1 g"), (ug, style.INK, "microgravity")):
        s = sweep(g, p)
        ax.plot(CA_GRID, s["phi"] * 100, color=colour, label=label,
                linestyle="-" if colour == style.INK else "--")

    ops = pd.read_csv(os.path.join(TABLES, "T05_operating_points.tsv"), sep="\t")
    ops = ops.set_index("arm")

    # Flight and ground control sit almost exactly on top of one another, because the
    # sealed canister draws BOTH arms down to the same place. That coincidence is the
    # point of the whole figure, so they are labelled once, as a pair.
    for arm, colour in (("VENTED", style.BLUE), ("GC", style.GREY), ("FLT", style.VERMILION)):
        r = ops.loc[arm]
        ax.plot(r["Ca_umol_mol"], r["phi_percent"], "o", color=colour, markersize=6.5,
                markeredgecolor="white", markeredgewidth=0.9, zorder=5)

    ax.annotate("flight AND ground control\nboth sit here, in the sealed\ncanister",
                xy=(ops.loc["FLT", "Ca_umol_mol"], ops.loc["FLT", "phi_percent"]),
                xytext=(168, 47), fontsize=6.7, color=style.VERMILION, ha="left",
                arrowprops=dict(arrowstyle="-|>", color=style.VERMILION, linewidth=0.8))
    ax.annotate("vented hardware\nholds ambient", xy=(ops.loc["VENTED", "Ca_umol_mol"],
                                                      ops.loc["VENTED", "phi_percent"]),
                xytext=(-14, 20), textcoords="offset points",
                fontsize=6.7, color=style.BLUE, ha="right",
                arrowprops=dict(arrowstyle="-|>", color=style.BLUE, linewidth=0.8))

    ax.set_xlabel("canister CO$_2$ ($\\mu$mol mol$^{-1}$)")
    ax.set_ylabel("oxygenation fraction $\\varphi$ (%)")
    ax.legend(loc="upper right")


def panel_b(ax, gbl, p) -> None:
    style.panel(ax, "b", "The gravity gap shrinks, it does not amplify")
    for scale in SCALES:
        earth, ug = gbl[scale]
        d = (sweep(ug, p)["phi"].values - sweep(earth, p)["phi"].values) * 100
        ax.plot(CA_GRID, d, color=SCALE_COLOUR[scale], label=scale)
    ax.axhline(0, color=style.GREY, linewidth=0.8)
    ax.set_xlabel("canister CO$_2$ ($\\mu$mol mol$^{-1}$)")
    ax.set_ylabel("$\\varphi$ in $\\mu$g $-$ $\\varphi$ at 1 g\n(percentage points)")
    ax.legend(loc="upper left", title="canopy scale", title_fontsize=6.8,
              bbox_to_anchor=(0.0, 1.02))
    ax.set_ylim(-0.12, 1.95)
    ax.annotate("as CO$_2$ runs out, so does the flux\nthrough the boundary layer —\n"
                "the two gravities converge",
                xy=(78, 0.06), xytext=(175, 0.42), fontsize=6.6, color=style.INK,
                arrowprops=dict(arrowstyle="-|>", color=style.INK, linewidth=0.7))
    ax.text(0.98, 0.02, "we expected this to widen", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.6, color=style.VERMILION, style="italic")


def panel_c(ax, gbl, p) -> None:
    style.panel(ax, "c", "The penalty that does survive")
    for scale in SCALES:
        earth, ug = gbl[scale]
        a_e, a_u = sweep(earth, p)["A"].values, sweep(ug, p)["A"].values
        with np.errstate(divide="ignore", invalid="ignore"):
            pen = np.where(a_e > 0.05, (a_u - a_e) / a_e * 100, np.nan)
        ax.plot(CA_GRID, pen, color=SCALE_COLOUR[scale], label=scale)
    ax.set_xlabel("canister CO$_2$ ($\\mu$mol mol$^{-1}$)")
    ax.set_ylabel("assimilation in $\\mu$g\nrelative to 1 g (%)")
    ax.legend(loc="center right", title="canopy scale", title_fontsize=6.8)
    ax.set_ylim(-10.4, -3.6)
    ax.text(0.98, 0.03, "roughly flat across the whole range —\n"
                        "a persistent 5–9 % deficit",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.6, color=style.INK)


def panel_d(ax) -> None:
    style.panel(ax, "d", "What actually moves the metabolism")
    prov = pd.read_csv(os.path.join(TABLES, "T06_predicted_compounds.tsv"),
                       sep="\t", comment="#")
    grav = prov["log2FC_FLT_vs_GC"].abs().max()
    hard = prov["log2FC_BRIC_vs_VENTED"].abs().max()
    ax.bar(["microgravity\n(flight vs ground)", "hardware\n(sealed vs vented)"],
           [grav, hard], color=[style.BLUE, style.VERMILION], width=0.5)
    for i, v in enumerate([grav, hard]):
        ax.text(i, v + hard * 0.03, f"{v:.2f}", ha="center", fontsize=8.5,
                fontweight="bold")
    ax.set_ylabel("largest predicted |log$_2$FC|\non a metabolite pool")
    ax.set_ylim(0, hard * 1.30)
    ax.text(0.5, hard * 0.66, f"{hard / grav:.0f}$\\times$ larger — but it hits\n"
                              "flight and ground alike, so it\ncancels out of the contrast",
            ha="center", va="center", fontsize=6.7, color=style.INK,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor=style.GREY_LIGHT, linewidth=0.6))


def main() -> None:
    style.apply()
    p = LeafParams()
    gbl = gbl_by_scale()

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6))
    panel_a(axes[0, 0], gbl, p)
    panel_b(axes[0, 1], gbl, p)
    panel_c(axes[1, 0], gbl, p)
    panel_d(axes[1, 1])
    fig.subplots_adjust(hspace=0.46, wspace=0.34, top=0.90, bottom=0.09)

    fig.suptitle("A drawn-down enclosure, not microgravity, is what reaches Rubisco",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.985)
    style.predicted_note(fig, "Farquhar–von Caemmerer–Berry model with Bernacchi (2001, 2003) "
                              "kinetics, driven by LunarLeaf-CFD boundary-layer conductance.")
    style.save(fig, "fig02_fvcb_operating_points", FIGURES)


if __name__ == "__main__":
    main()
