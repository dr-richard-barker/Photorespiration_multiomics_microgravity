#!/usr/bin/env python3
"""Figure: ABRS (Ventilated) vs BRIC-LED (Sealed) Multi-Omics Comparison.

Four panels:
  A: Operating Points & Gas Exchange: Assimilation (A) and Photorespiration (Vo)
     from sealed drawdown (100 ppm) to ventilated ISS cabin CO2 (3500 ppm).
  B: ABRS Organ-Specific Transcriptome: Gene set shifts across Shoots, Roots,
     Hypocotyls, and Whole Plants.
  C: Extended Hardware Ladder: Starvation vs Photosynthesis separation across
     7 studies (illuminated vs dark; sealed vs tape vs vented).
  D: Pathway Significance Contrast: BRIC-LED vs ABRS pathway p-values,
     separating hardware artifacts (ER/UPR stress) from conserved spaceflight responses.
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "figures")))
import style  # noqa: E402
from fvcb import LeafParams, solve_operating_point  # noqa: E402
from paths import ensure  # noqa: E402

ABRS_TABLES = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "tables"
))
ABRS_FIGURES = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "figures"
))


def plot_panel_a(ax):
    """Panel A: Operating Points across CO2 concentrations."""
    ca_vals = np.logspace(np.log10(60), np.log10(4000), 100)
    p = LeafParams()
    A_ug, Vo_ug, phi_ug = [], [], []

    for ca in ca_vals:
        sol = solve_operating_point(0.291, p, Ca=ca)  # micro-g rosette
        A_ug.append(sol["A"])
        Vo_ug.append(sol["Vo"])
        phi_ug.append(sol["phi"] * 100)

    ax.plot(ca_vals, A_ug, color=style.BLUE, lw=2.2, label=r"Assimilation $A$")
    ax.plot(ca_vals, Vo_ug, color=style.VERMILION, lw=2.2, label=r"Photorespiration $V_o$")

    # Mark BRIC operating point (100 ppm)
    ax.axvline(100, color=style.VERMILION, ls="--", alpha=0.7, lw=1.5)
    ax.text(105, 22, "Sealed BRIC\n(~100 ppm)", color=style.VERMILION, fontsize=8.5, weight="bold")

    # Mark Ambient Earth (400 ppm)
    ax.axvline(400, color=style.GREY, ls=":", alpha=0.7, lw=1.2)
    ax.text(420, 22, "Earth Ambient\n(400 ppm)", color=style.GREY, fontsize=8)

    # Mark ABRS operating point (3500 ppm)
    ax.axvline(3500, color=style.BLUE, ls="--", alpha=0.7, lw=1.5)
    ax.text(1800, 10, "Ventilated ABRS\n(~3500 ppm ISS cabin)", color=style.BLUE, fontsize=8.5, weight="bold")

    ax.set_xscale("log")
    ax.set_xlabel(r"Canister $CO_2$ ($C_a$, $\mu$mol mol$^{-1}$)", fontsize=9.5)
    ax.set_ylabel(r"Flux ($\mu$mol m$^{-2}$ s$^{-1}$)", fontsize=9.5)
    ax.set_title("A. FvCB Photosynthesis & Photorespiration Operating Range", fontsize=10, weight="bold", loc="left")
    ax.legend(loc="center right", frameon=False, fontsize=8.5)
    ax.set_ylim(-1, 32)


def plot_panel_b(ax):
    """Panel B: ABRS Organ Comparison."""
    t13_path = os.path.join(ABRS_TABLES, "T13_abrs_organ_comparison.tsv")
    if not os.path.exists(t13_path):
        ax.text(0.5, 0.5, "T13 table missing", ha="center")
        return

    df = pd.read_csv(t13_path, sep="\t")
    target_sets = ["photosynthesis_apparatus", "photorespiration_core", "carbon_starvation_DIN", "cell_wall"]
    sub = df[df["set"].isin(target_sets)].copy()

    organs = ["shoot", "hypocotyl", "root", "whole_plant"]
    x = np.arange(len(organs))
    width = 0.18

    set_meta = {
        "photosynthesis_apparatus": ("Photosystems", style.BLUE),
        "photorespiration_core": ("Photorespiration Core", style.GREY),
        "carbon_starvation_DIN": ("DIN Starvation", style.VERMILION),
        "cell_wall": ("Cell Wall Remodeling", style.ORANGE),
    }

    for i, sname in enumerate(target_sets):
        label, color = set_meta[sname]
        shifts = []
        for o in organs:
            r = sub[(sub["organ"] == o) & (sub["set"] == sname)]
            shifts.append(r["shift"].iloc[0] if not r.empty else 0.0)
        ax.bar(x + i * width - 1.5 * width, shifts, width, label=label, color=color, alpha=0.85)

    ax.axhline(0, color=style.INK, lw=0.8, alpha=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(["Shoots\n(Leaves)", "Hypocotyls", "Roots", "Whole\nPlant"], fontsize=9)
    ax.set_ylabel(r"Median Shift ($\log_2$FC vs background)", fontsize=9.5)
    ax.set_title("B. ABRS Organ-Specific Transcriptome Remodeling", fontsize=10, weight="bold", loc="left")
    ax.legend(loc="upper left", frameon=False, fontsize=8)


def plot_panel_c(ax):
    """Panel C: Extended Hardware Ladder."""
    t12_path = os.path.join(ABRS_TABLES, "T12_abrs_extended_ladder_contrasts.tsv")
    if not os.path.exists(t12_path):
        ax.text(0.5, 0.5, "T12 table missing", ha="center")
        return

    df = pd.read_csv(t12_path, sep="\t")
    key = "starvation vs photosynthesis"
    sub = df[df["contrast"] == key].sort_values("separation").reset_index(drop=True)

    y = np.arange(len(sub))
    colors = [style.ORANGE if r["light"] == "light" else style.BLUE for _, r in sub.iterrows()]

    bars = ax.barh(y, sub["separation"], color=colors, height=0.6, alpha=0.85)
    ax.axvline(0, color=style.INK, lw=0.8, alpha=0.5)

    labels = []
    for _, r in sub.iterrows():
        hw = r["hardware"]
        acc = r["accession"]
        labels.append(f"{acc} ({hw})")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel(r"Starvation vs Photosynthesis Separation ($\log_2$FC)", fontsize=9.5)
    ax.set_title("C. Extended Hardware Ladder (Illumination & Ventilation)", fontsize=10, weight="bold", loc="left")

    # Add custom legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=style.ORANGE, alpha=0.85, label="Illuminated hardware"),
        Patch(facecolor=style.BLUE, alpha=0.85, label="Dark hardware"),
    ]
    ax.legend(handles=legend_elements, loc="lower right", frameon=False, fontsize=8)


def plot_panel_d(ax):
    """Panel D: Pathway Significance Contrast (BRIC-LED vs ABRS)."""
    t10_path = os.path.join(ABRS_TABLES, "T10_abrs_vs_bric_comparison.tsv")
    if not os.path.exists(t10_path):
        ax.text(0.5, 0.5, "T10 table missing", ha="center")
        return

    df = pd.read_csv(t10_path, sep="\t")
    # Pick representative key pathways
    selected = [
        "Protein processing in endoplasmic reticulum",
        "Photosynthesis",
        "Starch and sucrose metabolism",
        "Biotic Stress",
        "Raffinose metabolism",
        "Plant hormone signal transduction",
        "Circadian rhythm - plant",
    ]
    sub = df[df["pathway"].isin(selected)].copy()

    # Convert p-values to -log10
    bric_p = [-np.log10(max(float(p), 1e-10)) for p in sub["bric_p_combined"]]
    abrs_p = []
    for p in sub["abrs_p_gene"]:
        try:
            val = float(p)
            abrs_p.append(-np.log10(max(val, 1e-10)))
        except ValueError:
            abrs_p.append(0.0)

    y = np.arange(len(sub))
    height = 0.35

    ax.barh(y + height / 2, bric_p, height, color=style.VERMILION, alpha=0.85, label="Sealed BRIC-LED")
    ax.barh(y - height / 2, abrs_p, height, color=style.BLUE, alpha=0.85, label="Ventilated ABRS")

    # Threshold line at p = 0.05 (-log10 = 1.30)
    ax.axvline(-np.log10(0.05), color=style.GREY, ls=":", lw=1.2)
    ax.text(-np.log10(0.05) + 0.1, len(sub) - 0.5, "p = 0.05", color=style.GREY, fontsize=8)

    ylabels = [
        p.replace("Protein processing in endoplasmic reticulum", "Protein processing in\nendoplasmic reticulum")
         .replace("Plant hormone signal transduction", "Plant hormone signal\ntransduction")
         .replace("Starch and sucrose metabolism", "Starch & sucrose\nmetabolism")
        for p in sub["pathway"]
    ]
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel(r"$-\log_{10}(p\text{-value})$", fontsize=9.5)
    ax.set_title("D. Pathway Significance: Sealed Hardware vs Controlled Airflow", fontsize=10, weight="bold", loc="left")
    ax.legend(loc="lower right", frameon=False, fontsize=8)


def main():
    style.apply()
    ensure(ABRS_FIGURES)
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    fig.subplots_adjust(hspace=0.35, wspace=0.52)

    plot_panel_a(axes[0, 0])
    plot_panel_b(axes[0, 1])
    plot_panel_c(axes[1, 0])
    plot_panel_d(axes[1, 1])

    png_path = os.path.join(ABRS_FIGURES, "fig_abrs_comparison.png")
    pdf_path = os.path.join(ABRS_FIGURES, "fig_abrs_comparison.pdf")

    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, bbox_inches="tight")
    plt.close()

    print(f"Generated Figure:")
    print(f"  {png_path}")
    print(f"  {pdf_path}")


if __name__ == "__main__":
    main()
