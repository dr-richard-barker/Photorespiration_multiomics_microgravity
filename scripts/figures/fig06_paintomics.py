#!/usr/bin/env python3
"""Figure 6 — what PaintOmics made of the three layers together.

a  identifier mapping, per omic per database
b  pathway enrichment: the significant ones, and where photorespiration lands
c  metabolite hub analysis — predicted compounds ranked by REAL differentially
   expressed genes in their KEGG neighbourhood
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

# Mapping counts as reported by PaintOmics on job m1z16Qg3DK, recorded in methods/.
MAPPING = pd.DataFrame({
    "omic": ["Gene expression", "Gene expression", "Proteomics", "Proteomics",
             "Metabolomics", "Metabolomics"],
    "database": ["KEGG", "MapMan", "KEGG", "MapMan", "KEGG", "MapMan"],
    "mapped": [20456, 20983, 5002, 4998, 21, 21],
    "total": [20983, 20983, 5160, 5160, 21, 21],
})
HIGHLIGHT = {"Photosynthesis", "Photosynthesis - antenna proteins", "photosynthesis",
             "Starch and sucrose metabolism"}


def panel_a(ax):
    style.panel(ax, "a", "Identifier mapping")
    # Horizontal: three long omic names will not fit as x tick labels in a narrow panel.
    omics = ["Gene expression", "Proteomics", "Metabolomics"]
    h = 0.36
    y = np.arange(len(omics))[::-1]
    for i, (db, colour) in enumerate((("KEGG", style.BLUE), ("MapMan", style.GREEN))):
        sub = MAPPING[MAPPING["database"] == db].set_index("omic").loc[omics]
        pct = (sub["mapped"] / sub["total"] * 100).values
        ax.barh(y + (0.5 - i) * h, pct, height=h, color=colour, label=db)
        for yi, val in zip(y, pct):
            ax.text(val + 2, yi + (0.5 - i) * h, f"{val:.0f}%", va="center",
                    fontsize=6.2)
    ax.set_yticks(y)
    ax.set_yticklabels(["Gene\nexpression", "Proteomics", "Metabolomics"], fontsize=6.8)
    ax.set_xlabel("features mapped (%)")
    ax.set_xlim(0, 132)
    ax.legend(loc="lower center", ncol=2, fontsize=6.4,
              bbox_to_anchor=(0.5, -0.30))
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)


def panel_b(ax):
    style.panel(ax, "b", "Pathway enrichment, all three layers")
    sig = pd.read_csv(os.path.join(TABLES, "T08_paintomics_significant.tsv"),
                      sep="\t", comment="#")
    carbon = pd.read_csv(os.path.join(TABLES, "T09_paintomics_carbon.tsv"),
                         sep="\t", comment="#")
    photoresp = carbon[carbon["pathway"].str.contains("Glyoxylate")].iloc[0]

    # A curated set, not an arbitrary top-N: the strongest hits overall, plus every
    # pathway the model made a claim about, plus photorespiration itself. Ranking by
    # p alone would omit the photosynthesis pathways entirely (they sit 13th-26th).
    strongest = sig.nsmallest(6, "p_combined_fisher")
    claimed = sig[sig["pathway"].isin(HIGHLIGHT)]
    top = pd.concat([strongest, claimed, photoresp.to_frame().T]).drop_duplicates("pathway")
    top["p_combined_fisher"] = pd.to_numeric(top["p_combined_fisher"])
    top = top.sort_values("p_combined_fisher", ascending=False)
    y = np.arange(len(top))
    vals = -np.log10(top["p_combined_fisher"])
    cols = [style.VERMILION if n in HIGHLIGHT else
            (style.INK if "Glyoxylate" in n else style.GREY) for n in top["pathway"]]
    ax.barh(y, vals, color=cols, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels([n if len(n) < 27 else n[:25] + "…" for n in top["pathway"]],
                       fontsize=6.2)
    ax.set_xlabel("$-$log$_{10}$ combined $p$ (Fisher)")
    ax.axvline(-np.log10(0.05), color=style.INK, linestyle="--", linewidth=0.8)
    ax.text(-np.log10(0.05), len(top) - 0.35, "  $p$ = 0.05", fontsize=6.2,
            va="center", ha="left", color=style.INK)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, float(vals.max()) * 1.55)
    ax.annotate("ranks LAST of 231", xy=(vals.iloc[0], 0), xytext=(2.9, 1.15),
                fontsize=6.5, color=style.INK, va="center",
                arrowprops=dict(arrowstyle="-|>", color=style.INK, linewidth=0.7))
    ax.text(0.98, 0.60, "red = a pathway the\nmodel predicted would fall",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.4,
            color=style.VERMILION)


def panel_c(ax):
    style.panel(ax, "c", "Metabolite hubs")
    hub = pd.read_csv(os.path.join(TABLES, "T10_metabolite_hubs.tsv"), sep="\t",
                      comment="#")
    hub = hub.sort_values("FDR", ascending=False)
    y = np.arange(len(hub))
    vals = -np.log10(hub["FDR"])
    cols = [style.VERMILION if b == "carbon starvation" else
            (style.BLUE if b == "photosynthate" else style.GREEN) for b in hub["block"]]
    ax.barh(y, vals, color=cols, height=0.68)
    ax.set_yticks(y)
    ax.set_yticklabels([c if len(c) < 26 else c[:24] + "…" for c in hub["compound"]],
                       fontsize=6.2)
    ax.axvline(-np.log10(0.05), color=style.INK, linestyle="--", linewidth=0.8)
    ax.set_xlabel("$-$log$_{10}$ FDR")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, float(vals.max()) * 1.75)
    ax.text(0.99, 0.03, "not one of the eight C2\nphotorespiratory intermediates\n"
                        "is a significant hub —\nthough all eight mapped",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=6.4,
            color=style.INK)


def main():
    style.apply()
    fig, axes = plt.subplots(1, 3, figsize=(7.5, 4.1),
                             gridspec_kw={"width_ratios": [0.80, 1.30, 1.00],
                                          "wspace": 1.15})
    panel_a(axes[0]); panel_b(axes[1]); panel_c(axes[2])
    fig.subplots_adjust(top=0.83, bottom=0.24, left=0.085, right=0.985)
    fig.suptitle("PaintOmics: photosynthesis enriched, photorespiration last of 231",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.985)
    style.predicted_note(fig, "Job m1z16Qg3DK, organism ath, KEGG + MapMan, AI interpretation "
                              "off. Hub analysis counts REAL differentially expressed genes.")
    style.save(fig, "fig06_paintomics", FIGURES)


if __name__ == "__main__":
    main()
