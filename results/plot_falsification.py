#!/usr/bin/env python3
"""Figure: the model's prediction against the measured OSD-522 omics.

Panel a  what the model predicted, per gene set, and what the data did
Panel b  the two set-vs-set contrasts, which are immune to the global metabolic shift
Panel c  why the gravity effect is small and the hardware effect is not
"""

from __future__ import annotations

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)

PRED_COLOUR = "#8a8fa3"
OBS_COLOUR = {"down": "#2b6cb0", "up": "#c05621", "no change": "#9aa3bd"}

LABELS = {
    "ath00630": "Glyoxylate & dicarboxylate\n(KEGG ath00630)",
    "ath00710": "Carbon fixation\n(KEGG ath00710)",
    "ath00500": "Starch & sucrose\n(KEGG ath00500)",
    "photorespiration_core": "Photorespiration\nC2 enzymes",
    "carbon_starvation_DIN": "Carbon starvation\n(DIN family)",
    "fermentation": "Fermentation\n(PDC/ADH/LDH)",
    "hypoxia_responsive": "Hypoxia-responsive",
    "photosynthesis_apparatus": "Photosystem &\nlight harvesting",
    "rubisco": "Rubisco",
}
ORDER = ["carbon_starvation_DIN", "photorespiration_core", "fermentation",
         "hypoxia_responsive", "ath00630", "ath00710", "ath00500",
         "photosynthesis_apparatus", "rubisco"]


def main() -> None:
    df = pd.read_csv(os.path.join(HERE, "falsification_check.tsv"), sep="\t")
    tx = df[(df["layer"] == "transcriptome")].set_index("set")

    fig = plt.figure(figsize=(12.5, 5.4), dpi=200)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.35, 1.0, 1.0], wspace=0.42)

    # ---- panel a: per-set observed shift, annotated with the prediction ------------
    ax = fig.add_subplot(gs[0, 0])
    rows = [s for s in ORDER if s in tx.index and np.isfinite(tx.loc[s, "shift"])]
    y = np.arange(len(rows))
    shifts = [tx.loc[s, "shift"] for s in rows]
    colours = [OBS_COLOUR.get(tx.loc[s, "observed"], "#9aa3bd") for s in rows]

    ax.barh(y, shifts, color=colours, height=0.66)
    ax.axvline(0, color="#333", lw=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(s, s) for s in rows], fontsize=7.4)
    ax.invert_yaxis()
    ax.set_xlabel("measured shift vs all other genes (log$_2$FC)", fontsize=8.5)
    ax.set_title("a  Gene-set response, flight vs ground", loc="left",
                 fontweight="bold", fontsize=9.5)

    # Annotations sit in a fixed column to the right of every bar, so they can never
    # collide with a bar, the axis labels, or the neighbouring panel.
    lo, hi = min(shifts), max(shifts)
    span = hi - lo
    ann_x = hi + span * 0.10
    ax.set_xlim(lo - span * 0.10, ann_x + span * 0.72)
    for i, s in enumerate(rows):
        pred, verdict = tx.loc[s, "predicted"], tx.loc[s, "verdict"]
        mark, colour = ("✓", "#2f855a") if verdict == "MATCH" else (
            ("✗", "#c53030") if verdict == "MISMATCH" else ("", "#555"))
        ax.text(ann_x, i, f"{mark} predicted {pred}", va="center", ha="left",
                fontsize=6.5, color=colour)
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    ax.tick_params(labelsize=7.5)

    # ---- panel b: set-vs-set contrasts --------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    cons = df[df["predicted"] == "A above B"]
    labs, vals, cols = [], [], []
    for layer, marker in (("transcriptome", "transcript"), ("proteome", "protein")):
        sub = cons[cons["layer"] == layer]
        for _, r in sub.iterrows():
            if not np.isfinite(r["shift"]):
                continue
            short = r["set"].replace(" vs ", "\nvs ")
            labs.append(f"{short}\n({marker})")
            vals.append(r["shift"])
            cols.append("#2f855a" if r["verdict"] == "SUPPORTS" else "#9aa3bd")
    yy = np.arange(len(labs))
    ax.barh(yy, vals, color=cols, height=0.6)
    ax.axvline(0, color="#333", lw=0.9)
    ax.set_yticks(yy)
    ax.set_yticklabels(labs, fontsize=6.6)
    ax.invert_yaxis()
    ax.set_xlabel("median difference A − B (log$_2$FC)", fontsize=8.5)
    ax.set_title("b  Set-vs-set, model says A above B", loc="left",
                 fontweight="bold", fontsize=9.5)
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    ax.tick_params(labelsize=7.5)

    # ---- panel c: gravity effect vs hardware effect --------------------------------
    ax = fig.add_subplot(gs[0, 2])
    prov = pd.read_csv(os.path.join(REPO, "metabolome", "compound_provenance.tsv"),
                       sep="\t", comment="#")
    grav = prov["log2FC_FLT_vs_GC"].abs().max()
    hard = prov["log2FC_BRIC_vs_VENTED"].abs().max()
    ax.bar(["microgravity\n(flight vs ground)", "hardware\n(sealed vs vented)"],
           [grav, hard], color=["#2b6cb0", "#c05621"], width=0.55)
    ax.set_ylabel("largest predicted |log$_2$FC| on a metabolite pool", fontsize=8.5)
    ax.set_title("c  What actually moves the metabolism", loc="left",
                 fontweight="bold", fontsize=9.5)
    for i, v in enumerate([grav, hard]):
        ax.text(i, v + hard * 0.025, f"{v:.2f}", ha="center", fontsize=8.5,
                fontweight="bold")
    # Sits in the empty space above the small bar, clear of both bars.
    ax.text(0.12, hard * 0.60,
            f"{hard/grav:.0f}× larger —\nbut it hits flight\nand ground alike,\n"
            "so it cancels out\nof the contrast the\nexperiment can see",
            ha="center", va="center", fontsize=6.8, color="#444",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                      edgecolor="#d5d8e0", linewidth=0.6))
    ax.set_ylim(0, hard * 1.18)
    ax.grid(axis="y", linestyle=":", alpha=0.5)
    ax.tick_params(labelsize=7.8)

    fig.suptitle("OSD-522 (BRIC-LED-001): a gas-transport model predicts the spaceflight "
                 "shoot response",
                 fontsize=11, fontweight="bold", x=0.012, ha="left", y=0.985)
    fig.text(0.012, 0.015,
             "Model built from LunarLeaf-CFD gas transport + FvCB photosynthesis, without "
             "reference to the omics. Metabolite layer is predicted, not measured.",
             fontsize=7, color="#555")
    fig.tight_layout(rect=[0, 0.035, 1, 0.945])

    out = os.path.join(HERE, "F1_model_vs_measurement.png")
    fig.savefig(out)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
