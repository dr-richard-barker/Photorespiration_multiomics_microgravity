#!/usr/bin/env python3
"""Figure 5 — the model's blind prediction against the OSD-522 measurement.

a  gene-set responses, with what the model predicted before seeing any data
b  set-vs-set contrasts, which are immune to the global metabolic depression

The model was built from gas transport and photosynthesis alone. It never saw the omics.
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

ORDER = ["carbon_starvation_DIN", "photorespiration_core", "fermentation",
         "hypoxia_responsive", "ath00630", "ath00710", "ath00500",
         "photosynthesis_apparatus", "rubisco"]
OBS_COLOUR = {"down": style.BLUE, "up": style.VERMILION, "no change": style.GREY}


def panel_a(ax, tx):
    style.panel(ax, "a", "Gene-set response, flight vs ground")
    rows = [s for s in ORDER if s in tx.index and np.isfinite(tx.loc[s, "shift"])]
    y = np.arange(len(rows))
    shifts = [tx.loc[s, "shift"] for s in rows]
    ax.barh(y, shifts, color=[OBS_COLOUR.get(tx.loc[s, "observed"], style.GREY)
                              for s in rows], height=0.66)
    ax.axvline(0, color=style.INK, linewidth=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([style.PRETTY.get(s, s) for s in rows], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("measured shift vs all other genes (log$_2$FC)")
    lo, hi = min(shifts), max(shifts)
    span = hi - lo
    ann_x = hi + span * 0.10
    ax.set_xlim(lo - span * 0.10, ann_x + span * 0.95)
    for i, s in enumerate(rows):
        pred, verdict = tx.loc[s, "predicted"], tx.loc[s, "verdict"]
        # Unicode ticks and crosses are not in the sans stack matplotlib resolves here
        # and render as tofu, so the verdict is carried by word and colour instead.
        mark, colour = ("agrees:", style.GREEN) if verdict == "MATCH" else (
            ("DIFFERS:", style.VERMILION) if verdict == "MISMATCH" else ("", style.GREY))
        ax.text(ann_x, i, f"{mark} predicted {pred}", va="center", ha="left",
                fontsize=6.5, color=colour)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)


def panel_b(ax, df):
    style.panel(ax, "b", "Set vs set — the sharper test")
    cons = df[df["predicted"] == "A above B"].dropna(subset=["shift"])
    labs, vals, cols = [], [], []
    for _, r in cons.iterrows():
        labs.append(f"{r['set'].replace(' vs ', chr(10) + 'vs ')}\n({r['layer']})")
        vals.append(r["shift"])
        cols.append(style.GREEN if r["verdict"] == "SUPPORTS" else style.GREY)
    y = np.arange(len(labs))
    ax.barh(y, vals, color=cols, height=0.6)
    ax.axvline(0, color=style.INK, linewidth=0.9)
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=6.4)
    ax.invert_yaxis()
    ax.set_xlabel("median difference A $-$ B (log$_2$FC)")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.text(0.98, 0.02, "green = supports the model", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.4, color=style.GREEN)


def main():
    style.apply()
    df = pd.read_csv(os.path.join(TABLES, "T07_falsification_osd522.tsv"), sep="\t")
    tx = df[df["layer"] == "transcriptome"].set_index("set")
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6),
                             gridspec_kw={"width_ratios": [1.35, 1.0], "wspace": 0.55})
    panel_a(axes[0], tx)
    panel_b(axes[1], df)
    fig.subplots_adjust(top=0.84, bottom=0.20, left=0.155, right=0.985)
    fig.suptitle("A model that never saw the data called eight of nine gene sets",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.985)
    style.predicted_note(fig, "Predictions fixed before any omics were downloaded; "
                              "gene sets from KEGG pathway membership and descriptions.")
    style.save(fig, "fig05_blind_test", FIGURES)


if __name__ == "__main__":
    main()
