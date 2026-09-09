#!/usr/bin/env python3
"""Figure 4 — what OSD-522 actually measured.

a  sample PCA of the 12 RNA-seq libraries
b  transcriptome volcano
c  proteome volcano
d  transcript against protein, for the 4,414 genes measured on both layers
"""
from __future__ import annotations
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genesets, style  # noqa: E402
from paths import CACHE, FIGURES, TABLES  # noqa: E402

COUNTS = os.path.join(CACHE, "GLDS-522_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv")
ALPHA = 0.05


def panel_a(ax):
    style.panel(ax, "a", "Sample structure (RNA-seq)")
    counts = pd.read_csv(COUNTS, index_col=0)
    counts = counts.loc[counts.sum(axis=1) >= 10]
    lcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1).T
    lcpm = lcpm.loc[:, lcpm.var(axis=0).sort_values(ascending=False).index[:2000]]
    pcs = PCA(n_components=2).fit(lcpm - lcpm.mean(axis=0))
    xy = pcs.transform(lcpm - lcpm.mean(axis=0))
    flight = [s.startswith("FT") for s in lcpm.index]
    for mask, colour, label in ((flight, style.VERMILION, "Space Flight"),
                                ([not f for f in flight], style.GREY, "Ground Control")):
        ax.scatter(xy[mask, 0], xy[mask, 1], s=52, color=colour, label=label,
                   edgecolor="white", linewidth=0.8)
    ax.set_xlabel(f"PC1 ({pcs.explained_variance_ratio_[0]*100:.0f}%)")
    ax.set_ylabel(f"PC2 ({pcs.explained_variance_ratio_[1]*100:.0f}%)")
    ax.legend(loc="best")


def volcano(ax, letter, title, df, lfc, p, n_label):
    style.panel(ax, letter, title)
    d = df.dropna(subset=[lfc, p]).copy()
    d["y"] = -np.log10(d[p].clip(lower=1e-300))
    sig = d[p] < ALPHA
    ax.scatter(d.loc[~sig, lfc], d.loc[~sig, "y"], s=3, color=style.GREY_LIGHT,
               linewidths=0, rasterized=True)
    ax.scatter(d.loc[sig & (d[lfc] > 0), lfc], d.loc[sig & (d[lfc] > 0), "y"], s=4,
               color=style.VERMILION, linewidths=0, rasterized=True)
    ax.scatter(d.loc[sig & (d[lfc] < 0), lfc], d.loc[sig & (d[lfc] < 0), "y"], s=4,
               color=style.BLUE, linewidths=0, rasterized=True)
    ax.axhline(-np.log10(ALPHA), color=style.INK, linewidth=0.7, linestyle="--")
    ax.set_xlabel("log$_2$ fold change (flight / ground)")
    ax.set_ylabel(f"$-$log$_{{10}}$ {n_label}")
    ax.text(0.98, 0.97, f"{int(sig.sum()):,} at FDR < {ALPHA}\n"
                        f"{int((sig & (d[lfc] > 0)).sum()):,} up   "
                        f"{int((sig & (d[lfc] < 0)).sum()):,} down",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5)
    return d


def panel_d(ax):
    style.panel(ax, "d", "Transcript vs protein")
    tx = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"),
                     sep="\t", index_col=0)
    pr = pd.read_csv(os.path.join(TABLES, "T02_osd522_proteome.tsv"), sep="\t", index_col=0)
    conv = genesets.uniprot_to_agi()
    pr = pr.assign(agi=[conv.get(a) for a in pr.index]).dropna(subset=["agi"])
    pr = pr.groupby("agi")["log2fc"].mean()
    both = pd.DataFrame({"tx": tx["log2fc"], "pr": pr}).dropna()
    ax.scatter(both["tx"], both["pr"], s=4, color=style.GREY, alpha=0.45,
               linewidths=0, rasterized=True)
    r, p = stats.spearmanr(both["tx"], both["pr"])
    lim = 4.2
    ax.plot([-lim, lim], [-lim, lim], color=style.GREY_LIGHT, linewidth=0.8, zorder=0)
    ax.axhline(0, color=style.GREY_LIGHT, linewidth=0.7)
    ax.axvline(0, color=style.GREY_LIGHT, linewidth=0.7)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("transcript log$_2$FC")
    ax.set_ylabel("protein log$_2$FC")
    ax.text(0.03, 0.97, f"{len(both):,} genes on both layers\nSpearman $\\rho$ = {r:+.2f}"
                        f"  ($p$ = {p:.1e})",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.5)


def main():
    style.apply()
    fig, axes = plt.subplots(2, 2, figsize=(7.3, 5.6))
    panel_a(axes[0, 0])
    tx = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"),
                     sep="\t", index_col=0)
    volcano(axes[0, 1], "b", "Transcriptome", tx, "log2fc", "fdr", "FDR")
    pr = pd.read_csv(os.path.join(TABLES, "T02_osd522_proteome.tsv"), sep="\t", index_col=0)
    volcano(axes[1, 0], "c", "Proteome", pr, "log2fc", "adj_p", "adjusted $p$")
    panel_d(axes[1, 1])
    fig.subplots_adjust(hspace=0.46, wspace=0.36, top=0.90, bottom=0.10)
    fig.suptitle("OSD-522: Arabidopsis shoots, 6 flight vs 6 ground, in BRIC-LED hardware",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.985)
    style.predicted_note(fig, "Transcriptome computed here with PyDESeq2; proteome ratios as "
                              "deposited by the original study. Both measured.")
    style.save(fig, "fig04_osd522_omics", FIGURES)


if __name__ == "__main__":
    main()
