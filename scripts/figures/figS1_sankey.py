#!/usr/bin/env python3
"""Supplementary Figure S1 — where the loci go: transcript bin to protein bin to pathway.

The static counterpart of the interactive Sankey in the data explorer, drawn from the same
`export_site_data.assemble_pathways` and `build_sankey` the site uses, so the two cannot
show different membership for the same pathway.

a  the Sankey itself, ribbon width = locus count
b  how far the MapMan reconstruction agrees with PaintOmics' own feature count

Panel b is here because panel a depends on it. PaintOmics publishes enrichment but not the
gene lists behind it, and MapMan has no equivalent of the KEGG REST API, so the MapMan half
of column 3 is reconstructed. Showing the Sankey without showing how well that
reconstruction lands would be asking the reader to take it on trust.

Database is marked the same way as Figure 6b — a square in the identifier-mapping colours —
rather than by node colour, because colour here already means DEG direction.
"""
from __future__ import annotations
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patheffects as pe
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import style  # noqa: E402
import export_site_data as ex  # noqa: E402
from paths import FIGURES  # noqa: E402

TOP_N = 10
BIN_COLOUR = {"up": style.VERMILION, "down": style.BLUE,
              "not significant": style.GREY, "not measured": style.GREY_LIGHT}
DB_COLOUR = {"KEGG": style.BLUE, "MapMan": style.GREEN}
NODE_W = 0.022
GAP = 0.018


def ribbon(ax, x0, y0, h0, x1, y1, h1, colour):
    cx = (x0 + x1) / 2
    verts = [(x0, y0), (cx, y0), (cx, y1), (x1, y1),
             (x1, y1 + h1), (cx, y1 + h1), (cx, y0 + h0), (x0, y0 + h0), (x0, y0)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=colour, edgecolor="none",
                           alpha=0.38, zorder=1))


def lay_out(nodes, floor=0.006):
    """Node heights proportional to share, on top of a floor, fitting the column exactly.

    The floor is taken out of the free space before it is shared, not after; taken after,
    the column overruns its own box — which it did, in the browser, until it was measured.
    """
    total = sum(n["n"] for n in nodes) or 1
    free = max(0.0, 1.0 - GAP * (len(nodes) - 1) - floor * len(nodes))
    y, out = 0.0, []
    for n in nodes:
        h = floor + n["n"] / total * free
        out.append({**n, "y": y, "h": h, "in_used": 0.0, "out_used": 0.0})
        y += h + GAP
    return out


def panel_a(ax, sk):
    style.panel(ax, "a", f"Transcript bin to protein bin to pathway "
                         f"(top {TOP_N} of 28)")
    view = sk["views"]["core"]
    paths = [p for p in sk["pathways"] if not p["broad"]][:TOP_N]
    keep = {p["name"] for p in paths}
    meta = {p["name"]: p for p in paths}

    p2p = [l for l in sk["pr_to_path"] if l["to"] in keep]
    totals = {n: sum(l["n"] for l in p2p if l["to"] == n) for n in keep}

    cols = [
        lay_out([{"k": b, "label": b, "n": view["tx_bins"].get(b, 0)}
                 for b in sk["bins"] if view["tx_bins"].get(b, 0)]),
        lay_out([{"k": b, "label": b, "n": view["pr_bins"].get(b, 0)}
                 for b in sk["bins"] if view["pr_bins"].get(b, 0)]),
        lay_out([{"k": p["name"], "label": p["name"], "n": totals[p["name"]]}
                 for p in paths if totals[p["name"]]]),
    ]
    x = [0.0, 0.42, 0.84]
    find = lambda ci, k: next(n for n in cols[ci] if n["k"] == k)  # noqa: E731

    for links, ci in ((view["tx_to_pr"], 0), (p2p, 1)):
        for l in sorted(links, key=lambda d: -d["n"]):
            a, b = find(ci, l["from"]), find(ci + 1, l["to"])
            ha, hb = l["n"] / a["n"] * a["h"], l["n"] / b["n"] * b["h"]
            ribbon(ax, x[ci] + NODE_W, a["y"] + a["out_used"], ha,
                   x[ci + 1], b["y"] + b["in_used"], hb, BIN_COLOUR[l["from"]])
            a["out_used"] += ha
            b["in_used"] += hb

    for ci, nodes in enumerate(cols):
        for n in nodes:
            is_path = ci == 2
            ax.add_patch(Rectangle((x[ci], n["y"]), NODE_W, n["h"], zorder=3,
                                   facecolor=style.GREY if is_path
                                   else BIN_COLOUR[n["k"]], edgecolor="none"))
            label = n["label"] if len(n["label"]) < 30 else n["label"][:28] + "…"
            if is_path:
                ax.plot(x[ci] + NODE_W + 0.016, n["y"] + n["h"] / 2, marker="s",
                        markersize=2.6, color=DB_COLOUR[meta[n["k"]]["db"]], zorder=4)
                ax.text(x[ci] + NODE_W + 0.034, n["y"] + n["h"] / 2,
                        f"{label}  ({n['n']:,})", va="center", ha="left", fontsize=6.2)
            else:
                # Column 2 has ribbons on both sides, so its labels necessarily sit over
                # flow. A halo keeps them readable without moving them off their node.
                t = ax.text(x[ci] - 0.012, n["y"] + n["h"] / 2, f"{label}  ({n['n']:,})",
                            va="center", ha="right", fontsize=6.2, zorder=5)
                if ci == 1:
                    t.set_path_effects([pe.withStroke(linewidth=2.2, foreground="white")])

    for xi, title, ha in ((x[0] + NODE_W / 2, "Transcriptome", "center"),
                          (x[1] + NODE_W / 2, "Proteome", "center"),
                          (x[2], "Pathway", "left")):
        ax.text(xi, -0.035, title, ha=ha, va="bottom", fontsize=7, fontweight="bold",
                color=style.INK)

    ax.set_xlim(-0.30, 1.30)
    ax.set_ylim(1.02, -0.075)
    ax.axis("off")


def panel_b(ax, sk):
    # style.panel's 7pt pad would put the letter on top of the two-line subtitle below, so
    # this panel sets the same title with room for it. Note `ax.get_title()` reads the
    # CENTRE title and returns "" for a loc="left" one — reusing it here silently blanked
    # the panel letter.
    ax.set_title("b  MapMan membership: reconstruction vs PaintOmics", loc="left",
                 fontweight="bold", fontsize=9.5, pad=26)
    mm = [p for p in sk["pathways"] if p["db"] == "MapMan"]
    mm = sorted(mm, key=lambda p: p["measured"])
    y = np.arange(len(mm))
    ours = np.array([p["measured"] for p in mm], dtype=float)
    theirs = np.array([p["features"] for p in mm], dtype=float)
    agree = np.abs(ours - theirs) <= 3

    h = 0.36
    ax.barh(y + h / 2, ours, height=h, color=style.GREEN, label="reconstructed here")
    ax.barh(y - h / 2, theirs, height=h, color=style.GREY,
            label="PaintOmics' own feature count")
    ax.set_xscale("symlog", linthresh=10)
    ax.set_yticks(y)
    ax.set_yticklabels([p["name"] if len(p["name"]) < 30 else p["name"][:28] + "…"
                        for p in mm], fontsize=6.2)
    ax.set_xlabel("loci measured here and placed on the diagram")
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, float(max(ours.max(), theirs.max())) * 2.6)
    # Only the disagreements are annotated, with how far apart they are. Labelling the
    # nine that match crowded the panel and said nothing the bars did not already show.
    for yi, a, o, t, name in zip(y, agree, ours, theirs, [p["name"] for p in mm]):
        if not a:
            ax.text(max(o, t) * 1.25, yi, f"{o / t:.1f}×", va="center", fontsize=6.2,
                    color=style.VERMILION, fontweight="bold")
    for tick, a in zip(ax.get_yticklabels(), agree):
        if not a:
            tick.set_color(style.VERMILION)
    ax.legend(loc="lower right", fontsize=6.2, frameon=False)
    # As a subtitle, NOT set_title — that would overwrite the panel letter style.panel set.
    ax.text(0.0, 1.012, f"{int(agree.sum())} of {len(mm)} agree to within three features; "
                        f"the {int((~agree).sum())} that do not are the largest maps, and "
                        f"the reason is\nnot recoverable from PaintOmics' published sources",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=6.4, color=style.INK)


def main():
    style.apply()
    P = ex.assemble_pathways(verbose=False)
    sk = ex.build_sankey(P["tx"], P["pr"], P["pr_agi"], P["pathways"], P["pathway_db"],
                         P["broad"], P["sig_paths"], P["unresolved"])

    fig, axes = plt.subplots(2, 1, figsize=(7.5, 8.6),
                             gridspec_kw={"height_ratios": [1.5, 1.0], "hspace": 0.22})
    panel_a(axes[0], sk)
    panel_b(axes[1], sk)
    fig.subplots_adjust(top=0.90, bottom=0.075, left=0.215, right=0.965)
    fig.suptitle("Supplementary Figure S1 — transcript, protein and pathway together",
                 fontsize=11, fontweight="bold", x=0.008, ha="left", y=0.985)
    style.predicted_note(
        fig, f"Ribbon width is a locus count. Columns 1 and 2 conserve loci exactly; "
             f"column 3 counts memberships. Three whole-ontology MapMan maps are excluded "
             f"(each covers more than {sk['broad_map_share']:.0%} of measured loci). "
             f"Blue square = KEGG, green = MapMan, as in Fig. 6a.")
    style.save(fig, "figS1_sankey", os.path.join(FIGURES, "supplementary"))


if __name__ == "__main__":
    main()
