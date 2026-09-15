#!/usr/bin/env python3
"""Assert the site's data can actually do what the site's controls promise.

Written after a real bug: selecting a gene set in the omics dashboard highlighted nothing,
because the volcano export thinned out non-significant points and the sets are largely
non-significant. `photorespiration_core` and `hypoxia_responsive` lost every member — the
two sets whose *absence of change* is the paper's central claim.

A control that silently does nothing is worse than no control, so this fails loudly.
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DOCS_DATA  # noqa: E402


def load(name):
    with open(os.path.join(DOCS_DATA, name)) as fh:
        return json.load(fh)


def main() -> int:
    failures: list[str] = []

    volcanoes = {
        "transcriptome": load("volcano_transcriptome.json"),
        "proteome": load("volcano_proteome.json"),
    }
    plotted = {layer: {p[0] for p in v["points"]} for layer, v in volcanoes.items()}

    print("Site data check — every highlightable feature must be plotted\n")
    for index_name in ("gene_sets.json", "pathway_membership.json"):
        index = load(index_name)
        print(f"  {index_name}  ({len(index)} entries)")
        worst = 0
        for name, layers in index.items():
            for layer in ("transcriptome", "proteome"):
                ids = set(layers[layer])
                missing = ids - plotted[layer]
                if missing:
                    failures.append(
                        f"{index_name}:{name}[{layer}] — {len(missing)} of {len(ids)} "
                        f"measured features are not in the volcano export "
                        f"(e.g. {sorted(missing)[:3]})")
                worst = max(worst, len(missing))
        print(f"    max missing across all entries and layers: {worst}")

    # The two sets that were entirely invisible before the fix. Both were 0; the expected
    # values are PLOTTABLE counts, which for hypoxia_responsive is 2 of 3 measured because
    # one member has no adjusted p after DESeq2 independent filtering.
    sets = load("gene_sets.json")
    for name, expect_tx in (("photorespiration_core", 13), ("hypoxia_responsive", 2)):
        got = len(sets[name]["transcriptome"])
        measured = sets[name]["measured"]["transcriptome"]
        state = "ok  " if got == expect_tx else "FAIL"
        note = "" if got == measured else f"; {measured - got} measured but not testable"
        print(f"  {state} {name}: {got} of {measured} measured transcripts plottable"
              f"{note} (was 0 before the fix)")
        if got != expect_tx:
            failures.append(f"{name}: expected {expect_tx} plottable transcripts, got {got}")

    # Sankey conservation: columns 1 and 2 are genes and must balance, in BOTH views —
    # the checkbox switches between them, so a view that does not balance is a chart that
    # silently loses loci when the reader ticks a box.
    s = load("sankey.json")
    for name, v in s["views"].items():
        for side, key, bins in (("tx_to_pr out of", "from", v["tx_bins"]),
                                ("tx_to_pr into", "to", v["pr_bins"])):
            for b, total in bins.items():
                flowed = sum(l["n"] for l in v["tx_to_pr"] if l[key] == b)
                if flowed != total:
                    failures.append(f"sankey[{name}] {side} '{b}': ribbons sum to "
                                    f"{flowed}, node is {total}")
        for bins in ("tx_bins", "pr_bins"):
            if sum(v[bins].values()) != v["n_genes"]:
                failures.append(f"sankey[{name}] {bins} sums to {sum(v[bins].values())}, "
                                f"not the {v['n_genes']} genes in the view")
        print(f"  sankey[{name}]: {v['n_genes']:,} loci conserved across the transcript and "
              f"protein columns; {v['n_memberships']:,} memberships "
              f"({v['n_memberships'] - v['n_genes']:,} in more than one pathway)")

    # Every pathway the Sankey can draw must be sizeable, and every one it can highlight
    # must be in the membership index the volcano reads.
    memberships = load("pathway_membership.json")
    for p in s["pathways"]:
        if not p["broad"] and p["name"] not in memberships:
            failures.append(f"sankey pathway '{p['name']}' is clickable but has no entry "
                            f"in pathway_membership.json")
        if p["broad"] and p["name"] in memberships:
            failures.append(f"broad map '{p['name']}' should be out of the highlight index")
    kegg = sum(1 for p in s["pathways"] if p["db"] == "KEGG")
    mapman = len(s["pathways"]) - kegg
    agree = sum(1 for p in s["pathways"]
                if p["db"] == "MapMan" and abs(p["measured"] - p["features"]) <= 3)
    print(f"  sankey: {len(s['pathways'])} of {s['n_pathways_total']} significant pathways "
          f"have membership ({kegg} KEGG, {mapman} MapMan); MapMan reconstruction agrees "
          f"with PaintOmics' feature count for {agree} of {mapman}")

    cov = load("pathway_coverage.json")
    if cov["resolved"] != cov["total"]:
        failures.append(f"pathway coverage is {cov['resolved']}/{cov['total']}, not complete")
    print(f"  pathway coverage: {cov['resolved']}/{cov['total']} significant pathways have "
          f"membership ({cov['kegg']} KEGG via REST, {cov['mapman']} MapMan via vendored "
          f"diagrams); {len(cov['broad'])} broad maps held out of the highlight index")

    print()
    if failures:
        print(f"FAILED — {len(failures)} problem(s):")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
