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

    # Sankey conservation: columns 1 and 2 are genes and must balance.
    s = load("sankey.json")
    for side, key, bins in (("tx_to_pr out of", "from", s["tx_bins"]),
                            ("tx_to_pr into", "to", s["pr_bins"])):
        for b, total in bins.items():
            flowed = sum(l["n"] for l in s["tx_to_pr"] if l[key] == b)
            if flowed != total:
                failures.append(f"sankey {side} '{b}': ribbons sum to {flowed}, node is {total}")
    print(f"  sankey: {s['n_genes']} genes conserved across transcript and protein columns")
    print(f"  sankey: {s['n_memberships']} memberships across {len(s['pathways'])} pathways "
          f"({s['n_memberships'] - s['n_genes']} loci in more than one pathway)")

    cov = load("pathway_coverage.json")
    print(f"  pathway coverage: {cov['resolved']}/{cov['total']} significant pathways have "
          f"KEGG membership; {len(cov['unresolved'])} MapMan bins omitted")

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
