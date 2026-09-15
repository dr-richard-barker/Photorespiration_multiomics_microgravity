#!/usr/bin/env python3
"""Export compact JSON for the interactive site.

The site is static — no server, no build step — so every dataset it explores ships as JSON
under docs/data/. Volcano data is the only bulky part; it is thinned to what a scatter plot
can actually resolve (rounded to 3 decimals, non-significant points subsampled) so the page
stays quick on a phone.
"""
from __future__ import annotations
import json, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genesets  # noqa: E402
from paths import DATA, DOCS_DATA, STUDY_REGISTRY, TABLES, ensure  # noqa: E402


def volcano(df, lfc, p, idname, alpha, keep_ns=4000, seed=0, must_keep=frozenset()):
    """Thin the cloud for speed, but never drop a feature the page can highlight.

    The thinning is what broke gene-set highlighting: members that were not significant
    fell outside the random sample and so were not in the DOM to highlight. Worse, the
    sets whose *absence of change* is the finding — photorespiration_core, hypoxia — are
    entirely non-significant, so they lost every member. `must_keep` forces them back in.
    """
    d = df.dropna(subset=[lfc, p]).copy()
    sig = d[p] < alpha
    forced = d.index.isin(must_keep) & ~sig
    ns = d[~sig & ~forced]
    if len(ns) > keep_ns:
        ns = ns.sample(keep_ns, random_state=seed)
    out = pd.concat([d[sig], d[forced], ns])
    return {
        "id": idname,
        "alpha": alpha,
        "n_total": int(len(d)),
        "n_sig": int(sig.sum()),
        "n_plotted": int(len(out)),
        "points": [[str(i), round(float(r[lfc]), 3),
                    round(float(-np.log10(max(r[p], 1e-300))), 3),
                    int(r[p] < alpha)]
                   for i, r in out.iterrows()],
    }


BINS = ["up", "down", "not significant", "not measured"]


def _bin(gene, table, pcol, fcol, alpha=0.05):
    if gene not in table.index:
        return "not measured"
    row = table.loc[gene]
    if pd.isna(row[pcol]):
        return "not significant"
    if row[pcol] < alpha:
        return "up" if row[fcol] > 0 else "down"
    return "not significant"


def build_sankey(tx, pr, pr_agi, pathways, sig_paths, unresolved):
    """Transcript DEG bin -> protein DEG bin -> pathway, sized by locus count.

    Columns 1 and 2 conserve genes exactly: every gene sits in one transcript bin and one
    protein bin. Column 3 counts MEMBERSHIPS, not genes, because a locus can belong to
    several pathways — the page says so, and the counts are integers rather than
    fractional splits because locus counts are what was asked for.
    """
    # Collapse the proteome to one row per gene: mean fold change, strongest adjusted p.
    pr_g = pr.assign(agi=[pr_agi.get(u) for u in pr.index]).dropna(subset=["agi"])
    pr_g = pr_g.groupby("agi").agg(log2fc=("log2fc", "mean"), adj_p=("adj_p", "min"))

    genes = sorted(set().union(*pathways.values())) if pathways else []
    rec = {g: (_bin(g, tx, "fdr", "log2fc"), _bin(g, pr_g, "adj_p", "log2fc"))
           for g in genes}

    tx_pr = {}
    for t, p_ in rec.values():
        tx_pr[(t, p_)] = tx_pr.get((t, p_), 0) + 1

    pr_path = {}
    for name, members in pathways.items():
        for g in members:
            if g not in rec:
                continue
            key = (rec[g][1], name)
            pr_path[key] = pr_path.get(key, 0) + 1

    order = {n: i for i, n in enumerate(sig_paths.pathway)}
    ranked = sorted(pathways, key=lambda n: order.get(n, 999))

    return {
        "bins": BINS,
        "pathways": [{"name": n, "size": len(pathways[n]),
                      "p": float(sig_paths.set_index("pathway")
                                 .loc[n, "p_combined_fisher"])}
                     for n in ranked],
        "tx_bins": {b: sum(1 for t, _ in rec.values() if t == b) for b in BINS},
        "pr_bins": {b: sum(1 for _, p_ in rec.values() if p_ == b) for b in BINS},
        "tx_to_pr": [{"from": t, "to": p_, "n": n} for (t, p_), n in sorted(tx_pr.items())],
        "pr_to_path": [{"from": p_, "to": n, "n": c}
                       for (p_, n), c in sorted(pr_path.items())],
        "n_genes": len(rec),
        "n_memberships": sum(len(m) for m in pathways.values()),
        "unresolved": sorted(unresolved),
        "n_pathways_total": int(len(sig_paths)),
    }


def main() -> int:
    ensure(DOCS_DATA)
    written = []

    def dump(name, obj):
        path = os.path.join(DOCS_DATA, name)
        with open(path, "w") as fh:
            json.dump(obj, fh, separators=(",", ":"))
        written.append((name, os.path.getsize(path)))

    # --- the CFD sweep the model explorer offers as presets ---------------------------
    cfd = pd.read_csv(os.path.join(DATA, "lunarleaf_gbl_sweep.csv"))
    dump("cfd_sweep.json", cfd.to_dict(orient="records"))

    # --- enclosure behaviour ----------------------------------------------------------
    ll = os.path.join(DATA, "lunarleaf")
    t7 = pd.read_csv(os.path.join(ll, "T7_hardware_timeseries.csv"))
    # VEGGIE and the open reference carry NaN in dishmean_co2 (they vent, so there is no
    # enclosure mean to report). json.dump would emit a bare NaN, which JSON.parse rejects,
    # so they become null and the page skips them.
    def clean(vals):
        return [None if pd.isna(v) else round(float(v), 5) for v in vals]
    dump("enclosure_timeseries.json",
         {case: {"step": g["step"].tolist(), "co2": clean(g["dishmean_co2"])}
          for case, g in t7.groupby("case")})
    dump("carbon_retention.json",
         pd.read_csv(os.path.join(ll, "T11_photosynthesis_feedback.csv"))
         .to_dict(orient="records"))

    # --- gene sets and pathways, per layer -------------------------------------------
    # Both membership indexes are precomputed per identifier space so the browser never
    # needs a UniProt->AGI map. hypoxia_responsive has 0 measured proteins: that is real
    # biology, and the page reports "0 of 3 measured" rather than silently doing nothing.
    sets = genesets.build()
    tx = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"),
                     sep="\t", index_col=0)
    pr = pd.read_csv(os.path.join(TABLES, "T02_osd522_proteome.tsv"), sep="\t", index_col=0)
    conv = genesets.uniprot_to_agi()
    pr_agi = {u: conv.get(u) for u in pr.index}

    sig_paths = pd.read_csv(os.path.join(TABLES, "T08_paintomics_significant.tsv"),
                            sep="\t", comment="#")
    resolved, unresolved = genesets.paintomics_pathway_ids(sig_paths.pathway, sig_paths.db)
    pathways = {name: genesets.pathway_genes(pid) for name, pid in resolved.items()}
    print(f"  pathways resolved to KEGG: {len(resolved)}/{len(sig_paths)} "
          f"({len(unresolved)} MapMan bins have no offline membership)")

    # A volcano needs a y-value, so a feature with no adjusted p cannot be plotted at all.
    # DESeq2's independent filtering leaves 1,658 transcripts with a fold change but no FDR.
    # The index therefore lists only PLOTTABLE features, and carries the measured counts
    # alongside so the page can say "2 of 3 measured" instead of quietly under-delivering.
    tx_plottable = set(tx.index[tx["fdr"].notna()])
    pr_plottable = set(pr.index[pr["adj_p"].notna()])
    tx_measured, pr_measured = set(tx.index), set(pr.index)

    def per_layer(members):
        agi = sorted(members & tx_plottable)
        upr_all = {u for u, a in pr_agi.items() if a in members and u in pr_measured}
        upr = sorted(upr_all & pr_plottable)
        return {
            "transcriptome": agi,
            "proteome": upr,
            "total": len(members),
            "measured": {"transcriptome": len(members & tx_measured),
                         "proteome": len(upr_all)},
        }

    dump("gene_sets.json", {n: per_layer(m) for n, m in sets.items()})
    dump("pathway_membership.json", {n: per_layer(m) for n, m in pathways.items()})
    dump("pathway_coverage.json",
         {"resolved": len(resolved), "total": int(len(sig_paths)),
          "unresolved": sorted(unresolved)})

    # Every gene the page can highlight must survive the volcano thinning.
    must_keep = set().union(*sets.values()) | (set().union(*pathways.values())
                                               if pathways else set())

    # --- measured omics, thinned for the browser but never below must_keep ------------
    dump("volcano_transcriptome.json",
         volcano(tx, "log2fc", "fdr", "AGI", 0.05, must_keep=must_keep))
    prot_keep = {u for u, a in pr_agi.items() if a in must_keep}
    dump("volcano_proteome.json",
         volcano(pr, "log2fc", "adj_p", "UniProt", 0.05, must_keep=prot_keep))

    # --- sankey: transcript bin -> protein bin -> pathway -----------------------------
    dump("sankey.json", build_sankey(tx, pr, pr_agi, pathways, sig_paths, unresolved))

    # --- results tables the site tabulates --------------------------------------------
    for name, path, kw in (
        ("enrichment.json", os.path.join(TABLES, "T08_paintomics_significant.tsv"),
         dict(sep="\t", comment="#")),
        ("carbon_pathways.json", os.path.join(TABLES, "T09_paintomics_carbon.tsv"),
         dict(sep="\t", comment="#")),
        ("metabolite_hubs.json", os.path.join(TABLES, "T10_metabolite_hubs.tsv"),
         dict(sep="\t", comment="#")),
        ("predicted_compounds.json", os.path.join(TABLES, "T06_predicted_compounds.tsv"),
         dict(sep="\t", comment="#")),
        ("ladder.json", os.path.join(TABLES, "T11_hardware_ladder.tsv"), dict(sep="\t")),
        ("ladder_contrasts.json", os.path.join(TABLES, "T12_ladder_contrasts.tsv"),
         dict(sep="\t")),
        ("falsification.json", os.path.join(TABLES, "T07_falsification_osd522.tsv"),
         dict(sep="\t")),
    ):
        df = pd.read_csv(path, **kw)
        dump(name, json.loads(df.to_json(orient="records")))

    # --- the study registry, including the exclusions and their reasons ---------------
    reg = [ln.rstrip("\n").split("\t") for ln in open(STUDY_REGISTRY)
           if not ln.startswith("#") and ln.strip()]
    dump("studies.json", [dict(zip(reg[0], r)) for r in reg[1:]])

    total = sum(s for _n, s in written)
    for name, size in written:
        print(f"  {name:34s} {size:>9,} bytes")
    print(f"  {'TOTAL':34s} {total:>9,} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
