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
import mapman  # noqa: E402
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


# Three MapMan diagrams are whole-ontology maps: `Overview` lists all 36 top-level bins and
# so covers 99.8 % of the measured loci, and `Regulation overview` and `Biotic Stress` cover
# 45 % and 44 %. The next largest is 17 %, so the break is a real one rather than a chosen
# line. Such a map cannot say where anything goes — its Sankey ribbons are the marginal
# distribution redrawn, and highlighting it on the volcano marks half the plot. They are
# therefore off by default, in both places, with the reason stated on the page and a
# checkbox on the Sankey. This is a display rule, not a filter on the result: the nodes
# ship with their real counts either way.
BROAD_MAP_SHARE = 0.25


def build_sankey(tx, pr, pr_agi, pathways, pathway_db, broad, sig_paths, unresolved):
    """Transcript DEG bin -> protein DEG bin -> pathway, sized by locus count.

    Columns 1 and 2 conserve genes exactly: every gene sits in one transcript bin and one
    protein bin. Column 3 counts MEMBERSHIPS, not genes, because a locus can belong to
    several pathways — the page says so, and the counts are integers rather than
    fractional splits because locus counts are what was asked for.

    Pathway nodes carry their database (KEGG membership from the REST API, MapMan from the
    vendored diagrams) and PaintOmics' own feature count, so a reader can see where the
    reconstruction and the server disagree instead of being handed one number.
    """
    # Collapse the proteome to one row per gene: mean fold change, strongest adjusted p.
    pr_g = pr.assign(agi=[pr_agi.get(u) for u in pr.index]).dropna(subset=["agi"])
    pr_g = pr_g.groupby("agi").agg(log2fc=("log2fc", "mean"), adj_p=("adj_p", "min"))

    rec = {g: (_bin(g, tx, "fdr", "log2fc"), _bin(g, pr_g, "adj_p", "log2fc"))
           for g in sorted(set().union(*pathways.values()) if pathways else set())}

    def view(names):
        """Columns 1 and 2 for one set of pathways. Genes are counted once, however many
        of the pathways they belong to, so the two columns balance."""
        genes = set().union(*(pathways[n] for n in names)) if names else set()
        tx_pr = {}
        for g in genes:
            tx_pr[rec[g]] = tx_pr.get(rec[g], 0) + 1
        return {
            "tx_bins": {b: sum(1 for g in genes if rec[g][0] == b) for b in BINS},
            "pr_bins": {b: sum(1 for g in genes if rec[g][1] == b) for b in BINS},
            "tx_to_pr": [{"from": t, "to": p_, "n": n}
                         for (t, p_), n in sorted(tx_pr.items())],
            "n_genes": len(genes),
            "n_memberships": sum(len(pathways[n]) for n in names),
        }

    pr_path = {}
    for name, members in pathways.items():
        for g in members:
            key = (rec[g][1], name)
            pr_path[key] = pr_path.get(key, 0) + 1

    order = {n: i for i, n in enumerate(sig_paths.pathway)}
    ranked = sorted(pathways, key=lambda n: order.get(n, 999))
    t8 = sig_paths.set_index("pathway")
    measured = set(tx.index) | {a for a in pr_agi.values() if a}

    return {
        "bins": BINS,
        "pathways": [{"name": n, "size": len(pathways[n]),
                      "db": pathway_db[n],
                      "measured": len(pathways[n] & measured),
                      "features": int(t8.loc[n, "features"]),
                      "broad": n in broad,
                      "p": float(t8.loc[n, "p_combined_fisher"])}
                     for n in ranked],
        # Columns 1 and 2 depend on which pathways are on screen, so both answers ship and
        # the checkbox switches between them rather than the page recomputing from a
        # per-gene table it would otherwise have to download.
        "views": {"core": view([n for n in ranked if n not in broad]),
                  "all": view(ranked)},
        "pr_to_path": [{"from": p_, "to": n, "n": c}
                       for (p_, n), c in sorted(pr_path.items())],
        "unresolved": sorted(unresolved),
        "n_pathways_total": int(len(sig_paths)),
        "n_measured": len(measured),
        "broad_map_share": BROAD_MAP_SHARE,
    }


def assemble_pathways(verbose: bool = True):
    """Everything the Sankey needs: the measured layers and both databases' membership.

    Shared by the site export and the supplementary Sankey figure, so the two cannot
    drift into showing different memberships for the same pathway.

    Two databases, two membership sources, and neither may guess which rows are its own:
    MapMan's lowercase "photosynthesis" bin collides with KEGG's "Photosynthesis", so both
    resolvers take the `db` column and only look up their own rows.
    """
    tx = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"),
                     sep="\t", index_col=0)
    pr = pd.read_csv(os.path.join(TABLES, "T02_osd522_proteome.tsv"), sep="\t", index_col=0)
    conv = genesets.uniprot_to_agi()
    pr_agi = {u: conv.get(u) for u in pr.index}

    sig_paths = pd.read_csv(os.path.join(TABLES, "T08_paintomics_significant.tsv"),
                            sep="\t", comment="#")
    kegg_ids, kegg_missing = genesets.paintomics_pathway_ids(sig_paths.pathway, sig_paths.db)
    mm_sets, mm_missing = mapman.paintomics_mapman_sets(sig_paths.pathway, sig_paths.db)
    pathways = {name: genesets.pathway_genes(pid) for name, pid in kegg_ids.items()}
    pathways.update(mm_sets)
    pathway_db = {n: "KEGG" for n in kegg_ids} | {n: "MapMan" for n in mm_sets}
    unresolved = [n for n in kegg_missing if n not in mm_sets] + mm_missing

    universe = set(tx.index) | {a for a in pr_agi.values() if a}
    broad = {n for n, m in pathways.items()
             if len(m & universe) > BROAD_MAP_SHARE * len(universe)}

    if verbose:
        print(f"  pathways with membership: {len(pathways)}/{len(sig_paths)} "
              f"({len(kegg_ids)} KEGG via REST, {len(mm_sets)} MapMan via vendored diagrams"
              + (f"; {len(unresolved)} unresolved" if unresolved else "") + ")")
    return dict(tx=tx, pr=pr, pr_agi=pr_agi, sig_paths=sig_paths, pathways=pathways,
                pathway_db=pathway_db, broad=broad, unresolved=unresolved,
                universe=universe, n_kegg=len(kegg_ids), n_mapman=len(mm_sets))


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
    # BOTH CO2 columns ship. The enclosure mean is undefined for a vented case — there is
    # no closed volume to average — so VEGGIE and the open reference are all-NaN there and
    # the chart used to draw three lines while its legend named five. The leaf-surface mean
    # is defined for every case, and is the quantity the FvCB coupling actually consumes,
    # so it is what the chart shows by default. json.dump would emit a bare NaN, which
    # JSON.parse rejects, so missing values become null.
    def clean(vals):
        return [None if pd.isna(v) else round(float(v), 5) for v in vals]
    dump("enclosure_timeseries.json",
         {case: {"step": g["step"].tolist(),
                 "surface": clean(g["surf_co2_mean"]),
                 "enclosure": clean(g["dishmean_co2"])}
          for case, g in t7.groupby("case")})
    dump("carbon_retention.json",
         pd.read_csv(os.path.join(ll, "T11_photosynthesis_feedback.csv"))
         .to_dict(orient="records"))

    # --- gene sets and pathways, per layer -------------------------------------------
    # Both membership indexes are precomputed per identifier space so the browser never
    # needs a UniProt->AGI map. hypoxia_responsive has 0 measured proteins: that is real
    # biology, and the page reports "0 of 3 measured" rather than silently doing nothing.
    sets = genesets.build()
    P = assemble_pathways()
    tx, pr, pr_agi = P["tx"], P["pr"], P["pr_agi"]
    sig_paths, pathways = P["sig_paths"], P["pathways"]
    pathway_db, broad, unresolved = P["pathway_db"], P["broad"], P["unresolved"]

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

    # Broad maps are excluded from the highlight index (see BROAD_MAP_SHARE): an overlay
    # that marks 44 % of the points marks nothing, and forcing 21,000 loci through the
    # volcano thinning would defeat the thinning entirely.
    highlightable = {n: m for n, m in pathways.items() if n not in broad}

    dump("gene_sets.json", {n: per_layer(m) for n, m in sets.items()})
    dump("pathway_membership.json", {n: per_layer(m) for n, m in highlightable.items()})
    dump("pathway_coverage.json",
         {"resolved": len(pathways), "total": int(len(sig_paths)),
          "kegg": P["n_kegg"], "mapman": P["n_mapman"],
          "unresolved": sorted(unresolved), "broad": sorted(broad),
          "broad_share": BROAD_MAP_SHARE})

    # Every gene the page can highlight must survive the volcano thinning.
    must_keep = set().union(*sets.values()) | (set().union(*highlightable.values())
                                               if highlightable else set())

    # --- measured omics, thinned for the browser but never below must_keep ------------
    dump("volcano_transcriptome.json",
         volcano(tx, "log2fc", "fdr", "AGI", 0.05, must_keep=must_keep))
    prot_keep = {u for u, a in pr_agi.items() if a in must_keep}
    dump("volcano_proteome.json",
         volcano(pr, "log2fc", "adj_p", "UniProt", 0.05, must_keep=prot_keep))

    # --- sankey: transcript bin -> protein bin -> pathway -----------------------------
    dump("sankey.json",
         build_sankey(tx, pr, pr_agi, pathways, pathway_db, broad, sig_paths,
                      unresolved))

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
