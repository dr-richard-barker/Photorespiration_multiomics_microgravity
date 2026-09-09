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


def volcano(df, lfc, p, idname, alpha, keep_ns=4000, seed=0):
    d = df.dropna(subset=[lfc, p]).copy()
    sig = d[p] < alpha
    ns = d[~sig]
    if len(ns) > keep_ns:
        ns = ns.sample(keep_ns, random_state=seed)
    out = pd.concat([d[sig], ns])
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

    # --- measured omics, thinned for the browser --------------------------------------
    tx = pd.read_csv(os.path.join(TABLES, "T01_osd522_transcriptome.tsv"),
                     sep="\t", index_col=0)
    dump("volcano_transcriptome.json", volcano(tx, "log2fc", "fdr", "AGI", 0.05))
    pr = pd.read_csv(os.path.join(TABLES, "T02_osd522_proteome.tsv"),
                     sep="\t", index_col=0)
    dump("volcano_proteome.json", volcano(pr, "log2fc", "adj_p", "UniProt", 0.05))

    # --- gene-set membership, so the volcano can highlight a set ----------------------
    sets = genesets.build()
    measured = set(tx.index)
    dump("gene_sets.json",
         {name: sorted(m & measured) for name, m in sets.items()})

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
