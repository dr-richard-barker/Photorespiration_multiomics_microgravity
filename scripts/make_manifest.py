#!/usr/bin/env python3
"""Write results/tables/MANIFEST.tsv — one row per derived table.

FAIR data needs to say, for every file, what it holds and which script produced it. This
generates that index and FAILS if any table in results/tables/ is undescribed, so a new
output cannot quietly appear without provenance.
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import TABLES  # noqa: E402

DESCRIBED = {
 "T01_osd522_transcriptome.tsv": ("scripts/01_dge_rnaseq.py",
   "OSD-522 shoot transcriptome, flight vs ground: log2FC, p, FDR per AGI locus",
   "NASA OSDR OSD-522 RSEM counts"),
 "T02_osd522_proteome.tsv": ("scripts/02_prep_proteomics.py",
   "OSD-522 shoot proteome, flight vs ground: deposited S/G log2 ratio and adjusted p",
   "NASA OSDR OSD-522 shoot SOL + MEM proteomics reports"),
 "T04_photorespiration_vs_gravity.csv": ("scripts/fvcb.py",
   "FvCB operating point across the CFD gravity/scale sweep: Cc, phi, A, Rp, Rp/A",
   "data/lunarleaf_gbl_sweep.csv"),
 "T05_operating_points.tsv": ("scripts/03_predict_metabolome.py",
   "The four modelled arms (flight, ground, sealed, vented): g_bl, Ca, Cc, phi, A, Vo, Vc",
   "LunarLeaf-CFD T13 + fvcb.py"),
 "T06_predicted_compounds.tsv": ("scripts/03_predict_metabolome.py",
   "21 PREDICTED metabolite pools with KEGG id, driver flux, tier and both contrasts. "
   "MODEL OUTPUT, NEVER MEASURED.",
   "fvcb.py fluxes; KEGG ids from rest.kegg.jp"),
 "T07_falsification_osd522.tsv": ("scripts/04_falsification_check.py",
   "Blind prediction vs measurement on OSD-522: gene-set shifts, Mann-Whitney p, verdicts",
   "T01, T02, KEGG gene sets"),
 "T08_paintomics_significant.tsv": ("scripts/05_submit_paintomics.py (job m1z16Qg3DK)",
   "All 31 pathways significant at combined Fisher p < 0.05, with per-omic p-values",
   "PaintOmics AI"),
 "T09_paintomics_carbon.tsv": ("scripts/05_submit_paintomics.py (job m1z16Qg3DK)",
   "The carbon pathways the model made claims about, significant or not",
   "PaintOmics AI"),
 "T10_metabolite_hubs.tsv": ("scripts/05_submit_paintomics.py (job m1z16Qg3DK)",
   "PaintOmics metabolite hub analysis: predicted compounds ranked by REAL DE genes nearby",
   "PaintOmics AI"),
 "T11_hardware_ladder.tsv": ("scripts/07_hardware_ladder.py",
   "Six OSDR studies x ten gene sets: median shift, Mann-Whitney p, model prediction",
   "data/study_registry.tsv + GeneLab DGE tables"),
 "T12_ladder_contrasts.tsv": ("scripts/07_hardware_ladder.py",
   "Set-vs-set separations per study — the illumination test and the enclosure test",
   "data/study_registry.tsv + GeneLab DGE tables"),
}


def main() -> int:
    present = sorted(f for f in os.listdir(TABLES)
                     if f.endswith((".tsv", ".csv")) and f != "MANIFEST.tsv")
    missing = [f for f in present if f not in DESCRIBED]
    stale = [f for f in DESCRIBED if f not in present]

    out = os.path.join(TABLES, "MANIFEST.tsv")
    with open(out, "w") as fh:
        fh.write("# Every derived table in this directory, what it holds, and what made it.\n")
        fh.write("# Regenerate with scripts/make_manifest.py, which fails on an "
                 "undescribed file.\n")
        fh.write("file\tbytes\tgenerating_script\tdescription\tsource_data\n")
        for f in present:
            script, desc, src = DESCRIBED.get(f, ("?", "UNDESCRIBED", "?"))
            size = os.path.getsize(os.path.join(TABLES, f))
            fh.write(f"{f}\t{size}\t{script}\t{desc}\t{src}\n")

    print(f"wrote {out}  ({len(present)} tables)")
    if stale:
        print(f"  note: described but not present: {', '.join(stale)}")
    if missing:
        print(f"  FAIL: undescribed tables: {', '.join(missing)}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
