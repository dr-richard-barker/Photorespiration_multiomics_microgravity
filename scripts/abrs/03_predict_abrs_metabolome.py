#!/usr/bin/env python3
"""Predict the metabolite layer for ABRS (TAGES / OSD-7 / OSD-16).

In the ABRS flight hardware:
  - Active airflow ventilation sweeps the leaf boundary layer.
  - ISS cabin atmosphere maintains CO2 at ~3,000-4,000 ppm (nominal 3,500 ppm).
  - Ground controls (OES chamber) faithfully replicated this 48-hour delayed playback.
  - Active catalytic ethylene scrubbers prevented volatile accumulation.

This script solves the FvCB operating points under ABRS environmental conditions
and translates the fluxes into predicted log2 fold changes for the 21 canonical
metabolite pools defined in the previous analysis (scripts/03_predict_metabolome.py).

Contrasts produced:
  1. FLT_vs_GC (ABRS): Flight vs Ground Control inside ventilated ABRS hardware
     holding ISS cabin CO2 (~3,500 ppm).
  2. ABRS_vs_BRIC: The cross-hardware atmospheric contrast (ABRS 3,500 ppm vs
     sealed BRIC-LED 100 ppm drawdown). Demonstrates the massive metabolic shift
     induced by enclosure ventilation.

Outputs:
  data/abrs/paintomics_upload/abrs_metabolomics_values.tab
  data/abrs/paintomics_upload/abrs_metabolomics_relevant.tab
  results/abrs/tables/T06_abrs_predicted_compounds.tsv
  results/abrs/tables/T05_abrs_operating_points.tsv
"""

from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))
from fvcb import LeafParams, load_cfd_sweep, solve_operating_point  # noqa: E402
from paths import ensure  # noqa: E402

ABRS_UPLOAD = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "data", "abrs", "paintomics_upload"
))
ABRS_TABLES = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "tables"
))

# Exact 21 compounds from previous analysis (C00xxx KEGG IDs)
COMPOUNDS = [
    # C2 photorespiratory cycle (tracks Vo)
    ("C00988", "2-Phosphoglycolate", "Vo", 1, "First product of Rubisco oxygenation"),
    ("C00160", "Glycolate", "Vo", 1, "2-PG dephosphorylated by PGLP1; exported to peroxisome"),
    ("C00048", "Glyoxylate", "Vo", 1, "Glycolate oxidised by GOX in the peroxisome"),
    ("C00037", "Glycine", "Vo", 1, "Glyoxylate transaminated; mitochondrial C2 cycle input"),
    ("C00065", "L-Serine", "Vo", 1, "Product of GDC/SHMT reaction, 2 Gly -> 1 Ser"),
    ("C00168", "Hydroxypyruvate", "Vo", 1, "Serine transaminated in peroxisome"),
    ("C00258", "D-Glycerate", "Vo", 1, "Hydroxypyruvate reduced by HPR1; returns to chloroplast"),
    ("C00014", "Ammonia", "Vo", 1, "Photorespiratory NH3 release at glycine decarboxylase"),

    # Net photosynthate (tracks A)
    ("C00197", "3-Phospho-D-glycerate", "A", 1, "First stable product of carboxylation"),
    ("C00089", "Sucrose", "A", 1, "Principal export product of net assimilation"),
    ("C00031", "D-Glucose", "A", 1, "Hexose pool fed by sucrose cleavage"),
    ("C00095", "D-Fructose", "A", 1, "Hexose pool fed by sucrose cleavage"),
    ("C00208", "Maltose", "A", 1, "Product of nocturnal starch turnover"),
    ("C00369", "Starch", "A", 1, "Transitory reserve; fills from surplus assimilate"),

    # Tier 2: one step beyond modelled flux
    ("C01182", "D-Ribulose 1,5-bisphosphate", "RuBP", 2,
     "Substrate upstream of Rubisco; moves inversely with log2(Cc)"),
    ("C00152", "L-Asparagine", "starvation", 2,
     "Canonical carbon-starvation marker (ASN1/DIN6), inversely proportional to A"),
    ("C00183", "L-Valine", "starvation", 2, "BCAA released by protein turnover during starvation"),
    ("C00123", "L-Leucine", "starvation", 2, "BCAA released by protein turnover during starvation"),
    ("C00407", "L-Isoleucine", "starvation", 2, "BCAA released by protein turnover during starvation"),
    ("C00064", "L-Glutamine", "starvation", 2, "N re-assimilation pool; rises as protein turnover outpaces growth"),
    ("C00025", "L-Glutamate", "starvation", 2, "Amino donor for transamination and N remobilisation hub"),
]


def log2_ratio(num: float, den: float) -> float:
    if not (num > 0 and den > 0):
        return float("nan")
    return math.log2(num / den)


def compute_abrs_operating_points(scale: str = "rosette",
                                  abrs_ca: float = 3500.0,
                                  bric_ca: float = 100.0) -> dict:
    """Solve FvCB operating points for ABRS vs BRIC conditions."""
    rows = load_cfd_sweep()
    at = {}
    for r in rows:
        if r["scale"] == scale:
            key = "earth" if r["gravity_g"] > 5.0 else ("ug" if r["gravity_g"] < 0.5 else None)
            if key and key not in at:
                at[key] = r

    p = LeafParams()
    solve = lambda r, ca: solve_operating_point(
        r["g_bl"], p, Ca=ca, O_excess=r.get("o2_excess_ppm", 0.0) or 0.0
    )

    # Operating points
    gc_abrs = solve(at["earth"], abrs_ca)
    flt_abrs = solve(at["ug"], abrs_ca)
    flt_bric = solve(at["ug"], bric_ca)

    return {
        "GC_ABRS": gc_abrs,
        "FLT_ABRS": flt_abrs,
        "FLT_BRIC": flt_bric,
        "_gbl": {"earth": at["earth"]["g_bl"], "ug": at["ug"]["g_bl"]},
    }


def compute_drivers(flt: dict, gc: dict, starvation_gain: float = 1.0) -> dict:
    d_vo = log2_ratio(flt["Vo"], gc["Vo"])
    d_a = log2_ratio(flt["A"], gc["A"])
    d_rubp = -log2_ratio(flt["Cc"], gc["Cc"])
    d_starv = -d_a * starvation_gain
    return {"Vo": d_vo, "A": d_a, "RuBP": d_rubp, "starvation": d_starv}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scale", default="rosette", choices=["leaf", "rosette", "canopy"])
    ap.add_argument("--abrs-ca", type=float, default=3500.0,
                    help="ISS cabin CO2 for ABRS in ppm (default: 3500)")
    ap.add_argument("--threshold", type=float, default=0.01,
                    help="minimum |log2FC| for relevant compounds")
    args = ap.parse_args()

    ensure(ABRS_UPLOAD, ABRS_TABLES)

    print("FvCB Photosynthesis & Metabolome Prediction for ABRS Hardware")
    print(f"  Scale: {args.scale}, ABRS Ca: {args.abrs_ca} ppm (ISS cabin)\n")

    pts = compute_abrs_operating_points(args.scale, abrs_ca=args.abrs_ca)

    gc = pts["GC_ABRS"]
    flt = pts["FLT_ABRS"]
    bric = pts["FLT_BRIC"]

    print("Operating Points:")
    print(f"  ABRS Ground Control (1g, {args.abrs_ca} ppm):")
    print(f"    Cc={gc['Cc']:.1f} ppm, A={gc['A']:.2f}, Vo={gc['Vo']:.2f}, phi={gc['phi']*100:.1f}%, limiting={gc['limiting']}")
    print(f"  ABRS Flight (micro-g, {args.abrs_ca} ppm):")
    print(f"    Cc={flt['Cc']:.1f} ppm, A={flt['A']:.2f}, Vo={flt['Vo']:.2f}, phi={flt['phi']*100:.1f}%, limiting={flt['limiting']}")
    print(f"  BRIC-LED Flight (micro-g, 100 ppm):")
    print(f"    Cc={bric['Cc']:.1f} ppm, A={bric['A']:.2f}, Vo={bric['Vo']:.2f}, phi={bric['phi']*100:.1f}%, limiting={bric['limiting']}")

    # Drivers for ABRS FLT vs GC
    d_abrs = compute_drivers(flt, gc)
    print("\nABRS FLT vs GC Drivers:")
    for k, v in d_abrs.items():
        print(f"  log2FC[{k:<10}] = {v:+.4f}")

    # Drivers for Cross-Hardware (ABRS FLT vs BRIC FLT)
    d_hw = compute_drivers(flt, bric)
    print("\nCross-Hardware (ABRS vs BRIC in flight) Drivers:")
    for k, v in d_hw.items():
        print(f"  log2FC[{k:<10}] = {v:+.4f}")

    # Write PaintOmics files for ABRS FLT vs GC
    val_path = os.path.join(ABRS_UPLOAD, "abrs_metabolomics_values.tab")
    with open(val_path, "w") as fh:
        fh.write("#compound\tFLT_vs_GC\n")
        for _cid, cname, driver, _tier, _basis in COMPOUNDS:
            fh.write(f"{cname}\t{d_abrs[driver]:.6f}\n")
    print(f"\nwrote {val_path}")

    # Relevant compounds: compounds with non-negligible predicted shift
    relevant = [c for c in COMPOUNDS if abs(d_abrs[c[2]]) >= args.threshold]
    rel_path = os.path.join(ABRS_UPLOAD, "abrs_metabolomics_relevant.tab")
    with open(rel_path, "w") as fh:
        for _cid, cname, *_ in relevant:
            fh.write(f"{cname}\n")
    print(f"wrote {rel_path} ({len(relevant)} compounds at |log2FC| >= {args.threshold})")

    # Write full provenance table
    t06_path = os.path.join(ABRS_TABLES, "T06_abrs_predicted_compounds.tsv")
    with open(t06_path, "w") as fh:
        fh.write("# Predicted metabolite pools for ABRS (TAGES / OSD-7).\n")
        fh.write("# MODEL OUTPUT under ABRS boundary conditions (ISS cabin CO2 ~3500 ppm, active ventilation).\n")
        fh.write("kegg_id\tcompound\tdriver\ttier\tlog2FC_FLT_vs_GC_ABRS\tlog2FC_ABRS_vs_BRIC\tstatus\tbasis\n")
        for cid, cname, driver, tier, basis in COMPOUNDS:
            fh.write(f"{cid}\t{cname}\t{driver}\t{tier}\t{d_abrs[driver]:+.6f}\t{d_hw[driver]:+.6f}\tPREDICTED\t{basis}\n")
    print(f"wrote {t06_path}")

    # Write operating points table
    t05_path = os.path.join(ABRS_TABLES, "T05_abrs_operating_points.tsv")
    with open(t05_path, "w") as fh:
        fh.write("condition\thardware\tgravity\tCa_ppm\tCc_ppm\tA_umol_m2_s\tVo_umol_m2_s\tphi_pct\tlimiting\n")
        fh.write(f"GC_ABRS\tABRS\t1g\t{args.abrs_ca}\t{gc['Cc']:.1f}\t{gc['A']:.2f}\t{gc['Vo']:.2f}\t{gc['phi']*100:.1f}\t{gc['limiting']}\n")
        fh.write(f"FLT_ABRS\tABRS\tmicro-g\t{args.abrs_ca}\t{flt['Cc']:.1f}\t{flt['A']:.2f}\t{flt['Vo']:.2f}\t{flt['phi']*100:.1f}\t{flt['limiting']}\n")
        fh.write(f"FLT_BRIC\tBRIC-LED\tmicro-g\t100.0\t{bric['Cc']:.1f}\t{bric['A']:.2f}\t{bric['Vo']:.2f}\t{bric['phi']*100:.1f}\t{bric['limiting']}\n")
    print(f"wrote {t05_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
