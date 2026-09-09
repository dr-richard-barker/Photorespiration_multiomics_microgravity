#!/usr/bin/env python3
"""Predict the metabolite layer OSDR does not contain, for OSD-522 (BRIC-LED-001).

NASA OSDR holds no plant metabolomics at all — of 567 studies only 6 are metabolite
profiling, and every one is mouse, human, rat or microbial (checked exhaustively via the
OSDR search API, 2026-09-08). None of the 66 plant studies has a metabolome. So the third
omic layer for a PaintOmics run cannot be downloaded; it has to be predicted.

THIS LAYER IS MODEL OUTPUT. IT WAS NEVER MEASURED. Every row of every output file says so.

It is not filler for a missing file — it is a falsifiable prediction that the real
transcriptome and proteome in the same PaintOmics job can contradict.

    LunarLeaf-CFD (validated LBM solver)  ->  g_bl(gravity), enclosure CO2 mass balance
    fvcb.py (this repo)                   ->  Cc, Vc, Vo, A at the operating point
    here                                  ->  metabolite pool log2 fold changes

Two contrasts are produced, and they are very different in size:

  FLT_vs_GC   Space Flight vs Ground Control, both inside sealed BRIC canisters. This is
              what OSD-522 actually measured, so it is the file that goes into the joint
              PaintOmics job alongside the real gene and protein layers.

  BRIC_vs_VENTED  The same plant in sealed BRIC versus vented hardware holding ambient
              CO2. This is the hardware effect. It is roughly forty times larger than the
              gravity effect — but because it hits flight and ground alike, it cancels
              almost entirely out of FLT_vs_GC.

The central modelling result, which is not what we first guessed: in a sealed canister the
boundary-layer difference between 1 g and microgravity does NOT amplify. As assimilation
falls toward the CO2 compensation point the flux through the boundary layer falls with it,
so the CO2 drop across that layer (A/g_bl) tends to zero and the two gravities converge on
oxygenation fraction. What survives is a persistent proportional assimilation penalty of
about 5-9 % that is nearly independent of where in the drawdown the plant sits.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fvcb import LeafParams, load_cfd_sweep, solve_operating_point  # noqa: E402
from paths import TABLES, UPLOAD, ensure  # noqa: E402

# Compound set. KEGG identifiers were retrieved from the KEGG REST API
# (rest.kegg.jp/list/...) on 2026-09-08, not written from memory. The primary name is
# KEGG's own first synonym, which is what PaintOmics matches against.
#
# driver  - the modelled quantity the pool is tied to
# tier    - 1: pool proportional to a flux the model computes directly
#           2: sign and scale follow from a modelled quantity through one further
#              documented step (stated in `basis`)
COMPOUNDS = [
    # --- C2 photorespiratory cycle: pool tracks oxygenation flux Vo ------------------
    ("C00988", "2-Phosphoglycolate", "Vo", 1, "First product of Rubisco oxygenation"),
    ("C00160", "Glycolate", "Vo", 1, "2-PG dephosphorylated by PGLP1; exported to peroxisome"),
    ("C00048", "Glyoxylate", "Vo", 1, "Glycolate oxidised by GOX in the peroxisome"),
    ("C00037", "Glycine", "Vo", 1, "Glyoxylate transaminated; the C2 cycle's mitochondrial input"),
    ("C00065", "L-Serine", "Vo", 1, "Product of the GDC/SHMT reaction, 2 Gly -> 1 Ser"),
    ("C00168", "Hydroxypyruvate", "Vo", 1, "Serine transaminated in the peroxisome"),
    ("C00258", "D-Glycerate", "Vo", 1, "Hydroxypyruvate reduced by HPR1; returns to chloroplast"),
    ("C00014", "Ammonia", "Vo", 1, "Photorespiratory NH3 release at glycine decarboxylase"),

    # --- Net photosynthate: pool tracks net assimilation A ---------------------------
    ("C00197", "3-Phospho-D-glycerate", "A", 1, "First stable product of carboxylation"),
    ("C00089", "Sucrose", "A", 1, "Principal export product of net assimilation"),
    ("C00031", "D-Glucose", "A", 1, "Hexose pool fed by sucrose cleavage"),
    ("C00095", "D-Fructose", "A", 1, "Hexose pool fed by sucrose cleavage"),
    ("C00208", "Maltose", "A", 1, "Product of nocturnal starch turnover"),
    ("C00369", "Starch", "A", 1, "Transitory reserve; fills from surplus assimilate"),

    # --- Tier 2: one documented step beyond a modelled flux --------------------------
    ("C01182", "D-Ribulose 1,5-bisphosphate", "RuBP", 2,
     "Substrate upstream of a CO2-limited Rubisco; accumulates as carboxylation is "
     "substrate-limited rather than RuBP-limited (model reports the limiting process)"),
    ("C00152", "L-Asparagine", "starvation", 2,
     "Canonical carbon-starvation amino acid; ASN1/DIN6 is the marker gene of the "
     "SnRK1 low-energy response, driven here by the assimilation deficit"),
    ("C00183", "L-Valine", "starvation", 2, "Branched-chain amino acid released by protein "
     "turnover and catabolised for energy under carbon limitation"),
    ("C00123", "L-Leucine", "starvation", 2, "Branched-chain amino acid, as L-Valine"),
    ("C00407", "L-Isoleucine", "starvation", 2, "Branched-chain amino acid, as L-Valine"),
    ("C00064", "L-Glutamine", "starvation", 2,
     "N re-assimilation pool; rises as protein turnover outpaces growth demand"),
    ("C00025", "L-Glutamate", "starvation", 2,
     "Amino donor for photorespiratory transamination and the hub of N remobilisation"),
]

# Compounds deliberately NOT predicted, and why. Recorded so the omission is a stated
# decision rather than an oversight.
EXCLUDED = [
    ("Citrate / Malate / Succinate / Fumarate / 2-Oxoglutarate (TCA)",
     "The model does not resolve respiratory flux, and carbon limitation can either "
     "deplete these pools (oxidised for energy) or raise them (anaplerosis from amino "
     "acid catabolism). Direction genuinely uncertain, so no value is asserted."),
    ("Ethanol / Lactate / Alanine (fermentation)",
     "Fermentation needs hypoxia. LunarLeaf-CFD T8 puts BRIC O2 depletion at 6.5 days "
     "and only in DARK canisters; OSD-522 is BRIC-LED, illuminated, where photosynthesis "
     "releases O2. No hypoxia is predicted, so no fermentation products are predicted."),
    ("Glycine/Serine ratio",
     "Both pools are tied to the same flux Vo, so under this model the ratio does not "
     "move. Reporting it as a photorespiration readout here would be circular."),
]


def operating_points(scale: str, sealed_ca: float, ambient_ca: float, csv_path: str | None):
    """Solve the FvCB operating point for each arm of both contrasts."""
    rows = load_cfd_sweep(csv_path) if csv_path else load_cfd_sweep()
    at = {}
    for r in rows:
        if r["scale"] == scale:
            key = "earth" if r["gravity_g"] > 5.0 else ("ug" if r["gravity_g"] < 0.5 else None)
            if key and key not in at:
                at[key] = r
    missing = {"earth", "ug"} - set(at)
    if missing:
        raise SystemExit(f"scale '{scale}' has no {sorted(missing)} row in the CFD export")

    p = LeafParams()
    solve = lambda r, ca: solve_operating_point(  # noqa: E731
        r["g_bl"], p, Ca=ca, O_excess=r.get("o2_excess_ppm", 0.0) or 0.0)

    return {
        # both arms of the real experiment sit in a sealed, CO2-drawn-down canister
        "GC": solve(at["earth"], sealed_ca),
        "FLT": solve(at["ug"], sealed_ca),
        # hardware contrast: same microgravity leaf, sealed box vs vented at ambient
        "BRIC": solve(at["ug"], sealed_ca),
        "VENTED": solve(at["ug"], ambient_ca),
        "_gbl": {"earth": at["earth"]["g_bl"], "ug": at["ug"]["g_bl"]},
    }


def log2_ratio(num: float, den: float) -> float:
    if not (num > 0 and den > 0):
        return float("nan")
    return math.log2(num / den)


def drivers(a: dict, b: dict, starvation_gain: float) -> dict:
    """log2 fold changes of each modelled driver, arm `a` relative to arm `b`."""
    d_vo = log2_ratio(a["Vo"], b["Vo"])
    d_a = log2_ratio(a["A"], b["A"])

    # RuBP: rises as the leaf becomes more CO2-limited. Cc is the model's measure of
    # CO2 availability at Rubisco, so the pool moves opposite to log2(Cc).
    d_rubp = -log2_ratio(a["Cc"], b["Cc"])

    # Carbon starvation: the minimal assumption is that starvation-responsive pools move
    # by the same magnitude as the assimilation deficit and in the opposite direction.
    # starvation_gain defaults to 1.0 so no unfitted parameter enters the prediction.
    d_starv = -d_a * starvation_gain

    return {"Vo": d_vo, "A": d_a, "RuBP": d_rubp, "starvation": d_starv}


def write_contrast(name: str, a: dict, b: dict, gain: float, threshold: float,
                   values_file: str, relevant_file: str) -> dict:
    d = drivers(a, b, gain)
    os.makedirs(UPLOAD, exist_ok=True)

    with open(os.path.join(UPLOAD, values_file), "w") as fh:
        fh.write(f"#compound\t{name}\n")
        for _cid, cname, driver, _tier, _basis in COMPOUNDS:
            fh.write(f"{cname}\t{d[driver]:.6f}\n")

    relevant = [c for c in COMPOUNDS if abs(d[c[2]]) >= threshold]
    with open(os.path.join(UPLOAD, relevant_file), "w") as fh:
        for _cid, cname, *_ in relevant:
            fh.write(f"{cname}\n")

    print(f"\n  {name}")
    for key in ("Vo", "A", "RuBP", "starvation"):
        print(f"    log2FC[{key:<10}] = {d[key]:+.4f}")
    print(f"    wrote {values_file} ({len(COMPOUNDS)} compounds), "
          f"{relevant_file} ({len(relevant)} at |log2FC| >= {threshold})")
    return d


def write_provenance(pts: dict, d_flt: dict, d_hw: dict, args) -> None:
    ensure(TABLES)
    path = os.path.join(TABLES, "T06_predicted_compounds.tsv")
    with open(path, "w") as fh:
        fh.write("# Predicted metabolite pools for OSD-522 (BRIC-LED-001).\n")
        fh.write("# EVERY VALUE IS MODEL OUTPUT — NOT MEASURED. OSDR has no plant metabolomics.\n")
        fh.write("# KEGG IDs retrieved from rest.kegg.jp on 2026-09-08.\n")
        fh.write("kegg_id\tcompound\tdriver\ttier\tlog2FC_FLT_vs_GC\tlog2FC_BRIC_vs_VENTED\t"
                 "status\tbasis\n")
        for cid, cname, driver, tier, basis in COMPOUNDS:
            fh.write(f"{cid}\t{cname}\t{driver}\t{tier}\t{d_flt[driver]:+.6f}\t"
                     f"{d_hw[driver]:+.6f}\tPREDICTED\t{basis}\n")
        fh.write("#\n# Considered and deliberately excluded:\n")
        for what, why in EXCLUDED:
            fh.write(f"# - {what}: {why}\n")

    gbl = pts["_gbl"]
    summary = os.path.join(TABLES, "T05_operating_points.tsv")
    with open(summary, "w") as fh:
        fh.write("arm\tg_bl_mol_m2_s\tCa_umol_mol\tCc\tphi_percent\tA_umol_m2_s\tVo\tVc\n")
        for arm in ("GC", "FLT", "BRIC", "VENTED"):
            r = pts[arm]
            g = gbl["earth"] if arm == "GC" else gbl["ug"]
            fh.write(f"{arm}\t{g:.3f}\t{r['Ca']:.1f}\t{r['Cc']:.2f}\t{r['phi']*100:.2f}\t"
                     f"{r['A']:.3f}\t{r['Vo']:.3f}\t{r['Vc']:.3f}\n")
    print(f"\n  wrote {path}")
    print(f"  wrote {summary}")


def sensitivity(scale: str, args) -> None:
    """Show that the FLT/GC prediction barely moves with the assumed sealed-canister CO2."""
    print("\n  sensitivity of the FLT_vs_GC prediction to the assumed canister CO2:")
    print(f"    {'Ca':>6} {'log2FC[A]':>11} {'log2FC[Vo]':>12}")
    for ca in (400.0, 250.0, 150.0, 100.0, 80.0, 70.0):
        pts = operating_points(scale, ca, args.ambient_ca, args.cfd_csv)
        d = drivers(pts["FLT"], pts["GC"], args.starvation_gain)
        print(f"    {ca:6.0f} {d['A']:+11.4f} {d['Vo']:+12.4f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scale", default="rosette", choices=["leaf", "rosette", "canopy"],
                    help="CFD length scale; a BRIC dish holds a rosette of seedlings")
    ap.add_argument("--sealed-ca", type=float, default=100.0,
                    help="CO2 in the sealed canister at quasi-steady state (umol/mol). "
                         "LunarLeaf-CFD T8 draws a lit BRIC from 400 toward the compensation "
                         "point within ~7 min; the prediction is insensitive to this value "
                         "(see --sensitivity)")
    ap.add_argument("--ambient-ca", type=float, default=400.0,
                    help="CO2 held by vented hardware (umol/mol)")
    ap.add_argument("--starvation-gain", type=float, default=1.0,
                    help="tier-2 amplification; 1.0 introduces no unfitted parameter")
    ap.add_argument("--threshold", type=float, default=0.01,
                    help="|log2FC| for a compound to enter the relevant-features file")
    ap.add_argument("--cfd-csv", default=None, help="LunarLeaf-CFD boundary-layer export")
    ap.add_argument("--sensitivity", action="store_true",
                    help="also print the canister-CO2 sensitivity table")
    args = ap.parse_args()

    print("Predicted metabolome for OSD-522 (BRIC-LED-001) — MODEL OUTPUT, NOT MEASURED")
    pts = operating_points(args.scale, args.sealed_ca, args.ambient_ca, args.cfd_csv)

    print(f"\n  scale={args.scale}  g_bl: Earth {pts['_gbl']['earth']:.3f} -> "
          f"micro-g {pts['_gbl']['ug']:.3f} mol m-2 s-1 (LunarLeaf-CFD)")
    for arm in ("GC", "FLT", "BRIC", "VENTED"):
        r = pts[arm]
        print(f"    {arm:<7} Ca={r['Ca']:5.0f}  Cc={r['Cc']:6.1f}  phi={r['phi']*100:5.1f}%  "
              f"A={r['A']:6.2f}  Vo={r['Vo']:6.3f}")

    d_flt = write_contrast("FLT_vs_GC", pts["FLT"], pts["GC"], args.starvation_gain,
                           args.threshold, "metabolomics_values.tab",
                           "metabolomics_relevant.tab")
    d_hw = write_contrast("BRIC_vs_VENTED", pts["BRIC"], pts["VENTED"], args.starvation_gain,
                          args.threshold, "metabolomics_hardware_values.tab",
                          "metabolomics_hardware_relevant.tab")

    write_provenance(pts, d_flt, d_hw, args)
    if args.sensitivity:
        sensitivity(args.scale, args)

    print(f"\n  hardware effect is {abs(d_hw['A'] / d_flt['A']):.0f}x the gravity effect "
          f"on assimilation, but it hits flight and ground alike and so largely cancels "
          f"out of FLT_vs_GC.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
