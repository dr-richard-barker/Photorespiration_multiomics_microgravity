#!/usr/bin/env python3
"""How much of the photosynthesis result depends on the irradiance the model is run at.

Every operating point, deficit and metabolite prediction in this repository is solved at
the fvcb.py default, Q = 1000 umol m-2 s-1 (saturating light). OSD-522 grew its plants under
60. This re-solves the same quantities at both irradiances and writes them side by side,
so the manuscript's low-light numbers come from a file rather than a one-off calculation.

  operating points  the three arms of Table 1 (ground sealed, flight sealed, ventilated),
                    rosette scale, with the same g_bl and O2 excess as
                    03_predict_metabolome.py uses
  compensation      where a lit sealed canister would actually stop (A = 0); Table 1
                    instead places the sealed arms at an assumed 100 ppm
  deficit           the microgravity assimilation deficit of Fig. 2c, (A_ug - A_1g)/A_1g,
                    over the same CO2 grid and the same "A_1g > 0.05" floor, per scale

Writes results/tables/T13_light_sensitivity.tsv.
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fvcb import LeafParams, load_cfd_sweep, solve_operating_point  # noqa: E402
from paths import TABLES, ensure  # noqa: E402

IRRADIANCES = (1000.0, 60.0)    # fvcb.py default; OSD-522 growth light
SEALED_CA, AMBIENT_CA = 100.0, 400.0   # as 03_predict_metabolome.py
CA_GRID = np.linspace(60, 420, 145)    # as fig02_fvcb_operating_points.py
SCALES = ["leaf", "rosette", "canopy"]


def cfd_rows() -> dict[str, dict[str, dict]]:
    """scale -> {"earth": row, "ug": row}, first matching row as in fig02 and 03."""
    out: dict[str, dict[str, dict]] = {}
    for r in load_cfd_sweep():
        key = "earth" if r["gravity_g"] > 5.0 else ("ug" if r["gravity_g"] < 0.5 else None)
        if key:
            out.setdefault(r["scale"], {}).setdefault(key, r)
    return out


def main() -> int:
    p = LeafParams()
    rows = cfd_rows()
    ensure(TABLES)
    path = os.path.join(TABLES, "T13_light_sensitivity.tsv")
    with open(path, "w") as fh:
        fh.write("# Irradiance sensitivity. Q = 1000 is what every other table uses; Q = 60 is "
                 "the OSD-522 growth light.\n")
        fh.write("# Derived by scripts/09_light_sensitivity.py from fvcb.py and the vendored "
                 "LunarLeaf-CFD T13.\n")
        fh.write("quantity\tQ_umol_m2_s\tscale\tCa_umol_mol\tphi_percent\tA_umol_m2_s\t"
                 "deficit_min_percent\tdeficit_max_percent\n")
        for Q in IRRADIANCES:
            ros = rows["rosette"]
            for arm, r, ca in (("operating_point:GC_sealed", ros["earth"], SEALED_CA),
                               ("operating_point:FLT_sealed", ros["ug"], SEALED_CA),
                               ("operating_point:VENTED", ros["ug"], AMBIENT_CA)):
                op = solve_operating_point(r["g_bl"], p, Ca=ca, Q=Q,
                                           O_excess=r.get("o2_excess_ppm", 0.0) or 0.0)
                fh.write(f"{arm}\t{Q:.0f}\trosette\t{ca:.0f}\t{op['phi'] * 100:.2f}\t"
                         f"{op['A']:.3f}\t\t\n")
            # Table 1 puts the sealed arms at an ASSUMED 100 ppm; this is where a lit canister
            # would actually stop, A = 0, found by bisection on Ca.
            g, lo, hi = ros["ug"]["g_bl"], 1.0, AMBIENT_CA
            for _ in range(60):
                mid = (lo + hi) / 2
                lo, hi = (mid, hi) if solve_operating_point(g, p, Ca=mid, Q=Q)["A"] < 0 \
                    else (lo, mid)
            op = solve_operating_point(g, p, Ca=hi, Q=Q)
            fh.write(f"compensation_point:FLT\t{Q:.0f}\trosette\t{hi:.1f}\t"
                     f"{op['phi'] * 100:.2f}\t{op['A']:.3f}\t\t\n")
            for scale in SCALES:
                e, u = rows[scale]["earth"]["g_bl"], rows[scale]["ug"]["g_bl"]
                a_e = np.array([solve_operating_point(e, p, Ca=float(c), Q=Q)["A"]
                                for c in CA_GRID])
                a_u = np.array([solve_operating_point(u, p, Ca=float(c), Q=Q)["A"]
                                for c in CA_GRID])
                ok = a_e > 0.05
                d = np.abs((a_u[ok] - a_e[ok]) / a_e[ok] * 100)
                fh.write(f"deficit\t{Q:.0f}\t{scale}\t{CA_GRID[ok].min():.0f}-"
                         f"{CA_GRID[ok].max():.0f}\t\t\t{d.min():.2f}\t{d.max():.2f}\n")
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
