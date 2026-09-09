"""
FvCB oxygenation calculator — coupling LunarLeaf-CFD g_bl(g) to photorespiration.

This is the first component of the subcellular photorespiration multi-omics model.
It takes the gravity-dependent leaf boundary-layer conductance g_bl(g) produced by
the LunarLeaf-CFD solver and pushes it through the conductance chain

    Ca --[g_bl]--> Cs --[g_s]--> Ci --[g_m]--> Cc   (CO2 reaching Rubisco)

into a Farquhar-von Caemmerer-Berry (FvCB) C3 model *with the Rubisco oxygenation
term made explicit*, so the output is not just net assimilation A but the
photorespiratory quantities:

    Vo/Vc      oxygenation : carboxylation ratio at Rubisco   (= 2 Gamma* / Cc)
    phi        oxygenation fraction  Vo / (Vc + Vo)
    Rp         photorespiratory CO2 release  = 0.5 Vo  (umol m-2 s-1)
    Rp / A     photorespiration relative to net assimilation

The scientific point: as gravity -> 0 the buoyant convection that sweeps the leaf
boundary layer collapses, g_bl falls (LunarLeaf-CFD: leaf 1.0 -> 0.49 mol m-2 s-1
from 1 g to microgravity), the 1/g_bl resistance term grows, Cc drops, Gamma*/Cc
rises, and Rubisco tips toward oxygenation -> more photorespiration. The effect
amplifies leaf -> rosette -> canopy exactly as g_bl falls further at denser scales.

Kinetic constants and their temperature responses are the in-vivo values of
Bernacchi et al. (2001) Plant Cell Environ. 24:253-259. Photosynthetic capacities
and conductances are mid-range Arabidopsis literature values, exposed as parameters
(they are placeholders until the stomatal/mesophyll coupling is added — see README).

Units: CO2 and O2 as mole fractions (umol mol-1); conductances mol m-2 s-1;
fluxes umol m-2 s-1; temperature degC. In mole-fraction units a drawdown across a
conductance is simply  dC [umol mol-1] = A [umol m-2 s-1] / g [mol m-2 s-1].
"""

from __future__ import annotations

import csv
import math
import os
from dataclasses import dataclass

R_GAS = 8.314  # J mol-1 K-1

# --- Rubisco kinetics at 25 C (Bernacchi et al. 2001, in vivo, mole-fraction basis) ---
KC25 = 404.9          # umol mol-1   Michaelis constant for CO2
KO25 = 278.4e3        # umol mol-1   Michaelis constant for O2 (278.4 mmol mol-1)
GAMMA_STAR25 = 42.75  # umol mol-1   CO2 compensation point w/o day respiration, at O_REF
O_REF = 210.0e3       # umol mol-1   O2 mole fraction the Gamma*25 above is defined at (21%)

# Arrhenius activation energies, J mol-1 (Bernacchi et al. 2001, 2003)
EA_KC = 79.43e3
EA_KO = 36.38e3
EA_GAMMA = 37.83e3
EA_VCMAX = 65.33e3
EA_JMAX = 43.9e3   # Bernacchi et al. 2003 (simple Arrhenius; peaked form deferred)
EA_RD = 46.39e3

# Rubisco specificity implied by the above: Sc/o = 0.5 * O_REF / Gamma*25 (mol/mol).
# Gamma* at any local O2 is then 0.5 * O_local / Sc/o  == Gamma*25 * O_local / O_REF.
SC_O_25 = 0.5 * O_REF / GAMMA_STAR25


def arrhenius(k25: float, Ea: float, Tleaf_C: float) -> float:
    """Temperature scaling of a rate/affinity constant from its 25 C value."""
    Tk = Tleaf_C + 273.15
    return k25 * math.exp(Ea * (Tk - 298.15) / (298.15 * R_GAS * Tk))


@dataclass
class LeafParams:
    """Photosynthetic capacities (25 C) and internal conductances to CO2."""
    Vcmax25: float = 90.0   # umol m-2 s-1   maximum carboxylation (Arabidopsis mid-range)
    Jmax25: float = 153.0   # umol m-2 s-1   max electron transport (~1.7 x Vcmax)
    Rd25: float = 1.35      # umol m-2 s-1   day (mitochondrial) respiration
    g_s: float = 0.20       # mol m-2 s-1    stomatal conductance to CO2
    g_m: float = 0.30       # mol m-2 s-1    mesophyll conductance to CO2
    theta: float = 0.70     # -             non-rectangular hyperbola curvature
    alpha: float = 0.30     # mol e- / mol photon   effective quantum yield of J


def electron_transport(Q: float, Jmax: float, theta: float, alpha: float) -> float:
    """Non-rectangular hyperbola J(Q): smaller root of theta*J^2 - (aQ+Jmax)J + aQ*Jmax."""
    aQ = alpha * Q
    b = aQ + Jmax
    disc = b * b - 4.0 * theta * aQ * Jmax
    return (b - math.sqrt(max(disc, 0.0))) / (2.0 * theta)


def _gamma_star(Tleaf_C: float, O_local: float) -> float:
    """Photorespiratory CO2 compensation point at leaf temperature and local O2."""
    return arrhenius(GAMMA_STAR25, EA_GAMMA, Tleaf_C) * (O_local / O_REF)


def assimilation_demand(Cc: float, O_local: float, Tleaf_C: float, Q: float,
                        p: LeafParams) -> dict:
    """Net A and its Rubisco- vs RuBP-limited components at a given Cc (the FvCB demand)."""
    Kc = arrhenius(KC25, EA_KC, Tleaf_C)
    Ko = arrhenius(KO25, EA_KO, Tleaf_C)
    Gs = _gamma_star(Tleaf_C, O_local)
    Vcmax = arrhenius(p.Vcmax25, EA_VCMAX, Tleaf_C)
    Jmax = arrhenius(p.Jmax25, EA_JMAX, Tleaf_C)
    Rd = arrhenius(p.Rd25, EA_RD, Tleaf_C)
    J = electron_transport(Q, Jmax, p.theta, p.alpha)

    drive = Cc - Gs
    Ac = drive * Vcmax / (Cc + Kc * (1.0 + O_local / Ko)) - Rd   # Rubisco-limited
    Aj = drive * J / (4.0 * Cc + 8.0 * Gs) - Rd                  # RuBP-regen-limited
    A = min(Ac, Aj)
    return {"A": A, "Ac": Ac, "Aj": Aj, "Rd": Rd, "Gamma_star": Gs, "J": J}


def photorespiration(Cc: float, A: float, Rd: float, Gamma_star: float) -> dict:
    """Rubisco carboxylation/oxygenation partitioning at the operating point.

    A = Vc - 0.5 Vo - Rd  and  0.5 Vo = Vc * Gamma*/Cc  ->  Vc = (A + Rd)/(1 - Gamma*/Cc).
    """
    frac = Gamma_star / Cc                    # = 0.5 * (Vo/Vc)
    Vo_over_Vc = 2.0 * frac
    Vc = (A + Rd) / (1.0 - frac) if Cc > Gamma_star else float("nan")
    Vo = Vc * Vo_over_Vc
    Rp = 0.5 * Vo                             # photorespiratory CO2 release
    phi = Vo / (Vc + Vo) if Vc > 0 else float("nan")
    return {"Vo_over_Vc": Vo_over_Vc, "phi": phi, "Vc": Vc, "Vo": Vo, "Rp": Rp}


def solve_operating_point(g_bl: float, p: LeafParams, Ca: float = 400.0,
                          O_excess: float = 0.0, Tleaf_C: float = 25.0,
                          Q: float = 1000.0, tol: float = 1e-6,
                          max_iter: int = 200) -> dict:
    """Solve the supply=demand fixed point for Cc given the CFD boundary-layer g_bl.

    Supply (diffusion):  Cc = Ca - A * (1/g_bl + 1/g_s + 1/g_m)
    Demand (FvCB):       A  = min(Ac, Aj)(Cc)
    O2 at Rubisco: O_REF + O_excess (the CFD surface O2 build-up; small at leaf scale
    but it raises Gamma* consistently, so both CFD gradients feed Rubisco).
    """
    r_tot = 1.0 / g_bl + 1.0 / p.g_s + 1.0 / p.g_m
    O_local = O_REF + O_excess
    Cc = 0.7 * Ca  # initial guess
    A = 0.0
    for _ in range(max_iter):
        dem = assimilation_demand(Cc, O_local, Tleaf_C, Q, p)
        A = dem["A"]
        Cc_new = Ca - A * r_tot
        Cc_new = max(Cc_new, dem["Gamma_star"] + 1e-6)  # keep Cc physical
        Cc_relaxed = 0.5 * (Cc + Cc_new)                # damped update for stability
        if abs(Cc_relaxed - Cc) < tol:
            Cc = Cc_relaxed
            break
        Cc = Cc_relaxed
    dem = assimilation_demand(Cc, O_local, Tleaf_C, Q, p)
    pr = photorespiration(Cc, dem["A"], dem["Rd"], dem["Gamma_star"])
    limiting = "Rubisco" if dem["Ac"] <= dem["Aj"] else "RuBP-regen"
    return {
        "g_bl": g_bl, "Ca": Ca, "Cc": Cc, "Ci": Ca - dem["A"] * (1.0 / g_bl + 1.0 / p.g_s),
        "A": dem["A"], "limiting": limiting, "Gamma_star": dem["Gamma_star"],
        "r_bl": 1.0 / g_bl, "r_tot": r_tot, **pr,
        "Rp_over_A": pr["Rp"] / dem["A"] if dem["A"] > 0 else float("nan"),
    }


# --- LunarLeaf-CFD boundary-layer export ------------------------------------------
# The gravity/geometry sweep is read from a CSV the CFD produces (its export_cfd.ts
# writes results/tables/T13_boundary_layer.csv; a snapshot is vendored under data/).
# Required columns: scenario, scale, gravity_g, g_bl_mol_m2_s, o2_excess_ppm.
# Any extra columns (delta_mm, Sherwood, dC_CO2_mean) are carried through for display.
DEFAULT_CFD_CSV = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "lunarleaf_gbl_sweep.csv"
)
_REQUIRED_COLS = {"scenario", "scale", "gravity_g", "g_bl_mol_m2_s", "o2_excess_ppm"}


def _g_label(gravity_g: float) -> str:
    """Human label for a gravity level (m/s^2) — falls back to the numeric value."""
    known = {9.81: "Earth", 3.71: "Mars", 1.62: "Moon", 0.0: "micro-g"}
    for g, name in known.items():
        if abs(gravity_g - g) < 0.05:
            return name
    return f"{gravity_g:.2f} m/s2"


def load_cfd_sweep(path: str = DEFAULT_CFD_CSV) -> list[dict]:
    """Parse a LunarLeaf-CFD boundary-layer export CSV into typed rows."""
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = _REQUIRED_COLS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing required column(s): {sorted(missing)}")
        rows = []
        for raw in reader:
            gravity_g = float(raw["gravity_g"])
            rows.append({
                "scenario": raw["scenario"],
                "scale": raw["scale"],
                "gravity_g": gravity_g,
                "g_ratio": gravity_g / 9.81,
                "g_label": _g_label(gravity_g),
                "g_bl": float(raw["g_bl_mol_m2_s"]),
                "o2_excess_ppm": float(raw["o2_excess_ppm"]),
                "delta_mm": float(raw["delta_mm"]) if raw.get("delta_mm") else None,
                "Sherwood": float(raw["Sherwood"]) if raw.get("Sherwood") else None,
            })
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def run_sweep(p: LeafParams | None = None, cfd_csv: str = DEFAULT_CFD_CSV, **kw) -> list[dict]:
    """Solve the FvCB operating point for every row of a CFD boundary-layer export."""
    p = p or LeafParams()
    rows = []
    for cfd in load_cfd_sweep(cfd_csv):
        r = solve_operating_point(cfd["g_bl"], p, O_excess=cfd["o2_excess_ppm"], **kw)
        r.update({k: cfd[k] for k in ("scenario", "scale", "g_label", "g_ratio",
                                      "o2_excess_ppm", "delta_mm", "Sherwood")})
        rows.append(r)
    return rows


def _fmt_table(rows: list[dict]) -> str:
    head = (f"{'scale':8} {'g':8} {'g_bl':>6} {'Cc':>7} {'Gamma*':>7} "
            f"{'Vo/Vc':>6} {'phi%':>6} {'A':>6} {'Rp':>6} {'Rp/A%':>6} {'lim':>10}")
    lines = [head, "-" * len(head)]
    for r in rows:
        lines.append(
            f"{r['scale']:8} {r['g_label']:8} {r['g_bl']:6.3f} {r['Cc']:7.1f} "
            f"{r['Gamma_star']:7.1f} {r['Vo_over_Vc']:6.3f} {100*r['phi']:6.1f} "
            f"{r['A']:6.2f} {r['Rp']:6.2f} {100*r['Rp_over_A']:6.1f} {r['limiting']:>10}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    import argparse
    import sys

    ap = argparse.ArgumentParser(description="FvCB oxygenation calculator driven by a "
                                             "LunarLeaf-CFD boundary-layer export CSV.")
    ap.add_argument("--csv", default=DEFAULT_CFD_CSV,
                    help="LunarLeaf-CFD boundary-layer export (default: data/lunarleaf_gbl_sweep.csv)")
    default_out = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "results", "tables", "T04_photorespiration_vs_gravity.csv")
    ap.add_argument("--out", default=default_out,
                    help="output CSV path")
    args = ap.parse_args()

    params = LeafParams()
    rows = run_sweep(params, cfd_csv=args.csv)
    print("FvCB oxygenation calculator — photorespiration vs gravity")
    print(f"(Ca=400 umol/mol, Tleaf=25C, Q=1000, g_s={params.g_s}, g_m={params.g_m} "
          f"mol m-2 s-1; g_bl read from {os.path.basename(args.csv)})\n")
    print(_fmt_table(rows))

    out = args.out
    cols = ["scenario", "scale", "g_label", "g_ratio", "g_bl", "Cc", "Ci", "Gamma_star",
            "Vo_over_Vc", "phi", "A", "Rp", "Rp_over_A", "limiting"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"\nwrote {out}")

    # --- self-checks (physical sanity) ---
    leaf = {r["g_label"]: r for r in rows if r["scale"] == "leaf"}
    ok = True
    def check(name, cond):
        global ok
        print(("  PASS " if cond else "  FAIL ") + name)
        ok = ok and cond
    print("\nself-checks:")
    check("Gamma* ~ 42.75 umol/mol at 25C/21% O2", abs(leaf["Earth"]["Gamma_star"] - 42.75) < 1.0)
    check("Earth leaf oxygenation fraction in 0.15-0.35", 0.15 < leaf["Earth"]["phi"] < 0.35)
    check("Cc falls Earth -> micro-g", leaf["micro-g"]["Cc"] < leaf["Earth"]["Cc"])
    check("oxygenation fraction rises Earth -> micro-g", leaf["micro-g"]["phi"] > leaf["Earth"]["phi"])
    check("Rp/A rises Earth -> micro-g", leaf["micro-g"]["Rp_over_A"] > leaf["Earth"]["Rp_over_A"])
    check("canopy micro-g most photorespiratory",
          max(rows, key=lambda r: r["Rp_over_A"])["scale"] == "canopy")
    sys.exit(0 if ok else 1)
