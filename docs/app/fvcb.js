/* fvcb.js — a faithful port of scripts/fvcb.py, so the model can be driven in a browser.
 *
 * The Python is the reference implementation. This file must agree with it to 1e-6 across a
 * parameter grid; scripts/check_js_parity.py enforces that, and any change here should be
 * made in both places or the check will fail.
 *
 * Kinetics are the in-vivo values of Bernacchi et al. (2001, 2003), mole-fraction basis.
 */

export const R_GAS = 8.314;          // J mol-1 K-1

export const KC25 = 404.9;           // umol mol-1
export const KO25 = 278.4e3;         // umol mol-1
export const GAMMA_STAR25 = 42.75;   // umol mol-1
export const O_REF = 210.0e3;        // umol mol-1  (21 % O2)

export const EA_KC = 79.43e3;
export const EA_KO = 36.38e3;
export const EA_GAMMA = 37.83e3;
export const EA_VCMAX = 65.33e3;
export const EA_JMAX = 43.9e3;
export const EA_RD = 46.39e3;

export const DEFAULT_PARAMS = {
  Vcmax25: 90.0,   // umol m-2 s-1
  Jmax25: 153.0,   // umol m-2 s-1
  Rd25: 1.35,      // umol m-2 s-1
  g_s: 0.20,       // mol m-2 s-1
  g_m: 0.30,       // mol m-2 s-1
  theta: 0.70,
  alpha: 0.30,
};

export function arrhenius(k25, Ea, TleafC) {
  const Tk = TleafC + 273.15;
  return k25 * Math.exp((Ea * (Tk - 298.15)) / (298.15 * R_GAS * Tk));
}

/* Non-rectangular hyperbola J(Q): the smaller root. */
export function electronTransport(Q, Jmax, theta, alpha) {
  const aQ = alpha * Q;
  const b = aQ + Jmax;
  const disc = b * b - 4.0 * theta * aQ * Jmax;
  return (b - Math.sqrt(Math.max(disc, 0.0))) / (2.0 * theta);
}

export function gammaStar(TleafC, Olocal) {
  return arrhenius(GAMMA_STAR25, EA_GAMMA, TleafC) * (Olocal / O_REF);
}

/* The FvCB demand function: net A at a given chloroplast CO2. */
export function assimilationDemand(Cc, Olocal, TleafC, Q, p) {
  const Kc = arrhenius(KC25, EA_KC, TleafC);
  const Ko = arrhenius(KO25, EA_KO, TleafC);
  const Gs = gammaStar(TleafC, Olocal);
  const Vcmax = arrhenius(p.Vcmax25, EA_VCMAX, TleafC);
  const Jmax = arrhenius(p.Jmax25, EA_JMAX, TleafC);
  const Rd = arrhenius(p.Rd25, EA_RD, TleafC);
  const J = electronTransport(Q, Jmax, p.theta, p.alpha);

  const drive = Cc - Gs;
  const Ac = (drive * Vcmax) / (Cc + Kc * (1.0 + Olocal / Ko)) - Rd;
  const Aj = (drive * J) / (4.0 * Cc + 8.0 * Gs) - Rd;
  return { A: Math.min(Ac, Aj), Ac, Aj, Rd, Gamma_star: Gs, J };
}

// Rubisco partitioning at the operating point.
//   A = Vc - 0.5 Vo - Rd  and  0.5 Vo = Vc (Gamma* / Cc)
//   so Vc = (A + Rd) / (1 - Gamma* / Cc).
// Line comments, not a block comment: "Gamma*/Cc" contains the block-comment terminator
// and silently ended the comment two lines early.
export function photorespiration(Cc, A, Rd, Gamma_star) {
  const frac = Gamma_star / Cc;
  const VoOverVc = 2.0 * frac;
  const Vc = Cc > Gamma_star ? (A + Rd) / (1.0 - frac) : NaN;
  const Vo = Vc * VoOverVc;
  return {
    Vo_over_Vc: VoOverVc,
    phi: Vc > 0 ? Vo / (Vc + Vo) : NaN,
    Vc, Vo,
    Rp: 0.5 * Vo,
  };
}

/* Solve supply = demand for Cc, given the CFD boundary-layer conductance.
 *   supply:  Cc = Ca - A * (1/g_bl + 1/g_s + 1/g_m)
 *   demand:  A  = min(Ac, Aj)(Cc)
 * Damped fixed point, identical to the Python. */
export function solveOperatingPoint(gBl, p = DEFAULT_PARAMS, opts = {}) {
  const { Ca = 400.0, O_excess = 0.0, TleafC = 25.0, Q = 1000.0,
          tol = 1e-6, maxIter = 200 } = opts;
  const rTot = 1.0 / gBl + 1.0 / p.g_s + 1.0 / p.g_m;
  const Olocal = O_REF + O_excess;

  let Cc = 0.7 * Ca;
  for (let i = 0; i < maxIter; i++) {
    const dem = assimilationDemand(Cc, Olocal, TleafC, Q, p);
    let CcNew = Ca - dem.A * rTot;
    CcNew = Math.max(CcNew, dem.Gamma_star + 1e-6);
    const relaxed = 0.5 * (Cc + CcNew);
    if (Math.abs(relaxed - Cc) < tol) { Cc = relaxed; break; }
    Cc = relaxed;
  }

  const dem = assimilationDemand(Cc, Olocal, TleafC, Q, p);
  const pr = photorespiration(Cc, dem.A, dem.Rd, dem.Gamma_star);
  return {
    g_bl: gBl, Ca, Cc,
    Ci: Ca - dem.A * (1.0 / gBl + 1.0 / p.g_s),
    A: dem.A,
    limiting: dem.Ac <= dem.Aj ? "Rubisco" : "RuBP-regen",
    Gamma_star: dem.Gamma_star,
    r_bl: 1.0 / gBl, r_tot: rTot,
    ...pr,
    Rp_over_A: dem.A > 0 ? pr.Rp / dem.A : NaN,
  };
}
