# Vendored LunarLeaf-CFD tables

Copied verbatim from the validated solver's `results/tables/` (T7, T10 and T13 regenerated;
see below) so this repository is self-contained and reproducible without that checkout.

| File | What it holds |
|---|---|
| `T1_measured_gas_exchange.csv` | measured Arabidopsis gas-exchange anchors the model is calibrated to |
| `T7_hardware_timeseries.csv` | enclosure CO₂ time series per hardware (BRIC light/dark, CARA light/dark, VEGGIE vented, open reference) — **all runs are microgravity** |
| `T8_enclosure_timescales.csv` | analytic sealed-enclosure timescales (a mass balance, so gravity-independent) |
| `T10_hardware_by_scale.csv` | surface gradients per hardware × leaf/rosette/canopy, with each run's peak Mach number (`ma_max`) |
| `T11_photosynthesis_feedback.csv` | 12 h carbon gain as % of Earth: VEGGIE 100, CARA 90, BRIC 1 |
| `T13_boundary_layer.csv` | `g_bl`, δ, Sherwood across gravity × scale — the FvCB coupling variable |

LunarLeaf-CFD is a D2Q9 lattice-Boltzmann solver that passed four validation gates
(lid-driven cavity, cylinder shedding, natural convection to 0.18 % on Nu, erfc diffusion),
builds cleanly in CI, and whose boundary-layer diagnostics were reproduced by driving the
deployed application.

## Which LunarLeaf-CFD version this is

The checkout lives at `~/Documents/LunarLeaf-CFD`. Each file here is byte-identical (same
git blob id) to `results/tables/` at these LunarLeaf commits:

| File | LunarLeaf commit |
|---|---|
| `T1`, `T8`, `T11` | identical to LunarLeaf HEAD (`fab6421`) |
| `T13_boundary_layer.csv`, and its copy `../lunarleaf_gbl_sweep.csv` | `e298785` (2026-10-07), unchanged through HEAD `fab6421` (2026-10-13) |
| `T10_hardware_by_scale.csv` | `e298785` (2026-10-07), unchanged through HEAD `fab6421` |
| `T7_hardware_timeseries.csv` | `e298785` (2026-10-07), unchanged through HEAD `fab6421` |

### T13 — re-vendored 2026-10-28

`e298785` re-ran the CFD exports to steady state (150 k steps; the 0 g canopy 900 k) after
the earlier 22–30 k-step snapshots were found to catch the 0 g scenes mid-relaxation, and
added a `gbl_drift_pct_last_5s` column. On 2026-10-28 the export was regenerated rather than
copied: `validation/export_cfd.ts` was run (`list | xargs … job`, then `merge`) on a clean
`git archive` of LunarLeaf HEAD `fab6421`, and its `T13_boundary_layer.csv` came out
byte-identical to the committed one (blob `399e184`), with every row's drift 0.00 %.
That file replaced both `T13_boundary_layer.csv` here and `../lunarleaf_gbl_sweep.csv`
(read by `scripts/fvcb.py` `load_cfd_sweep()`), and `scripts/run_all.sh` was re-run.

`g_bl` (mol m⁻² s⁻¹), previous snapshot (`87b81cf`) → current:

| Scale | Earth | microgravity |
|---|---|---|
| leaf | 1.000 → 0.997 | 0.494 → 0.443 |
| rosette | 0.540 → 0.540 | 0.291 → 0.272 |
| canopy | 0.148 → 0.133 | 0.109 → 0.058 |

Every `g_bl`-dependent number in the repository now comes from this export: the operating
points (`T04`, `T05`), the irradiance sensitivity and assimilation deficits (`T13`), and the
predicted metabolite layer (`T06`, `data/paintomics_upload/metabolomics_*values.tab`).
**The PaintOmics jobs (`m1z16Qg3DK`, `EVFahRGk7T`) were run before this update**, on the layer
predicted from the `87b81cf` snapshot; see `methods/PAINTOMICS_RESULTS.md`.

### T7 and T10 — re-vendored 2026-10-28

`e298785` also regenerated T7 and T10 (its pressure outlet for fan-driven scenes changes the
vented VEGGIE runs). On 2026-10-28 both were regenerated rather than copied, on a clean
`git archive` of LunarLeaf HEAD `fab6421`: `validation/hardware_sim.ts` for T7, and
`validation/scales.ts` (`list | xargs … job`, then `merge`) for T9/T10, each bundled with
esbuild as in LunarLeaf's `results/README.md` "Reproduce". T7 came out byte-identical to the
committed one (blob `2cda916`), T10 too (blob `dd70fd2`), and so did T9 (not vendored). The
LunarLeaf commits after `e298785` therefore do not change these tables.

T7 (previous snapshot `b3567f6` → current): only the `VEGGIE vented` `surf_co2_mean` series
changes, by at most 0.025, with a final value of −0.112 → −0.087. The BRIC, CARA and
open-reference series, and every `dishmean_co2` value, are unchanged. T7 feeds Fig. 1b
(`scripts/figures/fig01_hardware_atmosphere.py`) and the site's enclosure chart
(`docs/data/enclosure_timeseries.json`), not any `g_bl` number.

T10 (previous snapshot `32f27de` → current): only the VEGGIE rows change, and a peak-Mach
column `ma_max` is added. `surf_dC_H2O` / `surf_dC_CO2` for VEGGIE: leaf 0.078 / 0.101 →
0.068 / 0.087, rosette 0.236 / 0.322 → 0.213 / 0.294, canopy 0.252 / 0.314 → 0.286 / 0.406.
The canopy VEGGIE run peaks at Ma 0.32, above the 0.3 limit, and LunarLeaf's own `merge` step
flags it "treat with caution". T10 is not read by any script here.

The sibling repository `spaceflight-plant-hardware-cfd` is **not** used anywhere in this
analysis; see `methods/CFD_PROVENANCE_CONCERN.md`.
