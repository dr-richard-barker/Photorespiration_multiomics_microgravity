# Vendored LunarLeaf-CFD tables

Copied verbatim from the validated solver's `results/tables/` (T13 regenerated; see below)
so this repository is self-contained and reproducible without that checkout.

| File | What it holds |
|---|---|
| `T1_measured_gas_exchange.csv` | measured Arabidopsis gas-exchange anchors the model is calibrated to |
| `T7_hardware_timeseries.csv` | enclosure CO₂ time series per hardware (BRIC light/dark, CARA light/dark, VEGGIE vented, open reference) — **all runs are microgravity** |
| `T8_enclosure_timescales.csv` | analytic sealed-enclosure timescales (a mass balance, so gravity-independent) |
| `T10_hardware_by_scale.csv` | surface gradients per hardware × leaf/rosette/canopy |
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
| `T10_hardware_by_scale.csv` | `32f27de` (2026-07-23) — **older than `e298785`** |
| `T7_hardware_timeseries.csv` | `b3567f6` (2026-09-21, adds the CARA dark run) — **older than `e298785`** |

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

### T7 and T10 — still the older snapshot

`e298785` also regenerated T7 and T10; they have not been re-vendored here. Checked
2026-10-28: the vendored T7 differs from LunarLeaf HEAD only in the `VEGGIE vented`
`surf_co2_mean` series (the fan-outlet fix in `e298785`; maximum absolute difference 0.025,
final value −0.112 → −0.087); the BRIC, CARA and open-reference series are
identical. T7 feeds Fig. 1b and the site's enclosure chart, not any `g_bl` number. T10 is not
read by any script here.

The sibling repository `spaceflight-plant-hardware-cfd` is **not** used anywhere in this
analysis; see `methods/CFD_PROVENANCE_CONCERN.md`.
