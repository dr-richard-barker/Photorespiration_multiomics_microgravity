# Vendored LunarLeaf-CFD tables

Copied verbatim from the validated solver at
`~/Documents/cose-rollout/LunarLeaf-CFD/results/tables/` so this repository is
self-contained and reproducible without that checkout.

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

The checkout now lives at `~/Documents/LunarLeaf-CFD`. Checked 2026-10-27 by comparing
git blob ids, each file here is byte-identical to `results/tables/` at these LunarLeaf
commits:

| File | LunarLeaf commit |
|---|---|
| `T1`, `T8`, `T11` | identical to LunarLeaf HEAD |
| `T13_boundary_layer.csv` | `87b81cf` (2026-07-29) |
| `T10_hardware_by_scale.csv` | `32f27de` (2026-07-23) |
| `T7_hardware_timeseries.csv` | `b3567f6` (2026-09-21, adds the CARA dark run) |

**This snapshot is out of date.** LunarLeaf `e298785` (2026-10-07) re-ran the exports to
steady state and regenerated T7, T10 and T13. Microgravity `g_bl` fell: leaf 0.494 → 0.443,
rosette 0.291 → 0.272, canopy 0.109 → 0.058 mol m⁻² s⁻¹; on Earth leaf and rosette barely moved but canopy fell 0.148 → 0.133. Every
number in this repository that depends on `g_bl` — the 5 / 7 / 9 % microgravity
assimilation penalty, the operating-point table, the predicted metabolome — still uses the
older values. Re-vendoring the new tables and re-running `scripts/run_all.sh` would change
those numbers; that has not been done.

The sibling repository `spaceflight-plant-hardware-cfd` is **not** used anywhere in this
analysis; see `methods/CFD_PROVENANCE_CONCERN.md`.
