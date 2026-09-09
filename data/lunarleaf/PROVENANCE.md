# Vendored LunarLeaf-CFD tables

Copied verbatim from the validated solver at
`~/Documents/cose-rollout/LunarLeaf-CFD/results/tables/` so this repository is
self-contained and reproducible without that checkout.

| File | What it holds |
|---|---|
| `T1_measured_gas_exchange.csv` | measured Arabidopsis gas-exchange anchors the model is calibrated to |
| `T7_hardware_timeseries.csv` | enclosure CO₂ time series per hardware (BRIC light/dark, CARA tape, VEGGIE vented, open reference) — **all runs are microgravity** |
| `T8_enclosure_timescales.csv` | analytic sealed-enclosure timescales (a mass balance, so gravity-independent) |
| `T10_hardware_by_scale.csv` | surface gradients per hardware × leaf/rosette/canopy |
| `T11_photosynthesis_feedback.csv` | 12 h carbon gain as % of Earth: VEGGIE 100, CARA 90, BRIC 1 |
| `T13_boundary_layer.csv` | `g_bl`, δ, Sherwood across gravity × scale — the FvCB coupling variable |

LunarLeaf-CFD is a D2Q9 lattice-Boltzmann solver that passed four validation gates
(lid-driven cavity, cylinder shedding, natural convection to 0.18 % on Nu, erfc diffusion),
builds cleanly in CI, and whose boundary-layer diagnostics were reproduced by driving the
deployed application.

The sibling repository `spaceflight-plant-hardware-cfd` is **not** used anywhere in this
analysis; see `methods/CFD_PROVENANCE_CONCERN.md`.
