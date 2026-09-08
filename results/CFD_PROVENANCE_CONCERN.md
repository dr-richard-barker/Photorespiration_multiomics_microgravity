# Why this analysis does not use `spaceflight-plant-hardware-cfd`

Recorded 2026-09-08, during the OSD-522 PaintOmics work. **No action taken on that
repository** — this note exists so the decision made here is traceable, and so the concern
reaches someone who can check it.

## The concern

`~/Documents/spaceflight-plant-hardware-cfd` would have been the natural source of
hardware-resolved boundary-layer conductance: its `manuscript/tables/data/` holds
`T3_canopy_conductance_stagnation.csv` (`g_bl` per hardware, per gravity, per fan mode),
`T5_fan_failure_hypoxia.csv` and `T6_cara_bric_petridish_microclimates.csv`, covering
VEGGIE, APH, CHROMEX, CARA and BRIC.

Three things do not line up:

1. **No solver output exists anywhere in the repository.** No `runs/`, no
   `postProcessing/`, no solver logs, no residuals, no time directories. What is present is
   OpenFOAM *case templates* (`templates/templates_{aph,veggie,chromex}/` with
   `controlDict`, `fvSchemes`, `snappyHexMeshDict`) and STL geometry under `cases/` — the
   inputs to a run, not the outputs of one.

2. **The figure scripts do not plot solver results — they plot analytic expressions.**
   In `manuscript/figures/fig9_fan_failure.py`:

   ```python
   t = np.linspace(0, 20, 200)
   u_aph = 0.60 * np.exp(-t / 4.8)
   u_veg = 0.15 * np.exp(-t / 2.4)
   u_chr = 0.008 * np.exp(-t / 0.8)
   ```

   and, for the conductance collapse panel:

   ```python
   g_1g = 0.362 + (1.071 - 0.362) * np.exp(-t_min / 0.2)
   ```

   The constants (`0.362`) are the same values that appear in the T3/T5 tables. The tables
   and the figures are hard-coded together; neither is read from a simulation.

3. **The supplement nevertheless claims a convergence study.**
   `manuscript/npj_supplementary.tex:19`:

   > "Grid convergence index (GCI) analysis verified asymptotic convergence with Richardson
   > extrapolation error $<1.8\%$ across all five computational hardware domains."

   A GCI analysis requires at least three systematically refined meshes and their solutions.
   None are in the repository.

## What this does not establish

The runs may have been performed on a cluster with only the summary tables committed. That
would explain the absence of output while leaving the numbers sound. **This note is a flag,
not a verdict** — the check that settles it is whether the mesh-refinement solutions behind
the GCI claim still exist somewhere.

## What was done instead

This analysis uses **only** `cose-rollout/LunarLeaf-CFD`, which is genuinely validated: four
passed validation gates (lid-driven cavity, cylinder shedding, natural convection to 0.18 %
on Nu, erfc diffusion), a CI-verified build, and boundary-layer diagnostics reproduced by
driving the deployed application. Its `results/tables/T7, T8, T10, T13` supply everything
needed — enclosure CO₂ mass balance and `g_bl` by gravity and scale — so nothing was lost by
leaving the other repository out.

## Related

Two earlier findings in this portfolio have the same shape, which is why this one is
written down rather than assumed benign:

- `aeroleaf-cfd` presents generated noise as force traces; its README now says so.
- A fabricated npj DOI propagated through a manuscript template into several repositories.
