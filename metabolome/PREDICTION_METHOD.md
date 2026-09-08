# The predicted metabolite layer — what it is, and what it is not

**Every number in `metabolomics_*.tab` and `compound_provenance.tsv` is model output.
Nothing here was measured.** The files are labelled `PREDICTED` and this document exists so
that the label cannot be lost when the data is.

## Why a layer had to be predicted at all

NASA OSDR has no plant metabolomics. That is not an impression — it was checked
exhaustively on 2026-09-08 by pulling all **567** OSDR studies through the search API and
cross-tabulating organism against `Study Assay Measurement Type`:

- **6** studies are `metabolite profiling`: OSD-108, OSD-116, OSD-343, OSD-360, OSD-417,
  OSD-571 (a further handful carry metabolomics only in their protocol text: OSD-138,
  OSD-145, OSD-734, OSD-832, OSD-833, OSD-834).
- Every one is **mouse, human, rat or microbial**.
- **Zero of the 66 plant studies** has a metabolome, of any kind, from any hardware.

So a three-omic plant PaintOmics job cannot be assembled from OSDR alone. The third layer
is either predicted or absent.

## What the prediction is for

It is a **falsifiable hypothesis**, not a stand-in for missing data. It is generated
without ever looking at the OSD-522 omics, so the real transcriptome and proteome in the
same PaintOmics job can contradict it. `osdr/falsification_check.py` performs exactly that
test, and its result is reported whichever way it comes out.

## The chain

```
LunarLeaf-CFD          validated D2Q9 LBM solver — g_bl(gravity), enclosure CO2 mass balance
   |                   results/tables/T7, T8, T10, T13
   v
fvcb.py                Farquhar-von Caemmerer-Berry C3 model with the Rubisco oxygenation
   |                   term explicit; Bernacchi et al. (2001, 2003) in-vivo kinetics
   v
predict_metabolome.py  metabolite pool log2 fold changes
```

Only LunarLeaf-CFD is used for the gas transport. The sibling repository
`spaceflight-plant-hardware-cfd` is **not** used — see `results/CFD_PROVENANCE_CONCERN.md`.

## The physical setting

OSD-522 is BRIC-LED-001: Arabidopsis seedlings in a **sealed** canister, **illuminated**,
on SpaceX-13/14. LunarLeaf-CFD's enclosure mass balance (T8) says a lit BRIC canister draws
CO₂ from 400 ppm to near zero in about **7 minutes** — photosynthesis self-limits — and T11
puts BRIC's 12 h carbon gain at **1 %** of Earth's, against 90 % for CARA and 100 % for
vented VEGGIE.

The step that matters for interpreting OSD-522: **that drawdown is a sealed-box mass
balance, so it happens in the ground control too.** Both arms of the experiment sit near
the CO₂ compensation point. What differs between them is the boundary-layer conductance:
LunarLeaf-CFD T13 gives `g_bl` 0.540 → 0.291 mol m⁻² s⁻¹ (Earth → microgravity, rosette).

## The result we did not expect

We began from the hypothesis that CO₂ starvation would **amplify** the microgravity
boundary-layer penalty. The model says the opposite, and the model is right:

| Canister CO₂ (µmol/mol) | log2FC[A] flight/ground | log2FC[Vo] flight/ground |
|---:|---:|---:|
| 400 | −0.0959 | +0.0213 |
| 250 | −0.1049 | +0.0126 |
| 150 | −0.1110 | +0.0063 |
| 100 | −0.1140 | +0.0031 |
|  70 | −0.1158 | +0.0011 |

As assimilation falls toward the compensation point, the flux through the boundary layer
falls with it, so the CO₂ drop across that layer (`A/g_bl`) tends to **zero** and the two
gravities *converge* on oxygenation fraction. The gravity effect on photorespiratory
partitioning shrinks rather than amplifying.

What survives is a **persistent proportional assimilation penalty** — about 5 % at leaf
scale, 7 % at rosette, 9 % at canopy — that barely moves across the whole CO₂ range. That
robustness is the prediction's strength: it does not depend on knowing the canister's exact
steady-state CO₂, which nobody measured.

## Two contrasts, differing by a factor of 25

| | log2FC photosynthate | log2FC photorespiratory | what it is |
|---|---:|---:|---|
| `FLT_vs_GC` | −0.114 | +0.003 | what OSD-522 measured; goes into the joint PaintOmics job |
| `BRIC_vs_VENTED` | −2.808 | +0.242 | the hardware effect; a metabolome-only second job |

The hardware effect is **25× larger** than the gravity effect — but it hits flight and
ground alike, so it very nearly cancels out of the contrast the experiment can see. This is
the central claim of the analysis: the sealed canister dominates the plant's carbon
metabolism, and the experiment is structurally blind to it.

## How a flux becomes a predicted pool

One assumption, stated once and applied uniformly. For a pool `P` fed by flux `J` and
drained by first-order consumption of capacity `k`, quasi-steady state gives `P ≈ J/k`. Over
a flight-versus-ground contrast the consuming capacity is taken as unchanged, so

```
log2(P_flight / P_ground)  =  log2(J_flight / J_ground)
```

Pools therefore inherit the log2 ratio of whichever modelled flux feeds them.

**Tier 1** — pool tied directly to a flux the model computes:

- C2 photorespiratory intermediates (2-phosphoglycolate, glycolate, glyoxylate, glycine,
  serine, hydroxypyruvate, glycerate, ammonia) follow **Vo**.
- Net photosynthate (3-phosphoglycerate, sucrose, glucose, fructose, maltose, starch)
  follows **A**.

**Tier 2** — one further documented step from a modelled quantity:

- Ribulose-1,5-bisphosphate follows **−log2(Cc)**: it accumulates as Rubisco becomes
  substrate-limited rather than RuBP-limited.
- Carbon-starvation pools (asparagine, valine, leucine, isoleucine, glutamine, glutamate)
  follow **−log2(A)**. The `--starvation-gain` parameter that scales this defaults to
  **1.0**, so no unfitted parameter enters the prediction.

## What is deliberately not predicted

Stated so that the omissions read as decisions rather than oversights.

- **TCA intermediates** (citrate, malate, succinate, fumarate, 2-oxoglutarate) — the model
  does not resolve respiratory flux, and carbon limitation can deplete these pools by
  oxidation or raise them by anaplerosis. The direction is genuinely uncertain, so no value
  is asserted.
- **Fermentation products** (ethanol, lactate, alanine) — fermentation requires hypoxia.
  T8 puts BRIC O₂ depletion at 6.5 days and only in **dark** canisters. OSD-522 is
  illuminated, so photosynthesis is releasing O₂ and no hypoxia is predicted.
- **The glycine/serine ratio** — both pools are tied to the same flux `Vo`, so under this
  model the ratio cannot move. Presenting it as a photorespiration readout would be
  circular.

## Known limitations

- `Vcmax25`, `Jmax25`, `Rd25`, `g_s` and `g_m` in `fvcb.py` remain documented mid-range
  Arabidopsis literature values, not fits to spaceflight material.
- `g_s` is fixed; a real stomatal response to CO₂ and humidity would modify the conductance
  chain, and BRIC canisters are humid.
- Pools integrate flux over time, so measured changes could exceed these instantaneous
  quasi-steady-state values. The prediction is deliberately not scaled up to account for
  that, because the integration time is unconstrained.
- The 21-compound set is small by metabolomics standards. It is confined to compounds
  reachable from a modelled flux; breadth was traded for defensibility.
- **The prediction has only four distinct values.** Because every pool in a block inherits
  the same flux ratio, all eight C2 intermediates carry +0.0031, all six photosynthate
  pools −0.1140, and so on. That is a faithful statement of what the model knows — it
  resolves fluxes, not individual pool kinetics — but it has two consequences worth
  anticipating. PaintOmics' metabolite **class activity** test may behave degenerately on
  it, and the **hub analysis** will rank compounds by their gene neighbourhoods rather than
  by any variation in our values. Distinguishing pools within a block would need
  compound-specific turnover constants that this model does not contain and that would have
  to be invented to supply.

## Reproducing

```bash
python3 metabolome/predict_metabolome.py --sensitivity
```

KEGG compound identifiers were retrieved from `rest.kegg.jp` on 2026-09-08 and are recorded
per compound in `compound_provenance.tsv`. They were not written from memory.
