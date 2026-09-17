# What this analysis implies for future spaceflight plant experiments

Drawn from what the OSD-522 analysis actually produced, not from a wish list. Each item
states the finding it rests on, so the reasoning can be checked or rejected.

---

## The finding these recommendations rest on

A validated CFD solver and a Farquhar–von Caemmerer–Berry model, run without reference to
any omics data, predicted the OSD-522 shoot response — and the real transcriptome and
proteome then confirmed it:

| Prediction (model only) | Measured (OSD-522) | |
|---|---|---|
| carbon fixation genes down | median −0.082, p = 1.3e-4 | ✔ |
| starch/sucrose genes down | median −0.142, p = 3.5e-6 | ✔ |
| photosystem apparatus down | −0.215, p = 3.8e-6 (protein −0.92, p = 5.0e-4) | ✔ |
| Rubisco down | −0.381, p = 4.5e-3 | ✔ |
| carbon-starvation (DIN) markers up | +0.708, p = 2.6e-3 | ✔ |
| photorespiratory enzymes **not** induced | no shift, p = 0.09 | ✔ |
| photorespiration above Rubisco | +0.239, p = 6.8e-3 (protein +0.285, p = 0.049) | ✔ |
| no fermentation or hypoxia | no shift | ✔ |

And the modelling result that reframes it: in a **sealed** canister the enclosure effect on
carbon metabolism is **25× the microgravity effect** — but it acts on the flight and ground
arms alike, so it very nearly cancels out of the contrast the experiment is able to see.

**The plants were carbon-starved by their hardware, in both arms, and the experiment was
structurally blind to it.**

---

## 1. Put a CO₂ sensor inside the canister

The single highest-value change, and the cheapest.

LunarLeaf-CFD's enclosure mass balance draws a lit BRIC canister from 400 ppm to near zero
in about **7 minutes**, leaving the plants at the CO₂ compensation point for the rest of the
photoperiod, with 12 h carbon gain at 1 % of Earth's. No BRIC-class experiment has ever
logged the atmosphere its plants actually experienced. Every interpretation of those
datasets — including this one — currently rests on a modelled gas composition.

A logging CO₂/O₂ sensor in the canister headspace would convert the central assumption of
this analysis into a measurement. It needs no new hardware qualification path of its own
beyond the sensor, and it retrospectively strengthens or corrects an entire archive of
BRIC-derived spaceflight transcriptomes.

## 2. Measure a plant metabolome in space — any plant metabolome

There is not one. Of 567 OSDR studies, 6 are metabolite profiling, and all 6 are mouse,
human, rat or microbial; of 66 plant studies, none has a metabolome. This is why the third
layer of this PaintOmics job had to be predicted rather than downloaded.

The gap does not need a new flight to start closing. A **targeted GC-MS or LC-MS panel of
20–30 compounds** — the C2 photorespiratory intermediates, soluble sugars, starch, and the
starvation amino acids — run on **archived material from flights already flown** would be
the first plant spaceflight metabolome in existence. The 21-compound predicted set in
`results/tables/T06_predicted_compounds.tsv` is a ready-made target list, and each compound comes
with a quantitative prediction to test against.

## 3. Fly sealed and vented hardware side by side

This analysis separates hardware effects from microgravity effects **by model**. Nobody has
separated them **by measurement**.

The experiment: one genotype, one duration, one harvest protocol, flown simultaneously in
sealed BRIC-class canisters and in vented VEGGIE- or APH-class chambers, with transcriptome
and metabolome from both. The model predicts the two will diverge far more from each other
than either diverges from its own ground control — a 25× ratio is a large enough target to
survive considerable error in the estimate.

If that prediction holds, a substantial part of the published "spaceflight plant response"
is a hardware atmosphere response, and the field's comparative meta-analyses need to
stratify by enclosure before they pool.

## 4. Run the cheap ground falsification first

Most of item 3 is testable at 1 g, this year, without flying anything: grow the same
genotype in sealed BRIC-geometry canisters versus vented controls, on the ground, and
sequence them.

The model says the sealed-versus-vented contrast is overwhelmingly an enclosure effect, so
**most of the signature should reproduce on the ground**. If it does, that is a strong and
inexpensive result that reframes how the archive should be read. If it does not, the model
is missing something specific to flight, and finding out which part is itself the result.
This is the item to do first, because it is the one that can fail cheaply.

## 5. Reconsider the ground control

Because the sealed enclosure dominates and affects both arms, the flight-versus-ground
contrast subtracts away the largest thing happening to the plants. A ground control in
*identical sealed hardware* is the right control for isolating microgravity — and it is what
BRIC-LED-001 correctly used — but it is the wrong control for asking what the plants
experienced. Future designs would benefit from a **third arm**: vented hardware on the
ground, giving a physiologically unstressed baseline alongside the hardware-matched one.

With three arms, the hardware effect and the gravity effect become separately estimable
instead of confounded.

## 6. Do not use photorespiratory gene induction as a photorespiration readout

A methodological caution that falls directly out of the result. Photorespiratory
partitioning in the sealed canister is predicted to be very high — oxygenation fraction
rising from ~27 % to ~58 % as CO₂ is drawn down — while photorespiratory **gene expression**
does not move at all. Flux is set by substrate availability at Rubisco; the enzymes are
already abundant and constitutive.

Transcriptomics alone therefore cannot detect this state. That is a further argument for
item 2, and a caution against reading a flat photorespiratory transcript profile as evidence
that photorespiration is unchanged.

---

## What would refute the model

Stated so the argument is falsifiable rather than merely consistent:

- an in-canister CO₂ log showing the atmosphere stays near ambient — item 1 settles it;
- a measured metabolome showing photorespiratory intermediates strongly changed between
  flight and ground, where the model predicts +0.003 log2;
- a sealed-versus-vented contrast no larger than the flight-versus-ground contrast, against
  the predicted 25×.
