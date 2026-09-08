# Job 2 (hardware contrast) — a negative result, reported as one

**Job `EVFahRGk7T`** · organism `ath` · KEGG + MapMan · AI interpretation off
Input: `metabolomics_hardware_values.tab` only — 21 predicted compounds, `BRIC_vs_VENTED`.

## What came back

| | |
|---|---|
| KEGG pathways found | 57 |
| KEGG pathways significant | **0** |
| MapMan pathways found | **0** |
| Metabolite hub analysis | **not offered** |
| Metabolite class activity | not offered (needs replicates) |

## Why, and why it is not a failure of the submission

The job ran correctly — all 21 compounds registered and mapped. PaintOmics simply has
almost nothing to work with:

1. **The hub analysis needs a gene layer.** Its question is "which metabolites have
   differentially expressed *genes* concentrated around them in the network". With no
   transcriptome there are no DE neighbours to count, so the analysis does not appear.
2. **The enrichment has no contrast.** All 21 compounds are in the relevant-features file,
   because the model predicts every one of them changes between sealed and vented hardware.
   A relevant list identical to the input list cannot be over-represented against itself.
3. **MapMan is keyed on genes**, so a compound-only submission maps to none of its bins —
   the 0 is expected, not an error.

## What this means

The 25× hardware effect is a result of the **model**, and it is in
`metabolome/compound_provenance.tsv` and `metabolome/operating_points.tsv`
(photosynthate −2.81 log2FC, photorespiratory +0.24). PaintOmics adds nothing to it, and
this job should not be presented as independent support for it.

Testing the hardware contrast properly needs **measured omics from two hardware types**,
which is recommendation 3 in `FUTURE_EXPERIMENTS.md`. Until such data exists, a three-layer
job like job 1 is the only form in which PaintOmics can say anything about this question.
