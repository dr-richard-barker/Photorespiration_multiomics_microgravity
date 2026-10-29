# Metabolite class activity analysis — job oK2h66740Q

(Resubmission of `m1z16Qg3DK`; every class row below is identical in both job records.)

It ran despite no experimental-design file: with no replicates PaintOmics falls back to a
**binomial test of the relevant list against the threshold you declared** (p₀ = 0.05),
rather than the replicate-based class test.

| Class | Category | Members in relevant list | FDR (BH) |
|---|---|---:|---:|
| Amino acids | Peptides | 6/8 | **2.0e-6** |
| Monosaccharides | Carbohydrates | 2/2 | 0.0042 |
| Oligosaccharides | Carbohydrates | 2/2 | 0.0042 |
| Neurotransmitters | Hormones and transmitters | 1/2 | 0.1219 |
| Carboxylic acids | Organic acids | 0/1 | 1.0000 |

"1 of 5 classes at level 2 pass BH < 0.05 · 4 with fewer than 3 members (descriptive)."
Eight of the 21 compounds fall in no class — all matched a KEGG compound, but these eight are
not in KEGG BRITE: 2-phosphoglycolate, glycolate, hydroxypyruvate, D-glycerate, ammonia,
3-phospho-D-glycerate, RuBP and starch. PaintOmics' own `excluded.unclassified` list is longer
(10 names in the first run, 11 in this one) because it also carries capitalised duplicates —
Glyoxylate, Maltose and, this time, D-Glucose — of compounds it classifies under their
lower-case names. An earlier version of this note counted that list as 10 unclassified
compounds.

## Read this one with care

**This result is close to circular and should not be cited as evidence for the model.** The
"relevant list" it tests is our own prediction — the compounds the model says changed. So
the analysis is reporting that *our predicted set is concentrated in amino acids*, which is
a property of how the prediction was constructed (the tier-2 carbon-starvation block is six
amino acids), not an independent finding about the plants.

It becomes a real test only when the metabolite layer is measured. That is
recommendation 2 in `FUTURE_EXPERIMENTS.md`.

The genuinely informative metabolomics output from this job is the **hub analysis**
(`results/tables/T10_metabolite_hubs.tsv`), which is *not* circular: it ranks our compounds by the density
of **measured, real** differentially expressed genes around them in the KEGG network.
