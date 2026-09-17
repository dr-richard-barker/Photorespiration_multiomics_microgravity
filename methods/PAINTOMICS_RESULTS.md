# PaintOmics AI results — OSD-522 (BRIC-LED-001)

Two jobs were run on <https://paintomics.org> (Supercomputador Drago, CSIC), organism
`ath`, databases **KEGG + MapMan**, **AI interpretation off** (no data went to `llm.iiia.es`).

| Job | ID | Contents | Link |
|---|---|---|---|
| 1 | `m1z16Qg3DK` | transcriptome + proteome (real, OSDR) + metabolome (predicted) | <https://paintomics.org/?jobID=m1z16Qg3DK> |
| 2 | `EVFahRGk7T` | predicted metabolome only, sealed vs vented hardware | <https://paintomics.org/?jobID=EVFahRGk7T> |

## Identifier mapping — all three layers landed cleanly

| Omic | KEGG | MapMan |
|---|---|---|
| Gene expression (20,983) | 20,456 · 98% | 20,983 · 100% |
| Proteomics (5,160) | 5,002 · 97% | 4,998 · 97% |
| Metabolomics (21) | 21 · 100% | 21 · 100% |

AGI locus codes went straight into KEGG `ath` with no conversion, UniProt accessions were
converted by PaintOmics itself, and all 21 predicted compound **names** resolved — the
choice to name compounds by KEGG's own primary synonym paid off.

## The headline: the model's central prediction holds

The model said photorespiratory partitioning would be untouched by flight while
photosynthesis and carbon supply fell. The enrichment table says exactly that.

| Pathway | Combined p (Fisher) | Model predicted |
|---|---:|---|
| **Photosynthesis** | **0.00797** | down ✔ |
| **Starch and sucrose metabolism** | **0.01857** | down ✔ |
| **Photosynthesis – antenna proteins** | **0.02174** | down ✔ |
| **photosynthesis** (MapMan) | **0.03531** | down ✔ |
| Carbon fixation by Calvin cycle | 0.64351 | down ✘ (see below) |
| **Glyoxylate and dicarboxylate metabolism** (photorespiration) | **0.99367** | no change ✔ |

232 pathways tested; 32 significant at combined p < 0.05. **Photorespiration ranks 214th of
232** — in the bottom 8 %, and the least enriched of any pathway the model made a claim
about, which is the pathway it specifically predicted would not move. It is *not* last:
eighteen pathways rank below it, nine of them at p = 1.

*This corrects an earlier reading of this table.* "Last of 231" was written from the web
interface and was wrong in both numbers. The full job record, recovered afterwards with
`POST /pa_recover_job` and vendored at
`results/paintomics_raw/job1_m1z16Qg3DK_full.json.gz`, carries 232 pathway records, all
with matched features, and puts glyoxylate and dicarboxylate metabolism at p = 0.99367,
rank 214. The p-value itself was transcribed correctly; the rank and the total were not.

*The Calvin-cycle miss is a difference of test, not of direction.* PaintOmics tests whether a
pathway is **over-represented in the relevant-features list**; our own gene-set test
(`results/tables/T07_falsification_osd522.tsv`) asks whether the pathway is **shifted** against
background, and there KEGG ath00710 is significantly down (p = 1.3e-4). A pathway can be
coherently shifted without being enriched among the FDR < 0.05 genes.

## Metabolite hub analysis — the strongest independent result

This asks which of our compounds have **measured, real** differentially expressed genes
concentrated around them in the KEGG network. It is not circular: the metabolite values are
predicted, but the genes counted around them are the OSDR transcriptome.

7 of 15 compounds reach FDR < 0.05, and every one is a carbon-starvation or carbon-supply
node:

| Compound | KEGG | FDR | DE neighbours |
|---|---|---:|---:|
| L-Isoleucine | C00407 | 0.003 | 112 |
| D-Ribulose 1,5-bisphosphate | C01182 | 0.003 | 40 |
| L-Glutamate | C00025 | 0.003 | 215 |
| L-Leucine | C00123 | 0.003 | 112 |
| L-Valine | C00183 | 0.003 | 112 |
| 3-Phospho-D-glycerate | C00197 | 0.004 | 81 |
| L-Glutamine | C00064 | 0.029 | 279 |

**Not one of the eight C2 photorespiratory intermediates** — 2-phosphoglycolate, glycolate,
glyoxylate, glycine, serine, hydroxypyruvate, glycerate, ammonia — is a significant hub,
although all eight mapped. The branched-chain amino acids topping the list is the signature
of BCAA catabolism feeding the mitochondrial electron-transfer flavoprotein under carbon
limitation: the classic carbon-starvation energy route.

## What the model did not predict, and should not claim

The single strongest signal in the whole job is **Protein processing in the endoplasmic
reticulum, combined p = 8.6e-9** (genes 1.6e-8, proteins 0.024) — the unfolded protein
response. A gas-transport model has nothing to say about it. It is a real and well-known
spaceflight finding in BRIC hardware (OSD-321 / BRIC-22 was devoted to it), and it is a
reminder that the CO₂-starvation account explains part of this experiment, not all of it.

Also unpredicted and significant: Biotic Stress and Cellular response overview (MapMan),
glucosinolate biosynthesis, plant hormone signal transduction, JA synthesis, MAPK signalling,
circadian rhythm.

## The two weaker analyses, reported honestly

- **Metabolite class activity** (`methods/PAINTOMICS_CLASS_ACTIVITY.md`) returns amino acids at
  FDR 2.0e-6 — but the list it tests is our own prediction, so it largely restates how the
  prediction was built. Not evidence for the model.
- **Job 2, the hardware contrast** (`methods/PAINTOMICS_JOB2_OUTCOME.md`) returned 57 KEGG
  pathways, **0 significant**, no MapMan hits and no hub analysis. A metabolome-only
  submission cannot use PaintOmics' best tools. The 25× hardware effect stands as a model
  result; this job adds nothing to it.

## Files here

| File | Contents |
|---|---|
| `results/tables/T08_paintomics_significant.tsv` | all 31 significant pathways, with per-omic p-values |
| `results/tables/T09_paintomics_carbon.tsv` | the carbon pathways the model made claims about |
| `results/tables/T10_metabolite_hubs.tsv` | all 15 ranked hub compounds |
| `methods/PAINTOMICS_CLASS_ACTIVITY.md` | class test, with its circularity spelled out |
| `methods/PAINTOMICS_JOB2_OUTCOME.md` | why job 2 returned nothing |
| `job{1,2}_step1_response.json`, `job{1,2}_status.json` | raw server responses |

Reproduce with `python3 paintomics/submit_job.py --job 1`.
