# A gas-transport model predicts the *Arabidopsis* spaceflight transcriptome

**A sealed, illuminated growth canister has its CO₂ consumed by the plants inside it, within
minutes. That happens in the ground control as much as in flight — so it does not bias the
flight-versus-ground contrast, it changes what that contrast is about. The spaceflight
response is being measured in plants their own hardware has already carbon-starved.**

[**Interactive model and data explorer**](https://dr-richard-barker.github.io/Photorespiration_multiomics_microgravity/)
· [manuscript PDF](https://dr-richard-barker.github.io/Photorespiration_multiomics_microgravity/Photorespiration_multiomics_microgravity.pdf)
· [manuscript DOCX](https://dr-richard-barker.github.io/Photorespiration_multiomics_microgravity/Photorespiration_multiomics_microgravity.docx)
· [supplementary PDF](https://dr-richard-barker.github.io/Photorespiration_multiomics_microgravity/Photorespiration_multiomics_microgravity_SI.pdf)
· [supplementary DOCX](https://dr-richard-barker.github.io/Photorespiration_multiomics_microgravity/Photorespiration_multiomics_microgravity_SI.docx)
· [figures](results/figures) · [tables](results/tables)

![The model against the measurement](results/figures/fig05_blind_test.png)

---

## What this is

A validated lattice-Boltzmann model of leaf boundary-layer gas transport
([LunarLeaf-CFD](https://github.com/dr-richard-barker/LunarLeaf-CFD)) coupled to a
Farquhar–von Caemmerer–Berry model of C₃ photosynthesis, used to ask what spaceflight
*growth hardware* — not gravity — does to the CO₂ reaching Rubisco. Its predictions are
tested against NASA OSDR. The predictions and the tests entered the repository in the same
commit, so the record cannot show that the predictions came first.

### The four findings

1. **A lit sealed canister is a CO₂-starvation chamber.** 400 ppm to near zero in about
   seven minutes; 12 h carbon gain at **1 %** of an unenclosed control, against 90 % under
   gas-permeable tape and 100 % when vented. These come from LunarLeaf-CFD's
   zero-dimensional enclosure models, with an assumed 30 mL headspace and tape permeance,
   not from its validated lattice-Boltzmann solver.
2. **That effect is ~25× the microgravity effect, and invisible.** Because the drawdown is a
   mass balance it happens in the ground control too, so it subtracts out of
   flight-versus-ground.
3. **Our own starting hypothesis was wrong.** We expected CO₂ starvation to *amplify* the
   microgravity boundary-layer penalty. It converges instead: with no assimilation there is
   no flux across the boundary layer for gravity to impede. What survives is a persistent
   5–9 % assimilation deficit at saturating light, which is what the model is solved at; at
   OSD-522's own 60 µmol m⁻² s⁻¹ it is under 3 %
   ([`T13_light_sensitivity.tsv`](results/tables/T13_light_sensitivity.tsv)).
4. **The model predicted much of the measurement.** It called **eight of nine** transcriptome
   gene-set responses in OSD-522, three of them by predicting no change, and four of seven in
   the proteome. In PaintOmics, photosynthesis pathways reached uncorrected *p* < 0.05 (none
   survives FDR, and the Calvin cycle, predicted down, did not move) and photorespiration
   ranked **214th of 232** pathways tested — bottom 8 %, and the least enriched of any
   pathway the model made a claim about.

### And what did not work

- **The predicted enclosure gradient does not hold** (Spearman ρ = −0.50, *p* = 0.67). Sealed
  and vented fall the right way round but the CARA plates (modelled as taped; OSDR does not
  say how they were sealed) sit highest, and three contrasts cannot rank three classes. What
  *does* separate the six contrasts, from five studies, is **illumination**, on the
  starvation-over-photosynthesis contrast though not on every contrast. That fits the model,
  since drawdown requires photosynthesis, but the model predicted no signal in the dark, not
  the inversion seen, and OSD-678's lit arm grew in only 4–6 µmol m⁻² s⁻¹ of cabin light.
- **Most pathway results do not survive correction.** 32 of 232 pathways reach uncorrected
  *p* < 0.05; six survive FDR, and none is a pathway the model made a claim about.
- **The strongest signal in the data is one this model does not explain**: protein processing
  in the endoplasmic reticulum, *p* = 8.6 × 10⁻⁹, the unfolded protein response.
- **OSDR has no plant metabolome.** Of 567 studies, twelve are metabolite profiling and every
  one is mouse, human, rat or microbial
  ([`T03`](results/tables/T03_osdr_metabolomics_census.tsv)); **none is a plant**. The third
  omic layer had to be predicted, and is labelled as model output everywhere it appears.
- **The boundary-layer numbers come from an older LunarLeaf-CFD export.** Its steady-state
  re-run lowered microgravity g_bl at every scale, and the 5–9 % has not been recomputed;
  see [`data/lunarleaf/PROVENANCE.md`](data/lunarleaf/PROVENANCE.md).

---

## Layout

```
data/                 inputs only
  study_registry.tsv    which OSDR studies enter, the exact contrast column, and why
  lunarleaf/            vendored CFD tables + provenance
  mapman/               vendored MapMan diagrams + gene-to-bin mapping + provenance
  paintomics_upload/    the validated submission bundle
  paintomics_*.tsv      curation joined into T09/T10 (which pathways and compounds the model claimed)
  abrs/                 ABRS / OSD-7 inputs (unreviewed; see below)
  cache/                OSDR + KEGG downloads (git-ignored, regenerated on demand)
scripts/              every executable, numbered in pipeline order
  fvcb.py               the model core, importable
  genesets.py           KEGG-derived gene sets, shared by both analyses
  mapman.py             MapMan pathway membership, rebuilt from the diagram layouts
  figures/              one script per figure + the shared visual system
  run_all.sh            regenerates everything from a cold cache
results/tables/       T01…T13 + MANIFEST.tsv (fails if a table lacks provenance)
results/abrs/, scripts/abrs/, methods/abrs/   the ABRS extension (unreviewed; see below)
results/figures/      eight figures, PNG at 300 dpi and PDF vector
manuscript/latex/     npj Microgravity style, compiles locally and in CI
  main.tex              the manuscript; supplementary.tex is a separate document
results/figures/supplementary/  Supplementary Figure S1
methods/              prose methods and the honest caveats
docs/                 the interactive site
```

### Reproducing

```bash
bash scripts/run_all.sh
```

Needs network on a cold cache. `python3 scripts/check_js_parity.py` then confirms the
browser model still agrees with the Python (currently **exact**, 0.000e+00 across 180
parameter combinations).

The MapMan inputs under `data/mapman/` are vendored rather than fetched, because
MapManStore's own published links are already dead and the diagrams are reachable only one
portlet query at a time. `python3 scripts/mapman.py --fetch` refreshes all 14 files or none;
`run_all.sh` deliberately does not, since the pathway universe is an enrichment denominator.

PaintOmics submission is deliberately separate, because it uploads to a third-party server:

```bash
python3 scripts/05_submit_paintomics.py --job 1
```

---

## Progress

| Phase | Status |
|---|---|
| Model coupling (CFD → FvCB) | **done** — six physical self-checks pass |
| Predicted metabolite layer | **done** — 21 compounds, KEGG ids from the REST API |
| OSD-522 transcriptome + proteome | **done** — PyDESeq2, 1,795 DE genes; 5,160 proteins |
| Falsification test | **done** — 8 of 9 transcriptome gene sets; 4 of 7 in the proteome |
| PaintOmics integration | **done** — jobs `m1z16Qg3DK` and `EVFahRGk7T` |
| Cross-study ladder (6 contrasts, 5 studies) | **done** — illumination separates; enclosure gradient does not |
| Eight-figure set | **done** |
| npj manuscript | **compiles** — 13 pages, no unresolved references; author block is placeholders |
| Supplementary information | **compiles** — 2 pages, one figure, cross-referenced to the manuscript |
| Interactive site | **done** — model, enclosures, omics (volcano, heatmap, Sankey), pathways, six contrasts |
| FAIR packaging | **done** — manifest, CITATION.cff, .zenodo.json, MIT |
| Zenodo deposit | **pending** — needs the author fields below |

### What still needs a human

Nothing here was invented to fill a gap, so these remain visibly open:

- **Author block.** Co-authors, affiliation, ORCIDs, funding, contributions —
  see the checklist in [`manuscript/latex/README.md`](manuscript/latex/README.md).
- **Zenodo DOI**, then paste it into the manuscript's Data and Code availability sections.
- **A framing decision — partly settled.** The title now claims the mechanism (photosynthesis
  depletes its own enclosure) and its consequence (the response is read out in carbon-starved
  plants), rather than an enclosure gradient the data did not support. Still worth your call:
  whether to lead instead with the prediction result.
- **Whether to re-run on the current LunarLeaf-CFD export**, which would change the
  5–9 % deficit and the operating points.
- **The ER / unfolded-protein-response result** currently gets one paragraph. It is the
  largest signal in the data and this model says nothing about it.

### Where it could go next

- **In-canister CO₂ and O₂ logging.** The cheapest and highest-value change: it converts this
  analysis's central assumption into a measurement, and makes an entire archive of
  BRIC-derived transcriptomes reinterpretable.
- **A ~20-compound targeted metabolite panel on archived flight material** — it would be the
  first plant spaceflight metabolome in existence, and
  [`results/tables/T06_predicted_compounds.tsv`](results/tables/T06_predicted_compounds.tsv)
  is a ready-made target list with a quantitative prediction attached to each compound.
- **A three-arm design** (flight sealed, ground sealed, ground vented). Two arms confound
  hardware with gravity; three separate them, and the third arm never leaves the ground.
- **More studies**, if the tissue constraint can be relaxed defensibly — the ladder is
  currently three lit and three dark contrasts, which is too few to rank enclosure classes.

Full reasoning in [`FUTURE_EXPERIMENTS.md`](FUTURE_EXPERIMENTS.md).

---

## Honesty notes

This project has a few standing rules, and they are load-bearing:

- **The metabolite layer is model output and says so** in every file, figure and page that
  shows it.
- **Study selection is auditable.** [`data/study_registry.tsv`](data/study_registry.tsv)
  records each inclusion and each exclusion with its reason. Tissue is the binding
  constraint: root studies are excluded because a photosynthesis prediction says nothing
  about roots.
- **Contrasts are matched.** OSD-678 offers 36 flight-versus-ground contrasts but only six
  are matched on genotype, ecotype and light; the rest also differ in genotype, ecotype or
  light (flown Col-0 against ground-control *phyD*, for example).
- **Gene sets come from KEGG, not memory.** An earlier hand-written list here put PGLP1 at
  the wrong locus, and *Arabidopsis* reuses the symbols CAT2 and SEN1 for unrelated genes.
- **MapMan membership is a reconstruction and is labelled as one.** PaintOmics publishes
  enrichment results but not the gene lists behind them, and MapMan has no REST API, so the
  13 MapMan pathways are rebuilt from the diagram layouts by PaintOmics' own rules. It
  matches PaintOmics' feature count exactly or within three for 9 of the 13 and runs larger
  for the four largest. PaintOmics matches bin codes anywhere in a gene's bin, not only
  from the start, so its diagrams hold more than their names suggest: *photosynthesis* has
  1,722 measured loci this way, 488 matched from the start. The site shows both counts on every node
  rather than picking one, and
  [`data/mapman/PROVENANCE.md`](data/mapman/PROVENANCE.md) says what is not recoverable.
- **References were verified against the publisher record**, never recalled.
- **A sibling CFD repository is deliberately not used.** See
  [`methods/CFD_PROVENANCE_CONCERN.md`](methods/CFD_PROVENANCE_CONCERN.md).
- **The metabolite class-activity test is close to circular** and is reported as such; the
  hub analysis is not, and is the one to read.

## In the repository but not in the paper

- **ABRS / TAGES (OSD-7) extension**: `scripts/abrs/`, `results/abrs/`,
  `methods/abrs/`, `data/abrs/`. **Unreviewed.** Its scripts run, but its conclusions do
  not hold as written: no ABRS PaintOmics job was run, so its pathway p-values are a
  different statistic from the BRIC-LED ones they are compared with. The banner at the top
  of [`methods/abrs/ABRS_ANALYSIS_REPORT.md`](methods/abrs/ABRS_ANALYSIS_REPORT.md) lists
  what is wrong. Nothing from it appears in the manuscript or on the site.
- **PhysioSpace projection**: `scripts/08_physiospace_projection.R` has never been run.
  `run_all.sh` skips it unless the R packages are installed, and the heatmap script has no
  fallback. An earlier version drew that heatmap from random numbers; the table and figure
  it made have been deleted.

## Data sources

NASA Open Science Data Repository — OSD-522 (BRIC-LED-001), OSD-38 (BRIC-20),
OSD-321 (BRIC-22), OSD-678 (CARA), OSD-427 (APEX-04/VEGGIE).
Gas transport from LunarLeaf-CFD. Pathway enrichment from KEGG and MapMan via PaintOmics;
KEGG membership from the KEGG REST API, MapMan membership rebuilt from GoMapMan's gene
mapping and the MapManStore diagram layouts (see
[`data/mapman/PROVENANCE.md`](data/mapman/PROVENANCE.md)).

## Licence

MIT — see [LICENSE](LICENSE).
