# Vendored MapMan inputs

PaintOmics reports significant pathways from two databases. KEGG membership is fetched live
from the KEGG REST API (`scripts/genesets.py`). MapMan has no equivalent API, so the 13
MapMan pathways in `results/tables/T08_paintomics_significant.tsv` are reconstructed here
from the same inputs PaintOmics' own installer uses, by the same rules
(`scripts/mapman.py`).

| File | What it holds | Source |
|---|---|---|
| `gene-to-mapman_ath.tsv.gz` | Araport11 AGI locus → MapMan ontology bins, 33,341 loci over 1,384 bins | GoMapMan PaintOmics export, `protein_2018-05-25` |
| `diagrams/*.xml` (13) | the MapMan diagram layouts, each listing the ontology bins drawn on it | 3 from the GoMapMan 20-diagram tarball, 10 from the MapManStore archive |

Refresh with `python3 scripts/mapman.py --fetch`, which rewrites all 14 files or none.

## Why these are vendored rather than cached

Both upstream paths are fragile, and one is already broken. MapManStore's published static
links (`https://www.plabipd.de/img/mapman_36/*.zip`, linked from its own page) return 404;
the diagrams are only reachable through a Liferay portlet query, one request per diagram.
Vendoring keeps the repository archivable on Zenodo and reproducible offline.

Resource ids are taken from PaintOmics' `mapman_extra_diagrams.json`. MapMan **X4**-era
diagrams are deliberately not used: they carry a renumbered ontology that the 2018 gene
mapping does not match, so they would place the wrong genes without erroring. PaintOmics
excludes them for the same reason.

## How far the reconstruction agrees with PaintOmics

Compared against the `features` column of T08 — the number of our uploaded features
PaintOmics placed on each diagram — the reconstruction is **exact or within three features
for 9 of the 13** diagrams:

| Agreement | Diagrams |
|---|---|
| exact | Transport overview, AGPs, JA Synthesis, Prokaryotic Ribosome SSU 5S branch assembly, Sulphate Assimilation |
| within 3 | Large enzyme families overview, Raffinose metabolism, receptor like kinases, photosynthesis |
| **larger here** | Cellular response overview (3,617 vs 2,854), Regulation overview (9,560 vs 8,831), Biotic Stress (9,321 vs 4,435), Overview (21,056 vs 11,696) |

The four that disagree are the four largest maps. **Why they disagree is not recoverable
from PaintOmics' published sources, and this file does not guess.** What was tested:

- MapMan's bin inclusion rule is a prefix match, and PaintOmics implements it with an
  *unanchored* regular expression, so a diagram listing bin `5` also takes genes whose bin
  merely *contains* a segment `5` (`29.5.11.4.2`). Anchoring the match instead breaks two
  diagrams PaintOmics gets right — Raffinose metabolism collapses from 3,073 to 334 and
  photosynthesis from 1,722 to 488, against PaintOmics' 3,076 and 1,724 — so the unanchored
  rule is the one the server used, and is the one used here.
- Anchoring does **not** explain the four. Three of them barely move under it (Cellular
  response overview 3,617 → 3,614, Regulation overview 9,560 → 9,557, Biotic Stress
  9,321 → 8,678) and `Overview` is identical either way, because it lists all 36 top-level
  bins and so covers the whole mapping under any rule.
- Restricting membership to genes with a KEGG identifier does not explain it either
  (Overview 21,056 → 20,420 against PaintOmics' 11,696).

So the gap is real, unexplained, and roughly a factor of two on the largest maps.

Both counts ship to the site, and the interactive Sankey labels any diagram where they
disagree, rather than quietly presenting one of them as the truth.

## A MapMan "pathway" is a diagram, not an ontology bin

All 13 names are diagram titles; none is an ontology term, which is why matching them
against `ontology.obo` returns 0 of 13. Genes reach a diagram through the bins drawn on it.

One trap is worth recording because it fails silently: diagrams sometimes zero-pad a bin
segment (`14.01`, `13.1.5.3.02`) and the gene mapping never does (`14.1`). Compared
verbatim the two spellings never meet and the affected diagram imports with **no genes at
all**. `normalise_bin()` strips the padding; PaintOmics carries the same warning.
