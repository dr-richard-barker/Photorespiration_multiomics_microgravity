# PaintOmics submission — OSD-522 three-layer job

Everything in `upload/` is ready and validated (`python3 validate_upload.py`).
The site is <https://paintomics.org> ("PaintOmics AI"), running on Supercomputador Drago,
CSIC, Spain.

## Before you upload — what leaves this machine

| File | Origin | Sensitivity |
|---|---|---|
| `gene_expression_*.tab` | NASA OSDR OSD-522, already public | public data |
| `proteomics_*.tab` | NASA OSDR OSD-522, already public | public data |
| `metabolomics_*.tab` | this repository's model output | your unpublished prediction |

The gene and protein layers are derived from a public NASA repository. The metabolite layer
is **your own unpublished model output** — uploading it puts it on a third-party server.

**Keep the AI interpretation switched off** unless you decide otherwise. The site states
that enabling it sends your data to `llm.iiia.es` (IIIA-CSIC, EU) — a separate endpoint from
the pathway analysis itself.

## Job 1 — the three-layer job

**Step 1, organism and databases**

- Organism: type `Arabid` and choose **Arabidopsis thaliana (thale cress)** — KEGG code
  `ath`, confirmed present.
- Databases: **KEGG** (always on) and **MapMan**. MapMan is the plant-specific database and
  resolves photorespiration as its own bin, which KEGG does not — it folds it into the much
  broader `ath00630` glyoxylate and dicarboxylate map. Reactome and OmniPath are
  animal-oriented; leave them unticked.
- AI interpretation: **off**.

**Step 2, files** — drag `Gene expression`, `Proteomics` and `Metabolomics` into the
selected list, then:

| Slot | Data file | Relevant features file |
|---|---|---|
| Gene expression | `gene_expression_values.tab` | `gene_expression_relevant.tab` |
| Proteomics | `proteomics_values.tab` | `proteomics_relevant.tab` |
| Metabolomics | `metabolomics_values.tab` | `metabolomics_relevant.tab` |

Leave the metabolomics **experimental design** slot empty. It drives the metabolite class
activity test, which needs biological replicates; a prediction has none, and supplying
model-uncertainty draws as if they were replicates would misrepresent them.

**Step 3, identifier matching** — expect:

- genes as **AGI locus codes** (`AT1G01010`), which is exactly what KEGG `ath` uses, so
  these should map essentially one-to-one;
- proteins as **UniProt accessions**, which PaintOmics converts to Entrez itself;
- metabolites as **names**, which PaintOmics matches to KEGG compound IDs. Check its
  assignments against `../metabolome/compound_provenance.tsv`, which carries the KEGG ID we
  intend for each of the 21 compounds. Ambiguous names can be corrected on that screen.

## Job 2 — the hardware contrast (metabolome only)

Same organism and databases, one omic:

| Slot | Data file | Relevant features file |
|---|---|---|
| Metabolomics | `metabolomics_hardware_values.tab` | `metabolomics_hardware_relevant.tab` |

This is the sealed-BRIC versus vented-hardware comparison — the effect that is 25× the
gravity effect but cancels out of the flight-versus-ground contrast. It cannot be a fourth
column of job 1, because the gene and protein layers have no counterpart for it.

## What to look at in the results

Target maps, in order of expected value:

1. **MapMan photorespiration bin** — the sharpest test of the prediction.
2. **`ath00630`** glyoxylate and dicarboxylate metabolism (81 genes, 67 compounds).
3. **`ath00710`** carbon fixation, and **`ath00500`** starch and sucrose — where the model
   predicts a coordinated decline, and where the transcriptome already shows one.
4. **`ath01200`** carbon metabolism, for the overview.
5. **`ath00010`** glycolysis, and **`ath00250`/`ath00260`** amino acid metabolism, for the
   starvation response.

Also run the **metabolite hub analysis**, which asks which compounds have an unusual density
of significant genes around them. Because our metabolite layer is predicted and our gene
layer is measured, a hub that lights up is a place where the model and the real data agree
on the same neighbourhood — the most interesting kind of hit this job can produce.

## Expected outcome, written down before running

Recorded in advance so the run cannot be read retrospectively.
`osdr/falsification_check.py` has already tested the prediction directly against the data:

- carbon fixation, starch/sucrose, photosystem and Rubisco genes **down** — confirmed
  (p = 1.3e-4, 3.5e-6, 3.8e-6, 4.5e-3);
- DIN carbon-starvation markers **up** — confirmed (median +0.71, p = 2.6e-3);
- photorespiratory enzymes **not induced** — confirmed (p = 0.09, no shift), and sitting
  above Rubisco as predicted (+0.24, p = 6.8e-3; independently +0.29, p = 0.049 in the
  proteome);
- no fermentation or hypoxia signature — confirmed.

So PaintOmics should paint a coherent picture of a **carbon-starved, photosynthesis-
downregulated shoot with photorespiration untouched**. If it instead lights up
photorespiratory induction or a hypoxia response, the model is wrong and that is the finding.
