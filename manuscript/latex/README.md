# Manuscript — npj Microgravity style

```bash
cd manuscript/latex && latexmk -pdf main.tex && latexmk -pdf supplementary.tex
```

**In that order.** `supplementary.tex` uses `xr` to read `main.aux`, so that its caption can
say "Fig. 6a" without anyone typing the number. Build it alone from a clean checkout and the
reference renders as `??`, which latexmk reports as an undefined reference — loud rather than
a silently wrong figure number.

Builds locally on TeX Live 2026 (verified: 12 pages, no unresolved references or citations,
no overfull boxes). Also builds in CI — see `.github/workflows/build-manuscript.yml`.

## The two class patches, and why

`sn-jnl.cls` here is **already patched**, copied from a sibling repository. On TeX Live 2026
the stock Springer Nature class fails without them:

1. `\RequirePackage{manyfoot}` before the class's `\AtBeginDocument` block. The class calls
   `\SetFootnoteHook`, `\DeclareNewFootnote` and `\footinsA` but never loads `manyfoot`,
   giving "Undefined control sequence" at `\begin{document}`.
2. `\RequirePackage{xcolor}` for the same reason — the class calls `\definecolor` inside
   `\AtBeginDocument`.

`.latexmkrc` sets `$bibtex_use = 2` so latexmk actually runs BibTeX; without it the build
stops at a missing `.bbl`.

## Two BibTeX traps this file has already hit

- **An "at" sign starts an entry even inside a `%` comment.** A header comment that mentioned
  an at-entry literally aborted the parse of the first real reference.
- **A `%` INSIDE an entry is literal text**, not a comment, and fails with
  "expecting a , or }".

## TODO before submission

Every item below is a placeholder in `main.tex` or `references.bib`. Nothing was invented to
fill them, and the manuscript will read as unfinished until they are cleared.

**Authorship and front matter**
- [ ] Co-author list (`\author[1]{\fnm{[TODO}\sur{ co-author]}}`)
- [ ] Affiliation, city, country
- [ ] ORCID identifiers
- [ ] Acknowledgements
- [ ] Author contributions statement
- [ ] Funding statement

**Availability**
- [ ] Zenodo DOI for this repository, in Data availability and Code availability
- [ ] Repository URL in Code availability
- [ ] Archived release DOI for LunarLeaf-CFD (`lunarleaf` entry note)

**References**
- [ ] `osdr2025` — replace the corporate author with the full author list from the
      publisher record
- [ ] `kegg` — add the primary KEGG citation alongside the URL

**Scientific decisions the author should confirm**
- [ ] Whether to keep the hardware-confound framing as the headline, given the enclosure
      gradient did not hold (only illumination did)
- [ ] Whether the ER/unfolded-protein-response result, which this model does not explain,
      deserves more than the paragraph it currently has

## Figures

`figures/*.pdf` are copies, not symlinks, so the directory survives archiving on its own.
Regenerate with `bash scripts/run_all.sh` and re-copy:

```bash
cp results/figures/fig*.pdf results/figures/supplementary/fig*.pdf manuscript/latex/figures/
```

`scripts/run_all.sh` does this copy itself, so the compiled PDF cannot keep showing a figure
that has since been regenerated.
