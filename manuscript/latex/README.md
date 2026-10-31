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
- [x] Co-author list: none, sole author (Richard, 2026-10-30); placeholder removed from both files
- [x] Affiliation: independent researcher (Richard, 2026-10-30)
- [ ] City and country for the affiliation
- [x] ORCID: 0000-0001-5681-9857, checked on orcid.org and added to `.zenodo.json`. sn-jnl's
      `\orcid` needs `Orcidlogo.eps`, which the template did not ship, so enter it at submission
- [ ] Acknowledgements
- [x] Author contributions statement (sole-author form)
- [x] Competing interests: none, in the journal's wording (Richard, 2026-10-30, asked with
      CoSE's plant-growth hardware named)
- [x] AI-assistance statement: npj Microgravity's submission guidelines say LLM use should be
      documented in Methods, so it is the last Methods subsection (wording approved by Richard,
      2026-10-30)
- [x] Funding statement: none, in the journal's wording (Richard, 2026-10-30)

**Availability**
- [ ] Zenodo DOI for this repository, in Data availability and Code availability (Richard
      deposits; `.zenodo.json` is ready)
- [x] Repository URL in Code availability
- [ ] Archived release DOI for LunarLeaf-CFD (`lunarleaf` entry note; Richard deposits, and the
      repo has no GitHub release yet)

**References**
- [x] `osdr2025` — replace the corporate author with the full author list from the
      publisher record (done; verified against CrossRef)
- [x] `kegg` — add the primary KEGG citation alongside the URL (`kanehisa2000`; author
      list from PubMed, since CrossRef carries only the first author)

**Scientific decisions the author should confirm**
- [x] Keep the hardware-confound framing as the headline, although the enclosure gradient did
      not hold (only illumination did): kept as is (Richard, 2026-10-30)
- [x] The ER/unfolded-protein-response result, which this model does not explain, keeps its
      one paragraph (Richard, 2026-10-30)

## Figures

`figures/*.pdf` are copies, not symlinks, so the directory survives archiving on its own.
Regenerate with `bash scripts/run_all.sh` and re-copy:

```bash
cp results/figures/fig*.pdf results/figures/supplementary/fig*.pdf manuscript/latex/figures/
```

`scripts/run_all.sh` does this copy itself, so the compiled PDF cannot keep showing a figure
that has since been regenerated.
