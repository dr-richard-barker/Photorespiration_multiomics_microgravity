# ABRS / TAGES (OSD-7) extension — outputs

> **UNREVIEWED — not part of the manuscript or the site.** This ABRS / TAGES (OSD-7)
> extension was added in commit `edb9068` (2026-09-25). The 2026-10-27 `/abai check`
> found that its scripts re-run byte-identically, but that its conclusions do not hold
> as written, so it was taken out of the paper and the site on 2026-10-28 and kept here
> for later review. Do not cite any number below until these are resolved:
>
> - **No ABRS PaintOmics job was run.** The ABRS p-values come from a local one-sided
>   Fisher test on genes with unadjusted Welch p < 0.05 (no gene passes FDR in any
>   organ). They are set beside BRIC-LED's PaintOmics combined p-values as if the two
>   were the same statistic; they are not.
> - **The test has no direction**, so "induced" / "up-regulated" labels are unsupported.
> - **The ≈3,500 ppm ISS cabin CO₂ has no source.** Paul et al. 2013 report that
>   volatiles were scrubbed equally in flight and ground hardware, which undercuts the
>   "trapped ethylene" explanation for BRIC ER stress.
> - **Airflow is not modelled.** ABRS operating points reuse the BRIC LunarLeaf g_bl
>   rows; only C_a differs.
> - **OSD-7 is Affymetrix.** `data/study_registry.tsv` excludes it from the six-study
>   ladder because its log2 fold changes are not on the RNA-seq scale.
> - **OSD-16 proteomics was never analysed**; OSDR holds only raw iTRAQ spectra.
>
> Full findings: the 2026-10-27 check report (kept outside this repository).

Produced by `scripts/abrs/run_abrs.sh`; method notes in [`methods/abrs/ABRS_ANALYSIS_REPORT.md`](../../methods/abrs/ABRS_ANALYSIS_REPORT.md).
