#!/usr/bin/env python3
"""Differential expression for OSD-522 shoots: Space Flight vs Ground Control.

OSDR ships unnormalised RSEM counts for OSD-522 but no differential-expression table,
so we compute one. The design is a clean two-group comparison — 12 samples, all
"Plant Shoots", 6 Ground Control vs 6 Space Flight — read from the study runsheet
rather than assumed from sample-name prefixes.

Method: PyDESeq2 (the DESeq2 negative-binomial model — median-of-ratios size factors,
dispersion shrinkage, Wald test, Benjamini-Hochberg FDR) when it is importable. If it
is not, a documented fallback runs instead: DESeq2's own median-of-ratios normalisation
followed by a Welch t-test on log2 CPM with BH correction. The method that actually ran
is recorded in the output header and printed — it is never left ambiguous.

Outputs, in PaintOmics' input format:
  paintomics/upload/gene_expression_values.tab    #geneID <TAB> FLT_vs_GC   (log2 fold change)
  paintomics/upload/gene_expression_relevant.tab  one AGI per line, FDR < alpha
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import CACHE, UPLOAD  # noqa: E402

COUNTS = os.path.join(CACHE, "GLDS-522_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv")
RUNSHEET = os.path.join(CACHE, "GLDS-522_rna_seq_bulkRNASeq_v2_runsheet.csv")

FLIGHT, GROUND = "SpaceFlight", "GroundControl"
CONDITION = "FLT_vs_GC"


# --------------------------------------------------------------------------- load

def load_design() -> pd.DataFrame:
    """Sample -> condition, taken from the runsheet's factor-value columns."""
    rs = pd.read_csv(RUNSHEET)
    col_sample = next(c for c in rs.columns if c.strip().lower() == "sample name")
    col_flight = next(c for c in rs.columns if "spaceflight" in c.lower())
    col_part = next(c for c in rs.columns if "organism part" in c.lower())

    design = pd.DataFrame({
        "sample": rs[col_sample].astype(str),
        "organism_part": rs[col_part].astype(str),
        "raw_factor": rs[col_flight].astype(str),
    })
    mapping = {"space flight": FLIGHT, "ground control": GROUND}
    design["condition"] = design["raw_factor"].str.strip().str.lower().map(mapping)

    unmapped = design[design["condition"].isna()]
    if not unmapped.empty:
        raise SystemExit(f"unrecognised Spaceflight factor values: "
                         f"{sorted(unmapped['raw_factor'].unique())}")
    return design.set_index("sample")


def load_counts() -> pd.DataFrame:
    """Gene x sample count matrix, AGI-indexed."""
    counts = pd.read_csv(COUNTS, index_col=0)
    counts.index = counts.index.astype(str).str.strip()
    counts.index.name = "gene"
    return counts


def check_integrity(counts: pd.DataFrame, design: pd.DataFrame) -> None:
    """Assert the study is the shape we believe it is before trusting any result."""
    problems = []

    if list(counts.columns) != list(design.index):
        missing = set(counts.columns) ^ set(design.index)
        problems.append(f"counts columns and runsheet samples disagree: {sorted(missing)}")

    n_flt = int((design["condition"] == FLIGHT).sum())
    n_gc = int((design["condition"] == GROUND).sum())
    if (n_flt, n_gc) != (6, 6):
        problems.append(f"expected 6 flight + 6 ground, got {n_flt} + {n_gc}")

    parts = set(design["organism_part"].str.strip().str.lower())
    if parts != {"plant shoots"}:
        problems.append(f"expected shoots only, got organism parts {sorted(parts)}")

    agi = counts.index.str.match(r"^AT[1-5CM]G\d{5}$")
    frac = float(agi.mean())
    if frac < 0.95:
        problems.append(f"only {frac:.1%} of gene IDs look like AGI locus codes")

    if problems:
        raise SystemExit("data integrity check FAILED:\n  - " + "\n  - ".join(problems))

    print(f"  integrity OK: {counts.shape[0]:,} genes x {counts.shape[1]} samples, "
          f"{n_flt} flight vs {n_gc} ground, shoots only, {frac:.1%} AGI IDs")


# ----------------------------------------------------------------------- analysis

def run_pydeseq2(counts: pd.DataFrame, design: pd.DataFrame):
    """DESeq2 negative-binomial Wald test via PyDESeq2. Returns (results, method)."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    import pydeseq2

    metadata = design[["condition"]].copy()
    dds = DeseqDataSet(
        counts=counts.T.astype(int),          # PyDESeq2 wants samples x genes
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()
    stat = DeseqStats(dds, contrast=["condition", FLIGHT, GROUND], quiet=True)
    stat.summary()

    res = stat.results_df.rename(columns={"log2FoldChange": "log2fc", "padj": "fdr"})
    method = (f"PyDESeq2 {pydeseq2.__version__} — DESeq2 negative-binomial model, "
              f"median-of-ratios size factors, dispersion shrinkage, Wald test, BH FDR")
    return res[["log2fc", "pvalue", "fdr"]], method


def median_of_ratios(counts: pd.DataFrame) -> pd.Series:
    """DESeq2 size factors: median over genes of each sample's ratio to the geometric mean."""
    with np.errstate(divide="ignore"):
        log_counts = np.log(counts.replace(0, np.nan))
    log_gmean = log_counts.mean(axis=1)               # geometric mean per gene
    usable = np.isfinite(log_gmean)
    ratios = log_counts.loc[usable].sub(log_gmean[usable], axis=0)
    return np.exp(ratios.median(axis=0))


def bh_fdr(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg step-up adjusted p-values."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    vals = p[ok]
    n = vals.size
    order = np.argsort(vals)
    ranked = vals[order]
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]      # enforce monotonicity
    restored = np.empty(n)
    restored[order] = np.minimum(adj, 1.0)
    out[ok] = restored
    return out


def run_fallback(counts: pd.DataFrame, design: pd.DataFrame):
    """Median-of-ratios normalisation + Welch t-test on log2 CPM. Returns (results, method)."""
    from scipy import stats

    sf = median_of_ratios(counts)
    norm = counts.div(sf, axis=1)
    lcpm = np.log2(norm.div(norm.sum(axis=0), axis=1) * 1e6 + 1.0)

    flt = lcpm.loc[:, design.index[design["condition"] == FLIGHT]]
    gc = lcpm.loc[:, design.index[design["condition"] == GROUND]]

    t, p = stats.ttest_ind(flt, gc, axis=1, equal_var=False)
    res = pd.DataFrame({
        "log2fc": flt.mean(axis=1) - gc.mean(axis=1),
        "pvalue": p,
        "fdr": bh_fdr(p),
    }, index=counts.index)
    method = ("fallback — DESeq2 median-of-ratios size factors, Welch t-test on log2 CPM, "
              "BH FDR (PyDESeq2 unavailable)")
    return res, method


# ------------------------------------------------------------------------- output

def write_outputs(res: pd.DataFrame, method: str, alpha: float, min_count: int,
                  n_genes_in: int) -> None:
    os.makedirs(UPLOAD, exist_ok=True)

    values = res.dropna(subset=["log2fc"])
    values_path = os.path.join(UPLOAD, "gene_expression_values.tab")
    with open(values_path, "w") as fh:
        fh.write(f"#geneID\t{CONDITION}\n")
        for gene, row in values.iterrows():
            fh.write(f"{gene}\t{row['log2fc']:.6f}\n")

    sig = res[(res["fdr"] < alpha) & res["log2fc"].notna()]
    relevant_path = os.path.join(UPLOAD, "gene_expression_relevant.tab")
    with open(relevant_path, "w") as fh:
        for gene in sig.index:
            fh.write(f"{gene}\n")

    up = int((sig["log2fc"] > 0).sum())
    down = int((sig["log2fc"] < 0).sum())
    print(f"\n  method: {method}")
    print(f"  filtered {n_genes_in:,} -> {len(res):,} genes (row sum >= {min_count})")
    print(f"  wrote {values_path}  ({len(values):,} genes)")
    print(f"  wrote {relevant_path}  ({len(sig):,} genes at FDR < {alpha}; "
          f"{up:,} up in flight, {down:,} down)")

    # Provenance sits beside the upload files rather than inside them: PaintOmics'
    # parser expects the '#geneID' header line first and nothing else before the data.
    with open(os.path.join(UPLOAD, "gene_expression_METHOD.txt"), "w") as fh:
        fh.write(
            "OSD-522 (BRIC-LED-001) Arabidopsis shoots — Space Flight vs Ground Control\n"
            f"source:   NASA OSDR {os.path.basename(COUNTS)}\n"
            f"design:   6 Space Flight vs 6 Ground Control, all Plant Shoots (from runsheet)\n"
            f"filter:   genes with total count >= {min_count} across all samples\n"
            f"method:   {method}\n"
            f"contrast: {CONDITION} = log2(Space Flight / Ground Control)\n"
            f"relevant: FDR < {alpha} ({len(sig)} genes: {up} up, {down} down)\n"
        )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=0.05, help="FDR threshold (default 0.05)")
    ap.add_argument("--min-count", type=int, default=10,
                    help="drop genes whose total count across samples is below this")
    ap.add_argument("--force-fallback", action="store_true",
                    help="use the numpy/scipy fallback even if PyDESeq2 is available")
    args = ap.parse_args()

    for path in (COUNTS, RUNSHEET):
        if not os.path.exists(path):
            raise SystemExit(f"missing {path} — run osdr/fetch_osd522.py first")

    print("OSD-522 shoots — differential expression (Space Flight vs Ground Control)")
    design = load_design()
    counts = load_counts()
    check_integrity(counts, design)

    n_in = counts.shape[0]
    counts = counts.loc[counts.sum(axis=1) >= args.min_count]

    if args.force_fallback:
        res, method = run_fallback(counts, design)
    else:
        try:
            res, method = run_pydeseq2(counts, design)
        except ImportError as exc:
            print(f"  PyDESeq2 unavailable ({exc}); using fallback", file=sys.stderr)
            res, method = run_fallback(counts, design)

    write_outputs(res, method, args.alpha, args.min_count, n_in)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
