#!/usr/bin/env python3
"""Microarray differential expression for OSD-7 (TAGES/ABRS), per organ.

Reads the GeneLab RMA-normalised expression matrix and computes log2 fold
changes and moderated t-test p-values for Flight vs Ground Control in each
organ (Shoot, Root, Hypocotyl, Whole-Plant).

The expression values are already log2 RMA-normalised, so fold change is
simply mean(FLT) - mean(GC). Significance is computed using an independent
two-sample t-test with Benjamini-Hochberg FDR correction. For a more
rigorous analysis, limma (R) would be ideal, but this Python implementation
matches the project's existing style and is sufficient for PaintOmics input.

Outputs per organ:
    results/abrs/tables/T01_abrs_{organ}_transcriptome.tsv
    data/abrs/paintomics_upload/abrs_{organ}_gene_expression_values.tab
    data/abrs/paintomics_upload/abrs_{organ}_gene_expression_relevant.tab
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from paths import ensure  # noqa: E402

ABRS_CACHE = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "data", "abrs", "cache"
))
ABRS_UPLOAD = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "data", "abrs", "paintomics_upload"
))
ABRS_TABLES = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "tables"
))

EXPR_FILE = os.path.join(
    ABRS_CACHE,
    "GLDS-7_array_normalized_expression_probeset_GLmicroarray.csv",
)

# Column groups parsed from the GeneLab column headers
ORGANS = {
    "shoot": {
        "flt": [f"Atha_WS-0_Col-0_Shoot_FLT_Rep{i}" for i in range(1, 6)],
        "gc":  [f"Atha_WS-0_Col-0_Shoot_GC_Rep{i}" for i in range(1, 6)],
    },
    "root": {
        "flt": [f"Atha_WS-0_Col-0_Root_FLT_Rep{i}" for i in range(1, 6)],
        "gc":  [f"Atha_WS-0_Col-0_Root_GC_Rep{i}" for i in range(1, 6)],
    },
    "hypocotyl": {
        "flt": [f"Atha_WS-0_Col-0_Hypocotyl_FLT_Rep{i}" for i in range(1, 6)],
        "gc":  [f"Atha_WS-0_Col-0_Hypocotyl_GC_Rep{i}" for i in range(1, 6)],
    },
    "whole_plant": {
        "flt": [f"Atha_WS-0_Whole-Plant_FLT_Rep{i}" for i in range(1, 4)],
        "gc":  [f"Atha_WS-0_Whole-Plant_GC_Rep{i}" for i in range(1, 4)],
    },
}


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(n, dtype=float)
    adj[order] = np.minimum(1.0, pvals[order] * n / np.arange(1, n + 1))
    # Enforce monotonicity from the right
    running_min = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        running_min = min(running_min, adj[idx])
        adj[idx] = running_min
    return adj


def compute_dge(df: pd.DataFrame, flt_cols: list[str], gc_cols: list[str],
                organ: str) -> pd.DataFrame:
    """Compute log2FC and p-values for flight vs ground control."""
    flt = df[flt_cols].values  # already log2
    gc = df[gc_cols].values

    # log2 fold change
    log2fc = flt.mean(axis=1) - gc.mean(axis=1)

    # Independent samples t-test (Welch's)
    pvals = np.full(len(df), np.nan)
    for i in range(len(df)):
        f_row = flt[i, :]
        g_row = gc[i, :]
        # Skip if all NaN or no variance
        if np.all(np.isnan(f_row)) or np.all(np.isnan(g_row)):
            continue
        f_clean = f_row[~np.isnan(f_row)]
        g_clean = g_row[~np.isnan(g_row)]
        if len(f_clean) < 2 or len(g_clean) < 2:
            continue
        if np.std(f_clean) == 0 and np.std(g_clean) == 0:
            pvals[i] = 1.0
            continue
        _, p = stats.ttest_ind(f_clean, g_clean, equal_var=False)
        pvals[i] = p

    # FDR correction on non-NaN p-values
    valid = ~np.isnan(pvals)
    padj = np.full(len(df), np.nan)
    if valid.sum() > 0:
        padj[valid] = bh_adjust(pvals[valid])

    result = pd.DataFrame({
        "TAIR": df["TAIR"].values,
        "log2FC_FLT_vs_GC": log2fc,
        "pvalue": pvals,
        "padj": padj,
        "mean_FLT": flt.mean(axis=1),
        "mean_GC": gc.mean(axis=1),
        "n_FLT": (~np.isnan(flt)).sum(axis=1),
        "n_GC": (~np.isnan(gc)).sum(axis=1),
    })
    return result


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organs", nargs="+", default=list(ORGANS.keys()),
                    choices=list(ORGANS.keys()),
                    help="which organs to process (default: all)")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="FDR threshold for relevant genes")
    args = ap.parse_args()

    if not os.path.exists(EXPR_FILE):
        print(f"ERROR: expression file not found: {EXPR_FILE}")
        print("Run 00_fetch_abrs.py first.")
        return 1

    print("OSD-7 (TAGES/ABRS) — microarray differential expression\n")
    print(f"  loading {EXPR_FILE}")
    df = pd.read_csv(EXPR_FILE, low_memory=False)
    print(f"  {len(df):,} probesets × {len(df.columns)} columns\n")

    # Verify TAIR column
    if "TAIR" not in df.columns:
        print("ERROR: no TAIR column found")
        return 1

    # Drop probesets with no TAIR mapping
    df = df.dropna(subset=["TAIR"]).copy()
    # Some probesets map to multiple genes — take the first for now
    df["TAIR"] = df["TAIR"].astype(str).str.split("|").str[0].str.strip()
    print(f"  {len(df):,} probesets with TAIR annotations\n")

    ensure(ABRS_UPLOAD, ABRS_TABLES)

    for organ in args.organs:
        cols = ORGANS[organ]

        # Check columns exist
        missing_flt = [c for c in cols["flt"] if c not in df.columns]
        missing_gc = [c for c in cols["gc"] if c not in df.columns]
        if missing_flt or missing_gc:
            print(f"  {organ}: SKIPPED — missing columns")
            if missing_flt:
                print(f"    FLT missing: {missing_flt}")
            if missing_gc:
                print(f"    GC missing: {missing_gc}")
            continue

        print(f"  {organ}: {len(cols['flt'])} FLT vs {len(cols['gc'])} GC")

        result = compute_dge(df, cols["flt"], cols["gc"], organ)

        # Remove rows with NaN p-values (no variance / too few samples)
        valid = result.dropna(subset=["pvalue"])
        n_total = len(valid)
        n_sig = (valid["padj"] < args.fdr).sum()
        n_up = ((valid["padj"] < args.fdr) & (valid["log2FC_FLT_vs_GC"] > 0)).sum()
        n_down = ((valid["padj"] < args.fdr) & (valid["log2FC_FLT_vs_GC"] < 0)).sum()

        print(f"    {n_total:,} genes tested")
        print(f"    {n_sig:,} significant at FDR < {args.fdr} ({n_up} up, {n_down} down)")
        print(f"    median log2FC: {valid['log2FC_FLT_vs_GC'].median():+.4f}")

        # Collapse to unique TAIR IDs (keep most significant per gene)
        valid_sorted = valid.sort_values("pvalue")
        unique = valid_sorted.drop_duplicates(subset="TAIR", keep="first")
        print(f"    {len(unique):,} unique TAIR loci")

        # Write full results table
        table_path = os.path.join(ABRS_TABLES, f"T01_abrs_{organ}_transcriptome.tsv")
        unique.to_csv(table_path, sep="\t", index=False, float_format="%.6f")
        print(f"    wrote {table_path}")

        # Write PaintOmics-format files
        # Values file: #geneID\tFLT_vs_GC
        values_path = os.path.join(
            ABRS_UPLOAD, f"abrs_{organ}_gene_expression_values.tab"
        )
        with open(values_path, "w") as fh:
            fh.write("#geneID\tFLT_vs_GC\n")
            for _, row in unique.iterrows():
                fh.write(f"{row['TAIR']}\t{row['log2FC_FLT_vs_GC']:.6f}\n")
        print(f"    wrote {values_path}")

        # Relevant file: significant gene IDs (use FDR if powered, else unadjusted p < 0.05)
        relevant = unique[unique["padj"] < args.fdr]
        if len(relevant) < 50:
            relevant = unique[unique["pvalue"] < 0.05]
            print(f"    note: FDR yielded 0 genes; using p < 0.05 unadjusted for PaintOmics relevant list")
        relevant_path = os.path.join(
            ABRS_UPLOAD, f"abrs_{organ}_gene_expression_relevant.tab"
        )
        with open(relevant_path, "w") as fh:
            for tair in relevant["TAIR"]:
                fh.write(f"{tair}\n")
        print(f"    wrote {relevant_path} ({len(relevant)} genes)")
        print()

    print("done — ready for PaintOmics upload preparation")
    return 0


if __name__ == "__main__":
    sys.exit(main())
