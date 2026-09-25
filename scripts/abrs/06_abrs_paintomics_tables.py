#!/usr/bin/env python3
"""Derive ABRS PaintOmics tables and cross-hardware comparison.

Produces:
  results/abrs/tables/T08_abrs_paintomics_significant.tsv
    All pathways significant at Fisher combined p < 0.05 in ABRS shoot.
  results/abrs/tables/T09_abrs_paintomics_carbon.tsv
    Evaluation of the 11 canonical carbon-related pathways in ABRS.
  results/abrs/tables/T10_abrs_vs_bric_comparison.tsv
    Side-by-side comparison of pathway significance between sealed BRIC-LED
    and ventilated ABRS hardware, elucidating true microgravity effects
    versus hardware-induced confounds.
"""

from __future__ import annotations

import glob
import gzip
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))
import genesets  # noqa: E402
import mapman  # noqa: E402
from paths import DATA, REPO, TABLES, ensure  # noqa: E402

ABRS_TABLES = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "tables"
))
ABRS_RAW = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "paintomics_raw"
))
ABRS_UPLOAD = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "data", "abrs", "paintomics_upload"
))

BRIC_T08 = os.path.join(TABLES, "T08_paintomics_significant.tsv")
ORGANISM_SUFFIX = re.compile(r"\s*-\s*Arabidopsis thaliana.*$", re.I)


def kegg_clean_names() -> dict[str, str]:
    """Return pid -> exact KEGG pathway name (without organism suffix)."""
    out = {}
    for line in genesets.kegg("list/pathway/ath", "kegg_ath_pathway_names.tsv").splitlines():
        if "\t" not in line:
            continue
        pid, name = line.split("\t")
        clean_pid = pid.replace("path:", "").strip()
        clean_name = ORGANISM_SUFFIX.sub("", name).strip()
        out[clean_pid] = clean_name
    return out


def run_local_pathway_enrichment(tissue: str = "shoot") -> pd.DataFrame:
    """Compute Fisher exact enrichment across KEGG pathways and MapMan diagrams."""
    val_file = os.path.join(ABRS_UPLOAD, f"abrs_{tissue}_gene_expression_values.tab")
    rel_file = os.path.join(ABRS_UPLOAD, f"abrs_{tissue}_gene_expression_relevant.tab")

    vals = pd.read_csv(val_file, sep="\t", index_col=0)
    all_genes = set(vals.index)
    with open(rel_file) as fh:
        rel_genes = {line.strip() for line in fh if line.strip()}

    N = len(all_genes)
    K = len(rel_genes)
    genome_median = vals.iloc[:, 0].median()

    rows = []

    # 1. KEGG pathways
    pid_to_name = kegg_clean_names()
    for pid, clean_name in pid_to_name.items():
        try:
            pgenes = genesets.pathway_genes(pid)
        except Exception:
            continue
        p_universe = pgenes.intersection(all_genes)
        if len(p_universe) < 3:
            continue

        p_rel = p_universe.intersection(rel_genes)
        k_in = len(p_rel)
        k_out = len(p_universe) - k_in
        rest_in = K - k_in
        rest_out = (N - len(p_universe)) - rest_in

        _, p_val = stats.fisher_exact([[k_in, rest_in], [k_out, rest_out]], alternative="greater")

        lfc_p = vals.loc[vals.index.isin(p_universe)].iloc[:, 0]
        direction = "up" if lfc_p.median() > genome_median else "down"

        rows.append({
            "db": "K",
            "pathway_id": pid,
            "pathway": clean_name,
            "features": len(p_universe),
            "relevant_features": k_in,
            "median_log2FC": lfc_p.median(),
            "direction": direction,
            "p_gene": p_val,
        })

    # 2. MapMan diagrams
    mm_mapping = mapman.bin_to_genes()
    for dname in sorted(mapman.available()):
        dgenes = mapman.diagram_genes(dname, mm_mapping)
        p_universe = dgenes.intersection(all_genes)
        if len(p_universe) < 3:
            continue
        p_rel = p_universe.intersection(rel_genes)
        k_in = len(p_rel)
        k_out = len(p_universe) - k_in
        rest_in = K - k_in
        rest_out = (N - len(p_universe)) - rest_in

        _, p_val = stats.fisher_exact([[k_in, rest_in], [k_out, rest_out]], alternative="greater")
        lfc_p = vals.loc[vals.index.isin(p_universe)].iloc[:, 0]
        direction = "up" if lfc_p.median() > genome_median else "down"

        rows.append({
            "db": "M",
            "pathway_id": "-",
            "pathway": dname,
            "features": len(p_universe),
            "relevant_features": k_in,
            "median_log2FC": lfc_p.median(),
            "direction": direction,
            "p_gene": p_val,
        })

    df = pd.DataFrame(rows).sort_values("p_gene").reset_index(drop=True)
    return df


def main() -> int:
    ensure(ABRS_TABLES)

    print("Generating ABRS Pathway Enrichment and Cross-Hardware Comparison Tables\n")

    df_paths = run_local_pathway_enrichment("shoot")
    sig_paths = df_paths[df_paths["p_gene"] < 0.05].copy()

    sig_paths["p_gene_fmt"] = [f"{p:.4e}" if p < 1e-3 else f"{p:.5f}" for p in sig_paths["p_gene"]]
    sig_paths["median_lfc_fmt"] = [f"{lfc:+.4f}" for lfc in sig_paths["median_log2FC"]]

    t08_path = os.path.join(ABRS_TABLES, "T08_abrs_paintomics_significant.tsv")
    cols = ["db", "pathway", "features", "relevant_features", "median_lfc_fmt", "direction", "p_gene_fmt"]
    with open(t08_path, "w") as fh:
        fh.write("# ABRS (OSD-7 Shoot) Pathway Enrichment in Spaceflight\n")
        fh.write(f"# Total pathways tested: {len(df_paths)}, Significant at p < 0.05: {len(sig_paths)}\n")
        fh.write("db\tpathway\tfeatures\trelevant_features\tmedian_log2FC\tdirection\tp_gene\n")
        for _, r in sig_paths[cols].iterrows():
            fh.write("\t".join(str(x) for x in r.values) + "\n")
    print(f"  wrote {t08_path} ({len(sig_paths)} significant pathways out of {len(df_paths)} tested)")

    # T09: Carbon claims audit in ABRS
    claims = pd.read_csv(os.path.join(DATA, "paintomics_carbon_claims.tsv"), sep="\t")
    t09_rows = []
    for _, c in claims.iterrows():
        pname = c["pathway"]
        db = c["db"]
        match = df_paths[(df_paths["db"] == db) & (df_paths["pathway"].str.lower() == pname.lower())]
        if match.empty:
            match = df_paths[(df_paths["db"] == db) & (df_paths["pathway"].str.lower().str.contains(pname.lower()[:12]))]

        if not match.empty:
            m = match.iloc[0]
            rank = int(m.name) + 1
            t09_rows.append({
                "db": db,
                "pathway": pname,
                "predicted": c["predicted"],
                "features": m["features"],
                "relevant": m["relevant_features"],
                "median_log2FC": f"{m['median_log2FC']:+.4f}",
                "p_gene": f"{m['p_gene']:.4e}" if m['p_gene'] < 1e-3 else f"{m['p_gene']:.5f}",
                "abrs_rank": f"{rank}/{len(df_paths)}",
                "bric_outcome": c["outcome"],
                "abrs_outcome": "significant" if m["p_gene"] < 0.05 else "ns",
            })
        else:
            t09_rows.append({
                "db": db,
                "pathway": pname,
                "predicted": c["predicted"],
                "features": "-",
                "relevant": "-",
                "median_log2FC": "-",
                "p_gene": "-",
                "abrs_rank": "-",
                "bric_outcome": c["outcome"],
                "abrs_outcome": "not in local map",
            })

    t09_path = os.path.join(ABRS_TABLES, "T09_abrs_paintomics_carbon.tsv")
    pd.DataFrame(t09_rows).to_csv(t09_path, sep="\t", index=False)
    print(f"  wrote {t09_path}")

    # T10: Direct Side-by-Side Comparison of BRIC-LED vs ABRS
    if os.path.exists(BRIC_T08):
        bric_df = pd.read_csv(BRIC_T08, sep="\t", comment="#")
        comp_rows = []
        for _, b in bric_df.iterrows():
            bpname = b["pathway"]
            db = b["db"]
            match = df_paths[(df_paths["db"] == db) & (df_paths["pathway"].str.lower() == bpname.lower())]
            if match.empty:
                match = df_paths[(df_paths["db"] == db) & (df_paths["pathway"].str.lower().str.contains(bpname.lower()[:12]))]

            if not match.empty:
                m = match.iloc[0]
                abrs_p = m["p_gene"]
                abrs_dir = m["direction"]
                abrs_lfc = f"{m['median_log2FC']:+.4f}"
            else:
                abrs_p = np.nan
                abrs_dir = "-"
                abrs_lfc = "-"

            bric_p = float(b["p_combined_fisher"])
            if bric_p < 0.05 and (pd.isna(abrs_p) or abrs_p >= 0.05):
                interp = "BRIC-LED hardware-specific (CO2 starvation / sealed volatiles)"
            elif bric_p < 0.05 and abrs_p < 0.05:
                interp = "Conserved spaceflight response (Hardware-independent)"
            elif bric_p >= 0.05 and abrs_p < 0.05:
                interp = "ABRS-specific response"
            else:
                interp = "Non-significant in both"

            comp_rows.append({
                "db": db,
                "pathway": bpname,
                "bric_p_combined": f"{bric_p:.4e}" if bric_p < 1e-3 else f"{bric_p:.5f}",
                "abrs_p_gene": f"{abrs_p:.4e}" if not pd.isna(abrs_p) and abrs_p < 1e-3 else (f"{abrs_p:.5f}" if not pd.isna(abrs_p) else "-"),
                "abrs_median_log2FC": abrs_lfc,
                "abrs_direction": abrs_dir,
                "classification": interp,
            })

        t10_path = os.path.join(ABRS_TABLES, "T10_abrs_vs_bric_comparison.tsv")
        pd.DataFrame(comp_rows).to_csv(t10_path, sep="\t", index=False)
        print(f"  wrote {t10_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
