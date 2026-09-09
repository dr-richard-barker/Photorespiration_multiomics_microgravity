#!/usr/bin/env python3
"""OSD-522 shoot proteome -> PaintOmics format.

Unlike the RNA-seq, the OSD-522 proteomics tables ship already differential: each row
carries `S/G_Log2_ratio` (Spaceflight / Ground log2 ratio), `S/G_P-Value` and
`S/G_Adj_ P-Value` alongside a UniProt `Accession`. So no statistics are recomputed
here — the deposited values are used as they stand, and the work is reconciling the
two shoot fractions into one table.

Two fractions were measured on shoots, soluble (SOL) and membrane (MEM). A protein
detected in both gets one row: the fraction with the more significant adjusted p-value
wins, and every such decision is written to the log so the choice is auditable.

Identifiers stay as UniProt accessions — PaintOmics converts them itself (its own
worked example maps P17182 -> 13806). The AGI locus is carried into the log only, for
cross-referencing against the transcriptome.

Outputs:
  paintomics/upload/proteomics_values.tab    #proteinID <TAB> FLT_vs_GC  (log2 ratio)
  paintomics/upload/proteomics_relevant.tab  one accession per line, adj p < alpha
"""

from __future__ import annotations

import argparse
import os
import sys
import re

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import CACHE, UPLOAD  # noqa: E402

FRACTIONS = {
    "SOL": "GLDS-522_proteomics_GO_Shoot_SOL_Report_20220223_Proteins.csv",
    "MEM": "GLDS-522_proteomics_GO_Shoot_MEM_Report_20220223_Proteins.csv",
}
CONDITION = "FLT_vs_GC"

COL_ACC = "Accession"
COL_DESC = "Description"
COL_LFC = "S/G_Log2_ratio"
COL_P = "S/G_P-Value"
COL_ADJ = "S/G_Adj_P-Value"

# The two deposited fraction tables do not spell their headers identically — the soluble
# file has "S/G_Adj_ P-Value" (a space) and the membrane file "S/G_Adj._P-Value" (a dot),
# and the peptide-count columns differ the same way. Match on a normalised key instead of
# the literal string, so a header typo in either file does not silently drop a fraction.
def _norm(col: str) -> str:
    return re.sub(r"[^a-z0-9]", "", col.lower())


def resolve(df: pd.DataFrame, wanted: str) -> str | None:
    """Find the column whose normalised name matches `wanted`."""
    target = _norm(wanted)
    for col in df.columns:
        if _norm(col) == target:
            return col
    return None


def load_fraction(name: str, path: str) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    cols = {key: resolve(df, key) for key in (COL_ACC, COL_LFC, COL_ADJ, COL_P, COL_DESC)}
    missing = [k for k in (COL_ACC, COL_LFC, COL_ADJ) if cols[k] is None]
    if missing:
        raise SystemExit(f"{os.path.basename(path)}: missing columns {missing}\n"
                         f"  columns present: {list(df.columns)[:20]}")

    out = pd.DataFrame({
        "accession": df[cols[COL_ACC]].astype(str).str.strip(),
        "log2fc": pd.to_numeric(df[cols[COL_LFC]], errors="coerce"),
        "pvalue": pd.to_numeric(df[cols[COL_P]], errors="coerce")
                  if cols[COL_P] else np.nan,
        "adj_p": pd.to_numeric(df[cols[COL_ADJ]], errors="coerce"),
        "description": df[cols[COL_DESC]].astype(str) if cols[COL_DESC] else "",
        "fraction": name,
    })
    # 'GN=NAME' inside the UniProt-style description; kept for the log only.
    out["gene_name"] = out["description"].str.extract(r"GN=(\S+)", expand=False)
    out = out[out["accession"].str.len() > 0]
    out = out.dropna(subset=["log2fc"])
    print(f"  {name}: {len(df):,} rows -> {len(out):,} with a usable log2 ratio, "
          f"{out['accession'].nunique():,} unique accessions")
    return out


def reconcile(frames: list[pd.DataFrame], log_path: str) -> pd.DataFrame:
    """One row per accession; on collision keep the more significant adjusted p-value."""
    allrows = pd.concat(frames, ignore_index=True)

    # rank: smaller adjusted p wins; NaN adjusted p sorts last
    allrows["_rank"] = allrows["adj_p"].fillna(np.inf)
    allrows = allrows.sort_values(["accession", "_rank"], kind="mergesort")

    dup_acc = allrows["accession"][allrows["accession"].duplicated()].unique()
    kept = allrows.drop_duplicates(subset="accession", keep="first").set_index("accession")

    with open(log_path, "w") as fh:
        fh.write("# OSD-522 shoot proteome — SOL/MEM reconciliation\n")
        fh.write("# proteins seen in both fractions; the row with the smaller adjusted "
                 "p-value was kept\n")
        fh.write("accession\tgene\tkept_fraction\tkept_log2fc\tkept_adj_p\tdropped\n")
        for acc in dup_acc:
            rows = allrows[allrows["accession"] == acc]
            k = rows.iloc[0]
            dropped = "; ".join(
                f"{r.fraction} log2fc={r.log2fc:.3f} adj_p={r.adj_p:.3g}"
                for r in rows.iloc[1:].itertuples()
            )
            fh.write(f"{acc}\t{k.gene_name}\t{k.fraction}\t{k.log2fc:.4f}\t"
                     f"{k.adj_p:.4g}\t{dropped}\n")

    print(f"  reconciled: {len(kept):,} unique proteins "
          f"({len(dup_acc):,} seen in both fractions, logged)")
    return kept


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=0.05,
                    help="adjusted p-value threshold for the relevant-features file")
    args = ap.parse_args()

    print("OSD-522 shoot proteome — Space Flight vs Ground Control (deposited ratios)")
    frames = []
    for name, filename in FRACTIONS.items():
        path = os.path.join(CACHE, filename)
        if not os.path.exists(path):
            raise SystemExit(f"missing {path} — run osdr/fetch_osd522.py first")
        frames.append(load_fraction(name, path))

    os.makedirs(UPLOAD, exist_ok=True)
    prot = reconcile(frames, os.path.join(UPLOAD, "proteomics_RECONCILIATION.log"))

    values_path = os.path.join(UPLOAD, "proteomics_values.tab")
    with open(values_path, "w") as fh:
        fh.write(f"#proteinID\t{CONDITION}\n")
        for acc, row in prot.iterrows():
            fh.write(f"{acc}\t{row['log2fc']:.6f}\n")

    sig = prot[prot["adj_p"] < args.alpha]
    relevant_path = os.path.join(UPLOAD, "proteomics_relevant.tab")
    with open(relevant_path, "w") as fh:
        for acc in sig.index:
            fh.write(f"{acc}\n")

    lfc = prot["log2fc"]
    up, down = int((sig["log2fc"] > 0).sum()), int((sig["log2fc"] < 0).sum())
    print(f"  log2 ratio: median {lfc.median():+.4f}, mean {lfc.mean():+.4f}, "
          f"sd {lfc.std():.3f}, range {lfc.min():.2f}..{lfc.max():.2f}")
    print(f"  wrote {values_path}  ({len(prot):,} proteins)")
    print(f"  wrote {relevant_path}  ({len(sig):,} at adj p < {args.alpha}; "
          f"{up:,} up in flight, {down:,} down)")

    with open(os.path.join(UPLOAD, "proteomics_METHOD.txt"), "w") as fh:
        fh.write(
            "OSD-522 (BRIC-LED-001) Arabidopsis shoots — proteome, Space Flight vs Ground Control\n"
            f"source:   NASA OSDR {', '.join(FRACTIONS.values())}\n"
            "values:   the deposited 'S/G_Log2_ratio' column, used as-is — no statistics recomputed\n"
            "fractions: shoot soluble + shoot membrane merged; on collision the row with the\n"
            "           smaller adjusted p-value was kept (see proteomics_RECONCILIATION.log)\n"
            "ids:      UniProt accessions; PaintOmics maps them to Entrez itself\n"
            f"contrast: {CONDITION} = log2(Space Flight / Ground Control)\n"
            f"relevant: deposited 'S/G_Adj_ P-Value' < {args.alpha} "
            f"({len(sig)} proteins: {up} up, {down} down)\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
