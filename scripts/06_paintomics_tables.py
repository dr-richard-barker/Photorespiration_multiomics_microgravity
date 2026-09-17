#!/usr/bin/env python3
"""Derive T08, T09 and T10 from the PaintOmics job's own saved result.

These three tables used to be transcribed by hand from the PaintOmics web interface, and
MANIFEST.tsv credited them to `05_submit_paintomics.py`, which only submits a job and
saves its acknowledgement. Nothing in the repository could regenerate or check them, and
`check_job_status` is one-shot — the saved `job1_status.json` is the
"not on the queue anymore" reply, with no results in it.

The full result was recovered afterwards from `POST /pa_recover_job` with the job id, and
is vendored at `results/paintomics_raw/job1_m1z16Qg3DK_full.json.gz`. Everything measured
in these tables now comes from that file. What cannot come from it — which pathways the
model made a prior claim about, and which metabolic block each hub compound belongs to —
is curation, and lives in two small files under `data/` that this script joins in and
fails loudly on if a row is missing.

Re-deriving also corrected two counts that had been wrong everywhere they appeared: the
job tested 232 pathways, not 231, and 32 reached combined Fisher p < 0.05, not 31. The
missing one was KEGG's global "Metabolic pathways" map (ath01100, p = 0.0119).
"""

from __future__ import annotations

import gzip
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DATA, REPO, TABLES, ensure  # noqa: E402

JOB_ID = "m1z16Qg3DK"
RAW = os.path.join(REPO, "results", "paintomics_raw", f"job1_{JOB_ID}_full.json.gz")
ALPHA = 0.05
OMIC = {"p_gene": "Gene expression", "p_protein": "Proteomics",
        "p_metabolite": "Metabolomics"}

HEADER = [
    f"# PaintOmics AI job {JOB_ID} — OSD-522 (BRIC-LED-001) Arabidopsis shoots, "
    f"Space Flight vs Ground Control",
    "# organism ath; databases KEGG + MapMan; AI interpretation off",
    "# 3 omics: gene expression (real, OSDR), proteomics (real, OSDR), "
    "metabolomics (PREDICTED, model output)",
]


def load_job() -> dict:
    if not os.path.exists(RAW):
        raise SystemExit(
            f"{RAW} is missing. Recover it with a POST to "
            f"https://paintomics.org/pa_recover_job with jobID={JOB_ID}.")
    with gzip.open(RAW, "rt") as fh:
        return json.load(fh)


def pathway_frame(job: dict) -> pd.DataFrame:
    rows = []
    for p in job["pathwaysInfo"]:
        g = p.get("globalOmicPvalues") or {}
        row = {
            "db": "K" if p.get("source") == "KEGG" else "M",
            "pathway": p["name"],
            "features": len(p.get("matchedGenes") or []),
            "unique_metabolites": len(p.get("matchedCompounds") or []),
            "p_combined_fisher": p["combinedSignificancePvalues"]["Fisher"][0],
        }
        for col, key in OMIC.items():
            row[col] = g.get(key)
        rows.append(row)
    return pd.DataFrame(rows).sort_values("p_combined_fisher").reset_index(drop=True)


def pvalue(v) -> str:
    """The format these tables have always used: five decimals down to 1e-3, then
    scientific with an unpadded exponent, which keeps five significant figures either
    way. An absent omic layer prints as a dash."""
    if pd.isna(v):
        return "-"
    if v >= 1e-3:
        return f"{v:.5f}"
    mant, exp = f"{v:.4e}".split("e")
    return f"{mant}e{int(exp)}"


def fmt(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in list(OMIC) + ["p_combined_fisher"]:
        out[c] = [pvalue(v) for v in out[c]]
    return out


def write(path: str, df: pd.DataFrame, extra: list[str]) -> None:
    with open(path, "w") as fh:
        for line in HEADER + extra:
            fh.write(line + "\n")
        df.to_csv(fh, sep="\t", index=False)
    print(f"  wrote {os.path.relpath(path, REPO)}  ({len(df)} rows)")


def main() -> int:
    ensure(TABLES)
    job = load_job()
    paths = pathway_frame(job)
    total = len(paths)
    sig = paths[paths["p_combined_fisher"] < ALPHA]

    cols = ["db", "pathway", "features", "unique_metabolites",
            "p_gene", "p_protein", "p_metabolite", "p_combined_fisher"]
    write(os.path.join(TABLES, "T08_paintomics_significant.tsv"),
          fmt(sig)[cols],
          [f"# {total} pathways tested, {len(sig)} significant at combined Fisher "
           f"p < {ALPHA}. Sorted by combined p.",
           "# Derived from results/paintomics_raw/ by scripts/06_paintomics_tables.py."])

    # --- T09: the pathways the model made a prior claim about, significant or not -----
    claims = pd.read_csv(os.path.join(DATA, "paintomics_carbon_claims.tsv"), sep="\t")
    known = set(paths["pathway"])
    missing = [p for p in claims["pathway"] if p not in known]
    if missing:
        raise SystemExit(f"carbon claims name pathways the job does not have: {missing}")
    carbon = paths[paths["pathway"].isin(set(claims["pathway"]))].merge(
        claims[["pathway", "predicted", "outcome"]], on="pathway", how="left")
    gly = carbon[carbon["pathway"].str.contains("Glyoxylate")]
    rank = int(paths.index[paths["pathway"].str.contains("Glyoxylate")][0]) + 1
    write(os.path.join(TABLES, "T09_paintomics_carbon.tsv"),
          fmt(carbon)[cols + ["predicted", "outcome"]],
          [f"# The carbon pathways the model made claims about, significant or not.",
           f"# Glyoxylate and dicarboxylate metabolism (photorespiration) ranks {rank} "
           f"of {total} by combined p,",
           f"# at p = {float(gly['p_combined_fisher'].iloc[0]):.5g}: in the bottom "
           f"{100 * (total - rank + 1) / total:.0f}% but NOT last — "
           f"{total - rank} pathways rank below it.",
           "# Derived from results/paintomics_raw/ by scripts/06_paintomics_tables.py."])

    # --- T10: metabolite hub analysis -------------------------------------------------
    hub = pd.DataFrame(job["hubAnalysisResult"].values())
    best = (hub.sort_values("pvalue_adjust").groupby("name", as_index=False).first()
            .rename(columns={"name": "kegg_id", "pvalue_adjust": "FDR",
                             "DEN": "DE_neighbours", "step": "best_at_step"}))
    blocks = pd.read_csv(os.path.join(DATA, "paintomics_hub_blocks.tsv"), sep="\t")
    missing = [c for c in best["kegg_id"] if c not in set(blocks["kegg_id"])]
    if missing:
        raise SystemExit(f"hub compounds with no curated block: {missing}")
    # Five compounds tie at the same adjusted p, so rank ties by how much evidence sits
    # behind them rather than leaving the order to whatever the dict happened to yield.
    t10 = (best.merge(blocks, on="kegg_id", how="left")
           .sort_values(["FDR", "DE_neighbours"], ascending=[True, False])
           [["compound", "kegg_id", "FDR", "DE_neighbours", "best_at_step", "block"]])
    t10["FDR"] = t10["FDR"].round(3)
    n_sig = int((t10["FDR"] < ALPHA).sum())
    write(os.path.join(TABLES, "T10_metabolite_hubs.tsv"), t10,
          [f"# PaintOmics metabolite hub analysis (KEGG network): which metabolites have "
           f"differentially",
           f"# expressed genes concentrated around them. {len(t10)} of the 21 predicted "
           f"compounds entered",
           f"# the analysis; {n_sig} reach FDR < {ALPHA}. The DE genes counted are REAL; "
           f"the compounds are model output.",
           "# Derived from results/paintomics_raw/ by scripts/06_paintomics_tables.py."])

    print(f"\n  {total} pathways tested, {len(sig)} significant; "
          f"photorespiration ranks {rank} of {total}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
