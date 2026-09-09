#!/usr/bin/env python3
"""Score every study in the hardware ladder on the same gene sets.

The model's claim is not about one experiment, it is an ORDERING. LunarLeaf-CFD T11 puts
12 h carbon gain, as a percentage of Earth's, at 1 % for a sealed BRIC canister, 90 % under
CARA's micropore tape and 100 % when vented. If the carbon-starvation signature in flown
plants is driven by the enclosure rather than by microgravity, it should be strongest in
sealed hardware, intermediate under tape, and weakest when vented — and it should be absent
in a DARK enclosure of any kind, because without photosynthesis there is no CO2 drawdown.

`data/study_registry.tsv` decides which studies enter and names the exact matched
flight-vs-ground contrast column for each. Nothing here picks a contrast on its own.

Outputs
    results/tables/T11_hardware_ladder.tsv    every study x gene set
    results/tables/T12_ladder_contrasts.tsv   the set-vs-set separations, per study
"""

from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genesets  # noqa: E402
from paths import CACHE, STUDY_REGISTRY, TABLES, UPLOAD, ensure  # noqa: E402

DOWNLOAD = ("https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download"
            "?source=datamanager&file={f}")
CLASS_ORDER = ["sealed", "tape", "vented"]


def load_registry() -> pd.DataFrame:
    rows = []
    with open(STUDY_REGISTRY) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.rstrip("\n").split("\t"))
    df = pd.DataFrame(rows[1:], columns=rows[0])
    df["carbon_pct_earth"] = pd.to_numeric(df["carbon_pct_earth"], errors="coerce")
    return df


def fetch_dge(accession: str, filename: str) -> pd.DataFrame:
    """Download (and cache) a GeneLab differential-expression table."""
    # Registry ids like OSD-678-light name a contrast within one study, not a new accession.
    osd = "-".join(accession.split("-")[:2])
    path = os.path.join(CACHE, f"dge_{osd}.csv")
    if not os.path.exists(path):
        ensure(CACHE)
        url = DOWNLOAD.format(osd=osd, f=urllib.parse.quote(filename))
        req = urllib.request.Request(url, headers={"User-Agent": "photoresp-multiomics/1.0"})
        print(f"    fetching {filename}")
        with urllib.request.urlopen(req, timeout=300) as resp:
            body = resp.read()
        with open(path, "wb") as fh:
            fh.write(body)
    return pd.read_csv(path, low_memory=False)


def series_for(row: pd.Series) -> pd.Series:
    """AGI-indexed log2 fold change for one registry row."""
    acc, contrast = row["accession"], row["contrast"]

    if acc == "OSD-522":
        # Computed here rather than shipped by GeneLab; already in PaintOmics format.
        v = pd.read_csv(os.path.join(UPLOAD, "gene_expression_values.tab"),
                        sep="\t", index_col=0)
        return v.iloc[:, 0]

    df = fetch_dge(acc, row["dge_file"])
    col = f"Log2fc_({contrast.split('Log2fc_(')[1]}" if contrast.startswith("Log2fc_") \
        else contrast
    if col not in df.columns:
        raise SystemExit(f"{acc}: contrast column not present: {col}")

    idcol = df.columns[0]
    if idcol.upper() not in ("TAIR", "ENSEMBL", "GENE"):
        print(f"    note: {acc} identifier column is {idcol!r}")
    s = pd.Series(pd.to_numeric(df[col], errors="coerce").values,
                  index=df[idcol].astype(str).str.strip())
    s = s[~s.index.duplicated(keep="first")].dropna()
    return s


def score(values: pd.Series, members: set[str], name: str) -> dict:
    inside = values[values.index.isin(members)].dropna()
    outside = values[~values.index.isin(members)].dropna()
    if len(inside) < 3:
        return {"set": name, "n": len(inside), "median": np.nan, "shift": np.nan,
                "p": np.nan}
    _u, p = stats.mannwhitneyu(inside, outside, alternative="two-sided")
    return {"set": name, "n": len(inside), "median": float(inside.median()),
            "shift": float(inside.median() - outside.median()), "p": float(p)}


def separations(values: pd.Series, sets: dict[str, set[str]]) -> list[dict]:
    """Set-vs-set separations — immune to the global metabolic shift."""
    out = []
    for label, a_name, b_name, _why in genesets.CONTRASTS:
        a = values[values.index.isin(sets.get(a_name, set()))].dropna()
        b = values[values.index.isin(sets.get(b_name, set()))].dropna()
        if len(a) < 3 or len(b) < 3:
            out.append({"contrast": label, "n_a": len(a), "n_b": len(b),
                        "separation": np.nan, "p": np.nan})
            continue
        _u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        out.append({"contrast": label, "n_a": len(a), "n_b": len(b),
                    "separation": float(a.median() - b.median()), "p": float(p)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()

    reg = load_registry()
    included = reg[reg["include"] == "yes"]
    excluded = reg[reg["include"] != "yes"]

    print("Hardware ladder — do the enclosures order the response?")
    print(f"  registry: {len(included)} included, {len(excluded)} excluded")
    for _, r in excluded.iterrows():
        print(f"    excluded {r['accession']:<14} {r['reason'][:74]}")

    print("\n  building gene sets from KEGG")
    sets = genesets.build()

    set_rows, sep_rows = [], []
    for _, r in included.iterrows():
        print(f"\n  {r['accession']}  [{r['enclosure_class']}, {r['light']}] "
              f"{r['hardware']} — {r['tissue']}")
        try:
            values = series_for(r)
        except Exception as exc:  # noqa: BLE001
            print(f"    FAILED: {exc}")
            continue
        print(f"    {len(values):,} genes, median log2FC {values.median():+.4f}")

        for name, members in sets.items():
            s = score(values, members, name)
            s.update(accession=r["accession"], enclosure_class=r["enclosure_class"],
                     carbon_pct_earth=r["carbon_pct_earth"], light=r["light"],
                     hardware=r["hardware"], tissue=r["tissue"],
                     predicted=genesets.PREDICTION.get(name))
            set_rows.append(s)

        for sep in separations(values, sets):
            sep.update(accession=r["accession"], enclosure_class=r["enclosure_class"],
                       carbon_pct_earth=r["carbon_pct_earth"], light=r["light"],
                       hardware=r["hardware"])
            sep_rows.append(sep)

    ensure(TABLES)
    cols = ["accession", "enclosure_class", "carbon_pct_earth", "light", "hardware",
            "tissue", "set", "n", "median", "shift", "p", "predicted"]
    t11 = pd.DataFrame(set_rows)[cols]
    t11.to_csv(os.path.join(TABLES, "T11_hardware_ladder.tsv"), sep="\t", index=False)

    scols = ["accession", "enclosure_class", "carbon_pct_earth", "light", "hardware",
             "contrast", "n_a", "n_b", "separation", "p"]
    t12 = pd.DataFrame(sep_rows)[scols]
    t12.to_csv(os.path.join(TABLES, "T12_ladder_contrasts.tsv"), sep="\t", index=False)

    # ---- the two ordered questions, reported whichever way they fall -------------
    key = "starvation vs photosynthesis"
    sub = t12[t12["contrast"] == key].copy()
    print(f"\n  '{key}' separation, by illumination and enclosure")
    print(f"  {'study':<16} {'class':<8} {'light':<6} {'sep':>9} {'p':>10}")
    print("  " + "-" * 55)
    for _, r in sub.sort_values(["light", "carbon_pct_earth"]).iterrows():
        p = "n/a" if not np.isfinite(r["p"]) else f"{r['p']:.2e}"
        sep = "n/a" if not np.isfinite(r["separation"]) else f"{r['separation']:+.4f}"
        print(f"  {r['accession']:<16} {r['enclosure_class']:<8} {r['light']:<6} "
              f"{sep:>9} {p:>10}")

    lit = sub[sub["light"] == "light"].dropna(subset=["separation"])
    dark = sub[sub["light"] == "dark"].dropna(subset=["separation"])

    # TEST 1 — illumination. The model says the carbon-starvation signature needs
    # photosynthesis to have drawn the enclosure down, so it should appear only in the light.
    print("\n  TEST 1 — illumination")
    print(f"    lit  ({len(lit)}): {', '.join(f'{v:+.2f}' for v in lit['separation'])}")
    print(f"    dark ({len(dark)}): {', '.join(f'{v:+.2f}' for v in dark['separation'])}")
    if len(lit) and len(dark):
        if lit["separation"].min() > dark["separation"].max():
            print("    SEPARATES PERFECTLY: every lit study sits above every dark study.")
        else:
            print("    does not separate cleanly.")
        if len(lit) >= 3 and len(dark) >= 3:
            _u, pmw = stats.mannwhitneyu(lit["separation"], dark["separation"],
                                         alternative="greater")
            print(f"    one-sided Mann-Whitney (lit > dark): p = {pmw:.4f} "
                  f"(n={len(lit)} vs {len(dark)}; 0.05 is unreachable below n=4 per group)")

    # TEST 2 — enclosure gradient among lit studies only.
    print("\n  TEST 2 — enclosure gradient, lit studies only")
    if len(lit) >= 3:
        rho, prho = stats.spearmanr(lit["carbon_pct_earth"], lit["separation"])
        for _, r in lit.sort_values("carbon_pct_earth").iterrows():
            print(f"    {r['enclosure_class']:<8} carbon {r['carbon_pct_earth']:>3.0f}% "
                  f"-> separation {r['separation']:+.3f}  ({r['accession']})")
        print(f"    Spearman rho = {rho:+.3f}, p = {prho:.3f} "
              f"(model predicts negative rho)")
    else:
        print("    too few lit studies to rank")

    print(f"\n  wrote {TABLES}/T11_hardware_ladder.tsv and T12_ladder_contrasts.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
