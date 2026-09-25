#!/usr/bin/env python3
"""Extended Hardware Ladder: Integrating ABRS (OSD-7) with previous 6-study ladder.

The previous analysis evaluated 6 studies across sealed BRIC, tape CARA, and
vented VEGGIE hardware. This script extends the ladder by incorporating OSD-7
(ABRS / TAGES), providing the highest-quality ventilated hardware arm with
controlled airflow and CO2 regulation.

Contrasts evaluated across all studies:
  1. Starvation vs Photosynthesis (DIN starvation markers vs photosynthesis apparatus)
  2. Photorespiration vs Carbon Fixation (ath00630 vs ath00710)
  3. Photorespiration Core vs Rubisco (curated enzyme sets)

Outputs:
  results/abrs/tables/T11_abrs_extended_hardware_ladder.tsv
  results/abrs/tables/T12_abrs_extended_ladder_contrasts.tsv
  results/abrs/tables/T13_abrs_organ_comparison.tsv
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))
import genesets  # noqa: E402
from paths import TABLES, ensure  # noqa: E402

ABRS_TABLES = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "tables"
))

PREV_T11 = os.path.join(TABLES, "T11_hardware_ladder.tsv")
PREV_T12 = os.path.join(TABLES, "T12_ladder_contrasts.tsv")


def score(values: pd.Series, members: set[str], name: str) -> dict:
    inside = values[values.index.isin(members)].dropna()
    outside = values[~values.index.isin(members)].dropna()
    if len(inside) < 3:
        return {"set": name, "n": len(inside), "median": np.nan, "shift": np.nan, "p": np.nan}
    _u, p = stats.mannwhitneyu(inside, outside, alternative="two-sided")
    return {"set": name, "n": len(inside), "median": float(inside.median()),
            "shift": float(inside.median() - outside.median()), "p": float(p)}


def separations(values: pd.Series, sets: dict[str, set[str]]) -> list[dict]:
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
    ensure(ABRS_TABLES)

    print("Extended Hardware Ladder: Integrating ABRS (OSD-7)\n")

    # 1. Load gene sets
    sets = genesets.build()
    sets["ath04141"] = genesets.pathway_genes("ath04141")
    ann = genesets.annotation()
    cw = ann[ann["description"].str.contains("xyloglucan|pectin|expansin|cellulose synthase", case=False, regex=True, na=False)]
    sets["cell_wall"] = set(cw.index)

    # 2. Score ABRS organs
    abrs_organs = [
        ("OSD-7-shoot", "shoot", "vented", 100, "light", "ABRS", "shoots"),
        ("OSD-7-root", "root", "vented", 100, "light", "ABRS", "roots"),
        ("OSD-7-hypocotyl", "hypocotyl", "vented", 100, "light", "ABRS", "hypocotyls"),
        ("OSD-7-whole_plant", "whole_plant", "vented", 100, "light", "ABRS", "whole plant"),
    ]

    new_set_rows = []
    new_sep_rows = []
    organ_records = []

    for acc, organ, enc, carbon, light, hw, tissue in abrs_organs:
        t_path = os.path.join(ABRS_TABLES, f"T01_abrs_{organ}_transcriptome.tsv")
        if not os.path.exists(t_path):
            continue
        df = pd.read_csv(t_path, sep="\t")
        vals = pd.Series(df["log2FC_FLT_vs_GC"].values, index=df["TAIR"].values).dropna()

        # Score sets
        for sname, members in sets.items():
            sc = score(vals, members, sname)
            sc.update(accession=acc, enclosure_class=enc, carbon_pct_earth=carbon,
                      light=light, hardware=hw, tissue=tissue,
                      predicted=genesets.PREDICTION.get(sname))
            new_set_rows.append(sc)
            organ_records.append({"accession": acc, "organ": organ, "set": sname, **sc})

        # Separations
        for sep in separations(vals, sets):
            sep.update(accession=acc, enclosure_class=enc, carbon_pct_earth=carbon,
                       light=light, hardware=hw)
            new_sep_rows.append(sep)

    # 3. Combine with previous ladder tables
    if os.path.exists(PREV_T11):
        prev_t11 = pd.read_csv(PREV_T11, sep="\t")
        # Add OSD-7-shoot to the primary hardware ladder
        shoot_set_rows = [r for r in new_set_rows if r["accession"] == "OSD-7-shoot"]
        combined_t11 = pd.concat([prev_t11, pd.DataFrame(shoot_set_rows)], ignore_index=True)
        ext_t11_path = os.path.join(ABRS_TABLES, "T11_abrs_extended_hardware_ladder.tsv")
        combined_t11.to_csv(ext_t11_path, sep="\t", index=False)
        print(f"  wrote {ext_t11_path} ({len(combined_t11)} rows)")

    if os.path.exists(PREV_T12):
        prev_t12 = pd.read_csv(PREV_T12, sep="\t")
        shoot_sep_rows = [r for r in new_sep_rows if r["accession"] == "OSD-7-shoot"]
        combined_t12 = pd.concat([prev_t12, pd.DataFrame(shoot_sep_rows)], ignore_index=True)
        ext_t12_path = os.path.join(ABRS_TABLES, "T12_abrs_extended_ladder_contrasts.tsv")
        combined_t12.to_csv(ext_t12_path, sep="\t", index=False)
        print(f"  wrote {ext_t12_path} ({len(combined_t12)} rows)")

        # Summary of illumination and enclosure tests on extended ladder
        key = "starvation vs photosynthesis"
        sub = combined_t12[combined_t12["contrast"] == key].sort_values(["light", "carbon_pct_earth"])
        print(f"\nExtended Ladder '{key}' Contrast:")
        print(f"  {'accession':<18} {'class':<8} {'light':<6} {'hardware':<12} {'sep':>9} {'p':>10}")
        print("  " + "-" * 70)
        for _, r in sub.iterrows():
            print(f"  {r['accession']:<18} {r['enclosure_class']:<8} {r['light']:<6} {r['hardware']:<12} {r['separation']:>+9.4f} {r['p']:>10.2e}")

        lit = sub[sub["light"] == "light"].dropna(subset=["separation"])
        dark = sub[sub["light"] == "dark"].dropna(subset=["separation"])

        print("\nIllumination Test (Extended Ladder):")
        print(f"  Lit studies ({len(lit)}): min={lit['separation'].min():+.4f}, max={lit['separation'].max():+.4f}")
        print(f"  Dark studies ({len(dark)}): min={dark['separation'].min():+.4f}, max={dark['separation'].max():+.4f}")
        if lit["separation"].min() > dark["separation"].max():
            print("  -> PERFECT SEPARATION CONFIRMED: Every illuminated study sits above every dark study.")

    # 4. Write organ comparison table
    t13_path = os.path.join(ABRS_TABLES, "T13_abrs_organ_comparison.tsv")
    pd.DataFrame(organ_records).to_csv(t13_path, sep="\t", index=False)
    print(f"\n  wrote {t13_path} ({len(organ_records)} organ-set rows)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
