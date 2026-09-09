#!/usr/bin/env python3
"""Test the model's prediction against the real OSD-522 omics — before any PaintOmics run.

metabolome/predict_metabolome.py makes a specific, falsifiable claim about Space Flight
vs Ground Control inside a sealed BRIC-LED canister:

  * photorespiratory flux Vo is essentially UNCHANGED   (log2FC ~ +0.003)
  * net assimilation A is DOWN about 7-8 %              (log2FC ~ -0.114)
  * therefore carbon-starvation responses should be UP, and photorespiratory
    gene/protein expression should NOT be induced

This script scores the measured transcriptome and proteome against those directions.
It can refute the model, and a refutation is a result — nothing here is tuned to agree.

Gene sets come from KEGG's own pathway membership and gene symbols, fetched live from
rest.kegg.jp and cached. They are not typed from memory: an earlier hand-written list in
this session had PGLP1 at the wrong locus, and symbols such as CAT2 and SEN1 are ambiguous
in Arabidopsis (two different genes share each). Curated markers are therefore resolved by
symbol AND disambiguated on KEGG's description text, with every resolution printed.

Statistics: each set's log2 fold changes are compared with all remaining measured features
by a two-sided Mann-Whitney U test. That asks whether the set is shifted relative to the
rest of the experiment, which is the relevant question for a directional prediction.
"""

from __future__ import annotations

import argparse
import os
import sys
import re
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import genesets  # noqa: E402
from paths import TABLES, UPLOAD, ensure  # noqa: E402

PATHWAY_SETS = genesets.PATHWAY_SETS
MARKER_SETS = genesets.MARKER_SETS
PREDICTION = genesets.PREDICTION
CONTRASTS = genesets.CONTRASTS






def score(values: pd.Series, members: set[str], label: str, predicted: str | None) -> dict:
    """Mann-Whitney U of the set against everything else measured."""
    inside = values[values.index.isin(members)].dropna()
    outside = values[~values.index.isin(members)].dropna()
    if len(inside) < 3:
        return {"set": label, "n": len(inside), "median": np.nan, "p": np.nan,
                "observed": "too few", "predicted": predicted, "verdict": "n/a"}

    u, p = stats.mannwhitneyu(inside, outside, alternative="two-sided")
    med = float(inside.median())
    shift = med - float(outside.median())

    if p >= 0.05:
        observed = "no change"
    else:
        observed = "up" if shift > 0 else "down"

    if predicted is None:
        verdict = "—"
    elif predicted == "none":
        verdict = "MATCH" if observed == "no change" else "MISMATCH"
    else:
        verdict = "MATCH" if observed == predicted else "MISMATCH"

    return {"set": label, "n": len(inside), "median": med, "shift": shift,
            "p": float(p), "observed": observed, "predicted": predicted, "verdict": verdict}


def report(title: str, values: pd.Series, gene_sets: dict[str, set[str]]) -> pd.DataFrame:
    print(f"\n{title}")
    print(f"  {len(values):,} features; background median log2FC "
          f"{values.median():+.4f}")
    rows = [score(values, m, name, PREDICTION.get(name)) for name, m in gene_sets.items()]
    df = pd.DataFrame(rows)
    print(f"\n  {'set':<26} {'n':>5} {'median':>8} {'shift':>8} {'p':>10} "
          f"{'observed':>10} {'predicted':>10} {'verdict':>9}")
    print("  " + "-" * 94)
    for r in rows:
        p = "n/a" if not np.isfinite(r.get("p", np.nan)) else f"{r['p']:.2e}"
        med = "n/a" if not np.isfinite(r["median"]) else f"{r['median']:+.4f}"
        sh = "n/a" if not np.isfinite(r.get("shift", np.nan)) else f"{r['shift']:+.4f}"
        print(f"  {r['set']:<26} {r['n']:>5} {med:>8} {sh:>8} {p:>10} "
              f"{str(r['observed']):>10} {str(r['predicted']):>10} {r['verdict']:>9}")

    contrast_rows = run_contrasts(values, gene_sets)
    return pd.concat([df, pd.DataFrame(contrast_rows)], ignore_index=True)


def run_contrasts(values: pd.Series, gene_sets: dict[str, set[str]]) -> list[dict]:
    """Set-vs-set tests, which are immune to the global metabolic shift."""
    print(f"\n  set-vs-set contrasts (A minus B; the model says A should sit above B)")
    print(f"  {'contrast':<38} {'nA':>4} {'nB':>4} {'medA':>8} {'medB':>8} "
          f"{'A-B':>8} {'p':>10} {'verdict':>9}")
    print("  " + "-" * 94)
    out = []
    for label, a_name, b_name, _rationale in CONTRASTS:
        a = values[values.index.isin(gene_sets.get(a_name, set()))].dropna()
        b = values[values.index.isin(gene_sets.get(b_name, set()))].dropna()
        if len(a) < 3 or len(b) < 3:
            print(f"  {label:<38} {len(a):>4} {len(b):>4} "
                  f"{'':>8} {'':>8} {'':>8} {'':>10} {'too few':>9}")
            continue
        u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        diff = float(a.median() - b.median())
        # The model claims A sits above B, so support needs diff > 0 AND a real difference.
        verdict = "SUPPORTS" if (diff > 0 and p < 0.05) else (
            "no signal" if p >= 0.05 else "OPPOSES")
        print(f"  {label:<38} {len(a):>4} {len(b):>4} {a.median():+8.4f} "
              f"{b.median():+8.4f} {diff:+8.4f} {p:10.2e} {verdict:>9}")
        out.append({"set": label, "n": len(a), "median": float(a.median()),
                    "shift": diff, "p": float(p), "observed": f"{diff:+.4f}",
                    "predicted": "A above B", "verdict": verdict})
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--quiet-resolution", action="store_true",
                    help="do not print per-symbol disambiguation")
    args = ap.parse_args()

    print("Falsification check — model prediction vs measured OSD-522 omics")
    print("\n  building gene sets from KEGG (cached in osdr/cache/)")

    gene_sets = genesets.build(verbose=not args.quiet_resolution)

    # --- transcriptome -------------------------------------------------------------
    genes = pd.read_csv(os.path.join(UPLOAD, "gene_expression_values.tab"),
                        sep="\t", index_col=0).iloc[:, 0]
    tx = report("TRANSCRIPTOME — OSD-522 shoots, Space Flight vs Ground Control", genes,
                gene_sets)

    # --- proteome: map UniProt accessions to AGI via KEGG's own conversion ----------
    conv = genesets.uniprot_to_agi()

    prot = pd.read_csv(os.path.join(UPLOAD, "proteomics_values.tab"), sep="\t", index_col=0)
    prot_agi = prot.copy()
    prot_agi["agi"] = [conv.get(a) for a in prot_agi.index]
    mapped = prot_agi.dropna(subset=["agi"])
    print(f"\n  proteome: {len(mapped):,}/{len(prot):,} accessions mapped to AGI via KEGG")
    # Several accessions can map to one locus; average their log2 ratios so each gene
    # contributes once, matching how the transcriptome is keyed.
    value_col = prot.columns[0]
    per_gene = mapped.groupby("agi")[value_col].mean()
    pr = report("PROTEOME — OSD-522 shoots, Space Flight vs Ground Control", per_gene,
                gene_sets)

    ensure(TABLES)
    out = os.path.join(TABLES, "T07_falsification_osd522.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    combined = pd.concat([tx.assign(layer="transcriptome"), pr.assign(layer="proteome")])
    combined.to_csv(out, sep="\t", index=False)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
