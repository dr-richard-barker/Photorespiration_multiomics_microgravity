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
import re
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
CACHE = os.path.join(HERE, "cache")
UPLOAD = os.path.join(REPO, "paintomics", "upload")

KEGG = "https://rest.kegg.jp/{endpoint}"

# KEGG pathway sets — membership is authoritative, no symbol guessing involved.
PATHWAY_SETS = {
    "ath00630": "Glyoxylate & dicarboxylate metabolism (photorespiration)",
    "ath00710": "Carbon fixation in photosynthetic organisms",
    "ath00500": "Starch and sucrose metabolism",
    "ath00010": "Glycolysis / gluconeogenesis",
}

# Curated marker sets, matched on KEGG's own description text rather than on gene symbols.
# Symbols are unreliable here: Arabidopsis reuses CAT2 and SEN1 for unrelated genes, and
# several canonical markers (PDC1, the hypoxia-responsive family) carry no symbol at all in
# KEGG and are only findable by description. Every member that a regex selects is printed,
# so the sets can be audited rather than trusted.
MARKER_SETS = {
    # The C2 cycle proper. Enzyme families, so a description regex captures all isoforms.
    "photorespiration_core": r"phosphoglycolate phosphatase|glycolate oxidase|"
                             r"glycine decarboxylase|transhydroxymethyltransferase|"
                             r"hydroxypyruvate reductase|glycerate kinase|"
                             r"glutamate:glyoxylate aminotransferase",
    # DIN = DARK INDUCED, the canonical Arabidopsis sugar/carbon-starvation marker family.
    # A coherent, literature-defined family, unlike a hand-mixed starvation list.
    "carbon_starvation_DIN": r"\bDIN\d+;",
    # Fermentative entry point. PDC1 (AT4G33070) has no symbol in KEGG, only this description.
    "fermentation": r"pyruvate decarboxylase|alcohol dehydrogenase 1;|lactate dehydrogenase",
    "hypoxia_responsive": r"[Hh]ypoxia-responsive",
    # Light harvesting and photosystem subunits — the photosynthetic apparatus.
    "photosynthesis_apparatus": r"photosystem I{1,2} subunit|light harvesting complex|"
                                r"chlorophyll A/B binding",
    "rubisco": r"ribulose bisphosphate carboxylase",
}

# Direct set-vs-set contrasts. Every metabolic pathway set turns out to shift down together
# in this experiment (a global depression of metabolic transcripts in flight), which makes
# each set's comparison against "all other genes" partly a test of whether the gene is
# metabolic at all. These pairwise contrasts ask the sharper question the model actually
# makes a claim about: does photorespiration move differently from carbon fixation?
CONTRASTS = [
    ("photorespiration vs carbon fixation", "ath00630", "ath00710",
     "model: Vo flat while A falls, so photorespiration should sit ABOVE carbon fixation"),
    ("photorespiration_core vs rubisco", "photorespiration_core", "rubisco",
     "same claim, on curated enzyme sets rather than whole KEGG maps"),
    ("starvation vs photosynthesis", "carbon_starvation_DIN", "photosynthesis_apparatus",
     "model: carbon deficit, so starvation markers should sit ABOVE the apparatus"),
]

# What the model predicts each set should do in FLT vs GC. "none" = explicitly no change,
# which is a prediction the data can break just as much as a directional one.
PREDICTION = {
    "ath00630": "none",
    "photorespiration_core": "none",
    "carbon_starvation_DIN": "up",
    "fermentation": "none",          # lit canister: photosynthesis releases O2, no hypoxia
    "hypoxia_responsive": "none",    # same claim, independent gene set
    "ath00710": "down",
    "photosynthesis_apparatus": "down",
    "rubisco": "down",
    "ath00500": "down",
    "ath00010": None,      # not predicted; reported for context only
}


def kegg(endpoint: str, cache_name: str) -> str:
    path = os.path.join(CACHE, cache_name)
    if os.path.exists(path):
        return open(path).read()
    os.makedirs(CACHE, exist_ok=True)
    url = KEGG.format(endpoint=endpoint)
    req = urllib.request.Request(url, headers={"User-Agent": "photoresp-multiomics/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode()
    with open(path, "w") as fh:
        fh.write(text)
    return text


def load_annotation() -> pd.DataFrame:
    """AGI -> (symbol, description) from KEGG's own gene list."""
    rows = []
    for line in kegg("list/ath", "kegg_ath_genes.tsv").splitlines():
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        agi = parts[0].replace("ath:", "")
        desc = parts[3]
        symbol = desc.split(";")[0].strip() if ";" in desc else ""
        rows.append((agi, symbol, desc))
    return pd.DataFrame(rows, columns=["agi", "symbol", "description"]).set_index("agi")


def pathway_genes(pid: str) -> set[str]:
    text = kegg(f"link/ath/path:{pid}", f"kegg_{pid}_genes.tsv")
    return {ln.split("\t")[1].replace("ath:", "") for ln in text.splitlines() if "\t" in ln}


def resolve_markers(ann: pd.DataFrame, verbose: bool) -> dict[str, set[str]]:
    """Select each marker set by regex over KEGG's description, listing what it caught."""
    sets: dict[str, set[str]] = {}
    for set_name, pattern in MARKER_SETS.items():
        hits = ann[ann["description"].str.contains(pattern, case=False, regex=True, na=False)]
        sets[set_name] = set(hits.index)
        print(f"    {set_name}: {len(hits)} genes")
        if verbose:
            for agi, row in hits.iterrows():
                print(f"        {agi}  {row['description'][:72]}")
    return sets


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
    ann = load_annotation()
    print(f"    KEGG ath annotation: {len(ann):,} genes")

    gene_sets: dict[str, set[str]] = {}
    for pid, desc in PATHWAY_SETS.items():
        gene_sets[pid] = pathway_genes(pid)
        print(f"    {pid}: {len(gene_sets[pid])} genes — {desc}")
    gene_sets.update(resolve_markers(ann, verbose=not args.quiet_resolution))

    # --- transcriptome -------------------------------------------------------------
    genes = pd.read_csv(os.path.join(UPLOAD, "gene_expression_values.tab"),
                        sep="\t", index_col=0).iloc[:, 0]
    tx = report("TRANSCRIPTOME — OSD-522 shoots, Space Flight vs Ground Control", genes,
                gene_sets)

    # --- proteome: map UniProt accessions to AGI via KEGG's own conversion ----------
    conv = {}
    for line in kegg("conv/uniprot/ath", "kegg_ath_uniprot.tsv").splitlines():
        if "\t" not in line:
            continue
        agi, up = line.split("\t")
        conv[up.replace("up:", "")] = agi.replace("ath:", "")

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

    out = os.path.join(REPO, "results", "falsification_check.tsv")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    combined = pd.concat([tx.assign(layer="transcriptome"), pr.assign(layer="proteome")])
    combined.to_csv(out, sep="\t", index=False)
    print(f"\n  wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
