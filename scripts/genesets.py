#!/usr/bin/env python3
"""KEGG-derived Arabidopsis gene sets, shared by the single-study and cross-study analyses.

Both `04_falsification_check.py` (OSD-522 alone) and `07_hardware_ladder.py` (the enclosure
comparison) must score against *identical* sets, or the ladder is not comparable with the
result it is meant to extend. So the sets live here, once.

Membership is fetched live from `rest.kegg.jp` and cached. It is never typed from memory:
an earlier hand-written list in this project put PGLP1 at the wrong locus, and Arabidopsis
reuses the symbols CAT2 and SEN1 for two unrelated genes each. Pathway sets come from KEGG's
own pathway membership; curated sets are selected by regular expression over KEGG's
description text, because several canonical markers (PDC1, the hypoxia-responsive family)
carry no gene symbol in KEGG at all and are only findable that way.
"""

from __future__ import annotations

import os
import re
import sys
import urllib.request

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import CACHE, ensure  # noqa: E402

KEGG = "https://rest.kegg.jp/{endpoint}"

PATHWAY_SETS = {
    "ath00630": "Glyoxylate & dicarboxylate metabolism (photorespiration)",
    "ath00710": "Carbon fixation in photosynthetic organisms",
    "ath00500": "Starch and sucrose metabolism",
    "ath00010": "Glycolysis / gluconeogenesis",
}

MARKER_SETS = {
    "photorespiration_core": r"phosphoglycolate phosphatase|glycolate oxidase|"
                             r"glycine decarboxylase|transhydroxymethyltransferase|"
                             r"hydroxypyruvate reductase|glycerate kinase|"
                             r"glutamate:glyoxylate aminotransferase",
    # DIN = DARK INDUCED, the canonical Arabidopsis sugar/carbon-starvation marker family.
    "carbon_starvation_DIN": r"\bDIN\d+;",
    "fermentation": r"pyruvate decarboxylase|alcohol dehydrogenase 1;|lactate dehydrogenase",
    "hypoxia_responsive": r"[Hh]ypoxia-responsive",
    "photosynthesis_apparatus": r"photosystem I{1,2} subunit|light harvesting complex|"
                                r"chlorophyll A/B binding",
    "rubisco": r"ribulose bisphosphate carboxylase",
}

# What the CFD/FvCB model predicts each set should do in a flight-vs-ground contrast inside
# a SEALED, ILLUMINATED enclosure. "none" is a prediction the data can break just as much as
# a directional one. None means "not predicted; reported for context only".
PREDICTION = {
    "ath00630": "none",
    "photorespiration_core": "none",
    "carbon_starvation_DIN": "up",
    "fermentation": "none",
    "hypoxia_responsive": "none",
    "ath00710": "down",
    "photosynthesis_apparatus": "down",
    "rubisco": "down",
    "ath00500": "down",
    "ath00010": None,
}

# Set-vs-set contrasts. Every metabolic set shifts down together in these experiments (a
# global depression of metabolic transcripts in flight), which makes each set's comparison
# against "all other genes" partly a test of whether a gene is metabolic at all. These
# pairwise contrasts ask the sharper question the model actually makes a claim about.
CONTRASTS = [
    ("photorespiration vs carbon fixation", "ath00630", "ath00710",
     "model: Vo flat while A falls, so photorespiration should sit ABOVE carbon fixation"),
    ("photorespiration_core vs rubisco", "photorespiration_core", "rubisco",
     "same claim, on curated enzyme sets rather than whole KEGG maps"),
    ("starvation vs photosynthesis", "carbon_starvation_DIN", "photosynthesis_apparatus",
     "model: carbon deficit, so starvation markers should sit ABOVE the apparatus"),
]


def kegg(endpoint: str, cache_name: str) -> str:
    path = os.path.join(CACHE, cache_name)
    if os.path.exists(path):
        return open(path).read()
    ensure(CACHE)
    req = urllib.request.Request(KEGG.format(endpoint=endpoint),
                                 headers={"User-Agent": "photoresp-multiomics/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        text = resp.read().decode()
    with open(path, "w") as fh:
        fh.write(text)
    return text


def annotation() -> pd.DataFrame:
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


def uniprot_to_agi() -> dict[str, str]:
    """UniProt accession -> AGI locus, from KEGG's own conversion table."""
    conv = {}
    for line in kegg("conv/uniprot/ath", "kegg_ath_uniprot.tsv").splitlines():
        if "\t" not in line:
            continue
        agi, up = line.split("\t")
        conv[up.replace("up:", "")] = agi.replace("ath:", "")
    return conv


def build(verbose: bool = False) -> dict[str, set[str]]:
    """Every gene set, keyed by name."""
    sets: dict[str, set[str]] = {}
    for pid in PATHWAY_SETS:
        sets[pid] = pathway_genes(pid)
        if verbose:
            print(f"    {pid}: {len(sets[pid])} genes — {PATHWAY_SETS[pid]}")

    ann = annotation()
    for name, pattern in MARKER_SETS.items():
        hits = ann[ann["description"].str.contains(pattern, case=False, regex=True, na=False)]
        sets[name] = set(hits.index)
        if verbose:
            print(f"    {name}: {len(hits)} genes")
    return sets




# --- PaintOmics significant pathways -> KEGG ath ids --------------------------------
# The enrichment table (T08) names pathways but does not carry their KEGG ids, and the
# Sankey needs real gene membership. KEGG's own `list/pathway/ath` supplies the mapping.
#
# Two traps, both hit while writing this:
#
#   1. KEGG appends " - Arabidopsis thaliana (thale cress)" to every name, and splitting
#      on the FIRST " - " truncates "Photosynthesis - antenna proteins" to
#      "Photosynthesis". Strip only the TRAILING organism suffix. With that fixed, all 18
#      KEGG-side significant pathways resolve.
#   2. MapMan bin names can collide with KEGG pathway names — MapMan's lowercase
#      "photosynthesis" bin matched KEGG's "Photosynthesis" (ath00195) and would have been
#      silently double-counted alongside the genuine KEGG row. Callers must pass KEGG-side
#      names ONLY, which is why `db` is a required argument rather than a convention.
#
# MapMan pathways have no KEGG id and are resolved separately, by `scripts/mapman.py`,
# from vendored diagram layouts. They come back here as unresolved on purpose: this
# function answers only for KEGG.

ORGANISM_SUFFIX = re.compile(r"\s*-\s*Arabidopsis thaliana.*$", re.I)


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def kegg_pathway_name_index() -> dict[str, str]:
    """Normalised KEGG pathway name -> ath##### id."""
    out = {}
    for line in kegg("list/pathway/ath", "kegg_ath_pathway_names.tsv").splitlines():
        if "\t" not in line:
            continue
        pid, name = line.split("\t")
        out[_norm(ORGANISM_SUFFIX.sub("", name).strip())] = pid.replace("path:", "")
    return out


def paintomics_pathway_ids(names, db) -> tuple[dict[str, str], list[str]]:
    """Map PaintOmics pathway names to KEGG ath ids.

    `names` and `db` are parallel sequences; only rows whose `db` is "K" are looked up,
    because a MapMan bin name can collide with a KEGG pathway name (see above).

    Returns (resolved, unresolved). Unresolved entries are the MapMan rows, which
    `mapman.paintomics_mapman_sets` answers for, plus anything KEGG does not know for this
    organism — report those, never drop them quietly.
    """
    index = kegg_pathway_name_index()
    resolved, unresolved = {}, []
    for name, source in zip(names, db):
        pid = index.get(_norm(name)) if source == "K" else None
        if pid:
            resolved[name] = pid
        else:
            unresolved.append(name)
    return resolved, unresolved


if __name__ == "__main__":
    print("KEGG-derived gene sets")
    s = build(verbose=True)
    print(f"\n{len(s)} sets, {sum(len(v) for v in s.values())} memberships")
