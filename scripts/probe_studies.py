#!/usr/bin/env python3
"""Probe candidate OSDR studies and report what is actually in them.

Run once, by hand, to ground `data/study_registry.tsv` in real metadata rather than in
assumptions about what a study contains. It prints, per accession: title, hardware project,
organism, the organism-part and other factor values, the processed differential-expression
file if there is one, and every flight-vs-ground contrast column that file offers.

The registry is then written from this output with a human decision per study — this script
never writes the registry itself, because the include/exclude calls are judgements (root vs
shoot, genotype, platform) that should be recorded as such.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import os
import re
import sys
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import CACHE, ensure  # noqa: E402

SEARCH = "https://osdr.nasa.gov/osdr/data/search?term={term}&type=cgene&size=1000"
FILES = "https://osdr.nasa.gov/osdr/data/osd/files/{num}"
DOWNLOAD = ("https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download"
            "?source=datamanager&file={f}")

CANDIDATES = ["OSD-522", "OSD-38", "OSD-321", "OSD-678", "OSD-120",
              "OSD-218", "OSD-427", "OSD-193", "OSD-281", "OSD-16", "OSD-7"]

FLIGHT_WORDS = re.compile(r"space\s*flight|\bflt\b|\bflight\b", re.I)
GROUND_WORDS = re.compile(r"ground\s*control|\bgc\b|\bgnd\b", re.I)


def get(url: str, cache_name: str) -> bytes:
    ensure(CACHE)
    path = os.path.join(CACHE, cache_name)
    if os.path.exists(path):
        return open(path, "rb").read()
    req = urllib.request.Request(url, headers={"User-Agent": "photoresp-multiomics/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        body = resp.read()
    with open(path, "wb") as fh:
        fh.write(body)
    return body


def metadata() -> dict:
    """Accession -> study metadata, from one bulk search call."""
    raw = get(SEARCH.format(term=""), "osdr_all_studies.json")
    hits = json.loads(raw)["hits"]["hits"]
    return {h["_source"].get("Accession"): h["_source"] for h in hits}


def processed_files(acc: str) -> list[str]:
    num = acc.split("-")[1]
    raw = get(FILES.format(num=num), f"osdr_files_{acc}.json")
    out = []
    for _k, v in json.loads(raw).get("studies", {}).items():
        for f in v.get("study_files", []):
            out.append(f["file_name"])
    return out


def dge_file(files: list[str]) -> str | None:
    cands = [f for f in files if "differential_expression" in f and f.endswith(".csv")]
    if not cands:
        return None
    # prefer rRNA-removed when both exist — GeneLab's own preferred product
    rrna = [c for c in cands if "rRNArm" in c]
    return (rrna or cands)[0]


def contrasts(acc: str, filename: str) -> tuple[list[str], str]:
    """Every Log2fc_ column, plus the first column's name (the identifier column)."""
    url = DOWNLOAD.format(osd=acc, f=urllib.parse.quote(filename))
    raw = get(url, f"dge_{acc}.csv")
    text = raw.decode("utf-8", errors="replace")
    header = next(csv.reader(io.StringIO(text)))
    return [c for c in header if c.startswith("Log2fc_")], header[0]


def _terms(side: str) -> list[str]:
    return [t.strip() for t in side.split("&")]


def flight_vs_ground(cols: list[str]) -> list[str]:
    """Flight-vs-ground contrasts that are MATCHED on every other factor.

    A study like OSD-678 offers 36 contrasts with flight on the left and ground on the
    right, but most of them also swap genotype, ecotype or light treatment — comparing
    flown Col-0 against ground-control phyD confounds the very thing we are testing. Only
    contrasts whose non-flight terms are identical on both sides are returned.
    """
    keep = []
    for c in cols:
        m = re.match(r"Log2fc_\((.*)\)v\((.*)\)$", c)
        if not m:
            continue
        a, b = _terms(m.group(1)), _terms(m.group(2))
        if len(a) != len(b):
            continue
        a_rest = [t for t in a if not (FLIGHT_WORDS.search(t) or GROUND_WORDS.search(t))]
        b_rest = [t for t in b if not (FLIGHT_WORDS.search(t) or GROUND_WORDS.search(t))]
        a_flt = any(FLIGHT_WORDS.search(t) for t in a)
        b_gnd = any(GROUND_WORDS.search(t) for t in b)
        if a_flt and b_gnd and a_rest == b_rest:
            keep.append(c)
    return keep


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", nargs="*", default=CANDIDATES)
    args = ap.parse_args()

    meta = metadata()
    for acc in args.accessions:
        s = meta.get(acc, {})
        print("=" * 96)
        print(f"{acc}  {str(s.get('Study Title', ''))[:80]}")
        print(f"  project   : {s.get('Project Identifier', '')}")
        print(f"  organism  : {s.get('organism', '')}")
        print(f"  material  : {str(s.get('Material Type', ''))[:70]}")
        print(f"  factors   : {str(s.get('Study Factor Name', ''))[:70]}")
        print(f"  assay     : {str(s.get('Study Assay Measurement Type', ''))[:60]}"
              f" / {str(s.get('Study Assay Technology Type', ''))[:40]}")

        try:
            files = processed_files(acc)
        except Exception as exc:  # noqa: BLE001
            print(f"  FILES FAILED: {exc}")
            continue

        proteomes = [f for f in files if "proteom" in f.lower() and f.endswith(".csv")]
        dge = dge_file(files)
        print(f"  files     : {len(files)} total, {len(proteomes)} proteomics csv")
        if not dge:
            print("  DGE       : NONE — cannot use as a transcriptome arm")
            continue
        print(f"  DGE       : {dge}")

        try:
            cols, idcol = contrasts(acc, dge)
        except Exception as exc:  # noqa: BLE001
            print(f"  CONTRASTS FAILED: {exc}")
            continue

        fvg = flight_vs_ground(cols)
        print(f"  id column : {idcol}")
        print(f"  contrasts : {len(cols)} total, {len(fvg)} flight-vs-ground")
        for c in fvg[:10]:
            print(f"      {c}")
        if len(fvg) > 10:
            print(f"      … and {len(fvg) - 10} more")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
