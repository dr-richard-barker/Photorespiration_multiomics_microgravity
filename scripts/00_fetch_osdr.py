#!/usr/bin/env python3
"""Fetch the OSD-522 (BRIC-LED-001) processed files we need from NASA OSDR.

OSD-522 — "Integrative Transcriptomics and Proteomics Profiling of Arabidopsis thaliana
Elucidates Novel Mechanisms Underlying Spaceflight Adaptation" (SpaceX-13 / SpaceX-14,
BRIC-LED hardware, Arabidopsis seedling shoots, 6 Ground Control vs 6 Space Flight).

Only the four small processed files are pulled — not the 274-file raw archive:

  * RSEM unnormalised counts  (AGI-indexed gene x sample matrix)
  * the bulk RNA-seq runsheet (sample -> Organism Part / Spaceflight factor values)
  * shoot soluble proteome    (already differential: S/G log2 ratio + adjusted p)
  * shoot membrane proteome   (ditto)

Files are cached under osdr/cache/ and re-used on later runs; pass --force to re-download.
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.parse
import urllib.request

OSD = "OSD-522"
BASE = "https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download?source=datamanager&file={f}"

FILES = {
    "counts": "GLDS-522_rna_seq_RSEM_Unnormalized_Counts_GLbulkRNAseq.csv",
    "runsheet": "GLDS-522_rna_seq_bulkRNASeq_v2_runsheet.csv",
    "prot_shoot_sol": "GLDS-522_proteomics_GO_Shoot_SOL_Report_20220223_Proteins.csv",
    "prot_shoot_mem": "GLDS-522_proteomics_GO_Shoot_MEM_Report_20220223_Proteins.csv",
}

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import CACHE  # noqa: E402


def url_for(filename: str) -> str:
    return BASE.format(osd=OSD, f=urllib.parse.quote(filename))


def fetch(key: str, filename: str, force: bool = False) -> str:
    """Download one OSDR file into the cache; return its local path."""
    dest = os.path.join(CACHE, filename)
    if os.path.exists(dest) and not force:
        print(f"  cached  {key:16s} {filename} ({os.path.getsize(dest):,} bytes)")
        return dest

    os.makedirs(CACHE, exist_ok=True)
    url = url_for(filename)
    print(f"  fetch   {key:16s} {filename}")
    req = urllib.request.Request(url, headers={"User-Agent": "photoresp-multiomics/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        body = resp.read()
    if len(body) < 1024:
        raise RuntimeError(f"{filename}: response only {len(body)} bytes — likely an error page")
    with open(dest, "wb") as fh:
        fh.write(body)
    print(f"          -> {len(body):,} bytes")
    return dest


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true", help="re-download even if cached")
    args = ap.parse_args()

    print(f"OSD-522 (BRIC-LED-001) — fetching processed files into {CACHE}")
    for key, filename in FILES.items():
        try:
            fetch(key, filename, force=args.force)
        except Exception as exc:  # noqa: BLE001 — report and keep going
            print(f"  FAILED  {key}: {exc}", file=sys.stderr)
            return 1
    print("all files present")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
