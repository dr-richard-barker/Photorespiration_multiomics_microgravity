#!/usr/bin/env python3
"""Download processed OSD-7 (TAGES/ABRS) microarray data from NASA OSDR.

Downloads:
  1. RMA-normalized probeset expression matrix (GeneLab pipeline output)
  2. ISA-Tab metadata archive for sample-factor mapping

These are the GeneLab-processed files, not raw CEL files. The normalised
expression matrix is the starting point for differential expression analysis
in 01_dge_microarray.py.
"""

from __future__ import annotations

import os
import sys
import urllib.request
import zipfile

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from paths import ensure  # noqa: E402

ABRS_CACHE = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "data", "abrs", "cache"
)
ABRS_CACHE = os.path.normpath(ABRS_CACHE)

OSDR_DL = ("https://osdr.nasa.gov/geode-py/ws/studies/{osd}/download"
           "?source=datamanager&file={f}")
UA = {"User-Agent": "photoresp-multiomics-abrs/1.0"}

FILES = {
    # GeneLab RMA-normalized expression — the one we need for DGE
    "GLDS-7_array_normalized_expression_probeset_GLmicroarray.csv": {
        "osd": "OSD-7",
        "desc": "RMA-normalized probeset expression (20.7 MB)",
    },
    # ISA-Tab metadata — maps samples to factors (Flight/GC, organ)
    "OSD-7_metadata_OSD-7-ISA.zip": {
        "osd": "OSD-7",
        "desc": "ISA-Tab metadata archive (74 KB)",
    },
}


def download(filename: str, info: dict) -> str:
    """Download a single file from OSDR, returning the local path."""
    path = os.path.join(ABRS_CACHE, filename)
    if os.path.exists(path):
        sz = os.path.getsize(path)
        print(f"  cached  {filename} ({sz:,} bytes)")
        return path
    url = OSDR_DL.format(osd=info["osd"], f=filename)
    print(f"  fetch   {filename}")
    print(f"          {info['desc']}")
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=600) as resp:
        body = resp.read()
    ensure(ABRS_CACHE)
    with open(path, "wb") as fh:
        fh.write(body)
    print(f"          -> {len(body):,} bytes")
    return path


def extract_isa(zip_path: str) -> None:
    """Extract ISA-Tab files from the metadata zip."""
    dest = os.path.join(ABRS_CACHE, "isa")
    if os.path.isdir(dest):
        print(f"  ISA-Tab already extracted to {dest}")
        return
    os.makedirs(dest, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(dest)
    print(f"  extracted ISA-Tab to {dest}")
    for name in sorted(os.listdir(dest)):
        print(f"    {name}")


def main() -> int:
    print("OSD-7 (TAGES/ABRS) — downloading processed microarray data\n")
    paths = {}
    for fname, info in FILES.items():
        paths[fname] = download(fname, info)

    # Extract ISA-Tab for sample metadata
    isa_zip = paths.get("OSD-7_metadata_OSD-7-ISA.zip")
    if isa_zip:
        extract_isa(isa_zip)

    # Quick sanity check on the expression file
    expr_file = paths.get(
        "GLDS-7_array_normalized_expression_probeset_GLmicroarray.csv"
    )
    if expr_file:
        with open(expr_file) as fh:
            header = fh.readline().strip()
            n_cols = len(header.split(","))
            n_lines = sum(1 for _ in fh)
        print(f"\n  expression matrix: {n_lines:,} probesets × {n_cols} columns")

    print("\ndone — ready for 01_dge_microarray.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
