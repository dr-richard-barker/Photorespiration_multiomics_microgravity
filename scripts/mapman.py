#!/usr/bin/env python3
"""MapMan diagram membership for the 13 MapMan pathways PaintOmics called significant.

PaintOmics reports pathways from two databases. KEGG membership comes from the KEGG REST
API (`genesets.py`); MapMan membership has no REST API, so it is reconstructed here from
the same two inputs PaintOmics' own installer uses, and by the same rules.

**A MapMan "pathway" is a diagram, not an ontology bin.** All 13 names in T08 are diagram
titles, which is why matching them against `ontology.obo` returns 0 of 13. Each diagram is
an XML layout listing the ontology bins drawn on it; genes reach the diagram through the
gene-to-bin mapping.

Resolution follows `processMapManPathwaysData` in PaintOmics' `common_build_database.py`:

  1. Every `<Identifier id="...">` in the diagram is one bin.
  2. Bins are hierarchical and inclusion is recursive: bin `20.1` also takes `20.1.*`.
     The prefix pattern escapes the dots, or `18.4` would also match a bin `18X4`.
  3. Bin codes are normalised before comparison. Diagrams sometimes zero-pad a segment
     ("14.01", "13.1.5.3.02") and the gene mapping never does ("14.1"), so compared
     verbatim the two spellings never meet and the diagram imports with no genes at all.
     PaintOmics carries an explicit warning about this; the same normalisation is applied
     here. Only fully numeric segments are touched.

Inputs are vendored under `data/mapman/` — see its PROVENANCE.md. They are vendored rather
than cached because both upstream paths are fragile: MapManStore's own static
`img/mapman_36/*.zip` links are already dead (404), and the diagrams are only reachable
through a Liferay portlet query. `--fetch` refreshes them.
"""

from __future__ import annotations

import gzip
import itertools
import os
import re
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import DATA  # noqa: E402

MAPMAN = os.path.join(DATA, "mapman")
DIAGRAMS = os.path.join(MAPMAN, "diagrams")
GENE_TO_BIN = os.path.join(MAPMAN, "gene-to-mapman_ath.tsv.gz")


def normalise_bin(code: str) -> str:
    """Strip leading zeros from each numeric segment of a MapMan bin code."""
    if not code:
        return code
    segments = code.split(".")
    for i, seg in enumerate(segments):
        if seg.isdigit():
            segments[i] = seg.lstrip("0") or "0"
    return ".".join(segments)


def bin_to_genes() -> dict[str, set[str]]:
    """MapMan ontology bin -> AGI loci, from GoMapMan's Araport11 mapping."""
    out: dict[str, set[str]] = defaultdict(set)
    with gzip.open(GENE_TO_BIN, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            gene = parts[0].strip().upper()
            for code in parts[2].split():
                out[normalise_bin(code.strip())].add(gene)
    return dict(out)


def diagram_bins(name: str) -> set[str]:
    root = ET.parse(os.path.join(DIAGRAMS, f"{name}.xml")).getroot()
    return {normalise_bin(ident.get("id"))
            for area in root for ident in area if ident.get("id")}


def diagram_genes(name: str, mapping: dict[str, set[str]] | None = None) -> set[str]:
    """AGI loci placed on one MapMan diagram, recursive over the bin hierarchy."""
    mapping = bin_to_genes() if mapping is None else mapping
    genes: set[str] = set()
    for code in diagram_bins(name):
        pattern = re.compile(r"{0}(\.|\Z)".format(re.escape(code)))
        genes |= set(itertools.chain.from_iterable(
            v for k, v in mapping.items() if pattern.search(k)))
    return genes


def available() -> set[str]:
    return {f[:-4] for f in os.listdir(DIAGRAMS) if f.endswith(".xml")}


def paintomics_mapman_sets(names, db) -> tuple[dict[str, set[str]], list[str]]:
    """Membership for the MapMan rows of a PaintOmics enrichment table.

    `names` and `db` are parallel sequences; only rows whose `db` is "M" are looked up,
    the mirror of `genesets.paintomics_pathway_ids` — a MapMan bin name can collide with a
    KEGG pathway name (MapMan's lowercase "photosynthesis" against KEGG's "Photosynthesis"),
    so neither side may guess which database a row came from.

    Returns (resolved, unresolved); unresolved is anything with no vendored diagram.
    """
    have = available()
    mapping = bin_to_genes()
    resolved, unresolved = {}, []
    for name, source in zip(names, db):
        if source != "M":
            continue
        if name in have:
            resolved[name] = diagram_genes(name, mapping)
        else:
            unresolved.append(name)
    return resolved, unresolved


# --- refreshing the vendored inputs ---------------------------------------------------
# GoMapMan publishes the mapping PaintOmics uses; the diagrams come from MapManStore's
# portlet endpoint. Resource ids are from PaintOmics' own `mapman_extra_diagrams.json`
# except the three that GoMapMan's 20-diagram tarball already carries.

GOMAPMAN = ("https://gomapman.nib.si/api/GetFile/protein_2018-05-25%7Cpaintomics%7C"
            "gene-to-mapman_ath.tsv.gz")
MAPMAN_STORE = ("https://www.plabipd.de/portal/mapman"
                "?p_p_id=MapManDataDownload_WAR_MapManDataDownloadportlet"
                "&p_p_lifecycle=2&p_p_state=normal&p_p_mode=view"
                "&p_p_cacheability=cacheLevelPage&p_p_col_id=column-1&p_p_col_count=1"
                "&_MapManDataDownload_WAR_MapManDataDownloadportlet_Show=Pathways"
                "&_MapManDataDownload_WAR_MapManDataDownloadportlet_RessourceId={rid}"
                "&_MapManDataDownload_WAR_MapManDataDownloadportlet_Download=PathwayAnnotation")
GOMAPMAN_TARBALL = ("https://gomapman.nib.si/api/GetFile/protein_2018-05-25%7Cpaintomics%7C"
                    "mapman_pathways.tar.gz")

STORE_IDS = {
    "Cellular response overview": "76",
    "Large enzyme families overview": "108",
    "Transport overview": "52",
    "Overview": "103",
    "AGPs": "290",
    "Prokaryotic Ribosome SSU 5S branch assembly": "257",
    "Sulphate Assimilation": "90",
    "Raffinose metabolism": "139",
    "photosynthesis": "95",
    "Regulation overview": "74",
}
FROM_TARBALL = ["Biotic Stress", "JA Synthesis", "receptor like kinases"]


def fetch() -> int:
    """Re-download every vendored input. Writes only after all 14 downloads succeed."""
    import shutil
    import tarfile
    import tempfile
    import time
    import urllib.request

    def get(url, dest):
        req = urllib.request.Request(url, headers={"User-Agent": "photoresp-multiomics/1.0"})
        with urllib.request.urlopen(req, timeout=300) as resp, open(dest, "wb") as fh:
            shutil.copyfileobj(resp, fh)

    staging = tempfile.mkdtemp(prefix="mapman_")
    try:
        print(f"  gene-to-bin mapping from GoMapMan")
        get(GOMAPMAN, os.path.join(staging, "gene-to-mapman_ath.tsv.gz"))

        print(f"  {len(FROM_TARBALL)} diagrams from GoMapMan's 20-diagram tarball")
        tarpath = os.path.join(staging, "mapman_pathways.tar.gz")
        get(GOMAPMAN_TARBALL, tarpath)
        with tarfile.open(tarpath, "r:gz") as tar:
            for name in FROM_TARBALL:
                member = tar.extractfile(f"xml/{name}.xml")
                if member is None:
                    raise RuntimeError(f"{name} is not in the GoMapMan tarball")
                with open(os.path.join(staging, f"{name}.xml"), "wb") as fh:
                    fh.write(member.read())
        os.remove(tarpath)

        print(f"  {len(STORE_IDS)} diagrams from MapManStore")
        for name, rid in STORE_IDS.items():
            time.sleep(1.5)  # a Liferay portlet, not a CDN
            dest = os.path.join(staging, f"{name}.xml")
            get(MAPMAN_STORE.format(rid=rid), dest)
            with open(dest, errors="replace") as fh:
                if "<Image" not in fh.read(4096):
                    raise RuntimeError(f"{name} (id {rid}) did not return MapMan XML")

        # All or nothing: a partial refresh would leave the pathway universe — which is an
        # enrichment denominator — silently different from the one the results were built on.
        os.makedirs(DIAGRAMS, exist_ok=True)
        for fname in os.listdir(staging):
            dest = DIAGRAMS if fname.endswith(".xml") else MAPMAN
            shutil.move(os.path.join(staging, fname), os.path.join(dest, fname))
        print(f"  refreshed {len(STORE_IDS) + len(FROM_TARBALL)} diagrams + the mapping")
        return 0
    finally:
        shutil.rmtree(staging, ignore_errors=True)


if __name__ == "__main__":
    if "--fetch" in sys.argv:
        raise SystemExit(fetch())
    m = bin_to_genes()
    print(f"{len(m)} MapMan bins, "
          f"{len(set(itertools.chain.from_iterable(m.values())))} AGI loci\n")
    for name in sorted(available()):
        print(f"  {len(diagram_genes(name, m)):>6} loci  {name}")
