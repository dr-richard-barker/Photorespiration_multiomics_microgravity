#!/usr/bin/env python3
"""Count the metabolite-profiling studies in OSDR, and the organisms they cover.

The manuscript's claim that OSDR holds no plant metabolome rests on this census. It reads
the bulk study search that `probe_studies.py` caches (`data/cache/osdr_all_studies.json`,
fetched on first use) and writes every study whose assay measurement type mentions
metabolite profiling, with its organism, so the claim can be checked row by row.

A study counts as metabolite profiling when any value of "Study Assay Measurement Type"
contains "metabolite profiling" (case-insensitive). Writes
results/tables/T03_osdr_metabolomics_census.tsv.
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import TABLES, ensure  # noqa: E402
from probe_studies import SEARCH, get  # noqa: E402


def as_list(v) -> list[str]:
    if v is None:
        return []
    return [str(x) for x in v] if isinstance(v, list) else [str(v)]


def main() -> int:
    raw = get(SEARCH.format(term=""), "osdr_all_studies.json")
    hits = [h["_source"] for h in json.loads(raw)["hits"]["hits"]]
    metab = [s for s in hits
             if any("metabolite profiling" in m.lower()
                    for m in as_list(s.get("Study Assay Measurement Type")))]
    ensure(TABLES)
    path = os.path.join(TABLES, "T03_osdr_metabolomics_census.tsv")
    with open(path, "w") as fh:
        fh.write(f"# OSDR bulk study search (cached data/cache/osdr_all_studies.json): "
                 f"{len(hits)} studies, {len(metab)} with a metabolite-profiling assay.\n")
        fh.write("# Derived by scripts/osdr_census.py.\n")
        fh.write("accession\torganism\tmeasurement_types\ttitle\n")
        for s in sorted(metab, key=lambda s: int(str(s.get("Accession", "OSD-0"))[4:])):
            fh.write(f"{s.get('Accession')}\t{'; '.join(as_list(s.get('organism')))}\t"
                     f"{'; '.join(as_list(s.get('Study Assay Measurement Type')))}\t"
                     f"{str(s.get('Study Title', '')).strip()}\n")
    print(f"  wrote {path}  ({len(hits)} studies, {len(metab)} metabolite profiling)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
