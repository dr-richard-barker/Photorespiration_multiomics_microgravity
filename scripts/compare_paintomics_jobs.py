#!/usr/bin/env python3
"""Compare the first PaintOmics run with its resubmission, from the vendored job records.

Job 1 was first run (m1z16Qg3DK, 2026-09-08) on a metabolite layer predicted from an older
LunarLeaf-CFD export, and resubmitted (oK2h66740Q, 2026-10-29) on the current one. This
reports what moved between the two records, so the claim that nothing downstream did can be
checked rather than taken on trust. Job 2's first record was never saved past step 1, and the
server had deleted it by the time of the resubmission, so only the new job 2 is summarised.
"""

from __future__ import annotations

import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import PAINTOMICS_RAW  # noqa: E402

OLD = os.path.join(PAINTOMICS_RAW, "superseded", "job1_m1z16Qg3DK_full.json.gz")
NEW = os.path.join(PAINTOMICS_RAW, "job1_oK2h66740Q_full.json.gz")
JOB2 = os.path.join(PAINTOMICS_RAW, "job2_tnZU3q2Ie2_full.json.gz")
OMICS = ("Gene expression", "Proteomics", "Metabolomics")


def load(path: str) -> dict:
    with gzip.open(path, "rt") as fh:
        return json.load(fh)


def pathways(job: dict) -> dict:
    return {(p["name"], p.get("source")): p for p in job["pathwaysInfo"]}


def classes(job: dict) -> list:
    return sorted((lvl, c["name"], c["parent"], c["k"], c["n"], c["binomial"]["bh"])
                  for lvl, cs in job["classActivity"]["levels"].items() for c in cs)


def main() -> int:
    old, new = load(OLD), load(NEW)
    po, pn = pathways(old), pathways(new)
    print(f"job 1: {len(po)} pathways in {old['jobID']}, {len(pn)} in {new['jobID']}; "
          f"same set: {set(po) == set(pn)}")
    moved = {"combined": 0, **{o: 0 for o in OMICS}}
    for k in po.keys() & pn.keys():
        a, b = po[k], pn[k]
        if a["combinedSignificancePvalues"]["Fisher"] != b["combinedSignificancePvalues"]["Fisher"]:
            moved["combined"] += 1
        for o in OMICS:
            if (a.get("globalOmicPvalues") or {}).get(o) != (b.get("globalOmicPvalues") or {}).get(o):
                moved[o] += 1
    for what, n in moved.items():
        print(f"  {what:<16} p-values that differ: {n}")
    print(f"  hub analysis identical: {old['hubAnalysisResult'] == new['hubAnalysisResult']}")
    print(f"  class-activity classes identical: {classes(old) == classes(new)}")

    j2 = load(JOB2)
    ps = j2["pathwaysInfo"]
    by_db = {db: sum(1 for p in ps if p.get("source") == db) for db in ("KEGG", "MapMan")}
    sig = sum(1 for p in ps if p["combinedSignificancePvalues"]["Fisher"][0] < 0.05)
    print(f"job 2 ({j2['jobID']}): {by_db['KEGG']} KEGG and {by_db['MapMan']} MapMan pathways, "
          f"{sig} at combined p < 0.05; hub entries: {len(j2.get('hubAnalysisResult') or {})}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
