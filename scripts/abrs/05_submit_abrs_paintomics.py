#!/usr/bin/env python3
"""Submit the ABRS bundle to PaintOmics AI and poll to completion.

Follows the exact multipart protocol established in scripts/05_submit_paintomics.py.
Organism: Arabidopsis thaliana (ath)
Databases: KEGG, MapMan
aiConsent: false (protecting unpublished predicted metabolome)

Jobs supported:
  --tissue shoot       (Primary: 12d light-grown shoots, ventilated ABRS, FLT vs GC)
  --tissue root        (Organ-specific: 12d roots, ventilated ABRS, FLT vs GC)
  --tissue hypocotyl   (Organ-specific: 12d hypocotyls, ventilated ABRS, FLT vs GC)
  --tissue whole_plant (Organ-specific: intact 12d seedlings, ventilated ABRS, FLT vs GC)

Upon completion, automatically retrieves and gzips the full server state payload
from POST /pa_recover_job and saves it in results/abrs/paintomics_raw/.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time

import requests

BASE = "https://paintomics.org"
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..")))
from paths import ensure  # noqa: E402

ABRS_UPLOAD = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "data", "abrs", "paintomics_upload"
))
ABRS_RAW = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "results", "abrs", "paintomics_raw"
))

SPECIE = "ath"
DATABASES = ["KEGG", "MapMan"]


def build(tissue: str, description: str, ai_consent: bool, include_metabolome: bool):
    """Assemble multipart fields and file handles for /pa_step1."""
    data: list[tuple[str, str]] = [
        ("specie", SPECIE),
        ("jobDescription", description),
        ("experimentDesign", ""),
        ("aiConsent", "true" if ai_consent else "false"),
    ]
    data += [("databases[]", db) for db in DATABASES]

    omics = [
        {
            "omic_name": "Gene expression",
            "file": f"abrs_{tissue}_gene_expression_values.tab",
            "file_type": "Gene Expression file",
            "relevant_file": f"abrs_{tissue}_gene_expression_relevant.tab",
            "relevant_file_type": "Relevant Genes list",
            "match_type": "Gene",
            "enrichment": "genes",
        }
    ]

    if include_metabolome:
        omics.append({
            "omic_name": "Metabolomics",
            "file": "abrs_metabolomics_values.tab",
            "file_type": "Metabolomic quatification",
            "relevant_file": "abrs_metabolomics_relevant.tab",
            "relevant_file_type": "Relevant Compound list",
            "match_type": "Compound",
            "enrichment": "features",
        })

    files, handles = [], []
    for i, om in enumerate(omics):
        for key, value in (
            ("omic_name", om["omic_name"]),
            ("origin", "client"), ("filelocation", ""),
            ("file_type", om["file_type"]),
            ("relevant_origin", "client"), ("relevant_filelocation", ""),
            ("relevant_file_type", om["relevant_file_type"]),
            ("design_origin", ""), ("design_filelocation", ""),
            ("match_type", om["match_type"]),
            ("enrichment", om["enrichment"]),
        ):
            data.append((f"omic{i}_{key}", value))

        for field, fname in (("file", om["file"]), ("relevant_file", om["relevant_file"])):
            path = os.path.join(ABRS_UPLOAD, fname)
            if not os.path.exists(path):
                raise SystemExit(f"missing {path}")
            fh = open(path, "rb")
            handles.append(fh)
            files.append((f"omic{i}_{field}", (fname, fh, "text/plain")))

    return data, files, handles, omics


def poll(session: requests.Session, job_id: str, timeout: int) -> dict:
    """Follow check_job_status until job finishes."""
    started = time.time()
    last = None
    while time.time() - started < timeout:
        try:
            r = session.post(f"{BASE}/check_job_status/{job_id}", timeout=60)
            info = r.json()
        except Exception as exc:
            print(f"    poll error: {exc}")
            time.sleep(10)
            continue

        if info.get("success") and "geneBasedInputOmics" in info:
            g = info.get("geneBasedInputOmics") or []
            c = info.get("compoundBasedInputOmics") or []
            print(f"    registered omics: {len(g)} gene-based, {len(c)} compound-based")
            return info

        status = info.get("status") or info.get("state")
        message = info.get("message") or info.get("description") or ""
        if (status, message) != last:
            print(f"    [{int(time.time()-started):4d}s] {status}  {str(message)[:90]}")
            last = (status, message)

        if status in ("finished", "Finished", "done", "error", "Error", "failed"):
            return info
        time.sleep(8)

    print(f"    still running after {timeout}s — stopping poll")
    return {"status": "timeout"}


def recover_full_job(session: requests.Session, job_id: str, tissue: str) -> None:
    """Recover full server result payload via POST /pa_recover_job."""
    print(f"\n  Recovering full job payload from {BASE}/pa_recover_job...")
    try:
        r = session.post(f"{BASE}/pa_recover_job", data={"jobID": job_id}, timeout=120)
        if r.status_code == 200:
            payload = r.json()
            out_path = os.path.join(ABRS_RAW, f"job_abrs_{tissue}_{job_id}_full.json.gz")
            with gzip.open(out_path, "wt") as fh:
                json.dump(payload, fh)
            sz = os.path.getsize(out_path)
            print(f"  Successfully recovered and saved: {out_path} ({sz:,} bytes)")
        else:
            print(f"  Recovery HTTP error: {r.status_code}")
    except Exception as e:
        print(f"  Recovery failed: {e}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tissue", default="shoot",
                    choices=["shoot", "root", "hypocotyl", "whole_plant"],
                    help="organ/tissue to analyse (default: shoot)")
    ap.add_argument("--ai-consent", action="store_true",
                    help="enable AI interpretation (sends data to llm.iiia.es)")
    ap.add_argument("--timeout", type=int, default=900, help="seconds to poll")
    ap.add_argument("--no-metabolome", action="store_true",
                    help="exclude metabolome layer (e.g. for root analysis)")
    args = ap.parse_args()

    include_metabolome = not args.no_metabolome and args.tissue in ("shoot", "whole_plant")
    desc = f"OSD-7/TAGES ABRS Arabidopsis {args.tissue}, flight vs ground control"
    if include_metabolome:
        desc += " with FvCB-predicted metabolome"

    session = requests.Session()
    session.headers.update({"User-Agent": "photoresp-multiomics/1.0"})
    session.get(BASE, timeout=60)

    data, files, handles, omics = build(args.tissue, desc, args.ai_consent, include_metabolome)
    print(f"Submitting ABRS ({args.tissue}) to {BASE}/pa_step1")
    print(f"  Organism: {SPECIE}, Databases: {', '.join(DATABASES)}")
    print(f"  aiConsent: {'true' if args.ai_consent else 'false'}, Layers: {len(omics)}")

    for i, om in enumerate(omics):
        f_sz = os.path.getsize(os.path.join(ABRS_UPLOAD, om["file"]))
        r_sz = os.path.getsize(os.path.join(ABRS_UPLOAD, om["relevant_file"]))
        print(f"  omic{i} {om['omic_name']:<16} {om['file']} ({f_sz:,} B) + {om['relevant_file']} ({r_sz:,} B)")

    ensure(ABRS_RAW)

    try:
        r = session.post(f"{BASE}/pa_step1", data=data, files=files, timeout=900)
    finally:
        for fh in handles:
            fh.close()

    print(f"\n  HTTP {r.status_code}")
    try:
        payload = r.json()
    except ValueError:
        print(f"  non-JSON response:\n{r.text[:1500]}")
        return 1

    resp_file = os.path.join(ABRS_RAW, f"job_abrs_{args.tissue}_step1_response.json")
    with open(resp_file, "w") as fh:
        json.dump(payload, fh, indent=1)

    if not payload.get("success", True):
        print(f"  rejected: {json.dumps(payload)[:1200]}")
        return 1

    job_id = payload.get("jobID") or payload.get("jobId")
    print(f"  jobID: {job_id}")
    print(f"  web view: {BASE}/?jobID={job_id}")
    if not job_id:
        return 1

    print("\n  polling queue...")
    final = poll(session, job_id, args.timeout)
    status_file = os.path.join(ABRS_RAW, f"job_abrs_{args.tissue}_status.json")
    with open(status_file, "w") as fh:
        json.dump(final, fh, indent=1)

    # Recover full results
    recover_full_job(session, job_id, args.tissue)
    return 0


if __name__ == "__main__":
    sys.exit(main())
