#!/usr/bin/env python3
"""Submit the OSD-522 bundle to PaintOmics AI and poll the job to completion.

The browser route is not usable under automation: a file input cannot be driven without a
native file dialog, and the automation client blocks every cross-origin request the page
makes, so the page cannot be handed the files either. PaintOmics' own form posts a plain
multipart request to /pa_step1, so this script performs that same post directly.

The field names, file types, match types and enrichment types below were read off the live
ExtJS form (`form.getForm().getValues()`), not guessed.

  --job 1   three layers: genes + proteins + predicted metabolome  (FLT_vs_GC)
  --job 2   metabolome only: the sealed-vs-vented hardware contrast (BRIC_vs_VENTED)

`aiConsent` is sent as "false" by default: the AI interpretation forwards the data to a
separate LLM endpoint (llm.iiia.es), and the predicted metabolite layer is unpublished.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import requests

BASE = "https://paintomics.org"
HERE = os.path.dirname(os.path.abspath(__file__))
UPLOAD = os.path.join(HERE, "upload")
RESULTS = os.path.join(HERE, "results")

SPECIE = "ath"
DATABASES = ["KEGG", "MapMan"]

# One entry per omic slot, in the order PaintOmics numbers them.
JOB1 = [
    {"omic_name": "Gene expression", "file": "gene_expression_values.tab",
     "file_type": "Gene Expression file", "relevant_file": "gene_expression_relevant.tab",
     "relevant_file_type": "Relevant Genes list", "match_type": "Gene",
     "enrichment": "genes"},
    {"omic_name": "Proteomics", "file": "proteomics_values.tab",
     "file_type": "Proteomic quatification", "relevant_file": "proteomics_relevant.tab",
     "relevant_file_type": "Relevant proteins list", "match_type": "Gene",
     "enrichment": "features"},
    {"omic_name": "Metabolomics", "file": "metabolomics_values.tab",
     "file_type": "Metabolomic quatification", "relevant_file": "metabolomics_relevant.tab",
     "relevant_file_type": "Relevant Compound list", "match_type": "Compound",
     "enrichment": "features"},
]

JOB2 = [
    {"omic_name": "Metabolomics", "file": "metabolomics_hardware_values.tab",
     "file_type": "Metabolomic quatification",
     "relevant_file": "metabolomics_hardware_relevant.tab",
     "relevant_file_type": "Relevant Compound list", "match_type": "Compound",
     "enrichment": "features"},
]

DESCRIPTIONS = {
    1: "OSD-522 BRIC-LED-001 Arabidopsis shoots, flight vs ground: OSDR transcriptome + "
       "proteome with a CFD/FvCB-predicted metabolome",
    2: "Predicted metabolome only: sealed BRIC vs vented hardware (model output)",
}


def build(omics: list[dict], description: str, ai_consent: bool):
    """Assemble the multipart fields and file handles for /pa_step1."""
    data: list[tuple[str, str]] = [
        ("specie", SPECIE),
        ("jobDescription", description),
        ("experimentDesign", ""),
        ("aiConsent", "true" if ai_consent else "false"),
    ]
    data += [("databases[]", db) for db in DATABASES]

    files, handles = [], []
    for i, om in enumerate(omics):
        # `origin` says where the file comes from, and the server branches hard on it
        # (JobInformationManager.saveFiles): 'client' means the bytes are in this request,
        # 'mydata' means look up a previously stored file by name. The live form shows it
        # empty only because no file has been picked yet — choosing one sets 'client'.
        # Sending it empty is what made the first two submissions register zero omics.
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
            path = os.path.join(UPLOAD, fname)
            if not os.path.exists(path):
                raise SystemExit(f"missing {path}")
            fh = open(path, "rb")
            handles.append(fh)
            files.append((f"omic{i}_{field}", (fname, fh, "text/plain")))
    return data, files, handles


def poll(session: requests.Session, job_id: str, timeout: int) -> dict:
    """Follow check_job_status until the job leaves the queue."""
    started = time.time()
    last = None
    while time.time() - started < timeout:
        try:
            r = session.post(f"{BASE}/check_job_status/{job_id}", timeout=60)
            info = r.json()
        except Exception as exc:  # noqa: BLE001 — transient polling errors are expected
            print(f"    poll error: {exc}")
            time.sleep(10)
            continue

        if info.get("success") and "geneBasedInputOmics" in info:
            g = info.get("geneBasedInputOmics") or []
            c = info.get("compoundBasedInputOmics") or []
            print(f"    registered omics: {len(g)} gene-based, {len(c)} compound-based")
            for o in g + c:
                print(f"      {o.get('omicName')}: {o.get('omicSummary')}")
            return info
        status = info.get("status") or info.get("state")
        message = info.get("message") or info.get("description") or ""
        if (status, message) != last:
            print(f"    [{int(time.time()-started):4d}s] {status}  {str(message)[:90]}")
            last = (status, message)
        if status in ("finished", "Finished", "done", "error", "Error", "failed"):
            return info
        time.sleep(8)
    print(f"    still running after {timeout}s — stopping the poll, the job continues "
          f"server-side")
    return {"status": "timeout"}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--job", type=int, choices=[1, 2], required=True)
    ap.add_argument("--ai-consent", action="store_true",
                    help="enable the AI interpretation (sends data to llm.iiia.es)")
    ap.add_argument("--timeout", type=int, default=900, help="seconds to poll")
    args = ap.parse_args()

    omics = JOB1 if args.job == 1 else JOB2
    session = requests.Session()
    session.headers.update({"User-Agent": "photoresp-multiomics/1.0"})
    session.get(BASE, timeout=60)   # establish a session cookie

    data, files, handles = build(omics, DESCRIPTIONS[args.job], args.ai_consent)
    print(f"Submitting job {args.job} to {BASE}/pa_step1")
    print(f"  organism {SPECIE}, databases {', '.join(DATABASES)}, "
          f"aiConsent={'true' if args.ai_consent else 'false'}")
    for i, om in enumerate(omics):
        a = os.path.getsize(os.path.join(UPLOAD, om["file"]))
        b = os.path.getsize(os.path.join(UPLOAD, om["relevant_file"]))
        print(f"  omic{i} {om['omic_name']:<16} {om['file']} ({a:,} B) + "
              f"{om['relevant_file']} ({b:,} B)")

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

    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, f"job{args.job}_step1_response.json"), "w") as fh:
        json.dump(payload, fh, indent=1)

    if not payload.get("success", True):
        print(f"  rejected: {json.dumps(payload)[:1200]}")
        return 1

    job_id = payload.get("jobID") or payload.get("jobId")
    print(f"  jobID {job_id}")
    print(f"  view at {BASE}/?jobID={job_id}")
    if not job_id:
        print(f"  no jobID in response: {json.dumps(payload)[:800]}")
        return 1

    print("\n  polling the queue")
    final = poll(session, job_id, args.timeout)
    with open(os.path.join(RESULTS, f"job{args.job}_status.json"), "w") as fh:
        json.dump(final, fh, indent=1)
    print(f"\n  wrote {RESULTS}/job{args.job}_*.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
