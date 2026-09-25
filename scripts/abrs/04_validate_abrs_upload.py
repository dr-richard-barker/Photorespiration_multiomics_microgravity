#!/usr/bin/env python3
"""Validate ABRS PaintOmics upload files against specification."""

from __future__ import annotations

import math
import os
import sys

ABRS_UPLOAD = os.path.normpath(os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "..", "data", "abrs", "paintomics_upload"
))

CHECKS = [
    ("abrs_shoot_gene_expression_values.tab", "abrs_shoot_gene_expression_relevant.tab", "#geneID"),
    ("abrs_root_gene_expression_values.tab", "abrs_root_gene_expression_relevant.tab", "#geneID"),
    ("abrs_hypocotyl_gene_expression_values.tab", "abrs_hypocotyl_gene_expression_relevant.tab", "#geneID"),
    ("abrs_whole_plant_gene_expression_values.tab", "abrs_whole_plant_gene_expression_relevant.tab", "#geneID"),
    ("abrs_metabolomics_values.tab", "abrs_metabolomics_relevant.tab", "#compound"),
]


def check(val_name: str, rel_name: str, expected_hdr: str) -> bool:
    val_path = os.path.join(ABRS_UPLOAD, val_name)
    rel_path = os.path.join(ABRS_UPLOAD, rel_name)

    if not os.path.exists(val_path):
        print(f"FAIL: missing {val_name}")
        return False
    if not os.path.exists(rel_path):
        print(f"FAIL: missing {rel_name}")
        return False

    with open(val_path) as fh:
        v_lines = fh.read().splitlines()
    with open(rel_path) as fh:
        r_lines = fh.read().splitlines()

    # Check header
    if not v_lines[0].startswith(expected_hdr):
        print(f"FAIL: {val_name} header should start with {expected_hdr}, got {v_lines[0][:30]}")
        return False

    v_ids = set()
    for ln in v_lines[1:]:
        if not ln.strip():
            continue
        parts = ln.split("\t")
        if len(parts) < 2:
            print(f"FAIL: line not tab-separated: {ln[:30]}")
            return False
        fid, val = parts[0], parts[1]
        try:
            fval = float(val)
            if not math.isfinite(fval):
                print(f"FAIL: non-finite float in {val_name}: {val}")
                return False
        except ValueError:
            print(f"FAIL: cannot parse float in {val_name}: {val}")
            return False
        v_ids.add(fid)

    r_ids = set()
    for ln in r_lines:
        if not ln.strip():
            continue
        r_ids.add(ln.strip())

    missing_in_vals = r_ids - v_ids
    if missing_in_vals:
        print(f"FAIL: {len(missing_in_vals)} relevant IDs not in values file for {val_name}")
        return False

    print(f"PASS: {val_name:<42} ({len(v_ids):>5} features) + {rel_name:<42} ({len(r_ids):>5} relevant)")
    return True


def main() -> int:
    print("Validating ABRS PaintOmics Upload Bundle:\n")
    all_ok = True
    for v, r, hdr in CHECKS:
        ok = check(v, r, hdr)
        all_ok = all_ok and ok
    print(f"\nOverall: {'ALL CHECKS PASSED' if all_ok else 'SOME CHECKS FAILED'}")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
