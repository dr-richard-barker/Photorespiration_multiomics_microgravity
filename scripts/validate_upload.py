#!/usr/bin/env python3
"""Check the upload bundle against the structure of PaintOmics' own example files.

PaintOmics AI will offer to convert a non-conforming file in the browser, but a file it
has to guess at is a file we have lost control of. These checks are the ones the example
data implies: a '#'-prefixed header naming the ID column and one column per condition,
tab separators, no header at all in the relevant-features files, unique identifiers, and
no missing or non-finite values.

Run with --examples <dir> to compare against a downloaded PaintOmics example bundle.
"""

from __future__ import annotations

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paths import UPLOAD  # noqa: E402

# (values file, relevant file, expected header token, expected id pattern description)
BUNDLE = [
    ("gene_expression_values.tab", "gene_expression_relevant.tab", "#geneID", "AGI locus"),
    ("proteomics_values.tab", "proteomics_relevant.tab", "#proteinID", "UniProt accession"),
    ("metabolomics_values.tab", "metabolomics_relevant.tab", "#compound", "compound name"),
    ("metabolomics_hardware_values.tab", "metabolomics_hardware_relevant.tab", "#compound",
     "compound name"),
]

problems: list[str] = []
notes: list[str] = []


def fail(msg: str) -> None:
    problems.append(msg)


def check_values(path: str, expect_header: str) -> tuple[set[str], list[str]]:
    name = os.path.basename(path)
    if not os.path.exists(path):
        fail(f"{name}: missing")
        return set(), []

    with open(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        fail(f"{name}: empty")
        return set(), []

    header = lines[0]
    if not header.startswith("#"):
        fail(f"{name}: first line must start with '#', got {header[:40]!r}")
    if not header.startswith(expect_header):
        fail(f"{name}: header should begin {expect_header!r}, got {header.split(chr(9))[0]!r}")
    if "\t" not in header:
        fail(f"{name}: header has no tab — is it really tab-separated?")

    conditions = header.split("\t")[1:]
    ids, ncol = [], len(conditions) + 1
    for i, line in enumerate(lines[1:], start=2):
        if not line.strip():
            fail(f"{name}: blank line at {i}")
            continue
        parts = line.split("\t")
        if len(parts) != ncol:
            fail(f"{name}: line {i} has {len(parts)} fields, header has {ncol}")
            continue
        ids.append(parts[0])
        for j, val in enumerate(parts[1:]):
            try:
                x = float(val)
            except ValueError:
                fail(f"{name}: line {i} column {j+2} is not numeric: {val!r}")
                continue
            if math.isnan(x) or math.isinf(x):
                fail(f"{name}: line {i} column {j+2} is {val} — PaintOmics needs finite values")

    dupes = len(ids) - len(set(ids))
    if dupes:
        fail(f"{name}: {dupes} duplicate identifiers")
    print(f"  {name:<38} {len(ids):>6,} rows  {len(conditions)} condition(s): "
          f"{', '.join(conditions)}")
    return set(ids), conditions


def check_relevant(path: str, universe: set[str]) -> None:
    name = os.path.basename(path)
    if not os.path.exists(path):
        fail(f"{name}: missing")
        return
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not lines:
        fail(f"{name}: empty — PaintOmics needs at least one relevant feature")
        return
    if lines[0].startswith("#"):
        fail(f"{name}: relevant-features files carry no header, but line 1 starts with '#'")
    if any("\t" in ln for ln in lines):
        fail(f"{name}: contains a tab — this file is one identifier per line")

    missing = [x for x in lines if x not in universe]
    if missing:
        fail(f"{name}: {len(missing)} identifiers absent from the values file, "
             f"e.g. {missing[:3]}")
    dupes = len(lines) - len(set(lines))
    if dupes:
        fail(f"{name}: {dupes} duplicate identifiers")
    print(f"  {name:<38} {len(lines):>6,} relevant")


def compare_examples(example_dir: str) -> None:
    """Confirm our layout matches a downloaded PaintOmics example bundle byte-structurally."""
    ref = None
    for root, _dirs, files in os.walk(example_dir):
        if "gene_expression_values.tab" in files:
            ref = os.path.join(root, "gene_expression_values.tab")
            break
    if not ref:
        notes.append(f"no gene_expression_values.tab found under {example_dir} — skipped")
        return

    with open(ref) as fh:
        ref_header = fh.readline().rstrip("\n")
    with open(os.path.join(UPLOAD, "gene_expression_values.tab")) as fh:
        our_header = fh.readline().rstrip("\n")

    ref_tok, our_tok = ref_header.split("\t")[0], our_header.split("\t")[0]
    if ref_tok != our_tok:
        fail(f"header token differs from PaintOmics example: {ref_tok!r} vs {our_tok!r}")
    else:
        print(f"  header token matches PaintOmics example: {ref_tok!r}  ({ref})")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--examples", help="directory of a downloaded PaintOmics example bundle")
    args = ap.parse_args()

    print(f"Validating {UPLOAD}\n")
    all_conditions = {}
    for values, relevant, header, _idkind in BUNDLE:
        ids, conditions = check_values(os.path.join(UPLOAD, values), header)
        check_relevant(os.path.join(UPLOAD, relevant), ids)
        all_conditions[values] = conditions

    # The three files of the joint job must agree on their condition column, or PaintOmics
    # cannot line the omics up. The hardware file is a separate, metabolome-only job.
    joint = ["gene_expression_values.tab", "proteomics_values.tab", "metabolomics_values.tab"]
    named = {f: all_conditions.get(f) for f in joint}
    distinct = {tuple(v) for v in named.values() if v}
    if len(distinct) > 1:
        fail(f"joint-job files disagree on conditions: {named}")
    else:
        print(f"\n  joint job condition columns agree: {list(distinct)[0] if distinct else '—'}")

    if args.examples:
        compare_examples(args.examples)

    print()
    for n in notes:
        print(f"  note: {n}")
    if problems:
        print(f"FAILED — {len(problems)} problem(s):")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
