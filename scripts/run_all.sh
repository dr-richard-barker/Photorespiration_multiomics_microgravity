#!/usr/bin/env bash
# Regenerate every table and figure from scratch.
#
# Network is needed on a cold cache: OSDR downloads and KEGG lookups land in data/cache/,
# which is git-ignored. Everything after that is deterministic.
#
#   bash scripts/run_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."

echo "== 1. fetch OSD-522 =========================================================="
python3 scripts/00_fetch_osdr.py

echo "== 2. transcriptome (PyDESeq2) ==============================================="
python3 scripts/01_dge_rnaseq.py

echo "== 3. proteome ==============================================================="
python3 scripts/02_prep_proteomics.py

echo "== 4. model: FvCB sweep + predicted metabolome ==============================="
python3 scripts/fvcb.py
python3 scripts/03_predict_metabolome.py --sensitivity

echo "== 5. blind test against OSD-522 ============================================="
python3 scripts/04_falsification_check.py --quiet-resolution

echo "== 6. cross-study hardware ladder ============================================"
python3 scripts/07_hardware_ladder.py

echo "== 7. validate the PaintOmics upload bundle =================================="
python3 scripts/validate_upload.py

echo "== 8. figures ================================================================"
for f in scripts/figures/fig*.py; do python3 "$f"; done
# Fails if any figure carries a legend entry it never draws. Figure 1b shipped for weeks
# with five legend entries over three lines, because the column it plotted is empty for
# the two vented cases and matplotlib draws an all-NaN series as nothing, silently.
python3 scripts/check_figures.py
# The manuscript reads its own copy, so regenerating a figure without this step leaves
# the compiled PDF showing the old one — the exact drift the pre-publish check looks for.
cp results/figures/fig*.pdf results/figures/supplementary/fig*.pdf manuscript/latex/figures/

echo "== 9. site data + guards ====================================================="
python3 scripts/export_site_data.py
python3 scripts/check_site_data.py
python3 scripts/check_js_syntax.py
python3 scripts/check_js_parity.py

echo "== 10. manifest =============================================================="
python3 scripts/make_manifest.py

echo
echo "done. Tables in results/tables/, figures in results/figures/."
echo "PaintOmics submission is a separate, network-and-consent step:"
echo "  python3 scripts/05_submit_paintomics.py --job 1"
