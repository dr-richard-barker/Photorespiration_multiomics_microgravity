#!/usr/bin/env bash
# run_abrs.sh — End-to-end pipeline for ABRS/TAGES (OSD-7 & OSD-16) analysis
# and integrated comparison with the sealed BRIC-LED analysis.

set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO}"

echo "================================================================="
echo " NASA TAGES (ABRS) Multi-Omics Analysis & Hardware Integration   "
echo "================================================================="
echo

echo "Step 1: Fetching OSD-7 GeneLab microarray data & ISA-Tab metadata..."
python3 scripts/abrs/00_fetch_abrs.py

echo
echo "Step 2: Computing organ-specific microarray differential expression..."
python3 scripts/abrs/01_dge_microarray.py

echo
echo "Step 3: Simulating FvCB operating points & predicting ABRS metabolome..."
python3 scripts/abrs/03_predict_abrs_metabolome.py

echo
echo "Step 4: Validating ABRS PaintOmics upload bundle..."
python3 scripts/abrs/04_validate_abrs_upload.py

echo
echo "Step 5: Computing pathway enrichment & cross-hardware comparison tables..."
python3 scripts/abrs/06_abrs_paintomics_tables.py

echo
echo "Step 6: Executing extended hardware ladder analysis..."
python3 scripts/abrs/07_abrs_hardware_comparison.py

echo
echo "Step 7: Generating publication-ready comparison figure..."
python3 scripts/abrs/fig_abrs_comparison.py

echo
echo "================================================================="
echo " All ABRS analysis steps completed successfully!                "
echo " Tables:  results/abrs/tables/                                   "
echo " Figures: results/abrs/figures/                                  "
echo " Uploads: data/abrs/paintomics_upload/                           "
echo "================================================================="
