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
# All of the above is at saturating light; this records what changes at OSD-522's 60.
python3 scripts/09_light_sensitivity.py
# The "no plant metabolome in OSDR" census.
python3 scripts/osdr_census.py

echo "== 5. blind test against OSD-522 ============================================="
python3 scripts/04_falsification_check.py --quiet-resolution

echo "== 6. cross-study hardware ladder ============================================"
python3 scripts/07_hardware_ladder.py

echo "== 7. validate the PaintOmics upload bundle =================================="
python3 scripts/validate_upload.py

echo "== 7b. PaintOmics result tables ==============================================="
# T08/T09/T10 come from the saved job record, not from the web interface. They used to be
# transcribed by hand, with the manifest crediting a script that only submits the job.
python3 scripts/06_paintomics_tables.py

echo "== 8. figures ================================================================"
# fig09_physiospace.py is drawn in step 8b, only when its R step has run.
for f in scripts/figures/fig*.py; do
    [ "$f" = scripts/figures/fig09_physiospace.py ] && continue
    python3 "$f"
done
# Fails if any figure carries a legend entry it never draws. Figure 1b shipped for weeks
# with five legend entries over three lines, because the column it plotted is empty for
# the two vented cases and matplotlib draws an all-NaN series as nothing, silently.
python3 scripts/check_figures.py
# The manuscript reads its own copy, so regenerating a figure without this step leaves
# the compiled PDF showing the old one — the exact drift the pre-publish check looks for.
cp results/figures/fig*.pdf results/figures/supplementary/fig*.pdf manuscript/latex/figures/

echo "== 8b. PhysioSpace projection + heatmap (optional) ==========================="
# Needs the PhysioSpaceMethods and PlantPhysioSpace R packages. Without them the step is
# skipped with a warning, and so is the heatmap: there is no fallback. An earlier version
# swallowed the R failure here and fig09_physiospace.py then drew the heatmap from random
# numbers. With the packages present, a failure in either script stops the run.
if command -v Rscript &> /dev/null && Rscript -e 'quit(status = !all(c("PhysioSpaceMethods", "PlantPhysioSpace") %in% rownames(installed.packages())))' 2>/dev/null; then
    Rscript scripts/08_physiospace_projection.R
    python3 scripts/figures/fig09_physiospace.py
else
    echo "WARNING: PhysioSpace SKIPPED (R packages not installed). No T14 scores and no PhysioSpace figure were made." >&2
fi

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

