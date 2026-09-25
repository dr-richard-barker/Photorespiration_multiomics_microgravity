# ABAI Quality Control Audit Report

**Date:** 2026-09-25 15:43:19  
**Audit Verdict:** **CLEAR**  
**Repository Digest:** `8488af5bb86f0a1f9c205750e108246b795277a1efc69a55b4cf902e329e6f21`  
**Tracked Files Audited:** 297  
**Git HEAD:** `edb9068e4b2bbffa7e977747fa02f6d0e7bd7065`

---

## 1. Document & Figure Artifacts
- **Manuscript PDF:** `manuscript/latex/main.pdf` (Found manuscript/latex/main.pdf (777,286 bytes))
- **Supplementary PDF:** `manuscript/latex/supplementary.pdf` (Found manuscript/latex/supplementary.pdf (220,227 bytes))
- **Manuscript Word (.docx):** `manuscript/latex/main.docx` (Found manuscript/latex/main.docx (337,453 bytes))
- **Supplementary Word (.docx):** `manuscript/latex/supplementary.docx` (Found manuscript/latex/supplementary.docx (68,804 bytes))
- **Figures:** All 10 figures (Figures 1–9 and Figure S1) present in `manuscript/latex/figures/` and confirmed non-empty.

## 2. LaTeX Compilation & Typography Audit
- **Main Log (`main.log`):**
  - Fatal errors: 0
  - Undefined references (`??`): 0
  - Undefined citations (`[?]`): 0
  - Severe horizontal overflow (>5pt): 0
- **Supplementary Log (`supplementary.log`):**
  - Fatal errors: 0
  - Undefined references: 0
  - Undefined citations: 0
  - Severe horizontal overflow (>5pt): 0

## 3. Literature & Citation Integrity
- **BibTeX Database:** `references.bib` containing 30 peer-reviewed entries.
- **In-text Citations:** 29 citations verified.
- **Unresolved / Broken Citations:** 0
- **Fabricated DOI Sweeps:** Passed (0 flagged).

## 4. Quantitative Claim Verification (Grounding)
| Claim / Metric | Table Source | Manuscript Grounding | Status |
| :--- | :--- | :--- | :--- |
| ER Stress BRIC-LED vs ABRS | T08 / T10 | $p = 8.65\times 10^{-9}$ vs $p = 0.9456$ | **PASS** |
| Biotic Stress Conserved | T08 / T10 | $p = 2.10\times 10^{-4}$ vs $p = 0.0151$ | **PASS** |
| Hardware Ladder Parity | T11 / T12 | $\Delta_{\text{VEGGIE}} = +0.4142$ vs $\Delta_{\text{ABRS}} = +0.4146$ | **PASS** |
| FvCB Assimilation Rate | T05 | $A = 28.85\ \mu\text{mol m}^{-2}\text{ s}^{-1}$ vs $28.87$ | **PASS** |
| Oxygenation Fraction $\varphi$ | T05 | $\varphi = 2.6\%$ (ABRS) vs $53.7\%$ (BRIC-LED) | **PASS** |

## 5. Subsystem QA Suites
- `scripts/check_figures.py`: **PASS** (all 9 core figures + S1 validated against source series)
- `scripts/check_js_parity.py`: **PASS** (180 points between `docs/app/fvcb.js` and `scripts/fvcb.py` agree to $< 10^{-6}$)
- `scripts/check_js_syntax.py`: **PASS** (all ES modules parse cleanly)
- `scripts/check_site_data.py`: **PASS** (all gene sets and pathway memberships verified)

---
*Generated automatically by `scripts/abai_qc_audit.py` under the ABAI Quality Control framework.*
