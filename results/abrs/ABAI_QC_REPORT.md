# ABAI Quality Control Audit Report

**Date:** 2026-10-02 19:43:45  
**Audit Verdict:** **CLEAR**  
**Repository Digest:** `32c7cef132917a91c1e7ec5cfe7f11635164ffcc00e80006b2d9dc31760e491e`  
**Tracked Files Audited:** 299  
**Git HEAD:** `aef7f97050e2469ca4ba587e119d8b96b91600b4`

---

## 1. Document & Figure Artifacts
- **Manuscript PDF:** `manuscript/latex/main.pdf` (Found manuscript/latex/main.pdf (775,549 bytes))
- **Supplementary PDF:** `manuscript/latex/supplementary.pdf` (Found manuscript/latex/supplementary.pdf (220,227 bytes))
- **Manuscript Word (.docx):** `manuscript/latex/main.docx` (Found manuscript/latex/main.docx (335,557 bytes))
- **Supplementary Word (.docx):** `manuscript/latex/supplementary.docx` (Found manuscript/latex/supplementary.docx (68,803 bytes))
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
