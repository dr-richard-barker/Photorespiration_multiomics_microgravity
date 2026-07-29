# Photorespiration multi-omics — microgravity

The subcellular photorespiration model that sits **inside the leaf**, downstream of the
gravity-dependent gas transport solved **outside the leaf** by
[LunarLeaf-CFD](https://github.com/dr-richard-barker/LunarLeaf-CFD), and upstream of the
transcriptomic / metabolic response reviewed in
[Hypoxia_vs_elevated_CO2_in_spaceflight](https://github.com/dr-richard-barker/Hypoxia_vs_elevated_CO2_in_spaceflight).

Photorespiration is the missing middle between those two: it is governed by the O₂/CO₂ ratio
*at Rubisco*, and that ratio is set by a chain of conductances whose top link — the leaf
boundary layer — is exactly what LunarLeaf-CFD computes as a function of gravity.

```
Bulk air Ca ──[ g_bl : gravity-dependent, from LunarLeaf-CFD ]──▶ leaf surface Cs
   Cs ──[ g_s stomata ]──▶ intercellular Ci ──[ g_m mesophyll ]──▶ chloroplast Cc, O₂
                          │
                          ▼   Rubisco carboxylation vs oxygenation (Vo/Vc = 2Γ*/Cc)
        CHLOROPLAST → PEROXISOME → MITOCHONDRION   (the C2 photorespiratory cycle)
                          │
                          ▼   transcriptomic + metabolic signature (OSD-38 / CO2_RNAseq)
```

## What is here now — `fvcb.py`

`fvcb.py` is the first component: an **FvCB oxygenation calculator**. It takes the
gravity-dependent boundary-layer conductance `g_bl(g)` reported by LunarLeaf-CFD, pushes it
through the `Ca → Cs → Ci → Cc` conductance chain, and solves the supply = demand operating
point of a Farquhar–von Caemmerer–Berry C3 model **with the Rubisco oxygenation term made
explicit** — so the output is the photorespiratory quantities, not just net assimilation:

| symbol | meaning |
|---|---|
| `Cc` | CO₂ mole fraction reaching Rubisco (chloroplast stroma) |
| `Vo/Vc` | oxygenation : carboxylation ratio = `2Γ*/Cc` |
| `phi` | oxygenation fraction `Vo/(Vc+Vo)` |
| `Rp` | photorespiratory CO₂ release = `0.5·Vo` (µmol m⁻² s⁻¹) |
| `Rp/A` | photorespiration relative to net assimilation |

Run it:

```bash
python3 fvcb.py       # prints the gravity-sweep table, writes photorespiration_vs_gravity.csv,
                      # and runs physical self-checks (Γ*, φ range, monotonic trends)
```

No dependencies beyond the Python standard library.

### Result (default parameters: Ca 400 µmol/mol, 25 °C, Q 1000, g_s 0.20, g_m 0.30 mol m⁻² s⁻¹)

| scale | g | g_bl | Cc | φ (%) | A | Rp/A (%) |
|---|---|---|---|---|---|---|
| leaf | Earth | 1.000 | 239 | 26.4 | 17.3 | 23.5 |
| leaf | Mars | 0.847 | 237 | 26.5 | 17.1 | 23.7 |
| leaf | Moon | 0.719 | 235 | 26.7 | 17.0 | 24.0 |
| leaf | micro-g | 0.494 | 229 | 27.2 | 16.5 | 24.8 |
| rosette | micro-g | 0.291 | 217 | 28.3 | 15.6 | 26.7 |
| canopy | micro-g | 0.109 | 180 | 32.2 | 12.6 | 34.4 |

**Reading it honestly:** at the **single-leaf** scale the boundary-layer effect on
photorespiration is real but *modest* (Rp/A 23.5 → 24.8 % from 1 g to µg), because the boundary
layer is a minority of the total CO₂ diffusion resistance (`1/g_bl` vs the larger `1/g_s + 1/g_m`).
The effect becomes **substantial in a dense microgreen canopy in microgravity** (`g_bl` collapses to
0.109, Cc falls to 180 µmol/mol, Rp/A rises to 34.4 %) — matching LunarLeaf-CFD's finding that the
gas-transport penalty amplifies leaf → rosette → canopy.

## Provenance and assumptions (no fabrication)

- **Rubisco kinetics + temperature responses** — in-vivo values of Bernacchi et al. (2001)
  *Plant Cell Environ.* 24:253–259 (Kc25 404.9, Ko25 278.4 mmol/mol, Γ\*25 42.75 µmol/mol) with
  their Arrhenius activation energies; Jmax response from Bernacchi et al. (2003).
- **`g_bl(g)`** — the converged boundary-layer conductance reported by LunarLeaf-CFD
  (README / `results/tables/T2`), anchored so the Earth single leaf = 1.0 mol m⁻² s⁻¹. The CFD
  surface O₂ excess (T5) is fed in too; at leaf scale it is small (a few µmol/mol on 210 000) so
  the coupling is Cc-dominated — stated rather than hidden.
- **Placeholders, not fits** — `Vcmax25`, `Jmax25`, `Rd25`, `g_s`, `g_m` are mid-range Arabidopsis
  literature values exposed as `LeafParams`. They are not fitted to spaceflight data yet.

## Next steps toward the multi-omics model

1. Ingest a LunarLeaf-CFD CSV export directly (arbitrary gravity/geometry sweep) instead of the
   hard-coded `CFD_GBL` table.
2. Make `g_s` CO₂-/humidity-responsive and couple `g_m` — currently fixed (the same limitation
   LunarLeaf-CFD flags for its own feedback loop).
3. Feed the predicted photorespiratory flux into a compartmentalised (chloroplast → peroxisome →
   mitochondrion) flux model constrained by the OSD-38 / CO2_RNAseq transcriptomics
   (scFEA/FLUXestimator), closing the loop to the observed spaceflight signature.
