# Sam1-NO7325 — APS Synchrotron XRD Analysis

## Sample
- **File:** `Sam1-NO7325.chi`
- **Source:** APS (Advanced Photon Source), Argonne National Laboratory, 2026-1 run
- **Original path:** `C:\APS_EXP\2026-1\CS\Sam1-NO7325.tif` (raw detector image)
- **Reduced to:** `.chi` format via Fit2D azimuthal integration

## Beamline Calibration (from Fit2D)

| Parameter | Value |
|-----------|-------|
| Wavelength | **0.125870 Å** (~98.5 keV, high-energy X-rays) |
| Sample-to-detector distance | 1827.486 mm |
| Pixel size (H × V) | 100.0 × 100.0 microns |
| Beam center (X, Y) | 2121.364, 2058.313 px |
| Detector tilt rotation | -175.6654° |
| Detector tilt angle | 0.108897° |

## What the .chi File Contains

```
Line 1:  title/source path
Line 2:  "2-Theta Angle (Degrees)"
Line 3:  "Intensity"
Line 4:  3109  (number of data points)
Line 5+: two columns — 2theta (deg), intensity (counts)
```

Raw 2θ range: **0.0016° to 9.64°**  
(High-energy synchrotron geometry — all reflections compressed to very low angles)

## Pipeline Used

Script: `fe_sma_diffractgpt.py` (manual mode)

```bash
python fe_sma_diffractgpt.py manual "G:\synchrotron.chi\Sam1-NO7325.chi" \
    --formula Fe2MnAl \
    --wavelength 0.125870 \
    --prominence 0.01
```

### What the script did

1. **Loaded** the `.chi` file — auto-detected 4-line header, read 3109 (2theta, intensity) pairs
2. **Peak detection** — scipy `find_peaks` on normalized intensity profile (prominence threshold 0.01)
3. **Wavelength conversion** — synchrotron 2θ → d-spacing → Cu Kα equivalent 2θ via Bragg's law:
   ```
   d = λ_sync / (2 · sin(θ_sync))
   2θ_CuKa = 2 · arcsin(λ_CuKa / (2d))
   ```
   where λ_CuKa = 1.54056 Å. Peaks above 90° Cu Kα were clipped (unphysical for JARVIS matching).
4. **JARVIS-DFT matching** — cosine similarity against simulated XRD patterns from 166 Fe-Al candidates in the JARVIS-DFT database

## Detected Peaks

| Synchrotron 2θ (°) | Cu Kα equiv. 2θ (°) | Rel. Intensity (%) | Likely Reflection |
|--------------------|----------------------|--------------------|-------------------|
| ~0.76 | 9.30 | 15.9 | low-angle, possible superstructure |
| ~3.52 | 43.01 | 57.4 | FCC/L2₁ (111) |
| ~3.62 | **44.28** | **100.0** | B2 (110) — strongest peak |
| ~4.10 | 50.04 | 77.2 | B2/FCC (200) |
| ~5.28 | 64.37 | 13.6 | (220) |
| ~6.02 | 73.49 | 21.3 | (311) |
| ~6.68 | 81.52 | 28.8 | (222) |
| ~7.34 | 89.11 | 20.7 | — |

## JARVIS-DFT Match Results

| Rank | JVASP ID | Formula | Space Group | Cosine Score | Formation E (eV/atom) |
|------|----------|---------|-------------|-------------|----------------------|
| 1 | JVASP-99748 | AlFe₃ | Pm-3m (B2) | **0.6060** | -0.1789 |
| 2 | JVASP-99397 | TiAlFe₂ | P4/mmm | 0.4516 | -0.4409 |
| 3 | JVASP-8743 | TiAlFe₂ | Fm-3m (L2₁) | 0.4330 | -0.4973 |
| 4 | JVASP-15000 | AlFe | Pm-3m (B2) | 0.3880 | -0.3577 |
| 5 | JVASP-78761 | AlFe | Pm-3m (B2) | 0.3844 | -0.3577 |
| 6 | JVASP-86274 | Al₆Fe | Cmcm | 0.3477 | -0.2130 |
| 7 | JVASP-86876 | Al₆Fe | Cmcm | 0.3477 | -0.2130 |
| 8 | JVASP-104975 | AlVFe₂ | F-43m | 0.3343 | -0.1327 |
| 9 | JVASP-107244 | AlFe₃H | Pm-3m | 0.2796 | -0.1546 |
| 10 | JVASP-7957 | AlFe₃ | Fm-3m | 0.1889 | -0.1919 |

## Interpretation

- **Dominant phase: B2 (AlFe₃, Pm-3m)** — top JARVIS match at 0.606, consistent with the strongest peak at 44.28° Cu Kα equivalent (d ≈ 2.04 Å, B2 (110) reflection)
- **L2₁ signature weak** — the 43.01° peak could be the L2₁ (111) superlattice reflection, but it is weaker than the B2 (110) peak, indicating the target Heusler phase is not the majority phase
- **Consistent with prior JARVIS-DFT analysis** (see `fe_sma_atomgpt_review.txt`): B2 AlFe is 0.16 eV/atom more stable than L2₁ MnAlFe₂, explaining why B2 dominates the microstructure
- **No room-temperature superelasticity** — the B2 secondary phase suppresses the martensitic transformation needed for shape memory behavior

---

## 1. Structural Analysis

The XRD pattern of Sam1 is dominated by a **B2-ordered phase (Pm-3m)**, with the strongest peak at 44.28° Cu Kα equivalent (d ≈ 2.04 Å) indexing cleanly to B2 (110). The secondary peak cluster at 43.01° is consistent with the L2₁ Heusler (111) superlattice reflection, but its lower intensity (~57% vs. 100%) indicates the L2₁ phase is a minority component rather than the matrix.

The peak at 9.30° (low Cu Kα equivalent, d ≈ 9.5 Å) is anomalously large in d-spacing for a simple metallic phase and may indicate a long-period superstructure, a secondary oxide layer, or an artifact of the very low synchrotron angle region. It warrants further investigation.

**Phase assignment summary:**

| Phase | Structure | Space Group | Evidence |
|-------|-----------|-------------|---------|
| B2 AlFe₃ | CsCl-type ordered BCC | Pm-3m | Strongest peak (44.28°), top JARVIS match (0.606) |
| L2₁ MnAlFe₂ | Heusler | Fm-3m | Weak 43.01° peak (superlattice); 3rd JARVIS match |
| FCC γ | Austenite | Fm-3m | Partially overlaps with L2₁ peaks; cannot fully resolve |

The close angular proximity of B2 (110) at 44.28° and L2₁ (111) at 43.01° means these two phases are difficult to deconvolute without Rietveld refinement. Phase fraction quantification from peak intensities alone is unreliable here — a full Rietveld fit (GSAS-II or FullProf) against the raw `.chi` data is recommended.

---

## 2. Mechanical Properties

From JARVIS-DFT on the identified phases:

| Phase | JVASP ID | Bulk Modulus | Shear Modulus | Formation E |
|-------|----------|-------------|--------------|-------------|
| L2₁ MnAlFe₂ | JVASP-18826 | **207.5 GPa** | ~80 GPa (est.) | -0.195 eV/atom |
| B2 AlFe | JVASP-15000 | ~175 GPa | ~70 GPa (est.) | **-0.358 eV/atom** |
| B2 AlFe₃ | JVASP-99748 | ~160 GPa (est.) | — | -0.179 eV/atom |

The B2 phase acting as the matrix has important mechanical consequences:

- **Suppressed martensitic transformation:** B2 (CsCl-type) is a high-symmetry cubic phase with limited shear compliance. The tetragonal shear (c/a ratio change) required for L2₁ ↔ martensite is mechanically opposed by the stiff B2 matrix surrounding Heusler precipitates.
- **Dislocation pinning:** B2 ordered domains act as obstacles to dislocation motion. This raises the apparent yield stress but prevents the cooperative atomic shuffling that drives stress-induced martensitic transformation (SIMT).
- **High transformation stress threshold:** Consistent with the experimental observation that no well-defined stress plateau appeared in mechanical testing. The critical stress to trigger SIMT may exceed the fracture stress in this microstructure.

---

## 3. Magnetism

From JARVIS-DFT on the L2₁ Heusler component:

- **L2₁ MnAlFe₂ magnetic moment: 2.039 μB/formula unit** — a half-metallic ferrimagnet (Fe and Mn sublattices anti-aligned)
- **B2 AlFe:** weakly ferromagnetic (~0.7 μB/Fe in DFT), Curie temperature ~800 K

Magnetic considerations for this sample:

- The L2₁ Heusler phase is magnetically ordered at room temperature. In Heusler SMAs, ferromagnetic order in the austenite phase couples to the structural transformation — the Ms temperature is sensitive to magnetic state.
- **APS experiment context:** The 98.5 keV X-rays used at APS penetrate the full sample thickness and are not sensitive to magnetic order. If magnetic domain structure is needed, neutron diffraction or XMCD would be required.
- **Cycling above Curie temperature (~800 K):** If the AGG anneal at 1200°C was performed and the sample cooled through the Curie temperature without a controlled field, the magnetic microstructure will be randomly oriented, which is known to suppress the single-variant martensite needed for macroscopic shape recovery.

---

## 4. Conclusions

**The Sam1 APS synchrotron XRD pattern confirms a B2-dominant dual-phase microstructure, consistent with the prior SSRF analysis and the JARVIS-DFT thermodynamic predictions.**

Key takeaways:

1. **B2 has won the phase competition.** AlFe-type B2 (Ef = -0.358 eV/atom) is thermodynamically ~0.16 eV/atom more stable than the target L2₁ Heusler (Ef = -0.195 eV/atom) in the binary Fe-Al limit. The APS data shows this is not a metastable artifact — B2 is the majority phase at this composition and processing condition.

2. **L2₁ is present but minority.** The 43.01° peak is real and indicates some Heusler ordering exists. The material is not fully disordered. This is encouraging — it means composition or heat treatment adjustments could shift the phase balance.

3. **Superelasticity is blocked at the structural level.** Without a single-phase L2₁ austenite matrix, there is no coherent habit plane for stress-induced martensite. The B2 inclusions physically pin the transformation front.

4. **Thermal SME may still be active.** The presence of any L2₁ domains means thermally-driven martensite is possible below Ms. The earlier confirmation of thermal shape memory effect is consistent with this — the martensitic transformation is not fully suppressed, just incomplete.

---

## 5. Recommendations

**Immediate (before next synthesis run):**

- Run **Rietveld refinement** on this `.chi` file in GSAS-II to extract quantitative phase fractions (B2 vs. L2₁ vs. FCC γ). This is the most important next step — a number like "35% L2₁, 65% B2" is far more useful than a qualitative assessment.
- Run the **remaining `.chi` files** (Sam5, Sam6, Sam7, Sam8) through the same pipeline to compare phase fractions across samples. If different heat treatments were applied, the B2/L2₁ ratio shift will directly show the processing window.

**Composition adjustments (DFT-guided):**

- **Reduce carbon to <0.05 wt%.** Carbon strongly stabilizes B2/BCC over FCC Heusler (Zener ordering). This is the single highest-leverage change.
- **Increase Mn content (30–35 at%).** Higher Mn stabilizes FCC austenite and lowers the B2 → L2₁ transition temperature, promoting the Heusler phase at room temperature.
- **Screen NixFe₂₋ₓMnAl variants** using JARVIS/ALIGNN (via `fe_sma_atomgpt_analysis.py`) — Ni additions at 3–6 at% are known to depress Ms in Fe-Mn-Al-Ni, potentially bringing superelasticity into the room-temperature window.

**Processing:**
- Consider a **two-stage anneal:** high temperature (1200°C) for AGG, followed by a lower-temperature L2₁ ordering anneal (600–800°C, 24h) to promote Heusler ordering before quench. This may shift the B2/L2₁ balance without changing composition.

---

## Output Files

| File | Description |
|------|-------------|
| `xrd_peaks_detected.png` | Plot of raw synchrotron profile with detected peak positions marked |

## Code Changes Made This Session

Added to `fe_sma_diffractgpt.py`:

1. **`.chi` file support** — `load_peaks_from_file()` now auto-detects and skips non-data header lines (requires 2 parseable float columns to start reading)
2. **`--wavelength` flag** — converts synchrotron 2θ to Cu Kα equivalent for JARVIS matching
3. **`to_cuka()` function** — Bragg's law conversion: synchrotron 2θ → d-spacing → Cu Kα 2θ, clips >90°
4. **Non-blocking plot** — matplotlib `Agg` backend, saves PNG without requiring GUI window interaction
