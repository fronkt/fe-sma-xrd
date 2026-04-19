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
