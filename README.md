# Fe-SMA XRD / DiffractGPT Pipeline

Computational materials science tools for the Fe-Mn-Al-Si-Ni-C shape memory alloy project (SMST 2026 / ISEF).

Uses **DiffractGPT** (AtomGPT.org) and **JARVIS-DFT** to identify crystal structures from XRD patterns and cross-validate experimental results with DFT thermodynamics.

## Scripts

| File | Purpose |
|------|---------|
| `fe_sma_diffractgpt.py` | XRD image/data -> DiffractGPT -> predicted POSCAR structure |
| `fe_sma_atomgpt_analysis.py` | JARVIS-DFT phase stability review of the Fe-SMA system |
| `fe_sma_atomgpt_review.txt` | Generated analysis report |
| `fe_sma_xrd_notes.txt` | Full usage notes and reference |

## Quick Start

```bash
pip install -r requirements.txt

# Demo (no files needed):
python fe_sma_diffractgpt.py demo

# With an XRD image:
python fe_sma_diffractgpt.py image path/to/xrd.jpg --formula Fe2MnAl

# With a peaks CSV (columns: 2theta intensity):
python fe_sma_diffractgpt.py manual peaks.csv --formula Fe2MnAl

# JARVIS-DFT phase analysis:
python fe_sma_atomgpt_analysis.py
```

## AGAPI Key (AtomGPT.org)

Set your API key to enable cloud DiffractGPT structure prediction:

```bash
export AGAPI_KEY=sk-your-key-here   # Linux/macOS
set AGAPI_KEY=sk-your-key-here      # Windows CMD
```

Without a key the script falls back to local JARVIS-DFT cosine-similarity XRD matching.

## Key Findings (JARVIS-DFT)

| Phase | Formation E | Bulk Mod | Notes |
|-------|------------|---------|-------|
| L2₁ MnAlFe₂ (target γ-phase) | −0.195 eV/atom | 207 GPa | Heusler Fm-3m |
| B2 AlFe (observed secondary) | −0.358 eV/atom | 175 GPa | 0.16 eV/atom more stable |

The B2 secondary phase is thermodynamically favored over the target L2₁ Heusler — consistent with synchrotron XRD showing dual-phase microstructure and explaining the absence of room-temperature superelasticity.

## References

- [DiffractGPT paper](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c03137)
- [AtomGPT paper](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c01126)
- [JARVIS-DFT database](https://jarvis.nist.gov/)
- [AtomGPT.org API](https://atomgpt.org)
