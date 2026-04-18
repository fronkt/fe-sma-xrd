# Fe-SMA XRD / DiffractGPT Pipeline

Computational materials science tools for the Fe-Mn-Al-Si-Ni-C shape memory alloy project (SMST 2026 / ISEF).

Uses **DiffractGPT** (AtomGPT.org), **JARVIS-DFT**, **Materials Project**, and **DeepMind GNoME** to identify crystal structures from XRD patterns and cross-validate experimental results with DFT thermodynamics and AI-predicted stability.

## Scripts

| File | Purpose |
|------|---------|
| `fe_sma_diffractgpt.py` | XRD image/data → DiffractGPT → predicted POSCAR structure |
| `fe_sma_atomgpt_analysis.py` | JARVIS-DFT phase stability review of the Fe-SMA system |
| `gnome_optimade_query.py` | **[GNoME Option 1]** Query GNoME via OPTIMADE API for Fe-SMA structures |
| `gnome_screen.py` | **[GNoME Option 2]** Download + screen DeepMind's 520k-entry GNoME CSV locally |
| `gnome_predict.py` | **[GNoME Option 3]** CHGNet local energy/stability prediction (GNoME-compatible GNN) |
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

## GNoME Integration

Three complementary pathways to leverage DeepMind's 520,000-structure discovery dataset:

### Option 1 — OPTIMADE API (database search)
```bash
# Query GNoME for all structures containing Fe, Mn, Al
python gnome_optimade_query.py Fe Mn Al
# Outputs: gnome_optimade_results.json
```
Queries the community OPTIMADE mirror at `optimade-gnome.odbx.science`. Also wired into Claude Desktop as an MCP server — you can ask Claude to query it in natural language.

### Option 2 — Raw Data Screener (local data science)
```bash
# Download GNoME stable_materials_summary.csv (~150 MB) and screen for Fe-SMA system
python gnome_screen.py
# First run downloads the CSV; subsequent runs use the cache.
# Outputs: data/gnome/fe_sma_gnome_hits.csv
```
Screens all 520k+ GNoME entries for compositions within the Fe-Mn-Al-Si-Ni-C system, sorted by decomposition energy.

### Option 3 — Local GNN Prediction (CHGNet)
```bash
# Predict energy/stability of a structure using CHGNet (GNoME-compatible GNN)
python gnome_predict.py path/to/structure.POSCAR
python gnome_predict.py path/to/structure.cif

# Relax a structure (Python API):
from gnome_predict import relax_structure
result = relax_structure("my_structure.POSCAR", steps=200)
```
Uses **CHGNet v0.3** (PyTorch, Windows-native) as a local surrogate for GNoME stability prediction. Automatically cross-validates against JARVIS-DFT. Works on any POSCAR or CIF file; handles malformed POSCARs with prepended XRD data.

## References

- [DiffractGPT paper](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c03137)
- [AtomGPT paper](https://pubs.acs.org/doi/10.1021/acs.jpclett.4c01126)
- [JARVIS-DFT database](https://jarvis.nist.gov/)
- [AtomGPT.org API](https://atomgpt.org)
- [GNoME paper — Merchant et al., Nature 2023](https://www.nature.com/articles/s41586-023-06735-9)
- [GNoME dataset (Google Cloud)](https://storage.googleapis.com/gdm_materials_discovery/gnome_data/stable_materials_summary.csv)
- [GNoME OPTIMADE mirror](https://optimade-gnome.odbx.science)
- [CHGNet](https://chgnet.lbl.gov/)
- [MatGL](https://matgl.ai/)
