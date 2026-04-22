"""
DiffractGPT pipeline for Fe-SMA XRD analysis.

Two modes:
  1. image  — load an XRD image, detect peaks automatically, run DiffractGPT
  2. manual — supply a CSV/text file of (2theta, intensity) pairs directly

Usage:
  python fe_sma_diffractgpt.py image  path/to/xrd_image.jpg  --formula Fe2MnAl
  python fe_sma_diffractgpt.py manual path/to/peaks.csv      --formula Fe2MnAl
  python fe_sma_diffractgpt.py demo                          # runs built-in demo

Requirements:
  pip install numpy scipy matplotlib pillow jarvis-tools

AGAPI / AtomGPT.org (optional, for cloud DiffractGPT):
  Register at https://atomgpt.org to get an API key, then set:
    export AGAPI_KEY="sk-your-key-here"
  Without a key the script falls back to jarvis-tools XRD matching.
"""

import os
import sys
import argparse
import numpy as np
from pathlib import Path


# ── helpers ──────────────────────────────────────────────────────────────────

def load_xrd_from_image(image_path: str, row_frac: float = 0.5) -> tuple[np.ndarray, np.ndarray]:
    """
    Extract a 1-D intensity profile from an XRD image by averaging rows near
    row_frac (default: middle of the image).  Works best on 1-D line-detector
    or already-integrated synchrotron output images.
    """
    from PIL import Image as PILImage

    img = PILImage.open(image_path).convert("L")
    arr = np.array(img, dtype=float)

    # Average a band of rows (±5% of height) around row_frac
    h, w = arr.shape
    r0 = max(0, int((row_frac - 0.05) * h))
    r1 = min(h, int((row_frac + 0.05) * h))
    profile = arr[r0:r1, :].mean(axis=0)

    # Map pixel x-axis to a plausible 2θ range (10°–90°)
    two_theta = np.linspace(10, 90, w)
    return two_theta, profile


def detect_peaks(two_theta: np.ndarray, intensity: np.ndarray,
                 min_prominence: float = 0.05) -> tuple[np.ndarray, np.ndarray]:
    """Return (2theta_peaks, intensity_peaks) above prominence threshold."""
    from scipy.signal import find_peaks

    i_norm = intensity / intensity.max()
    prominence = min_prominence * i_norm.max()
    peaks, props = find_peaks(i_norm, prominence=prominence, distance=5)

    peak_2theta = two_theta[peaks]
    peak_intensity = (i_norm[peaks] / i_norm[peaks].max() * 100).round(1)
    return peak_2theta, peak_intensity


def load_peaks_from_file(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Load (2theta, intensity) from CSV, whitespace, or .chi file."""
    lines = Path(path).read_text().splitlines()

    # .chi: skip header lines until we hit a line that starts with a float
    skip = 0
    for line in lines:
        stripped = line.strip()
        if not stripped:
            skip += 1
            continue
        try:
            parts = stripped.split()
            if len(parts) >= 2:
                float(parts[0]); float(parts[1])
                break
        except (ValueError, IndexError):
            pass
        skip += 1

    data = np.loadtxt(path, skiprows=skip, usecols=(0, 1))
    return data[:, 0], data[:, 1]


CUKA = 1.54056  # Å


def to_cuka(two_theta: np.ndarray, wavelength: float) -> np.ndarray:
    """Convert synchrotron 2θ → Cu Kα 2θ via d-spacing (Bragg's law)."""
    theta_rad = np.deg2rad(two_theta / 2)
    d = wavelength / (2 * np.sin(theta_rad))
    sin_cuka = CUKA / (2 * d)
    valid = sin_cuka <= 1.0
    result = np.full_like(two_theta, np.nan)
    result[valid] = np.rad2deg(2 * np.arcsin(sin_cuka[valid]))
    return result


def format_for_diffractgpt(formula: str, two_theta: np.ndarray,
                            intensity: np.ndarray) -> str:
    """
    Format peaks in the DiffractGPT text format:
      Line 1: formula
      Lines 2+: 2theta intensity
    """
    lines = [formula]
    for t, i in zip(two_theta, intensity):
        lines.append(f"{t:.3f} {i:.1f}")
    return "\n".join(lines)


def run_xrd_match_jarvis(two_theta: np.ndarray, intensity: np.ndarray) -> list[dict]:
    """
    Match the XRD pattern against JARVIS-DFT using cosine similarity.
    Uses jarvis.analysis.diffraction.xrd.XRD to simulate reference patterns.
    """
    try:
        from jarvis.db.figshare import data as jarvis_data
        from jarvis.core.atoms import Atoms
        from jarvis.analysis.diffraction.xrd import XRD

        print("Loading JARVIS-DFT database for XRD matching (first run downloads ~41 MB)...")
        dft = jarvis_data("dft_3d")

        # Build experimental pattern vector on a shared Cu Kα 2θ grid
        grid = np.linspace(max(two_theta.min(), 10), min(two_theta.max(), 90), 500)
        sigma = 0.2  # Gaussian broadening (degrees)

        def smear(thetas, intensities):
            profile = np.zeros_like(grid)
            for t, i in zip(thetas, intensities):
                profile += i * np.exp(-0.5 * ((grid - t) / sigma) ** 2)
            return profile / (profile.max() + 1e-9)

        exp_profile = smear(two_theta, intensity)

        # Score a subset of Fe-Al-Mn entries for speed
        candidates = [m for m in dft if "Fe" in m["formula"] and "Al" in m["formula"]]
        print(f"Scoring {len(candidates)} Fe-Al candidates against experimental pattern...")

        results = []
        for mat in candidates[:60]:
            try:
                atoms = Atoms.from_dict(mat["atoms"])
                xrd = XRD()
                xrd.simulate(atoms)
                ref_profile = smear(xrd.two_theta_array, xrd.intensity_array)
                cos_sim = float(np.dot(exp_profile, ref_profile) /
                                (np.linalg.norm(exp_profile) * np.linalg.norm(ref_profile) + 1e-9))
                results.append({
                    "jid": mat["jid"],
                    "formula": mat["formula"],
                    "spg": mat.get("spg_symbol", "?"),
                    "score": round(cos_sim, 4),
                    "formation_energy": mat.get("formation_energy_peratom"),
                })
            except Exception:
                continue

        results.sort(key=lambda x: x["score"], reverse=True)
        return results[:10]

    except Exception as e:
        print(f"[XRD match error] {e}")
        return []


def run_diffractgpt_api(formula: str, peaks_text: str, api_key: str) -> str | None:
    """Call AtomGPT.org DiffractGPT endpoint if API key is available."""
    try:
        import httpx
        r = httpx.get(
            "https://atomgpt.org/diffractgpt/query",
            params={"formula": formula, "peaks": peaks_text, "APIKEY": api_key},
            timeout=60,
        )
        r.raise_for_status()
        # Response is plain text POSCAR (with comment header lines starting with #)
        text = r.text.strip()
        if not text:
            return None
        # Strip comment lines to isolate the POSCAR block
        poscar_lines = [l for l in text.splitlines() if not l.startswith("#")]
        return "\n".join(poscar_lines).strip() or None
    except Exception as e:
        print(f"[AGAPI DiffractGPT error] {e}")
        return None


# Reference d-spacing windows for Fe-SMA phases (Cu Kα, λ=1.54056 Å).
# Tuples: (d_min Å, d_max Å, short label, color)
_PHASE_REFS = [
    (2.07, 2.15, "FCC γ (111)",        "#e07b39"),
    (1.98, 2.07, "B2 (110)\nL2₁(220)", "#2ca02c"),
    (1.78, 1.87, "FCC γ (200)",        "#e07b39"),
    (1.42, 1.46, "B2 (200)\nL2₁(400)", "#2ca02c"),
    (1.27, 1.30, "FCC γ (220)\nL2₁(420)", "#9467bd"),
    (1.16, 1.19, "B2 (211)\nL2₁(422)", "#2ca02c"),
]


def label_phases(peak_2theta_cuka: np.ndarray) -> list[tuple[str, str]]:
    """Return (label, color) for each Cu Kα 2θ peak from d-spacing matching."""
    out = []
    for t in peak_2theta_cuka:
        d = CUKA / (2 * np.sin(np.deg2rad(t / 2)))
        lbl, col = f"d={d:.3f}Å", "gray"
        for d_min, d_max, name, color in _PHASE_REFS:
            if d_min <= d <= d_max:
                lbl, col = name, color
                break
        out.append((lbl, col))
    return out


def plot_xrd(two_theta: np.ndarray, intensity: np.ndarray,
             peak_2theta: np.ndarray, peak_intensity: np.ndarray,
             title: str = "XRD Profile with Detected Peaks",
             cuka_peaks: np.ndarray | None = None):
    """Plot XRD profile with detected peaks.

    Args:
        cuka_peaks: Cu Kα equivalent 2θ values used for phase labeling. If
                    provided, each peak is annotated with its phase assignment
                    derived from d-spacing matching. If None, only the angle is
                    shown (raw synchrotron mode).
    """
    import matplotlib
    matplotlib.use("Agg")  # non-interactive — saves file, no GUI window
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.plot(two_theta, intensity / intensity.max() * 100, lw=1.2, color="steelblue", label="Profile")
    ax.scatter(peak_2theta, peak_intensity, color="red", zorder=5, s=40, label=f"{len(peak_2theta)} peaks")

    phases = label_phases(cuka_peaks) if cuka_peaks is not None else None

    for idx, (t, i) in enumerate(zip(peak_2theta, peak_intensity)):
        if phases:
            lbl, col = phases[idx]
            ax.annotate(
                f"{t:.2f}°\n{lbl}",
                (t, i),
                textcoords="offset points",
                xytext=(0, 8),
                ha="center",
                fontsize=6.5,
                color=col,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=col, alpha=0.75, lw=0.6),
            )
        else:
            ax.annotate(f"{t:.1f}°", (t, i), textcoords="offset points", xytext=(0, 6),
                        ha="center", fontsize=7, color="darkred")

    if phases:
        from matplotlib.patches import Patch
        seen = {}
        for lbl, col in phases:
            seen.setdefault(lbl, col)
        legend_handles = [ax.lines[0], ax.collections[0]]
        for lbl, col in seen.items():
            legend_handles.append(Patch(facecolor=col, edgecolor=col, alpha=0.6, label=lbl))
        ax.legend(handles=legend_handles, fontsize=7, loc="upper right",
                  title="Phase key", title_fontsize=7)
    else:
        ax.legend()

    ax.set_xlabel("2θ (degrees)")
    ax.set_ylabel("Relative Intensity (%)")
    ax.set_title(title)
    plt.tight_layout()
    out = "xrd_peaks_detected.png"
    plt.savefig(out, dpi=150)
    plt.close()
    print(f"Saved peak plot -> {out}")


# ── demo data (FCC γ-phase peaks typical for Fe-Mn-Al-Si-Ni-C @ Cu Kα) ──────

DEMO_PEAKS = np.array([
    # 2theta,  relative intensity
    [43.8,  100.0],   # FCC (111)
    [50.9,   42.0],   # FCC (200)
    [74.7,   20.0],   # FCC (220)
    [43.2,   18.0],   # Fe3Al (110) / B2 overlap
    [31.5,    8.0],   # Fe3Al (110) superlattice
])


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="DiffractGPT Fe-SMA pipeline")
    parser.add_argument("mode", choices=["image", "manual", "demo"])
    parser.add_argument("input", nargs="?", help="Path to image or peaks file")
    parser.add_argument("--formula", default="Fe2MnAl",
                        help="Chemical formula for DiffractGPT (default: Fe2MnAl)")
    parser.add_argument("--row", type=float, default=0.5,
                        help="Fractional row to sample from image (0–1, default 0.5)")
    parser.add_argument("--prominence", type=float, default=0.05,
                        help="Min peak prominence fraction (default 0.05)")
    parser.add_argument("--wavelength", type=float, default=None,
                        help="X-ray wavelength in Angstroms (e.g. 0.12587 for APS 98.5 keV). "
                             "Omit for Cu Kα data (1.54056 Å). Converts synchrotron 2θ to Cu Kα for JARVIS.")
    args = parser.parse_args()

    api_key = os.environ.get("AGAPI_KEY", "")

    # ── STEP 1: get 2θ, intensity ────────────────────────────────────────────
    if args.mode == "demo":
        print("=== Running DEMO with synthetic Fe-Mn-Al-Si-Ni-C XRD peaks ===\n")
        peak_2theta = DEMO_PEAKS[:, 0]
        peak_intensity = DEMO_PEAKS[:, 1]
        args.formula = "Fe2MnAlSiNiC"

    elif args.mode == "image":
        assert args.input, "Provide path to XRD image"
        print(f"Loading image: {args.input}")
        two_theta, intensity = load_xrd_from_image(args.input, row_frac=args.row)
        peak_2theta, peak_intensity = detect_peaks(two_theta, intensity, args.prominence)
        plot_xrd(two_theta, intensity, peak_2theta, peak_intensity,
                 title=f"XRD — {Path(args.input).name}")

    elif args.mode == "manual":
        assert args.input, "Provide path to peaks CSV or .chi file (columns: 2theta intensity)"
        print(f"Loading peaks file: {args.input}")
        two_theta_raw, intensity_raw = load_peaks_from_file(args.input)

        # Peak detection on full profile
        peak_2theta_raw, peak_intensity_raw = detect_peaks(two_theta_raw, intensity_raw, args.prominence)

        # Convert synchrotron 2th -> Cu Ka 2th if wavelength provided
        if args.wavelength and args.wavelength != CUKA:
            print(f"Converting synchrotron 2th (lam={args.wavelength} A) -> Cu Ka 2th (lam={CUKA} A)...")
            peak_2theta_cuka = to_cuka(peak_2theta_raw, args.wavelength)
            valid = ~np.isnan(peak_2theta_cuka) & (peak_2theta_cuka <= 90.0)
            peak_2theta = peak_2theta_cuka[valid]
            peak_intensity = peak_intensity_raw[valid]
            print(f"  Synchrotron range: {peak_2theta_raw.min():.3f}-{peak_2theta_raw.max():.3f} deg")
            print(f"  Cu Ka equivalent:  {peak_2theta.min():.2f}-{peak_2theta.max():.2f} deg")
            # Plot synchrotron profile with Cu Kα-equivalent peak positions and phase labels
            plot_xrd(two_theta_raw, intensity_raw, peak_2theta_raw, peak_intensity_raw,
                     title=f"XRD — {Path(args.input).name} (synch. 2θ, peaks labeled in Cu Kα equiv.)",
                     cuka_peaks=peak_2theta)
        else:
            peak_2theta = peak_2theta_raw
            peak_intensity = peak_intensity_raw
            plot_xrd(two_theta_raw, intensity_raw, peak_2theta_raw, peak_intensity_raw,
                     title=f"XRD — {Path(args.input).name}",
                     cuka_peaks=peak_2theta_raw)

    print(f"\nDetected {len(peak_2theta)} peaks for formula {args.formula}:")
    for t, i in zip(peak_2theta, peak_intensity):
        print(f"  2th={t:.2f} deg  I={i:.1f}%")

    # ── STEP 2: format peaks ─────────────────────────────────────────────────
    peaks_text = format_for_diffractgpt(args.formula, peak_2theta, peak_intensity)
    print(f"\nDiffractGPT input text:\n{'-'*40}\n{peaks_text}\n{'-'*40}")

    # ── STEP 3: DiffractGPT (cloud) or JARVIS match (local) ─────────────────
    if api_key:
        print(f"\nRunning DiffractGPT via AtomGPT.org API...")
        poscar = run_diffractgpt_api(args.formula, peaks_text, api_key)
        if poscar:
            out_path = "diffractgpt_predicted.POSCAR"
            Path(out_path).write_text(poscar)
            print(f"Predicted structure saved → {out_path}")
            print(poscar[:500])
        else:
            print("DiffractGPT returned no structure. Falling back to JARVIS matching.")
            api_key = ""

    if not api_key:
        print("\nNo AGAPI key — running local JARVIS XRD database matching instead.")
        print("(Get a free API key at https://atomgpt.org to enable DiffractGPT structure prediction)\n")
        matches = run_xrd_match_jarvis(peak_2theta, peak_intensity)
        if matches:
            print(f"Top JARVIS matches (cosine similarity):")
            print(f"{'JID':<15} {'Formula':<18} {'SG':<12} {'Score':>6}  {'Ef (eV/at)':>12}")
            print("-" * 65)
            for m in matches:
                ef = f"{m['formation_energy']:.4f}" if m['formation_energy'] else "N/A"
                print(f"{m['jid']:<15} {m['formula']:<18} {m['spg']:<12} {m['score']:>6.4f}  {ef:>12}")
        else:
            print("No matches found.")

    print("\nDone.")


if __name__ == "__main__":
    main()
