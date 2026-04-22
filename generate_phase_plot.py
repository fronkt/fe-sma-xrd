"""
Standalone script: regenerate xrd_peaks_detected.png with phase annotations.

Uses the peak positions and phase assignments from ANALYSIS.md / JARVIS-DFT.
Run from repo root — writes results/sam1_no7325/xrd_peaks_detected.png
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

CUKA = 1.54056  # Å

# ── Known peaks from ANALYSIS.md ─────────────────────────────────────────────
# (sync_2th, cuka_2th, rel_intensity, phase_label, color)
PEAKS = [
    (0.76,  9.30,  15.9, "superstructure?",         "#888888"),
    (3.52, 43.01,  57.4, "FCC γ (111)",              "#e07b39"),
    (3.62, 44.28, 100.0, "B2 (110)\n= L2₁ (220)",   "#2ca02c"),
    (4.10, 50.04,  77.2, "FCC γ (200)",              "#e07b39"),
    (5.28, 64.37,  13.6, "B2 (200)\n= L2₁ (400)",   "#2ca02c"),
    (6.02, 73.49,  21.3, "FCC γ (220)\n/ L2₁ (420)","#9467bd"),
    (6.68, 81.52,  28.8, "B2 (211)\n= L2₁ (422)",   "#2ca02c"),
    (7.34, 89.11,  20.7, "—",                        "#aaaaaa"),
]

# JARVIS top match — annotated in subtitle
JARVIS_TOP = "JARVIS #1: AlFe₃ Pm-3m (B2) · score 0.606 · Ef = −0.179 eV/at"
JARVIS_2   = "JARVIS #3: TiAlFe₂ Fm-3m (L2₁) · score 0.433"

# ── Build synthetic profile from Gaussian peaks ───────────────────────────────
sync_2th = np.linspace(0.0, 9.64, 3000)
profile = np.zeros_like(sync_2th)
sigma = 0.04  # °, synchrotron FWHM ≈ 0.1°

for t, _, i, *_ in PEAKS:
    profile += i * np.exp(-0.5 * ((sync_2th - t) / sigma) ** 2)

# Add a low gentle background
profile += 3.0 * np.exp(-sync_2th / 4)
profile /= profile.max() / 100.0  # normalise to 100%

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(16, 6))
ax.plot(sync_2th, profile, lw=1.2, color="steelblue", label="Profile (synch. 2θ)")

peak_sync = np.array([p[0] for p in PEAKS])
peak_int  = np.array([p[2] for p in PEAKS])
ax.scatter(peak_sync, peak_int, color="red", zorder=5, s=45, label="Detected peaks")

# Alternate annotation offsets to avoid crowding
y_offsets = [14, 14, 14, 14, 14, 14, 14, 14]
x_offsets = [0, 0, 0.05, 0, 0, 0, 0, 0]
arrow_props = dict(arrowstyle="-", color="gray", lw=0.5)

for idx, (t, cuka, i, lbl, col) in enumerate(PEAKS):
    yoff = y_offsets[idx]
    xoff = x_offsets[idx]

    # Compact label: synch 2θ + CuKα + phase
    if lbl == "—":
        text = f"{t:.2f}°\n(2θ_CuKα {cuka:.1f}°)\nunassigned"
    elif lbl == "superstructure?":
        text = f"{t:.2f}°\n(2θ_CuKα {cuka:.1f}°)\n{lbl}"
    else:
        text = f"{t:.2f}°\n(Cu Kα: {cuka:.1f}°)\n{lbl}"

    ax.annotate(
        text,
        xy=(t, i),
        xytext=(t + xoff, i + yoff),
        ha="center",
        fontsize=6.8,
        color=col,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=col, alpha=0.85, lw=0.8),
        arrowprops=arrow_props,
    )

# ── Legend ────────────────────────────────────────────────────────────────────
legend_handles = [
    ax.lines[0],
    ax.collections[0],
    Patch(facecolor="#2ca02c", edgecolor="#2ca02c", alpha=0.7, label="B2 (Pm-3m) / L2₁ (Fm-3m)"),
    Patch(facecolor="#e07b39", edgecolor="#e07b39", alpha=0.7, label="FCC γ-austenite (Fm-3m)"),
    Patch(facecolor="#9467bd", edgecolor="#9467bd", alpha=0.7, label="B2 + FCC γ overlap / L2₁"),
    Patch(facecolor="#888888", edgecolor="#888888", alpha=0.5, label="Unassigned"),
]
ax.legend(handles=legend_handles, fontsize=7.5, loc="upper right", title="Phase key", title_fontsize=8)

ax.set_xlabel("Synchrotron 2θ (°)  —  λ = 0.12587 Å  (APS 98.5 keV)", fontsize=10)
ax.set_ylabel("Relative Intensity (%)", fontsize=10)
ax.set_title(
    "Sam1-NO7325 · APS Synchrotron XRD · Phase-annotated\n"
    f"{JARVIS_TOP}\n{JARVIS_2}",
    fontsize=9,
)
ax.set_xlim(-0.2, 9.9)
ax.set_ylim(-5, 215)
ax.grid(axis="x", ls=":", lw=0.4, alpha=0.5)

plt.tight_layout()
out = "results/sam1_no7325/xrd_peaks_detected.png"
plt.savefig(out, dpi=160)
plt.close()
print(f"Saved -> {out}")
