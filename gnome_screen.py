"""
Option 2: GNoME Raw Data Screener
Downloads stable_materials_summary.csv from DeepMind's GCS bucket and screens
for compositions relevant to the Fe-SMA project (Fe-Mn-Al-Si-Ni-C system).
"""

import os
import sys
import httpx
import pandas as pd
from pathlib import Path

GCS_BASE = "https://storage.googleapis.com/gdm_materials_discovery/gnome_data"
DATA_DIR = Path(__file__).parent / "data" / "gnome"

SUMMARY_CSV = "stable_materials_summary.csv"
R2SCAN_CSV = "stable_materials_r2scan.csv"


def download_file(filename: str, force: bool = False) -> Path:
    dest = DATA_DIR / filename
    if dest.exists() and not force:
        print(f"  Already cached: {dest}")
        return dest

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    url = f"{GCS_BASE}/{filename}"
    print(f"  Downloading {filename} from GCS... (~150 MB, please wait)")

    with httpx.Client(timeout=300, follow_redirects=True) as client:
        with client.stream("GET", url) as resp:
            resp.raise_for_status()
            total = int(resp.headers.get("content-length", 0))
            downloaded = 0
            with open(dest, "wb") as f:
                for chunk in resp.iter_bytes(chunk_size=1024 * 1024):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded / total * 100
                        print(f"\r  {pct:.1f}%  ({downloaded/1e6:.1f} MB / {total/1e6:.1f} MB)", end="", flush=True)
    print(f"\n  Saved -> {dest}")
    return dest


def load_summary(force_download: bool = False) -> pd.DataFrame:
    path = download_file(SUMMARY_CSV, force=force_download)
    print("  Loading CSV into DataFrame...")
    df = pd.read_csv(path, low_memory=False)
    print(f"  Loaded {len(df):,} entries, columns: {list(df.columns)}")
    return df


def screen_fe_sma(df: pd.DataFrame, target_elements: set | None = None) -> pd.DataFrame:
    """
    Filter GNoME dataset for entries whose composition is a subset of the
    Fe-SMA target system: Fe-Mn-Al-Si-Ni-C.
    """
    if target_elements is None:
        target_elements = {"Fe", "Mn", "Al", "Si", "Ni", "C"}

    # GNoME CSV has an 'Elements' column as a space-separated string
    def all_in_system(elements_str: str) -> bool:
        if pd.isna(elements_str):
            return False
        return set(str(elements_str).split()) <= target_elements

    mask = df["Elements"].apply(all_in_system)
    filtered = df[mask].copy()

    # Sort by decomposition energy (stability)
    if "Decomposition Energy (eV/atom)" in filtered.columns:
        filtered = filtered.sort_values("Decomposition Energy (eV/atom)")
    elif "Formation Energy (eV/atom)" in filtered.columns:
        filtered = filtered.sort_values("Formation Energy (eV/atom)")

    return filtered


def summarize(df: pd.DataFrame) -> None:
    print(f"\n{'='*70}")
    print(f"GNoME Fe-SMA screen results: {len(df):,} entries")
    print(f"{'='*70}")

    show_cols = [c for c in [
        "Composition", "Elements", "Decomposition Energy (eV/atom)",
        "Formation Energy (eV/atom)", "Space Group Number", "Volume per Atom (A^3)"
    ] if c in df.columns]

    if not show_cols:
        print(df.head(20).to_string())
    else:
        print(df[show_cols].head(30).to_string(index=False))

    # Binary subsystem counts
    if "Elements" in df.columns:
        n_elements = df["Elements"].apply(lambda x: len(str(x).split()) if pd.notna(x) else 0)
        print(f"\nBy number of elements:")
        print(n_elements.value_counts().sort_index().to_string())


if __name__ == "__main__":
    force = "--force" in sys.argv

    print("=== GNoME Raw Data Screener ===\n")
    df = load_summary(force_download=force)
    results = screen_fe_sma(df)
    summarize(results)

    out_path = DATA_DIR / "fe_sma_gnome_hits.csv"
    results.to_csv(out_path, index=False)
    print(f"\nSaved {len(results):,} hits -> {out_path}")
