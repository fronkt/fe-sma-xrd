"""
Option 1: GNoME via OPTIMADE API
Query the community-hosted GNoME OPTIMADE endpoint for structures.
Falls back to the Materials Project OPTIMADE endpoint if GNoME is down.
"""

import httpx
import json
from typing import Optional

GNOME_BASE = "https://optimade-gnome.odbx.science/v1"
MP_BASE = "https://optimade.materialsproject.org/v1"


def query_structures(
    elements: list[str],
    base_url: str = GNOME_BASE,
    max_results: int = 50,
) -> list[dict]:
    """
    Query OPTIMADE for structures containing all specified elements.
    Returns list of structure entries with id, formula, and attributes.
    """
    elements_filter = ",".join(f'"{e}"' for e in elements)
    params = {
        "filter": f"elements HAS ALL {elements_filter}",
        "page_limit": min(max_results, 100),
    }

    with httpx.Client(timeout=30) as client:
        resp = client.get(f"{base_url}/structures", params=params)
        resp.raise_for_status()
        data = resp.json()

    return data.get("data", [])


def query_gnome_fe_sma(elements: Optional[list[str]] = None) -> list[dict]:
    """Query GNoME for Fe-SMA relevant structures (Fe-Mn-Al system by default)."""
    if elements is None:
        elements = ["Fe", "Mn", "Al"]

    print(f"Querying GNoME OPTIMADE for elements: {elements}")
    try:
        results = query_structures(elements, base_url=GNOME_BASE)
        print(f"  GNoME endpoint: {len(results)} structures found")
        return results
    except Exception as e:
        print(f"  GNoME endpoint failed ({e}), falling back to MP OPTIMADE...")
        results = query_structures(elements, base_url=MP_BASE)
        print(f"  MP OPTIMADE fallback: {len(results)} structures found")
        return results


def print_summary(entries: list[dict]) -> None:
    print(f"\n{'ID':<30} {'Formula':<20} {'Sites'}")
    print("-" * 60)
    for e in entries:
        attrs = e.get("attributes", {})
        formula = attrs.get("chemical_formula_reduced") or attrs.get("chemical_formula_hill", "?")
        nsites = attrs.get("nsites", "?")
        print(f"{e.get('id', '?'):<30} {formula:<20} {nsites}")


if __name__ == "__main__":
    import sys
    elements = sys.argv[1:] if len(sys.argv) > 1 else ["Fe", "Mn", "Al"]
    entries = query_gnome_fe_sma(elements)
    print_summary(entries)

    out_path = "gnome_optimade_results.json"
    with open(out_path, "w") as f:
        json.dump(entries, f, indent=2)
    print(f"\nSaved {len(entries)} entries -> {out_path}")
