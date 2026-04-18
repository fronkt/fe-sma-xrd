"""
Option 3: GNoME-compatible Local Structure Prediction via CHGNet
Uses CHGNet (PyTorch-based, Windows-safe) as the local GNN model for energy
and stability prediction. CHGNet is trained on the same MP/GNoME DFT dataset.

Usage:
    python gnome_predict.py                         # runs on Fe2MnAl POSCAR
    python gnome_predict.py path/to/structure.cif   # custom structure
    python gnome_predict.py --batch data/gnome/fe_sma_gnome_hits.csv  # batch from GNoME screen
"""

import sys
from pathlib import Path


def _load_structure(structure_path: str):
    """Load structure, stripping any prepended non-POSCAR lines."""
    from pymatgen.core import Structure
    import tempfile, os
    path = Path(structure_path)
    if path.suffix.upper() in (".POSCAR", "") or "POSCAR" in path.name:
        lines = path.read_text(errors="replace").splitlines()
        # Find POSCAR start: first line that is a comment/title (not float pair)
        start = 0
        for i, line in enumerate(lines):
            parts = line.split()
            if len(parts) >= 1:
                try:
                    float(parts[0])
                    if len(parts) == 2:
                        continue  # looks like XRD data
                except ValueError:
                    start = i
                    break
        if start > 0:
            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".POSCAR", delete=False)
            tmp.write("\n".join(lines[start:]))
            tmp.close()
            s = Structure.from_file(tmp.name)
            os.unlink(tmp.name)
            return s
    return Structure.from_file(structure_path)


def predict_structure(structure_path: str) -> dict:
    """Run CHGNet single-point prediction on a structure file."""
    from chgnet.model import CHGNet

    print(f"  Loading structure: {structure_path}")
    structure = _load_structure(structure_path)
    print(f"  Formula: {structure.composition.reduced_formula}  |  Sites: {len(structure)}")

    print("  Loading CHGNet model...")
    model = CHGNet.load()

    print("  Running prediction...")
    prediction = model.predict_structure(structure)

    result = {
        "formula": structure.composition.reduced_formula,
        "energy_per_atom_eV": float(prediction["e"]),
        "forces_max_eV_A": float(abs(prediction["f"]).max()) if prediction.get("f") is not None else None,
        "stress_GPa": prediction.get("s"),
        "magmom_per_site": prediction.get("m"),
    }

    return result


def relax_structure(structure_path: str, steps: int = 100) -> dict:
    """Relax a structure using CHGNet as the force field."""
    from chgnet.model import CHGNet
    from chgnet.model.dynamics import StructOptimizer

    print(f"  Relaxing structure: {structure_path} (max {steps} steps)")
    structure = _load_structure(structure_path)
    model = CHGNet.load()
    relaxer = StructOptimizer(model=model)

    result = relaxer.relax(structure, steps=steps, fmax=0.05, verbose=True)
    final = result["final_structure"]
    traj = result["trajectory"]

    return {
        "formula": structure.composition.reduced_formula,
        "initial_energy_eV_atom": float(traj.energies[0]) / len(structure),
        "final_energy_eV_atom": float(traj.energies[-1]) / len(final),
        "steps_taken": len(traj.energies),
        "converged": len(traj.energies) < steps,
        "final_structure": final,
    }


def compare_with_jarvis(formula: str, predicted_energy: float) -> None:
    """Cross-validate CHGNet prediction against JARVIS-DFT if available."""
    try:
        from jarvis.db.figshare import data as jdata
        print("\n  Cross-validating with JARVIS-DFT...")
        dft_3d = jdata("dft_3d")
        hits = [e for e in dft_3d if e.get("formula") == formula and e.get("formation_energy_peratom") is not None]
        if hits:
            jarvis_e = hits[0]["formation_energy_peratom"]
            delta = predicted_energy - jarvis_e
            print(f"  JARVIS formation E: {jarvis_e:.4f} eV/atom")
            print(f"  CHGNet predicted E: {predicted_energy:.4f} eV/atom")
            print(f"  Δ = {delta:+.4f} eV/atom")
        else:
            print(f"  No JARVIS-DFT entry found for {formula}")
    except Exception as e:
        print(f"  JARVIS comparison skipped: {e}")


def batch_predict_from_csv(csv_path: str, max_entries: int = 20) -> None:
    """Predict energies for top entries from the GNoME screen CSV."""
    import pandas as pd
    from chgnet.model import CHGNet
    from pymatgen.core import Structure

    df = pd.read_csv(csv_path)
    model = CHGNet.load()
    results = []

    for _, row in df.head(max_entries).iterrows():
        formula = row.get("Composition", row.get("formula", "?"))
        print(f"  {formula}...", end=" ", flush=True)
        try:
            # GNoME CSV doesn't include structure coords — use pymatgen to build
            from pymatgen.core import Composition
            comp = Composition(formula)
            print(f"(formula only, no structure — skipping relaxation)")
            results.append({"formula": formula, "note": "no structure in CSV"})
        except Exception as e:
            results.append({"formula": formula, "error": str(e)})

    import json
    out = Path(csv_path).parent / "chgnet_batch_results.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nBatch results -> {out}")


if __name__ == "__main__":
    args = sys.argv[1:]

    if "--batch" in args:
        csv_idx = args.index("--batch") + 1
        csv_path = args[csv_idx] if csv_idx < len(args) else "data/gnome/fe_sma_gnome_hits.csv"
        batch_predict_from_csv(csv_path)
        sys.exit(0)

    # Default: predict on Fe2MnAl POSCAR from DiffractGPT
    default_poscar = Path(__file__).parent / "diffractgpt_predicted.POSCAR"
    structure_path = args[0] if args else str(default_poscar)

    if not Path(structure_path).exists():
        print(f"Structure file not found: {structure_path}")
        print("Run fe_sma_diffractgpt.py demo first to generate diffractgpt_predicted.POSCAR")
        sys.exit(1)

    print("=== CHGNet Structure Prediction (GNoME-compatible) ===\n")

    result = predict_structure(structure_path)
    print(f"\nResults:")
    print(f"  Formula:              {result['formula']}")
    print(f"  Energy/atom:          {result['energy_per_atom_eV']:.4f} eV/atom")
    if result["forces_max_eV_A"] is not None:
        print(f"  Max force:            {result['forces_max_eV_A']:.4f} eV/Å")

    compare_with_jarvis(result["formula"], result["energy_per_atom_eV"])

    print("\nTo relax the structure, run:")
    print(f"  >>> from gnome_predict import relax_structure")
    print(f"  >>> relax_structure('{structure_path}')")
