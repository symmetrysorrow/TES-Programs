"""Measure global steady heat balance from Elmer SaveScalars reactions.

This is deliberately based on HeatSolve's ``Calculate Loads`` plus the
SaveScalars ``boundary sum`` operator.  Elemental tetrahedron gradients are
not used for the acceptance decision.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


_FLOAT = re.compile(r"^[\s]*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)(?:\s+|$)")


def read_save_scalars(path: Path) -> tuple[list[str], list[float]]:
    names: list[str] = []
    values: list[float] = []
    in_names = False
    names_path = path.with_name(path.name + ".names")
    for raw in names_path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith("Variables in columns"):
            in_names = True
            continue
        if in_names:
            match = re.match(r"(\d+):\s*(.*)$", line)
            if match:
                names.append(match.group(2).strip())
                continue
            if line and not line.startswith("Variables"):
                in_names = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        tokens = raw.split()
        if len(tokens) == len(names):
            try:
                values = [float(token) for token in tokens]
            except ValueError:
                continue
    if not names or not values:
        raise ValueError(f"Could not parse SaveScalars data: {path}")
    if len(names) != len(values):
        raise ValueError(f"SaveScalars columns/data mismatch: {path}")
    return names, values


def measure(
    path: Path,
    *,
    hot_mask: str = "reaction_hot",
    bath_mask: str = "reaction_bath",
    hot_boundary_id: str = "2",
    bath_boundary_id: str = "1",
    delta_t: float = 0.01,
    tolerance: float = 1.0e-4,
) -> dict:
    names, values = read_save_scalars(path)

    def find(mask: str) -> tuple[str, float]:
        needle = f"over bc {mask}".lower()
        for name, value in zip(names, values):
            if needle in name.lower():
                return name, value
        raise ValueError(f"SaveScalars mask not found: {mask}")

    def find_operator(operator: str, mask: str) -> tuple[str, float] | None:
        needle = f"{operator}:".lower()
        mask_needle = f"over bc {mask}".lower()
        for name, value in zip(names, values):
            if name.lower().startswith(needle) and mask_needle in name.lower():
                return name, value
        return None

    hot_name, hot = find(hot_mask)
    bath_name, bath = find(bath_mask)
    hot_flux = find_operator("diffusive flux", hot_boundary_id)
    bath_flux = find_operator("diffusive flux", bath_boundary_id)
    scale = max(abs(hot), abs(bath), 1.0e-300)
    residual = hot + bath
    relative_residual = abs(residual) / scale
    return {
        "source": "Elmer HeatSolve Calculate Loads + SaveScalars boundary sum",
        "save_scalars_file": str(path),
        "hot_reaction_W": hot,
        "bath_reaction_W": bath,
        "hot_column": hot_name,
        "bath_column": bath_name,
        "global_energy_balance_residual_W": residual,
        "global_energy_balance_relative_residual": relative_residual,
        "delta_T_K": delta_t,
        "effective_conductance_W_per_K": abs(hot) / delta_t,
        "hot_diffusive_flux_W": hot_flux[1] if hot_flux else None,
        "bath_diffusive_flux_W": bath_flux[1] if bath_flux else None,
        "diffusive_flux_columns": {
            "hot": hot_flux[0] if hot_flux else None,
            "bath": bath_flux[0] if bath_flux else None,
        },
        "tolerance_relative": tolerance,
        "status": "PASS" if relative_residual <= tolerance else "FAIL",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save_scalars", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--hot-mask", default="reaction_hot")
    parser.add_argument("--bath-mask", default="reaction_bath")
    parser.add_argument("--hot-boundary-id", default="2")
    parser.add_argument("--bath-boundary-id", default="1")
    parser.add_argument("--delta-t", type=float, default=0.01)
    parser.add_argument("--tolerance", type=float, default=1.0e-4)
    args = parser.parse_args()
    record = measure(
        args.save_scalars,
        hot_mask=args.hot_mask,
        bath_mask=args.bath_mask,
        hot_boundary_id=args.hot_boundary_id,
        bath_boundary_id=args.bath_boundary_id,
        delta_t=args.delta_t,
        tolerance=args.tolerance,
    )
    text = json.dumps(record, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0 if record["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
