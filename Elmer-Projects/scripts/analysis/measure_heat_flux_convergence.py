"""Measure reconstructed interface fluxes on several result/mesh levels."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.evaluate_physical_parity import (
    PAIR_DEFS,
    body_boundary_flux_balance,
    heat_flux_consistency,
    mesh_data,
)


def parse_level(value: str) -> tuple[str, Path, Path]:
    label, mesh, result = value.split("=", 2)
    return label, Path(mesh), Path(result)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--level", action="append", required=True, help="label=mesh_dir=result_file")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    levels = []
    for raw in args.level:
        label, mesh, result = parse_level(raw)
        bodies, boundaries, nodes, elements, elements_by_body, faces, _ = mesh_data(mesh)
        flux = heat_flux_consistency(mesh, result, args.project)
        interfaces = {}
        for _, _, interface, *_ in PAIR_DEFS:
            entry = flux[interface]
            interfaces[interface] = {
                "left": entry["left"],
                "right": entry["right"],
                "absolute_imbalance_W": entry["absolute_imbalance_W"],
                "normalized_imbalance": entry["normalized_imbalance"],
                "local_flux_jump": entry["local_flux_jump"],
                "surface_mean_temperature_difference_K": entry["surface_mean_temperature_difference_K"],
                "status": entry["status"],
            }
        levels.append(
            {
                "label": label,
                "mesh": str(mesh.resolve()),
                "result": str(result.resolve()),
                "node_count": len(nodes),
                "tetrahedron_count": len(elements),
                "body_count": len(elements_by_body),
                "interface_facet_counts": {
                    interface: {
                        "left": len(faces[left]),
                        "right": len(faces[right]),
                    }
                    for left, right, interface, *_ in PAIR_DEFS
                },
                "interfaces": interfaces,
                "body_boundary_flux_balance": body_boundary_flux_balance(mesh, result, args.project),
            }
        )

    report = {
        "method": "piecewise-linear tetrahedral gradient q=-k grad(T) dot n_outward",
        "result_field_selection": "last Perm block in each Elmer ASCII result",
        "levels": levels,
        "interpretation": {
            "epsilon_Q": "absolute_imbalance_W / max(abs(left_flux_W), abs(right_flux_W), 1e-30)",
            "local_flux_jump": "paired conformal facet integrated flux sum; unavailable for node-disjoint Mortar sides",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
