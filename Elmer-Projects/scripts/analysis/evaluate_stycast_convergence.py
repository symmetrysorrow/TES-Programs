"""Measure mesh volume and contact-facet convergence for Stycast."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.evaluate_physical_parity import contact_areas, expected_geometry, mesh_data
def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--mesh", action="append", required=True, help="label=mesh directory")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    expected, expected_contact = expected_geometry(args.project)
    rows = []
    for item in args.mesh:
        label, raw_path = item.split("=", 1)
        mesh = Path(raw_path)
        bodies, _, nodes, _, elements_by_body, _, _ = mesh_data(mesh)
        reverse = {value: key for key, value in bodies.items()}
        stycast_id = next(key for key, value in reverse.items() if value == "Stycast")
        volume = 0.0
        from scripts.analysis.evaluate_physical_parity import tetra_volume
        for conn in elements_by_body[stycast_id]:
            volume += tetra_volume([nodes[node] for node in conn])
        areas = contact_areas(mesh)
        rows.append(
            {
                "label": label,
                "mesh": str(mesh.resolve()),
                "node_count": sum(1 for line in (mesh / "mesh.nodes").read_text(encoding="utf-8").splitlines() if line.strip()),
                "stycast_element_count": len(elements_by_body[stycast_id]),
                "stycast_volume_m3": volume,
                "stycast_analytic_volume_m3": expected["Stycast"],
                "stycast_relative_volume_error": (volume - expected["Stycast"]) / expected["Stycast"],
                "contact_areas": areas,
                "contact_area_relative_errors": {
                    name: {
                        side: (float(values[side]) - expected_contact[name]) / expected_contact[name]
                        for side in ("left_m2", "right_m2")
                    }
                    for name, values in areas.items()
                },
            }
        )
    order = {"coarse": 0, "medium": 1, "fine": 2}
    rows.sort(key=lambda row: order.get(row["label"], 99))
    report = {
        "project": str(args.project.resolve()),
        "definition": "498 um analytic Stycast cylinder; conformal shared-node meshes",
        "levels": rows,
        "status": "REPORT_ONLY",
        "note": "Production acceptance still requires a monotone convergence decision and a production-size run.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
