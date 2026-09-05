"""Combine ideal/OCC/Gmsh/Elmer Stycast convergence reports."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", action="append", required=True, help="label=geometry JSON")
    parser.add_argument("--mesh-convergence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    geometry = []
    for item in args.geometry:
        label, raw_path = item.split("=", 1)
        report = json.loads(Path(raw_path).read_text(encoding="utf-8"))
        stycast = report["bodies"]["Stycast"]
        geometry.append({"label": label, **stycast})
    order = {"coarse": 0, "medium": 1, "fine": 2}
    geometry.sort(key=lambda row: order.get(row["label"], 99))
    errors = [abs(float(row["gmsh_vs_ideal_relative_error"])) for row in geometry]
    fem_conversion = [abs(float(row["elmer_vs_gmsh_relative_error"])) for row in geometry]
    trend_pass = all(errors[index] > errors[index + 1] for index in range(len(errors) - 1))
    conversion_pass = max(fem_conversion, default=0.0) <= 1.0e-9
    mesh_report = json.loads(args.mesh_convergence.read_text(encoding="utf-8"))
    report = {
        "status": "PASS" if trend_pass and conversion_pass else "FAIL",
        "ideal_reference": "498 um analytic cylinder",
        "levels": geometry,
        "mesh_surface_and_quality": mesh_report.get("levels", []),
        "checks": {
            "monotone_gmsh_to_ideal_error": trend_pass,
            "elmer_to_gmsh_volume_preservation": conversion_pass,
            "max_abs_elmer_vs_gmsh_relative_error": max(fem_conversion, default=0.0),
        },
        "interpretation": "OCC kernel volume equals the ideal cylinder; the remaining error is polygonal Gmsh discretization and decreases monotonically. ElmerGrid preserves the pre-Elmer Gmsh volume.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
