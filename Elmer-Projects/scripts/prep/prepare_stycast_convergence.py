"""Create a conformal-only Stycast mesh-convergence project."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_physical_parity_conformal.json"
OUTPUT = ROOT / "elmer_project_stycast_convergence.json"


def main() -> int:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    base = copy.deepcopy(project["meshes"]["mesh_physical_parity_conformal"])
    levels = {
        "mesh_stycast_convergence_coarse": (45.0e-6, 90.0e-6),
        "mesh_stycast_convergence_medium": (30.0e-6, 60.0e-6),
        "mesh_stycast_convergence_fine": (22.5e-6, 45.0e-6),
    }
    for name, (h_min, h_max) in levels.items():
        entry = copy.deepcopy(base)
        entry["dir"] = name
        entry["recipe"]["mesh_overrides"] = {
            "mesh_min": h_min,
            "mesh_max": h_max,
            "mesh_min_mode": "fixed",
            "mesh_max_mode": "fixed",
        }
        args = list(entry["recipe"]["elmergrid_args"])
        args[args.index("-out") + 1] = name
        entry["recipe"]["elmergrid_args"] = args
        entry["notes"] = (
            "Conformal shared-node Stycast volume convergence level; "
            f"fixed target h_min={h_min:.12g} m, h_max={h_max:.12g} m."
        )
        project["meshes"][name] = entry
    project["cases"] = {}
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
