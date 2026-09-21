"""Prepare an isolated finer no-mortar conformal mesh probe.

The existing conformal fine mesh is retained.  This creates a separate
project JSON and mesh registry entry with a smaller global characteristic
length, so ``build_mesh.py`` can generate a reproducible refinement candidate
without overwriting an existing mesh.

Usage::

    python scripts/analysis/prepare_nomortar_refinement_probe.py
    python build_mesh.py mesh_singlepixel_conformal_gpu_refine20 --project \
        artifacts/comparison/nomortar_refinement_probe/project.json
"""
from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_singlepixel_conformal_gpu.json"
OUTDIR = ROOT / "artifacts/comparison/nomortar_refinement_probe"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh-name", default="mesh_singlepixel_conformal_gpu_refine20")
    parser.add_argument("--mesh-min-um", type=float, default=20.0)
    parser.add_argument("--mesh-max-um", type=float, default=40.0)
    args = parser.parse_args()
    if args.mesh_min_um <= 0 or args.mesh_max_um < args.mesh_min_um:
        parser.error("mesh sizes must be positive and max >= min")

    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    base = copy.deepcopy(project["meshes"]["mesh_singlepixel_conformal_gpu_fine"])
    base["dir"] = args.mesh_name
    base["recipe"]["mesh_overrides"] = {
        "mesh_min": args.mesh_min_um * 1e-6,
        "mesh_max": args.mesh_max_um * 1e-6,
        "mesh_min_mode": "fixed",
        "mesh_max_mode": "fixed",
    }
    base["recipe"]["elmergrid_args"][-1] = args.mesh_name
    base["notes"] = (
        f"Isolated no-mortar conformal refinement probe: global mesh sizes "
        f"{args.mesh_min_um:g}/{args.mesh_max_um:g} um; existing fine mesh is untouched."
    )
    project["meshes"][args.mesh_name] = base
    OUTDIR.mkdir(parents=True, exist_ok=True)
    output = OUTDIR / "project.json"
    output.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"project={output}")
    print(f"mesh={args.mesh_name}")
    print(f"build_command=python build_mesh.py {args.mesh_name} --project {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
