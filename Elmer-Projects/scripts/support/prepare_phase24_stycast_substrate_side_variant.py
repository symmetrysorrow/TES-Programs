#!/usr/bin/env python3
"""Build the paired TES-side/substrate-side Stycast diagnostic mesh.

The parent is the existing Phase24 density-probe project.  Both opt-in
refinement hooks are enabled: the first Stycast layer (TES side) and the
last Stycast layer (substrate side).  No production project or mesh is
modified intentionally; the generated mesh has its own directory and
provenance record.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/phase24_stycast_interface_discretization_rootcause/mesh_phase24_stycast_density_10um.project.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h-um", type=float, default=10.0)
    parser.add_argument("--artifact-dir", type=Path, default=ROOT / "artifacts/phase24_stycast_substrate_side_rootcause")
    parser.add_argument("--mesh-name", default="mesh_phase24_stycast_density_10um_plus_substrate_10um")
    args = parser.parse_args()
    if args.h_um <= 0.0:
        raise SystemExit("--h-um must be positive")

    source_project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_name = next(iter(source_project["meshes"]))
    entry = copy.deepcopy(source_project["meshes"][source_name])
    entry["dir"] = args.mesh_name
    entry["recipe"]["elmergrid_args"][-1] = args.mesh_name
    entry["notes"] = (
        "Diagnostic paired refinement: TES-side first Stycast layer and "
        "Stycast/substrate final layer at the same target h; substrate mesh "
        "and all physical parameters remain inherited."
    )
    project = copy.deepcopy(source_project)
    project["meshes"] = {args.mesh_name: entry}
    project["cases"] = {}
    args.artifact_dir.mkdir(parents=True, exist_ok=True)
    project_path = args.artifact_dir / f"{args.mesh_name}.project.json"
    project_path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    env = os.environ.copy()
    env["STYCAST_INTERFACE_REFINE_H"] = str(args.h_um * 1.0e-6)
    env["STYCAST_SUBSTRATE_INTERFACE_REFINE_H"] = str(args.h_um * 1.0e-6)
    subprocess.run(
        [sys.executable, str(ROOT / "build_mesh.py"), args.mesh_name, "--project", str(project_path)],
        cwd=ROOT,
        env=env,
        check=True,
    )

    mesh_path = ROOT / "work" / "meshes" / args.mesh_name
    provenance = {
        "case": "phase24_tes_plus_substrate_side_refined",
        "parent_mesh": str(ROOT / "work/meshes/mesh_phase24_stycast_density_10um"),
        "parent_project": str(SOURCE),
        "changed_entities": [
            "Stycast first layer at TES contact: target h=10 um",
            "Stycast final layer at SiO2_2 substrate contact: target h=10 um",
        ],
        "unchanged_entities": [
            "TES mesh",
            "substrate geometry and mesh recipe",
            "Stycast total thickness and geometry",
            "materials, TES law, circuit, bath BC, geometry dimensions",
        ],
        "environment": {
            "STYCAST_INTERFACE_REFINE_H": args.h_um * 1.0e-6,
            "STYCAST_SUBSTRATE_INTERFACE_REFINE_H": args.h_um * 1.0e-6,
        },
        "scope": "diagnostic only; production mesh is not overwritten",
    }
    (mesh_path / "CONTROL_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    manifest = {
        "case": "phase24_tes_plus_substrate_side_refined",
        "target_h_um": args.h_um,
        "mesh": str(mesh_path),
        "project": str(project_path),
        "parent_mesh": provenance["parent_mesh"],
        "changed_entities": provenance["changed_entities"],
        "unchanged_entities": provenance["unchanged_entities"],
    }
    (args.artifact_dir / "mesh_variant_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
