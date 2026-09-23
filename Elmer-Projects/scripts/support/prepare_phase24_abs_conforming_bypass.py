#!/usr/bin/env python3
"""Build the single-interface Stycast->abs conforming diagnostic mesh.

The source project and all physical parameters come from the validated
10-um TES-side Stycast refinement.  Only the Stycast/abs contact is made
periodic/shared; TES/membrane and TES/Stycast remain mortar interfaces.
"""
from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/phase24_stycast_interface_discretization_rootcause/mesh_phase24_stycast_density_10um.project.json"
OUT = ROOT / "artifacts/phase24_abs_connectivity_mortar_bypass"
MESH_NAME = "mesh_phase24_abs_conforming_bypass_stycast_abs"
MESH_DIR = ROOT / "work/meshes" / MESH_NAME
ELMERGRID = Path(r"C:\Program Files\Elmer 26.1-Release\bin\ElmerGrid.exe")


def main() -> int:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_mesh = source["meshes"]["mesh_phase24_stycast_density_10um"]
    project = copy.deepcopy(source)
    project["elmer_overrides"] = {
        **project.get("elmer_overrides", {}),
        "fragment_mortar_interfaces": True,
        "conformal_shared_interfaces": False,
        "conformal_mortar_interfaces": True,
        "conformal_contact_interfaces": ["stycast_abs"],
    }
    project["mesh"] = {
        **project.get("mesh", {}),
        **source_mesh["recipe"].get("mesh_overrides", {}),
    }
    project["geometry"] = project["geometries"][source_mesh["geometry"]]
    project["notes"] = (
        "Diagnostic only: same 10-um TES-side Stycast refinement and all physical "
        "parameters as the mortar baseline; only Stycast->abs uses shared nodes."
    )
    OUT.mkdir(parents=True, exist_ok=True)
    project_path = OUT / f"{MESH_NAME}.project.json"
    project_path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    sys.path.insert(0, str(ROOT))
    from scripts.support.reconcile_project import reconcile_project

    resolved = reconcile_project(project)
    resolved["geometry"] = resolved["geometries"][source_mesh["geometry"]]
    resolved_path = OUT / f"{MESH_NAME}.resolved.json"
    resolved_path.write_text(json.dumps(resolved, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    env = os.environ.copy()
    env["STYCAST_INTERFACE_REFINE_H"] = "1.0e-5"
    subprocess.run(
        [sys.executable, str(ROOT / "generate_project_geometry.py"), str(resolved_path)],
        cwd=ROOT,
        env=env,
        check=True,
    )
    if not ELMERGRID.is_file():
        raise FileNotFoundError(ELMERGRID)
    if MESH_DIR.exists():
        raise FileExistsError(f"refusing to overwrite existing diagnostic mesh: {MESH_DIR}")
    subprocess.run(
        [
            str(ELMERGRID), "14", "2", "gmsh/project.msh", "-merge", "1e-10",
            "-out", str(MESH_DIR),
        ],
        cwd=ROOT,
        check=True,
    )
    provenance = {
        "mesh": MESH_NAME,
        "source_project": str(SOURCE),
        "source_mesh": "mesh_phase24_stycast_density_10um",
        "interface_change": "Stycast -> abs only: mortar to periodic/shared-node continuity",
        "remaining_mortar_interfaces": ["Membrane_SiNx -> TES", "TES -> Stycast"],
        "fixed": [
            "geometry, materials, TES law, circuit constants, bath, power points",
            "mesh_min=45 um, mesh_max=90 um, Stycast layers=32",
            "STYCAST_INTERFACE_REFINE_H=10 um",
        ],
        "production_mesh_overwritten": False,
        "mesh_dir": str(MESH_DIR),
    }
    (OUT / "mortar_bypass_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
