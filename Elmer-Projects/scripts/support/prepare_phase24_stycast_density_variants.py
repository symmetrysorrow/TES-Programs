#!/usr/bin/env python3
"""Build isolated Stycast-contact density probes for the Phase24 audit.

The probes are diagnostic meshes only.  They inherit the already-qualified
Phase24 project, change only the opt-in ``STYCAST_INTERFACE_REFINE_H`` hook,
and write each mesh into a distinct ``work/meshes`` directory.  No production
project or production mesh is edited.
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
SOURCE = ROOT / "artifacts/phase24_gate4_5_nomortar/restart_refinement/mortar/project.json"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h-um", nargs="+", type=float, required=True)
    parser.add_argument("--prefix", default="mesh_phase24_stycast_density")
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=ROOT / "artifacts/phase24_stycast_interface_discretization_rootcause",
    )
    args = parser.parse_args()

    source_project = json.loads(SOURCE.read_text(encoding="utf-8"))
    base_entry = source_project["meshes"]["mesh_singlepixel_gpu_fine_stycast32_mortar"]
    args.artifact_dir.mkdir(parents=True, exist_ok=True)
    manifest: list[dict] = []

    for h_um in args.h_um:
        if h_um <= 0.0:
            raise SystemExit("--h-um values must be positive")
        tag = f"{h_um:g}".replace(".", "p")
        mesh_name = f"{args.prefix}_{tag}um"
        mesh_dir = mesh_name
        project = copy.deepcopy(source_project)
        entry = copy.deepcopy(base_entry)
        entry["dir"] = mesh_dir
        entry["recipe"]["elmergrid_args"][-1] = mesh_dir
        entry["notes"] = (
            f"Diagnostic Stycast-contact density probe: "
            f"STYCAST_INTERFACE_REFINE_H={h_um:g} um; "
            "TES-side mesh and all physical parameters remain inherited."
        )
        project["meshes"] = {mesh_name: entry}
        project["cases"] = {}
        project_path = args.artifact_dir / f"{mesh_name}.project.json"
        project_path.write_text(
            json.dumps(project, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )

        env = os.environ.copy()
        env["STYCAST_INTERFACE_REFINE_H"] = str(h_um * 1.0e-6)
        subprocess.run(
            [
                sys.executable,
                str(ROOT / "build_mesh.py"),
                mesh_name,
                "--project",
                str(project_path),
            ],
            cwd=ROOT,
            env=env,
            check=True,
        )
        mesh_path = ROOT / "work" / "meshes" / mesh_dir
        provenance = {
            "mesh": mesh_name,
            "mesh_dir": str(mesh_path),
            "source_project": str(SOURCE),
            "environment": {"STYCAST_INTERFACE_REFINE_H": h_um * 1.0e-6},
            "scope": "diagnostic only; TES mesh, outer geometry, materials, circuit, and bath are unchanged",
        }
        (mesh_path / "CONTROL_PROVENANCE.json").write_text(
            json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
        )
        manifest.append(
            {
                "case": mesh_name,
                "target_h_um": h_um,
                "mesh": str(mesh_path),
                "project": str(project_path),
            }
        )

    (args.artifact_dir / "density_mesh_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
