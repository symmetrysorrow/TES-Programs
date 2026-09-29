"""Build one-at-a-time trace-resistance controlled diagnostic meshes."""

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


VARIANTS = {
    "tes_membrane": {
        "mesh_name": "mesh_phase24_trace_tes_membrane_historical",
        "environment": {"TES_MEMBRANE_TRACE_REFINE_H": 10.6e-6},
        "changed": "TES/Membrane direct trace neighborhood only; TES/Stycast top trace remains parent resolution",
    },
    "membrane_substrate": {
        "mesh_name": "mesh_phase24_trace_membrane_substrate_historical",
        "environment": {"MEMBRANE_SUBSTRATE_TRACE_REFINE_H": 10.0e-6},
        "changed": "downstream Membrane-to-substrate trace neighborhood below Membrane_SiNx only; TES/Membrane trace remains parent resolution",
    },
}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--variant", choices=["all", *VARIANTS], default="all")
    parser.add_argument("--artifact-dir", type=Path, default=ROOT / "artifacts/phase24_trace_resistance_controlled")
    args = parser.parse_args()

    source_project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_name = next(iter(source_project["meshes"]))
    args.artifact_dir.mkdir(parents=True, exist_ok=True)
    selected = VARIANTS if args.variant == "all" else {args.variant: VARIANTS[args.variant]}
    manifest = []
    for key, spec in selected.items():
        project = copy.deepcopy(source_project)
        entry = copy.deepcopy(project["meshes"][source_name])
        entry["dir"] = spec["mesh_name"]
        entry["recipe"]["elmergrid_args"][-1] = spec["mesh_name"]
        entry["notes"] = (
            f"Trace-resistance controlled diagnostic: {spec['changed']}. "
            "Materials, geometry dimensions, TES law, circuit, bath, and production mesh are unchanged."
        )
        project["meshes"] = {spec["mesh_name"]: entry}
        project["cases"] = {}
        project_path = args.artifact_dir / f"{spec['mesh_name']}.project.json"
        project_path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

        env = os.environ.copy()
        # Preserve the accepted refined mortar parent: the Stycast-side
        # TES/Stycast trace stays at h=10 um in both one-at-a-time tests.
        env["STYCAST_INTERFACE_REFINE_H"] = "1.0e-5"
        env.update({name: str(value) for name, value in spec["environment"].items()})
        subprocess.run(
            [sys.executable, str(ROOT / "src" / "build_mesh.py"), spec["mesh_name"], "--project", str(project_path)],
            cwd=ROOT,
            env=env,
            check=True,
        )
        mesh_path = ROOT / "workspace/work" / "meshes" / spec["mesh_name"]
        provenance = {
            "case": key,
            "mesh": str(mesh_path),
            "parent_mesh": str(ROOT / "workspace/work/meshes/mesh_phase24_stycast_density_10um"),
            "parent_project": str(SOURCE),
            "environment": spec["environment"],
            "changed": spec["changed"],
            "unchanged": [
                "materials",
                "TES law and circuit",
                "bath temperature/BC",
                "geometry dimensions",
                "TES-Stycast coupling and parent trace resolution",
                "production mesh",
            ],
            "scope": "diagnostic only",
        }
        (mesh_path / "CONTROL_PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
        manifest.append({"variant": key, "mesh": str(mesh_path), "project": str(project_path), **provenance})

    (args.artifact_dir / "variant_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
