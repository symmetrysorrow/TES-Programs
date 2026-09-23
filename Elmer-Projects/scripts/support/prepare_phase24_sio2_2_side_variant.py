#!/usr/bin/env python3
"""Build a diagnostic variant with absorber target-side refinement.

The TES-side Stycast h=10 um refinement and the Stycast final-layer h=10 um
refinement remain enabled.  This adds only a shallow contact-footprint
refinement box inside the actual target-side ``abs`` body at the Stycast
contact; no production mesh or physics is changed.
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
SOURCE = ROOT / (
    "artifacts/phase24_stycast_substrate_side_rootcause/"
    "mesh_phase24_stycast_density_10um_plus_substrate_10um.project.json"
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h-um", type=float, default=10.0)
    parser.add_argument(
        "--artifact-dir",
        type=Path,
        default=ROOT / "artifacts/phase24_residual_candidates_substrate_target_side",
    )
    parser.add_argument(
        "--mesh-name",
        default="mesh_phase24_stycast_density_10um_plus_substrate_target_10um",
    )
    args = parser.parse_args()
    if args.h_um <= 0.0:
        raise SystemExit("--h-um must be positive")

    source_project = json.loads(SOURCE.read_text(encoding="utf-8"))
    source_name = next(iter(source_project["meshes"]))
    entry = copy.deepcopy(source_project["meshes"][source_name])
    entry["dir"] = args.mesh_name
    entry["recipe"]["elmergrid_args"][-1] = args.mesh_name
    entry["notes"] = (
        "Diagnostic target-side refinement: TES-side Stycast first layer and "
        "Stycast final layer remain h=10 um; the absorber target side at the "
        "Stycast contact is refined only over the contact footprint."
    )
    project = copy.deepcopy(source_project)
    project["meshes"] = {args.mesh_name: entry}
    project["cases"] = {}
    args.artifact_dir.mkdir(parents=True, exist_ok=True)
    project_path = args.artifact_dir / f"{args.mesh_name}.project.json"
    project_path.write_text(
        json.dumps(project, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    env = os.environ.copy()
    env["STYCAST_INTERFACE_REFINE_H"] = str(args.h_um * 1.0e-6)
    env["STYCAST_SUBSTRATE_INTERFACE_REFINE_H"] = str(args.h_um * 1.0e-6)
    env["SUBSTRATE_CONTACT_TARGET_REFINE_H"] = str(args.h_um * 1.0e-6)
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "build_mesh.py"),
            args.mesh_name,
            "--project",
            str(project_path),
        ],
        cwd=ROOT,
        env=env,
        check=True,
    )

    mesh_path = ROOT / "work" / "meshes" / args.mesh_name
    provenance = {
        "case": "phase24_tes_plus_stycast_plus_substrate_target_refined",
        "target_h_um": args.h_um,
        "mesh": str(mesh_path),
        "parent_mesh": str(
            ROOT / "work/meshes/mesh_phase24_stycast_density_10um_plus_substrate_10um"
        ),
        "parent_project": str(SOURCE),
        "changed_entities": [
            "abs target body over the Stycast contact footprint: target h=10 um",
        ],
        "fixed_diagnostic_refinements": [
            "Stycast first layer at TES contact: h=10 um",
            "Stycast final layer at SiO2_2 contact: h=10 um",
        ],
        "unchanged_entities": [
            "TES mesh",
            "Stycast total thickness and geometry",
            "substrate geometry and material outside the contact refinement box",
            "outer substrate geometry outside the contact footprint",
            "materials, TES law, circuit, bath BC, geometry dimensions",
        ],
        "environment": {
            "STYCAST_INTERFACE_REFINE_H": args.h_um * 1.0e-6,
            "STYCAST_SUBSTRATE_INTERFACE_REFINE_H": args.h_um * 1.0e-6,
            "SUBSTRATE_CONTACT_TARGET_REFINE_H": args.h_um * 1.0e-6,
        },
        "scope": "diagnostic only; production mesh is not overwritten",
    }
    (mesh_path / "CONTROL_PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    manifest = {
        "case": provenance["case"],
        "target_h_um": args.h_um,
        "mesh": str(mesh_path),
        "project": str(project_path),
        "parent_mesh": provenance["parent_mesh"],
        "changed_entities": provenance["changed_entities"],
        "fixed_diagnostic_refinements": provenance["fixed_diagnostic_refinements"],
        "unchanged_entities": provenance["unchanged_entities"],
    }
    (args.artifact_dir / "mesh_variant_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
