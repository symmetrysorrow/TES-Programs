#!/usr/bin/env python3
"""Build a TES-frozen, Stycast-side-only interface refinement mesh."""

from __future__ import annotations

import copy
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts/phase24_gate4_5_nomortar/restart_refinement/mortar/project.json"
OUT = ROOT / "artifacts/phase24_stycast_only_refinement"
PROJECT = OUT / "project.json"
MESH_NAME = "mesh_phase24_stycast_only_refined"
MESH_DIR = MESH_NAME
REFINE_H = 5.0e-6


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    entry = copy.deepcopy(project["meshes"]["mesh_singlepixel_gpu_fine_stycast32_mortar"])
    entry["dir"] = MESH_DIR
    entry["recipe"]["elmergrid_args"][-1] = MESH_DIR
    entry["notes"] = (
        "Diagnostic TES-frozen control: STYCAST_INTERFACE_REFINE_H=5 um "
        "refines only the first Stycast layer at the TES contact; TES mesh "
        "and outer Stycast-substrate interface remain unrefined."
    )
    project["meshes"] = {MESH_NAME: entry}
    project["cases"] = {}
    PROJECT.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    env = os.environ.copy()
    env["STYCAST_INTERFACE_REFINE_H"] = str(REFINE_H)
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "build_mesh.py"),
            MESH_NAME,
            "--project",
            str(PROJECT),
        ],
        cwd=ROOT,
        env=env,
        check=True,
    )

    control = {
        "mesh": MESH_NAME,
        "mesh_dir": str(ROOT / "work/meshes" / MESH_DIR),
        "source_project": str(SOURCE),
        "environment": {"STYCAST_INTERFACE_REFINE_H": REFINE_H},
        "scope": "first Stycast layer only; TES-side mesh and outer Stycast-substrate interface are not locally refined",
    }
    (ROOT / "work/meshes" / MESH_DIR / "CONTROL_PROVENANCE.json").write_text(
        json.dumps(control, indent=2) + "\n", encoding="utf-8"
    )
    print(PROJECT)
    print(ROOT / "work/meshes" / MESH_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
