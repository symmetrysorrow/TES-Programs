"""Prepare an isolated Phase24 mesh/project for native coupling diagnostics.

The production Phase24 mesh is left untouched.  This diagnostic copy uses the
same geometry and material parameters but explicitly disables the conformal
TES--membrane imprint, so that the two meshes expose the same nonconforming
TES-to-membrane mortar topology as the historical reference.
"""
from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "artifacts" / "phase24_gate4_5_nomortar" / "restart_refinement" / "mortar" / "project.json"
OUT = ROOT / "artifacts" / "phase24_native_coupling_diagnostics"
PROJECT = OUT / "project.json"
MESH_NAME = "mesh_phase24_native_coupling_diag"
MESH_DIR = MESH_NAME


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    entry = copy.deepcopy(project["meshes"]["mesh_singlepixel_gpu_fine_stycast32_mortar"])
    entry["dir"] = MESH_DIR
    entry["recipe"]["elmergrid_args"][-1] = MESH_DIR
    entry["recipe"]["elmer_overrides"] = {
        **entry["recipe"].get("elmer_overrides", {}),
        "fragment_mortar_interfaces": False,
        "conformal_shared_interfaces": False,
        "conformal_mortar_interfaces": False,
    }
    project["meshes"] = {MESH_NAME: entry}
    project["cases"] = {}
    PROJECT.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    subprocess.run(
        [sys.executable, str(ROOT / "build_mesh.py"), MESH_NAME, "--project", str(PROJECT)],
        cwd=ROOT,
        check=True,
    )
    print(PROJECT)
    print(ROOT / "work" / "meshes" / MESH_DIR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
