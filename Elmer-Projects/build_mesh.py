"""Build (or rebuild) a mesh from the `meshes` registry in elmer_project.json.

    python build_mesh.py <mesh_name>

Pipeline: generate_project_geometry.py (gmsh, requires the external
Thermal-and-Electoric-Sim repository) -> ElmerGrid -> <dir>/ plus a
PROVENANCE.json recording the recipe, input hashes and tool outputs.

Recipe `mesh_overrides` are applied on top of the project's `mesh` section,
and recipe `parameter_overrides` are applied on top of the project's
`parameter_expressions` section, via a resolved copy in generated/ - the
source elmer_project.json is never modified.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
MESH_ROOT = ROOT / "work" / "meshes"
ELMERGRID = r"C:\Program Files\Elmer 26.1-Release\bin\ElmerGrid.exe"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_provenance(mesh_name: str, entry: dict, project_json: Path, verified: bool) -> None:
    mesh_dir = MESH_ROOT / entry["dir"]
    mesh_files = {
        p.name: sha256(p)
        for p in sorted(mesh_dir.glob("mesh.*"))
        if p.is_file()
    }
    provenance = {
        "mesh": mesh_name,
        "written": datetime.now().isoformat(timespec="seconds"),
        "recipe": entry.get("recipe", {}),
        "notes": entry.get("notes", ""),
        "recipe_verified": verified,
        "project_json_sha256": sha256(project_json),
        "mesh_files_sha256": mesh_files,
    }
    (mesh_dir / "PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    print(f"wrote {entry['dir']}/PROVENANCE.json (recipe_verified={verified})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh")
    parser.add_argument(
        "--record-only",
        action="store_true",
        help="only write PROVENANCE.json for the existing directory (recipe_verified=false)",
    )
    parser.add_argument(
        "--project",
        default="elmer_project.json",
        help="path to the project JSON (default: elmer_project.json). Use an "
        "alternate file to try different parameters without touching the "
        "main config; give its meshes/cases distinct names so outputs don't "
        "collide with the main project's generated/ and mesh directories.",
    )
    args = parser.parse_args()

    project_path = Path(args.project)
    project_json = project_path if project_path.is_absolute() else ROOT / project_path
    sys.path.insert(0, str(ROOT))
    from scripts.support.reconcile_project import reconcile_project

    raw_project = json.loads(project_json.read_text(encoding="utf-8"))
    model = reconcile_project(raw_project)
    entry = model.get("meshes", {}).get(args.mesh)
    if entry is None:
        raise SystemExit(f"mesh '{args.mesh}' is not in the meshes registry")

    if args.record_only:
        write_provenance(args.mesh, entry, project_json, verified=False)
        return 0

    geometry_name = entry.get("geometry")
    if not geometry_name:
        raise SystemExit(f"mesh '{args.mesh}' has no 'geometry' key")
    geometries = model.get("geometries", {})
    if geometry_name not in geometries:
        raise SystemExit(f"geometry '{geometry_name}' is not in the geometries registry")

    recipe = entry["recipe"]
    overrides = recipe.get("mesh_overrides", {})
    parameter_overrides = recipe.get("parameter_overrides", {})

    # The `geometry` tree injected below must be resolved against the
    # *overridden* parameter_expressions (dimensions like dabs_dx/tes_pitch
    # can differ per mesh). `model`/`geometries` above were reconciled with
    # the project's global parameter_expressions only, so when overrides are
    # present, re-reconcile a second time with them merged in and pull the
    # geometry tree from that instead. (generate_project_geometry.py's own
    # internal reconcile_project() call only re-derives the `geometries`
    # *registry*, not the legacy singular `geometry` key it actually reads -
    # so the tree written into `resolved["geometry"]` must already be
    # correct on disk.)
    if parameter_overrides:
        override_raw = dict(raw_project)
        override_raw["parameter_expressions"] = {
            **raw_project.get("parameter_expressions", {}),
            **parameter_overrides,
        }
        geometries_for_geometry_key = reconcile_project(override_raw).get("geometries", {})
    else:
        geometries_for_geometry_key = geometries

    resolved = dict(model)
    resolved["mesh"] = {**model.get("mesh", {}), **overrides}
    resolved["parameter_expressions"] = {
        **model.get("parameter_expressions", {}),
        **parameter_overrides,
    }
    # generate_project_geometry.py and the vendored loader still read a single
    # top-level `geometry` tree; inject the selected registry entry under that
    # legacy key so they need no changes.
    resolved["geometry"] = geometries_for_geometry_key[geometry_name]
    resolved_path = ROOT / "generated" / "_mesh_build_input.json"
    resolved_path.parent.mkdir(exist_ok=True)
    resolved_path.write_text(json.dumps(resolved, indent=2) + "\n", encoding="utf-8")

    print(f"[1/2] geometry + gmsh mesh via {recipe['generator']} (external repo required)")
    subprocess.run(
        [sys.executable, str(ROOT / recipe["generator"]), str(resolved_path)],
        cwd=ROOT,
        check=True,
    )

    elmergrid_args = list(recipe["elmergrid_args"])
    if "-out" in elmergrid_args:
        out_index = elmergrid_args.index("-out") + 1
        elmergrid_args[out_index] = str(MESH_ROOT / elmergrid_args[out_index])
    MESH_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"[2/2] ElmerGrid {' '.join(elmergrid_args)}")
    subprocess.run([ELMERGRID, *elmergrid_args], cwd=ROOT, check=True)

    write_provenance(args.mesh, entry, project_json, verified=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
