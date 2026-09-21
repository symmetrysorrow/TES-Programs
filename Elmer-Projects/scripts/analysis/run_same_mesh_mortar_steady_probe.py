"""Run a same-mesh mortar-vs-no-mortar steady probe.

The probe reuses the validated mortar boundary blocks from the existing
``mesh_refined_3x`` reference SIF, changes only the mesh database to the
fine conformal mesh, and writes a separate result/log under artifacts and
the fine mesh directory.  It is intentionally a MUMPS steady solve so that
linear-iteration error is not part of this isolation test.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "generated/cases/case_tes_steady_singlepixel_conformal_gpu_fine.sif"
OUTDIR = ROOT / "artifacts/comparison/same_mesh_mortar_steady_probe"
MESH_NAME = "mesh_singlepixel_conformal_gpu_fine"
CASE = "case_tes_steady_singlepixel_conformal_same_mesh_mortar_probe"
RESULT = ROOT / "work/meshes" / MESH_NAME / f"{CASE}.result"
DEFAULT_SOLVER = Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver.exe")


def make_sif() -> str:
    text = SOURCE.read_text(encoding="utf-8")
    text = text.replace(
        'Mesh DB "work/meshes" "mesh_refined_3x"',
        f'Mesh DB "work/meshes" "{MESH_NAME}"',
    )
    text = text.replace("Linear System Direct Method = Umfpack", "Linear System Direct Method = MUMPS")
    text = text.replace("Apply Mortar BCs = False", "Apply Mortar BCs = True")
    text = re.sub(
        r'(^\s*TES Series File\s*=\s*String\s+")[^"]+("\s*$)',
        rf'\g<1>artifacts/comparison/same_mesh_mortar_steady_probe/{CASE}_series.csv\g<2>',
        text,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r'(^\s*TES Iteration Series File\s*=\s*String\s+")[^"]+("\s*$)',
        rf'\g<1>artifacts/comparison/same_mesh_mortar_steady_probe/{CASE}_iterations.csv\g<2>',
        text,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r"^\s*Solver Input File\s*=.*$",
        f"  Solver Input File = artifacts/comparison/same_mesh_mortar_steady_probe/{CASE}.sif",
        text,
        flags=re.MULTILINE,
    )
    text = re.sub(
        r"^\s*Output File\s*=.*$",
        f"  Output File = ../work/meshes/{MESH_NAME}/{CASE}.result",
        text,
        count=1,
        flags=re.MULTILINE,
    )
    text = text.replace(
        "case_tes_steady_singlepixel_conformal_mortar_reference",
        CASE,
    )
    # The conformal fine mesh carries the same named boundary IDs as the
    # mortar reference.  Add only the explicit contact blocks; the no-mortar
    # source already contains the bath BC and all body/material definitions.
    mortar_blocks = """
Boundary Condition 2
  Target Boundaries(1) = 1104
  Name = \"TES bottom mortar\"
  Mortar BC = 3
  Galerkin Projector = True
  Plane Projector = True
End

Boundary Condition 3
  Target Boundaries(1) = 1305
  Name = \"Membrane_SiNx top mortar\"
End

Boundary Condition 4
  Target Boundaries(1) = 1204
  Name = \"Stycast bottom mortar\"
  Mortar BC = 5
  Galerkin Projector = True
  Plane Projector = True
End

Boundary Condition 5
  Target Boundaries(1) = 1105
  Name = \"TES top mortar\"
End

Boundary Condition 6
  Target Boundaries(1) = 1205
  Name = \"Stycast top mortar\"
  Mortar BC = 7
  Galerkin Projector = True
  Plane Projector = True
End

Boundary Condition 7
  Target Boundaries(1) = 1004
  Name = \"Pb bottom mortar\"
End
"""
    boundary_end = re.search(r"Boundary Condition 1\b.*?^End\s*$", text, flags=re.MULTILINE | re.DOTALL)
    if boundary_end is None:
        raise ValueError("bath boundary block not found")
    text = text[:boundary_end.end()] + "\n" + mortar_blocks.strip("\n") + text[boundary_end.end():]
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=DEFAULT_SOLVER)
    parser.add_argument("--runtime-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sif_path = OUTDIR / f"{CASE}.sif"
    sif_path.write_text(make_sif(), encoding="utf-8", newline="\n")
    if args.dry_run:
        print(sif_path)
        return 0

    if not args.solver.is_file():
        raise FileNotFoundError(args.solver)
    if not (ROOT / "work/meshes" / MESH_NAME / "mesh.header").is_file():
        raise FileNotFoundError(ROOT / "work/meshes" / MESH_NAME / "mesh.header")
    env = os.environ.copy()
    path_parts = [str(args.runtime_bin), str(args.solver.parent), str(args.solver.parent.parent / "share" / "elmersolver" / "lib")]
    env["PATH"] = os.pathsep.join(path_parts + [env.get("PATH", "")])
    log_path = OUTDIR / "solver.log"
    with log_path.open("w", encoding="utf-8", newline="\n") as log:
        completed = subprocess.run(
            [str(args.solver), str(sif_path)],
            cwd=ROOT,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    (OUTDIR / "result_path.txt").write_text(str(RESULT) + "\n", encoding="utf-8")
    print(f"exit_code={completed.returncode}")
    print(f"sif={sif_path}")
    print(f"log={log_path}")
    print(f"result={RESULT}")
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
