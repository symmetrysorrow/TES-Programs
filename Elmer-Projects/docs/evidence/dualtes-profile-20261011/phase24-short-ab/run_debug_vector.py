#!/usr/bin/env python3
"""Reproduce Phase24 vector-assembly failures with an isolated Debug runtime."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
from datetime import datetime


HERE = Path(__file__).resolve().parent
EVIDENCE = HERE.parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(ranks: int, root: Path) -> int:
    root = root.resolve()
    label = f"dualtes_vector_debug_{ranks}rank"
    source = HERE / "dualtes_vector.sif"
    sif = source.read_text(encoding="utf-8")
    sif = re.sub(
        r'(TES(?: 2)? Series File = String ")[^"]+',
        lambda m: m.group(1) + f"{label}_{'2' if m.group(1).startswith('TES 2') else '1'}_series.csv",
        sif,
    )
    sif = re.sub(
        r'(TES(?: 2)? Iteration Series File = String ")[^"]+',
        lambda m: m.group(1) + f"{label}_{'2' if m.group(1).startswith('TES 2') else '1'}_iterations.csv",
        sif,
    )
    mesh = re.search(r'Mesh DB\s+"[^"]+"\s+"([^"]+)"', sif).group(1)
    for tes in (1, 2):
        original = HERE / f"phase24_short_ab_dualtes_vector_{tes}.state"
        target = root / "workspace/work/meshes" / mesh / f"{label}_{tes}.state"
        shutil.copy2(original, target)
        sif = re.sub(
            rf'(TES {tes} State File = String ")[^"]+',
            lambda m: m.group(1) + str(Path("workspace/work/meshes") / mesh / target.name),
            sif,
        )
    sif_path = HERE / f"{label}.sif"
    sif_path.write_text(sif, encoding="utf-8")

    temp = Path(os.environ["TEMP"])
    debug_runtime = temp / "phase24-debug-runtime"
    debug_build = temp / "phase24-debug-build"
    solver = debug_runtime / "bin/ElmerSolver_mpi.exe"
    mpi = Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe")
    env = os.environ.copy()
    path_parts = [
        root / "build/udf",
        debug_build / "fem/src/modules",
        debug_build / "fem/src",
        debug_build / "matc/src",
        debug_build / "fhutiter/src",
        debug_build / "fem/src/binio",
        debug_build / "mathlibs/src/arpack",
        debug_build / "mathlibs/src/parpack",
        debug_runtime / "bin",
        debug_runtime / "share/elmersolver/lib",
        Path(r"G:\Elmer-Migration\msys64\ucrt64\bin"),
        Path(r"G:\Elmer-Migration\msys64\usr\bin"),
        Path(r"G:\Elmer-Migration\elmer\bin"),
        mpi.parent,
    ]
    env["ELMER_HOME"] = str(debug_runtime)
    env["PATH"] = os.pathsep.join([str(p) for p in path_parts] + [env.get("PATH", "")])
    env["OMP_NUM_THREADS"] = "1"
    command = [str(mpi), "-n", str(ranks), str(solver), str(sif_path)]
    started = datetime.now().astimezone().isoformat(timespec="seconds")
    log_path = HERE / f"{label}.solver.log"
    with log_path.open("w", encoding="utf-8") as log:
        proc = subprocess.run(command, cwd=root, env=env, stdout=log,
                              stderr=subprocess.STDOUT, check=False)
    output_hashes = {}
    for tes in (1, 2):
        for kind in ("series", "iterations"):
            generated = root / f"{label}_{tes}_{kind}.csv"
            if generated.is_file():
                shutil.copy2(generated, HERE / generated.name)
                output_hashes[generated.name] = sha256(HERE / generated.name)
        state = root / "workspace/work/meshes" / mesh / f"{label}_{tes}.state"
        if state.is_file():
            shutil.copy2(state, HERE / state.name)
            output_hashes[state.name] = sha256(HERE / state.name)
    manifest_path = HERE / "debug_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {"runs": []}
    manifest["runs"] = [r for r in manifest["runs"] if r["name"] != label]
    manifest["runs"].append({
        "name": label,
        "ranks": ranks,
        "started": started,
        "finished": datetime.now().astimezone().isoformat(timespec="seconds"),
        "command": command,
        "exit_code": proc.returncode,
        "source_sif_sha256": sha256(source),
        "debug_sif_sha256": sha256(sif_path),
        "solver_sha256": sha256(solver),
        "heat_solve_debug_sha256": sha256(debug_runtime / "share/elmersolver/lib/HeatSolve.dll"),
        "libelmersolver_debug_sha256": sha256(debug_build / "fem/src/libelmersolver.dll"),
        "initial_state_sha256": {
            f"{label}_{tes}.state": sha256(HERE / f"phase24_short_ab_dualtes_vector_{tes}.state")
            for tes in (1, 2)
        },
        "generated_output_sha256": output_hashes,
        "log": log_path.name,
    })
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return proc.returncode


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranks", type=int, required=True, choices=(1, 6))
    parser.add_argument("--repo-root", type=Path, required=True,
                        help="Runtime workspace root containing workspace/work/meshes")
    args = parser.parse_args()
    raise SystemExit(run(args.ranks, args.repo_root))
