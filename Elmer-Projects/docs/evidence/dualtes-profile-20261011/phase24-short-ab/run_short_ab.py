#!/usr/bin/env python3
"""Run bounded Phase24 short A/B diagnostics using the preserved 4-step inputs."""
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
VARIANTS = {
    "dualtes_hypre_report": (EVIDENCE / "dualtes-short" / "case.sif", False, False),
    "dualtes_vector": (EVIDENCE / "dualtes-short" / "case.sif", True, False),
    "dualtes_vector_static": (EVIDENCE / "dualtes-short" / "case.sif", True, True),
    "singlepixel_vector": (EVIDENCE / "singlepixel-short" / "case.sif", True, False),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source: Path, name: str, vector: bool, static_reuse: bool) -> Path:
    sif = source.read_text(encoding="utf-8")
    mesh = re.search(r'Mesh DB\s+"([^"]+)"\s+"([^"]+)"', sif)
    if not mesh:
        raise RuntimeError(f"No mesh database in {source}")
    mesh_name = mesh.group(2)
    restart = re.search(r'Restart File\s*=\s*"?([^"\r\n]+)', sif)
    if not restart:
        raise RuntimeError(f"No restart file in {source}")
    restart_base = Path(restart.group(1).strip()).name.removesuffix(".result")
    # Keep the mesh and restart paths relative, as the production UDF expects.
    # Give each run unique files in the runtime root and mesh directory.
    sif = re.sub(r'(TES(?: [12])? Series File = String \")[^\"]+',
                 lambda m: m.group(1) + f"phase24_short_ab_{name}_{'2' if m.group(1).startswith('TES 2') else '1'}_series.csv", sif)
    sif = re.sub(r'(TES(?: [12])? Iteration Series File = String \")[^\"]+',
                 lambda m: m.group(1) + f"phase24_short_ab_{name}_{'2' if m.group(1).startswith('TES 2') else '1'}_iterations.csv", sif)
    def copy_state(match: re.Match[str]) -> str:
        tes = "2" if match.group(1).startswith("TES 2") else "1"
        source_state = REPO / "workspace/work/meshes" / mesh_name / f"{restart_base}_{tes}.state"
        if not source_state.is_file() and tes == "1":
            source_state = REPO / "workspace/work/meshes" / mesh_name / f"{restart_base}.state"
        if not source_state.is_file():
            raise FileNotFoundError(f"Missing captured initial TES state: {source_state}")
        target_state = REPO / "workspace/work/meshes" / mesh_name / f"phase24_short_ab_{name}_{tes}.state"
        shutil.copy2(source_state, target_state)
        return match.group(1) + str(Path("workspace/work/meshes") / mesh_name / target_state.name)

    sif = re.sub(r'(TES(?: [12])? State File = String \")([^\"]+)', copy_state, sif)
    opts = [f"  Phase24 Vector Assembly = Logical {'True' if vector else 'False'}"]
    opts.append(f"  Phase24 Static Matrix Reuse = Logical {'True' if static_reuse else 'False'}")
    sif = sif.replace("  Phase24 Wall Profiling = Logical True", "  Phase24 Wall Profiling = Logical True\n" + "\n".join(opts))
    target = HERE / f"{name}.sif"
    target.write_text(sif, encoding="utf-8")
    return target


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True,
                        help="Existing runtime workspace containing workspace/work/meshes and build/udf")
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--runtime-bin", type=Path, required=True)
    parser.add_argument("--toolchain-bin", type=Path, required=True)
    parser.add_argument("--mpi-bin", type=Path, default=Path(r"C:\Program Files\Microsoft MPI\Bin"))
    parser.add_argument("--only", choices=tuple(VARIANTS), action="append")
    args = parser.parse_args()
    global REPO
    REPO = args.repo_root.resolve()
    solver, runtime, toolchain = (p.resolve() for p in (args.solver, args.runtime_bin, args.toolchain_bin))
    mpi_bin = args.mpi_bin.resolve()
    for path in (solver, runtime, toolchain, mpi_bin, REPO / "workspace/work/meshes"):
        if not path.exists():
            raise FileNotFoundError(path)
    env = os.environ.copy()
    prefix = solver.parent.parent
    env["ELMER_HOME"] = str(prefix)
    env["PATH"] = os.pathsep.join(map(str, (solver.parent, prefix / "share/elmersolver/lib",
                                             REPO / "build/udf", runtime, toolchain, mpi_bin,
                                             env.get("PATH", ""))))
    chosen = args.only or list(VARIANTS)
    if "dualtes_hypre_report" in chosen:
        env["PHASE24_HYPRE_REPORT"] = "1"
    manifest_path = HERE / "manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        manifest = {"created": datetime.now().astimezone().isoformat(timespec="seconds"),
                    "solver": str(solver), "solver_sha256": sha256(solver),
                    "heat_solver_sha256": sha256(prefix / "share/elmersolver/lib/HeatSolve.dll"),
                    "mpi_ranks": 6, "variants": []}
    for name in chosen:
        source, vector, static_reuse = VARIANTS[name]
        sif = prepare(source, name, vector, static_reuse)
        sif_text = sif.read_text(encoding="utf-8")
        initial_states = {Path(p).name: sha256(REPO / p.replace("\\", "/"))
                          for p in re.findall(r'TES(?: [12])? State File = String "([^"]+)', sif_text)}
        started = datetime.now().astimezone().isoformat(timespec="seconds")
        command = [str(mpi_bin / "mpiexec.exe"), "-n", "6", str(solver), str(sif)]
        with (HERE / f"{name}.solver.log").open("w", encoding="utf-8") as log:
            proc = subprocess.run(command, cwd=REPO, env=env, stdout=log,
                                  stderr=subprocess.STDOUT, check=False)
        generated = [REPO / f"phase24_short_ab_{name}_{i}_{kind}.csv"
                     for i in (1, 2) for kind in ("series", "iterations")]
        generated.extend(REPO / "workspace/work/meshes" / re.search(
            r'Mesh DB\s+"[^"]+"\s+"([^"]+)"', sif_text).group(1) /
            f"phase24_short_ab_{name}_{i}.state" for i in (1, 2))
        for path in generated:
            if path.is_file():
                shutil.copy2(path, HERE / path.name)
        manifest["variants"] = [v for v in manifest["variants"] if v["name"] != name]
        manifest["variants"].append({"name": name, "source_sif": str(source),
                                      "source_sif_sha256": sha256(source), "sif": sif.name,
                                      "sif_sha256": sha256(sif), "vector_assembly": vector,
                                      "static_matrix_reuse": static_reuse,
                                      "initial_tes_state_sha256": initial_states,
                                      "generated_output_sha256": {p.name: sha256(HERE / p.name)
                                                                  for p in generated if (HERE / p.name).is_file()},
                                      "started": started,
                                      "finished": datetime.now().astimezone().isoformat(timespec="seconds"),
                                      "exit_code": proc.returncode,
                                      "log": f"{name}.solver.log"})
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        if proc.returncode:
            raise SystemExit(proc.returncode)


if __name__ == "__main__":
    main()
