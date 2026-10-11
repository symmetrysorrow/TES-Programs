#!/usr/bin/env python3
"""Re-run the isolated four-step DualTES Vector-OFF cache A/B pair.

The supplied runtime root must be an isolated copy containing the DualTES mesh,
its 6-way partition, restart shards, ElmerSolver_mpi.exe, UDF DLLs and Lua
scripts. This runner resets both TES circuit state files from this evidence
folder before each case so OFF and ON start from the same state.
"""
from __future__ import annotations
import argparse, datetime, hashlib, json, os, shutil, subprocess, time
from pathlib import Path

HERE = Path(__file__).resolve().parent
STATE_NAMES = [f"phase24_short_ab_dualtes_vector_{i}.state" for i in (1, 2)]
SIF_NAMES = [f"dualtes_field_{label}_6rank.sif" for label in ("off", "on")]
MESH_TAG = "mesh_toml_95333801e1a3"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runtime-root", type=Path, required=True)
    ap.add_argument("--build-root", type=Path, required=True)
    ap.add_argument("--mpi-exe", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--udf-dir", type=Path)
    ap.add_argument("--elmer-extra-bin", type=Path)
    args = ap.parse_args()
    runtime = args.runtime_root.resolve()
    build = args.build_root.resolve()
    output = args.output_dir.resolve()
    mesh = runtime / "workspace" / "work" / "meshes" / MESH_TAG
    if not (runtime / "bin" / "ElmerSolver_mpi.exe").is_file():
        raise SystemExit(f"missing isolated solver executable under {runtime / 'bin'}")
    if not (mesh / "partitioning.6").is_dir() or not (mesh / "case_toml_2cc534365acf_steady.result.0").is_file():
        raise SystemExit(f"isolated runtime is missing 6-rank mesh/restart data under {mesh}")
    if not (runtime / MESH_TAG).is_dir():
        raise SystemExit(f"restart/output mesh alias missing: {runtime / MESH_TAG}")
    output.mkdir(parents=True, exist_ok=True)
    udf = args.udf_dir.resolve() if args.udf_dir else runtime / "build" / "udf"
    path_parts = [
        udf, build / "fem" / "src" / "modules", build / "fem" / "src",
        build / "matc" / "src", build / "fhutiter" / "src",
        build / "fem" / "src" / "binio", build / "mathlibs" / "src" / "arpack",
        build / "mathlibs" / "src" / "parpack", runtime / "bin",
        runtime / "share" / "elmersolver" / "lib",
        Path(r"G:\Elmer-Migration\msys64\ucrt64\bin"),
        Path(r"G:\Elmer-Migration\msys64\usr\bin"),
    ]
    if args.elmer_extra_bin:
        path_parts.append(args.elmer_extra_bin.resolve())
    path_parts.append(args.mpi_exe.resolve().parent)
    env = os.environ.copy()
    env.update(ELMER_HOME=str(runtime), OMP_NUM_THREADS="1", OMP_DYNAMIC="FALSE")
    env["PATH"] = os.pathsep.join([str(x) for x in path_parts] + [env.get("PATH", "")])
    records = []
    for label in ("off", "on"):
        for state_name in STATE_NAMES:
            shutil.copy2(HERE / state_name, mesh / state_name)
        sif_source = HERE / f"dualtes_field_{label}_6rank.sif"
        sif_runtime = runtime / f"dualtes_matc_cache_{label}.sif"
        shutil.copy2(sif_source, sif_runtime)
        log = output / f"dualtes_field_{label}_6rank.solver.log"
        start_wall = time.perf_counter()
        started = datetime.datetime.now().astimezone().isoformat(timespec="seconds")
        with log.open("w", encoding="utf-8") as stream:
            proc = subprocess.run(
                [str(args.mpi_exe.resolve()), "-n", "6",
                 str(runtime / "bin" / "ElmerSolver_mpi.exe"), sif_runtime.name],
                cwd=runtime, env=env, stdout=stream, stderr=subprocess.STDOUT,
                check=False,
            )
        elapsed = time.perf_counter() - start_wall
        outputs = []
        for tes in (1, 2):
            for kind in ("series", "iterations"):
                name = f"dualtes_field_{label}_{tes}_{kind}.csv"
                src = runtime / name
                if src.exists():
                    dst = output / name
                    shutil.copy2(src, dst)
                    outputs.append({"file": name, "sha256": sha256(dst)})
        result = runtime / MESH_TAG / f"dualtes_field_{label}.result"
        for src in ([result] if result.exists() else []) + sorted(result.parent.glob(result.name + ".*")):
            dst = output / src.name
            shutil.copy2(src, dst)
            outputs.append({"file": src.name, "sha256": sha256(dst), "bytes": dst.stat().st_size})
        record = {
            "label": label, "returncode": proc.returncode,
            "started": started, "elapsed_wall_seconds": elapsed,
            "log": log.name, "sif_sha256": sha256(sif_source), "outputs": outputs,
        }
        records.append(record)
        print(json.dumps(record, ensure_ascii=False))
    (output / "run_manifest.json").write_text(
        json.dumps({"ranks": 6, "steps": 4, "vector_assembly": False,
                    "build": "isolated Debug -O0; see parent manifest", "runs": records},
                   indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0 if all(x["returncode"] == 0 for x in records) else 1

if __name__ == "__main__":
    raise SystemExit(main())
