"""Run an Elmer case defined in elmer_project.json (redesign plan, Phase 2).

    python run.py <case_name> [--dry-run] [--skip-sync] [--force-deps]

- Regenerates generated/cases/ first (unless --skip-sync).
- Resolves the restart chain: if the case has `restart_from` and the
  dependency's `.result` file is missing, the dependency is run first
  (recursively). `--force-deps` reruns dependencies even when present.
- Runs ElmerSolver from the repository root, captures the log, and collects
  the case outputs (VTU/EP from the mesh directory, the series CSV from the
  root, the solver log) into results/<case>/. The `.result` file stays in the
  mesh directory because it is the restart interface between cases.
- Writes results/<case>/manifest.json with input hashes and the run outcome.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ELMERSOLVER = r"C:\Program Files\Elmer 26.1-Release\bin\ElmerSolver.exe"


def load_model(project_path: Path) -> dict:
    sys.path.insert(0, str(ROOT))
    from scripts.support.reconcile_project import reconcile_project

    return reconcile_project(json.loads(project_path.read_text(encoding="utf-8")))


def find_case_insensitive(path: Path) -> Path | None:
    """Return an existing path matching *path*, tolerating a different case
    on the final component. Windows/NTFS is case-preserving on create but
    case-insensitive on lookup; the dual-TES UDF's OPEN(..., STATUS='REPLACE')
    reuses whichever case an earlier run created for a series CSV, so the
    file that lands on disk does not always match the case written into the
    SIF (e.g. '..._l_series.csv' on disk vs '..._L_series.csv' in Constants).
    Returns None if no match is found."""
    if path.exists():
        return path
    if not path.parent.exists():
        return None
    target = path.name.lower()
    for candidate in path.parent.iterdir():
        if candidate.name.lower() == target:
            return candidate
    return None


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def runtime_environment(
    elmer_solver: str, runtime_bin: str | None
) -> tuple[dict[str, str], Path, Path]:
    """Return a loader environment pinned to one Elmer installation."""
    solver = Path(elmer_solver).resolve()
    if not solver.is_file():
        raise FileNotFoundError(f"ElmerSolver not found: {solver}")
    prefix = solver.parent.parent
    parts = [str(solver.parent), str(prefix / "share" / "elmersolver" / "lib")]
    if runtime_bin:
        runtime = Path(runtime_bin).resolve()
        if not runtime.is_dir():
            raise FileNotFoundError(f"runtime DLL directory not found: {runtime}")
        parts.append(str(runtime))
    env = os.environ.copy()
    env["ELMER_HOME"] = str(prefix)
    env["PATH"] = os.pathsep.join([*parts, env.get("PATH", "")])
    return env, solver, prefix


def write_runtime_sif(
    source: Path, out_dir: Path, udf_dll: str | None
) -> tuple[Path, bool]:
    """Create an execution-only SIF with an explicitly pinned UDF DLL."""
    if udf_dll is None:
        return source, False
    dll = Path(udf_dll).resolve()
    if not dll.is_file():
        raise FileNotFoundError(f"UDF DLL not found: {dll}")
    token = '"tes_transient_heat_source_t0"'
    text = source.read_text(encoding="utf-8")
    target_functions = ("TESTransientHeatSource", "AbsorberWindowPulseHeatSource")
    uses_target = any(name in text for name in target_functions)
    if not uses_target:
        return source, False
    if token not in text:
        raise ValueError(f"{source}: expected Procedure library token {token} was not found")
    out_dir.mkdir(parents=True, exist_ok=True)
    runtime_sif = out_dir / "runtime.sif"
    runtime_sif.write_text(
        text.replace(token, f'"{dll.as_posix()}"'),
        encoding="utf-8",
        newline="\n",
    )
    return runtime_sif, True


def runtime_artifacts(
    solver: Path, prefix: Path, udf_dll: str | None
) -> dict[str, str]:
    """Best-effort file hashes for loader-provenance in a run manifest."""
    paths = {"solver": solver, "libelmersolver": solver.parent / "libelmersolver.dll"}
    heat = prefix / "share" / "elmersolver" / "lib" / "HeatSolve.dll"
    if heat.exists():
        paths["HeatSolve"] = heat
    if udf_dll:
        paths["udf"] = Path(udf_dll).resolve()
    return {name: sha256(path) for name, path in paths.items() if path.is_file()}


def mpi_launcher(env: dict[str, str]) -> str:
    """Resolve mpiexec against the pinned child environment, not the parent."""
    launcher = shutil.which("mpiexec", path=env["PATH"])
    if launcher is None:
        raise FileNotFoundError("mpiexec was not found in the selected runtime PATH")
    return str(Path(launcher).resolve())


def restart_chain(model: dict, case_name: str) -> list[str]:
    """Dependency-first list of cases to consider for *case_name*."""
    chain: list[str] = []
    seen: set[str] = set()

    def visit(name: str) -> None:
        if name in seen:
            raise ValueError(f"restart_from cycle involving '{name}'")
        seen.add(name)
        spec = model["cases"].get(name)
        if spec is None:
            raise ValueError(f"unknown case '{name}'")
        dep = spec.get("restart_from")
        if dep:
            visit(dep)
        chain.append(name)

    visit(case_name)
    return chain


def mesh_dir_of(model: dict, case_name: str) -> Path:
    spec = model["cases"][case_name]
    return ROOT / "work" / "meshes" / model["meshes"][spec["mesh"]]["dir"]


def result_file_of(model: dict, case_name: str) -> Path | None:
    spec = model["cases"][case_name]
    if not spec.get("output_result"):
        return None
    serial = mesh_dir_of(model, case_name) / f"{case_name}.result"
    # MPI ResultOutput writes one restart file per rank, conventionally using
    # ``.result.0``, ``.result.1``, ... .  Returning rank 0 here makes the
    # dependency resolver treat a completed parallel steady run as a valid
    # restart source while the SIF still references the common base name.
    parallel_rank0 = mesh_dir_of(model, case_name) / f"{case_name}.result.0"
    return serial if serial.exists() or not parallel_rank0.exists() else parallel_rank0


def preexisting_restart_paths(model: dict, case_name: str, mpi_procs: int) -> list[Path]:
    """Resolve the complete restart interface required by an explicit case."""
    spec = model["cases"][case_name]
    if not spec.get("preexisting_restart"):
        return []
    base = spec.get("restart_file_base")
    if not base:
        raise ValueError(f"{case_name}: preexisting_restart requires restart_file_base")
    state = spec.get("state_file")
    if not state:
        raise ValueError(f"{case_name}: preexisting_restart requires state_file")
    mesh = mesh_dir_of(model, case_name)
    results = (
        [mesh / f"{base}.result"]
        if mpi_procs == 1
        else [mesh / f"{base}.result.{rank}" for rank in range(mpi_procs)]
    )
    return [*results, ROOT / state]


def validate_preexisting_restart(
    model: dict, case_name: str, mpi_procs: int
) -> list[Path]:
    needed = preexisting_restart_paths(model, case_name, mpi_procs)
    missing = [str(path) for path in needed if not path.is_file()]
    if missing:
        raise FileNotFoundError("preexisting restart missing: " + ", ".join(missing))
    return needed


def run_case(
    model: dict,
    case_name: str,
    project_path: Path,
    elmer_solver: str,
    mpi_procs: int,
    udf_dll: str | None = None,
    runtime_bin: str | None = None,
) -> int:
    spec = model["cases"][case_name]
    sif = ROOT / "generated" / "cases" / f"{case_name}.sif"
    mesh_dir = mesh_dir_of(model, case_name)
    out_dir = ROOT / "results" / case_name
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "solver.log"
    runtime_sif, udf_applied = write_runtime_sif(sif, out_dir, udf_dll)
    env, solver_path, prefix = runtime_environment(elmer_solver, runtime_bin)

    inputs = {
        project_path.name: sha256(project_path),
        sif.name: sha256(sif),
        "mesh.header": sha256(mesh_dir / "mesh.header"),
    }
    preexisting_inputs = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in validate_preexisting_restart(model, case_name, mpi_procs)
    }
    root_dlls = ("tes_heat_source_t0.dll",)
    if udf_dll is None:
        root_dlls = ("tes_transient_heat_source_t0.dll", *root_dlls)
    for dll in root_dlls:
        if (ROOT / dll).exists():
            inputs[dll] = sha256(ROOT / dll)

    started = datetime.now().isoformat(timespec="seconds")
    display_sif = (
        runtime_sif.relative_to(ROOT)
        if runtime_sif.is_relative_to(ROOT)
        else runtime_sif
    )
    print(f"[{case_name}] ElmerSolver {display_sif} (log: {log_path.relative_to(ROOT)})")
    with log_path.open("w", encoding="utf-8") as log:
        command = [str(solver_path), str(runtime_sif)]
        if mpi_procs > 1:
            command = [mpi_launcher(env), "-n", str(mpi_procs), *command]
        proc = subprocess.run(
            command,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
        )
    finished = datetime.now().isoformat(timespec="seconds")

    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    convergence = re.findall(r"ComputeChange:.*", log_text)[-4:]
    fatal = re.findall(r"(?:ERROR|Fatal).*", log_text)[:3]

    collected: list[str] = []
    for pattern in (f"{case_name}_t*.vtu", f"{case_name}*.pvtu", f"{case_name}.ep"):
        for src in mesh_dir.glob(pattern):
            shutil.move(str(src), out_dir / src.name)
            collected.append(src.name)
    series = spec.get("series_file")
    if series:
        found = find_case_insensitive(ROOT / series)
        if found is not None:
            # Single-instance cases (single-pixel): one series CSV, exact name.
            shutil.move(str(found), out_dir / series)
            collected.append(series)
        else:
            # Dual-TES cases: the builder writes one series CSV per circuit
            # instance with an 'L'/'R' side tag inserted before the
            # '_series.csv' suffix (scripts/support/build_cases.py,
            # _side_series_file) instead of the single base name.
            sys.path.insert(0, str(ROOT))
            from scripts.support.build_cases import _side_series_file

            for side in ("L", "R"):
                side_name = _side_series_file(series, side)
                found_side = find_case_insensitive(ROOT / side_name)
                if found_side is not None:
                    shutil.move(str(found_side), out_dir / side_name)
                    collected.append(side_name)
    iteration_series = spec.get("iteration_series_file")
    if iteration_series:
        found_iteration = find_case_insensitive(ROOT / iteration_series)
        if found_iteration is not None:
            shutil.move(str(found_iteration), out_dir / iteration_series)
            collected.append(iteration_series)
    result_file = result_file_of(model, case_name)

    manifest = {
        "case": case_name,
        "started": started,
        "finished": finished,
        "exit_code": proc.returncode,
        "inputs_sha256": inputs,
        "preexisting_restart_inputs_sha256": preexisting_inputs,
        "runtime_sif": str(runtime_sif.relative_to(ROOT)) if runtime_sif != sif else None,
        "runtime_sif_sha256": sha256(runtime_sif),
        "udf_applied": udf_applied,
        "udf_dll": str(Path(udf_dll).resolve()) if udf_applied else None,
        "udf_sha256": sha256(Path(udf_dll).resolve()) if udf_applied else None,
        "solver": str(solver_path),
        "elmer_prefix": str(prefix),
        "runtime_bin": str(Path(runtime_bin).resolve()) if runtime_bin else None,
        "runtime_artifacts_sha256": runtime_artifacts(
            solver_path, prefix, udf_dll if udf_applied else None
        ),
        "mesh": spec["mesh"],
        "restart_from": spec.get("restart_from"),
        "collected_outputs": sorted(collected),
        "result_file": str(result_file.relative_to(ROOT)) if result_file else None,
        "convergence_tail": convergence,
        "errors": fatal,
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    status = "OK" if proc.returncode == 0 and not fatal else "FAILED"
    print(f"[{case_name}] {status} (exit {proc.returncode}); outputs -> {out_dir.relative_to(ROOT)}")
    for line in convergence[-2:]:
        print(f"[{case_name}]   {line.strip()}")
    for line in fatal:
        print(f"[{case_name}]   {line.strip()}")
    return proc.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case")
    parser.add_argument("--dry-run", action="store_true", help="print the plan without running")
    parser.add_argument("--skip-sync", action="store_true", help="do not regenerate generated/cases first")
    parser.add_argument("--force-deps", action="store_true", help="rerun restart dependencies even if their .result exists")
    parser.add_argument(
        "--project",
        default="elmer_project.json",
        help="path to the project JSON (default: elmer_project.json). Use an "
        "alternate file to try different parameters without touching the "
        "main config; its case must use a distinct name so results/ and "
        "generated/cases/ outputs don't collide with the main project's.",
    )
    parser.add_argument(
        "--udf-dll",
        help="pin tes_transient_heat_source_t0 to this DLL in results/<case>/runtime.sif",
    )
    parser.add_argument(
        "--runtime-bin",
        help="runtime DLL directory appended after the selected Elmer install's bin/module paths",
    )
    parser.add_argument(
        "--elmer-solver",
        default=ELMERSOLVER,
        help="path to ElmerSolver executable; use this to select an alternate build such as the HYPRE/MPI solver",
    )
    parser.add_argument(
        "--mpi-procs",
        type=int,
        default=1,
        help="number of MPI ranks; requires a mesh/partitioning.N directory when N > 1",
    )
    args = parser.parse_args()

    project_path = Path(args.project)
    if not project_path.is_absolute():
        project_path = ROOT / project_path

    if not args.skip_sync:
        subprocess.run(
            [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project_path)],
            cwd=ROOT,
            check=True,
        )

    model = load_model(project_path)
    chain = restart_chain(model, args.case)
    validate_preexisting_restart(model, args.case, args.mpi_procs)
    target_preexisting = bool(model["cases"][args.case].get("preexisting_restart"))
    if target_preexisting and args.force_deps:
        raise ValueError("--force-deps cannot be used with preexisting_restart")

    plan: list[str] = []
    for name in chain[:-1]:
        result = result_file_of(model, name)
        if result is None:
            raise ValueError(f"'{name}' is a restart dependency but does not set output_result")
        if not target_preexisting and (args.force_deps or not result.exists()):
            plan.append(name)
        else:
            print(f"[{name}] restart field {result.relative_to(ROOT)} exists - skipping (use --force-deps to rerun)")
    plan.append(args.case)

    print("run plan: " + " -> ".join(plan))
    if args.dry_run:
        return 0

    for name in plan:
        code = run_case(
            model, name, project_path, args.elmer_solver, args.mpi_procs,
            args.udf_dll, args.runtime_bin,
        )
        if code != 0:
            print(f"aborting chain: {name} failed")
            return code
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
