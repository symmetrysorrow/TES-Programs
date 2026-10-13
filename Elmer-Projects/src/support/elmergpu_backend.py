"""Execute compiled project/model cases through the ElmerGPU WSL runtime."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

from support.scenario_config import ROOT
from support.reconcile_project import resolve_parameters
from support.vendored.dimensioned_expression import evaluate_dimensioned_expression

GPU_DISTRO = "ElmerGPU"
GPU_PREFIX = "/opt/elmer-gpu"

STEADY_FAST_KEYWORDS = (
    "Phase24 Prism Lumped Preconditioner = Logical True",
    "Phase24 Precond Matrix Reuse = Logical True",
    "BoomerAMG Aggressive Levels = Integer 1",
    "TES Inner Circuit Newton = Logical True",
    "Phase24 Coupled Temperature Tolerance = Real 1e-7",
)

TRANSIENT_FAST_KEYWORDS = (
    "Phase24 Prism Lumped Preconditioner = Logical True",
    "Phase24 Precond Matrix Reuse = Logical True",
    "Phase24 Coupled Temperature Tolerance = Real 1e-7",
    "TES Inner Circuit Predict From Temperature = Logical True",
    "TES Inner Circuit Newton = Logical True",
    "TES Inner Circuit Newton Slope Reuse = Logical True",
    "BoomerAMG Aggressive Levels = Integer 1",
)


def _merge_keywords(existing: object, defaults: tuple[str, ...], field: str) -> list[str]:
    if not isinstance(existing, list) or any(not isinstance(item, str) for item in existing):
        raise ValueError(f"{field} must be an array of SIF keyword strings")
    present = {item.partition("=")[0].strip().casefold() for item in existing}
    return [*existing, *(item for item in defaults if item.partition("=")[0].strip().casefold() not in present)]


def configure_fast_gpu_case(project_path: Path, case_name: str) -> None:
    """Apply the validated Phase24 fast preset to one generated case."""
    project = json.loads(project_path.read_text(encoding="utf-8"))
    try:
        case = project["cases"][case_name]
    except KeyError as exc:
        raise ValueError(f"generated project has no case {case_name!r}") from exc

    is_transient = case.get("template") == "pulse"
    defaults = TRANSIENT_FAST_KEYWORDS if is_transient else STEADY_FAST_KEYWORDS
    if is_transient:
        parameters = resolve_parameters(project.get("parameter_expressions", {}))
        pulse = case.get("pulse", {})
        try:
            start = evaluate_dimensioned_expression(str(pulse["start"]), variables=parameters).value
            duration = evaluate_dimensioned_expression(str(pulse["duration"]), variables=parameters).value
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{case_name}: ElmerGPU pulse start and duration must be dimensioned expressions") from exc
        event = f"Phase24 BDF1 Step Start Times(2) = Real {start:.17g} {start + duration:.17g}"
        defaults = (event, *defaults)
    case["extra_solver_keywords"] = _merge_keywords(
        case.get("extra_solver_keywords", []), defaults, "solver.extra_solver_keywords"
    )
    solver = case.setdefault("solver", {})
    solver.setdefault("boomer_amg_strong_threshold", 0.25)
    solver.setdefault("boomer_amg_relax_type", 8)
    solver.setdefault("boomer_amg_coarsen_type", 8)
    project_path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _wsl_path(path: Path) -> str:
    resolved = path.resolve()
    if not resolved.drive:
        return resolved.as_posix()
    return f"/mnt/{resolved.drive[0].lower()}{resolved.as_posix()[2:]}"


def ensure_partition(mesh_dir: Path, mpi_procs: int) -> None:
    if mpi_procs < 1:
        raise ValueError("MPI ranks must be at least 1")
    if mpi_procs == 1 or (mesh_dir / f"partitioning.{mpi_procs}").is_dir():
        return

    command = (
        f"export ELMER_HOME={GPU_PREFIX}/elmer "
        f"LD_LIBRARY_PATH={GPU_PREFIX}/elmer/lib/elmersolver "
        f"PATH={GPU_PREFIX}/elmer/bin:$PATH; "
        f"cd '{_wsl_path(mesh_dir.parent)}' && "
        f"ElmerGrid 2 2 '{mesh_dir.name}' -metiskway {mpi_procs}"
    )
    result = subprocess.run(
        ["wsl.exe", "-d", GPU_DISTRO, "-u", "root", "--exec", "bash", "-lc", command],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode or not (mesh_dir / f"partitioning.{mpi_procs}").is_dir():
        raise RuntimeError(
            f"ElmerGPU mesh partitioning failed for {mesh_dir.name} ({mpi_procs} ranks):\n"
            f"{result.stdout[-2000:]}{result.stderr[-2000:]}"
        )


def run_case(project_path: Path, case_name: str, mesh_name: str, mpi_procs: int) -> None:
    if os.name != "nt":
        raise RuntimeError("The ElmerGPU backend requires Windows with the ElmerGPU WSL distribution")
    project_path = project_path.resolve()
    mesh_dir = ROOT / "workspace" / "work" / "meshes" / mesh_name
    if not (mesh_dir / "mesh.header").is_file():
        raise FileNotFoundError(f"ElmerGPU mesh is missing: {mesh_dir}")

    subprocess.run(
        [sys.executable, str(ROOT / "src" / "sync_elmer_parameters.py"), str(project_path)],
        cwd=ROOT,
        check=True,
    )
    ensure_partition(mesh_dir, mpi_procs)

    from support.run_phase24_gpu_elmergpu import run_wsl

    record = run_wsl(project_path, case_name, mesh_name, threads=1, mpi=mpi_procs)
    print(f"ElmerGPU: {record['elapsed_s']:.0f} s ({record['launcher_log']})")
    if record["exit_code"]:
        raise RuntimeError(f"ElmerGPU case {case_name} failed; see {record['launcher_log']}")
