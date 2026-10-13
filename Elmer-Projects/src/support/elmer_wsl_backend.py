"""Execute compiled project/model cases through the Elmer WSL runtime."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

from support.scenario_config import ROOT
from support.reconcile_project import resolve_parameters
from support.vendored.dimensioned_expression import evaluate_dimensioned_expression

WSL_DISTRO = os.environ.get("TES_WSL_DISTRO", "ElmerGPU")
ELMER_PREFIX = "/opt/tes-elmer-cpu"

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


def configure_fast_case(project_path: Path, case_name: str) -> None:
    """Apply validated shared fast settings to one generated case.

    Phase24-prefixed solver keywords are retained because they are the stable
    SIF interface implemented by the Elmer build.
    """
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
            raise ValueError(f"{case_name}: Elmer WSL pulse start and duration must be dimensioned expressions") from exc
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


def ensure_partition(mesh_dir: Path, mpi_procs: int, *, dry_run: bool = False) -> None:
    if mpi_procs < 1:
        raise ValueError("MPI ranks must be at least 1")
    if mpi_procs == 1:
        return

    partition_dir = mesh_dir / f"partitioning.{mpi_procs}"
    expected = [partition_dir / f"part.{rank}.{suffix}" for rank in range(1, mpi_procs + 1)
                for suffix in ("nodes", "elements", "boundary")]
    if partition_dir.is_dir():
        missing = [path.name for path in expected if not path.is_file() or path.stat().st_size == 0]
        if missing:
            raise RuntimeError(
                f"incomplete {mpi_procs}-rank mesh partition at {partition_dir}; "
                f"missing or empty: {', '.join(missing)}"
            )
        return
    if dry_run:
        return

    command = (
        f"export ELMER_HOME={ELMER_PREFIX}/elmer "
        f"LD_LIBRARY_PATH={ELMER_PREFIX}/elmer/lib/elmersolver "
        f"PATH={ELMER_PREFIX}/elmer/bin:$PATH; "
        f"cd '{_wsl_path(mesh_dir.parent)}' && "
        f"ElmerGrid 2 2 '{mesh_dir.name}' -metiskway {mpi_procs}"
    )
    result = subprocess.run(
        ["wsl.exe", "-d", WSL_DISTRO, "-u", "root", "--exec", "bash", "-lc", command],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode or not (mesh_dir / f"partitioning.{mpi_procs}").is_dir():
        raise RuntimeError(
            f"Elmer WSL mesh partitioning failed for {mesh_dir.name} ({mpi_procs} ranks):\n"
            f"{result.stdout[-2000:]}{result.stderr[-2000:]}"
        )
    missing = [path.name for path in expected if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(
            f"ElmerGrid created an invalid {mpi_procs}-rank partition at {partition_dir}; "
            f"missing or empty: {', '.join(missing)}"
        )


def run_wsl(project: Path, case: str, mesh: str, threads: int = 1, elmer: str = "elmer",
            extra_env: list[str] | None = None, mpi: int = 1) -> dict:
    """Run the common Elmer project runner inside the ElmerGPU WSL distro."""
    (ROOT / mesh).mkdir(parents=True, exist_ok=True)
    bash = f"""set -euo pipefail
export PATH={ELMER_PREFIX}/venv/bin:{ELMER_PREFIX}/elmer/bin:/usr/sbin:/usr/bin:/sbin:/bin
export ELMER_HOME={ELMER_PREFIX}/{elmer}
export LD_LIBRARY_PATH={ELMER_PREFIX}/udf:{ELMER_PREFIX}/hypre/lib:{ELMER_PREFIX}/{elmer}/lib/elmersolver:/usr/lib/wsl/lib
export OMP_NUM_THREADS={threads}
{"".join(f"export {kv}" + chr(10) for kv in (extra_env or []))}
cd '{_wsl_path(ROOT)}'
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
{ELMER_PREFIX}/venv/bin/python src/run.py '{case}' --project '{_wsl_path(project)}' --skip-sync --mpi-procs {mpi} --backend elmer_wsl --elmer-solver {ELMER_PREFIX}/{elmer}/bin/{'ElmerSolver_mpi' if mpi > 1 else 'ElmerSolver'} --runtime-bin ''
"""
    log = project.parent / f"{case}.wsl_launcher.log"
    start = time.monotonic()
    with log.open("w", encoding="utf-8") as handle:
        result = subprocess.run(
            ["wsl.exe", "-d", WSL_DISTRO, "-u", "root", "--exec", "bash", "-lc", bash],
            cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT,
        )
    return {"case": case, "exit_code": result.returncode,
            "elapsed_s": time.monotonic() - start, "launcher_log": str(log)}


def run_case(project_path: Path, case_name: str, mesh_name: str, mpi_procs: int) -> None:
    if os.name != "nt":
        raise RuntimeError("The Elmer WSL backend requires Windows with the ElmerGPU WSL distribution")
    project_path = project_path.resolve()
    mesh_dir = ROOT / "workspace" / "work" / "meshes" / mesh_name
    if not (mesh_dir / "mesh.header").is_file():
        raise FileNotFoundError(f"Elmer WSL mesh is missing: {mesh_dir}")

    subprocess.run(
        [sys.executable, str(ROOT / "src" / "sync_elmer_parameters.py"), str(project_path)],
        cwd=ROOT,
        check=True,
    )
    ensure_partition(mesh_dir, mpi_procs)
    record = run_wsl(project_path, case_name, mesh_name, threads=1, mpi=mpi_procs)
    print(f"Elmer WSL: {record['elapsed_s']:.0f} s ({record['launcher_log']})")
    if record["exit_code"]:
        raise RuntimeError(f"Elmer WSL case {case_name} failed; see {record['launcher_log']}")
