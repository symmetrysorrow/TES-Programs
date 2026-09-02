"""Run the optimized SinglePixel mesh with an original or hybrid time grid.

This keeps the production-v2 spatial mesh and physics, but replaces its
provisional 10 us / 100 us / 1 ms tail schedule with the exact timestep
schedule used by ``case_tes_mpi_comsol_grid_full_uniform_continuous``.

One-command usage from the repository root::

    python scripts/prep/run_singlepixel_prod_v2_original_timegrid.py

The default linear solver is the MPI-safe HYPRE/BoomerAMG iterative solver.
Use ``--linear-system mumps`` only for a deliberate MUMPS retry.

The ``hybrid`` time grid keeps the validated fine rise-time grid through
100 us after the pulse and then joins the original COMSOL-grid tail.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PROJECT = ROOT / "elmer_project_comsol_timegrid.json"
OUTPUT_PROJECT = ROOT / "elmer_project_singlepixel_prod_v2_original_timegrid.json"

BASE_MESH = "mesh_singlepixel_prod_v2"
STEADY = "case_tes_steady_singlepixel_prod_v2_original_timegrid"
PULSE = "case_tes_pulse_singlepixel_prod_v2_original_timegrid"
HYBRID_STEADY = f"{STEADY}_hybrid"
HYBRID_PULSE = f"{PULSE}_hybrid"
KNOWN_SERIAL_STEADY = "case_tes_steady_singlepixel_prod_v2"
KNOWN_SERIAL_ITERATIONS = (
    ROOT
    / "results"
    / KNOWN_SERIAL_STEADY
    / f"{KNOWN_SERIAL_STEADY}_iterations.csv"
)

# Validated production-v2 rise-time grid from the earlier optimization run.
PULSE_PREFIX = [
    ["18[us]", 1],
    ["1[us]", 2],
    ["1[ns]", 1],
    ["10[ns]", 10],
    ["100[ns]", 9],
    ["1[us]", 9],
]
FINE_STEP = ["0.625[us]", 144]

DEFAULT_SOLVER = Path(
    r"D:\Github\TES-Programs\tools\elmer-hypre\install-phase13-step-commit\bin\ElmerSolver.exe"
)
DEFAULT_RUNTIME_BIN = Path(r"C:\msys64\ucrt64\bin")


def timestep_seconds(token: str) -> float:
    value, unit = token[:-1].split("[")
    return float(value) * {"s": 1.0, "ms": 1.0e-3, "us": 1.0e-6, "ns": 1.0e-9}[unit]


def hybrid_timesteps(original: list[list[object]]) -> list[list[object]]:
    """Use the validated rise grid, then preserve original tail timestamps."""
    early = [*PULSE_PREFIX, FINE_STEP]
    early_end = sum(
        timestep_seconds(str(token)) * int(count)
        for token, count in early
    )

    cumulative = 0.0
    for index, (token, count) in enumerate(original):
        group_end = cumulative + timestep_seconds(str(token)) * int(count)
        if group_end > early_end:
            bridge = group_end - early_end
            result = [*early]
            if bridge > 1.0e-18:
                result.append([f"{bridge:.17g}[s]", 1])
            # The bridge lands exactly on an original-grid sample.  Appending
            # subsequent original groups therefore preserves the COMSOL tail
            # endpoint and all of its later time points.
            result.extend(original[index + 1 :])
            return result
        cumulative = group_end
    raise ValueError("original time grid ends before the hybrid rise grid")


def truncate_timesteps(
    timesteps: list[list[object]], max_steps: int
) -> list[list[object]]:
    """Return the first *max_steps* while preserving grouped SIF stages."""
    if max_steps < 1:
        raise ValueError("max_steps must be >= 1")
    remaining = max_steps
    result: list[list[object]] = []
    for token, count in timesteps:
        take = min(int(count), remaining)
        if take:
            result.append([token, take])
            remaining -= take
        if remaining == 0:
            return result
    raise ValueError(
        f"requested {max_steps} smoke steps, but schedule has only "
        f"{max_steps - remaining}"
    )


def steady_state_from_iterations(path: Path) -> tuple[float, float, float, float, float]:
    """Recover the persisted circuit state from a converged steady log."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"steady iteration log is empty: {path}")
    last = rows[-1]
    return (
        float(last["tes_temperature_K"]),
        float(last["raw_current_A"]),
        float(last["tes_resistance_ohm"]),
        float(last["relaxed_power_W"]),
        float(last["previous_current_A"]),
    )


def write_known_steady_state(state_path: Path) -> None:
    """Create a fresh transient seed paired with the known serial result."""
    values = steady_state_from_iterations(KNOWN_SERIAL_ITERATIONS)
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(
        "".join(f"{value:24.16E}" for value in values) + "\n",
        encoding="ascii",
        newline="\n",
    )


def ensure_partition(mpi_procs: int) -> str:
    if mpi_procs <= 1:
        return BASE_MESH

    mesh_root = ROOT / "work" / "meshes"
    # Keep a separate partition directory for each MPI size.  Reusing the
    # 4-way directory for a different rank count would leave Elmer with the
    # wrong partitioning.N metadata.
    partitioned = (
        f"{BASE_MESH}_repart_x"
        if mpi_procs == 4
        else f"{BASE_MESH}_repart_x{mpi_procs}"
    )
    partition_dir = mesh_root / partitioned / f"partitioning.{mpi_procs}"
    if not partition_dir.is_dir():
        subprocess.run(
            [
                "ElmerGrid",
                "2",
                "2",
                BASE_MESH,
                "-partition",
                str(mpi_procs),
                "1",
                "1",
                "-out",
                partitioned,
            ],
            cwd=mesh_root,
            check=True,
        )
    # ElmerGrid's partition mode writes the partitioning directory but may
    # omit the serial mesh metadata at the output root.  build_cases.py and
    # ElmerSolver both use that metadata, so retain a complete mesh directory
    # without modifying the original mesh.
    source_dir = mesh_root / BASE_MESH
    target_dir = mesh_root / partitioned
    for name in ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary", "mesh.names", "entities.sif"):
        source = source_dir / name
        target = target_dir / name
        if source.is_file() and not target.exists():
            shutil.copy2(source, target)
    return partitioned


def make_project(
    mesh_dir: str,
    linear_system: str,
    reuse_steady: bool,
    reuse_known_serial_steady: bool,
    time_grid: str,
    smoke_steps: int | None = None,
    smoke_label: str = "amgx",
) -> str:
    steady_case = HYBRID_STEADY if time_grid == "hybrid" else STEADY
    pulse_case = HYBRID_PULSE if time_grid == "hybrid" else PULSE
    if smoke_steps is not None:
        pulse_case = f"{pulse_case}_{smoke_label}_smoke_{smoke_steps}step"
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    mesh_name = mesh_dir
    project["meshes"][mesh_name] = {
        "geometry": "single_pixel",
        "dir": mesh_dir,
        "notes": (
            "Optimized SinglePixel production-v2 mesh with the original "
            "COMSOL-grid transient schedule."
        ),
    }

    steady = {
        "template": "steady",
        "mesh": mesh_name,
        "heat_source": "circuit_inner",
        "initial_temperature": "T_0",
        "output_result": True,
        "vtu": False,
        "steady_state_max_iterations": 1,
        "output_intervals": 1,
        "solver": {
            "nonlinear_max_iterations": 120,
            "nonlinear_convergence_tolerance": 1e-8,
            "nonlinear_relaxation_factor": 1.0,
            "steady_state_convergence_tolerance": 1e-8,
            "linear_system": linear_system,
        },
        # The UDF opens this path relative to run.py's repository-root CWD.
        "state_file": f"work/meshes/{mesh_dir}/{steady_case}.state",
        "series_file": f"{steady_case}_series.csv",
        "iteration_series_file": f"{steady_case}_iterations.csv",
        "output_file_path": f"../work/meshes/{mesh_dir}/{steady_case}.result",
    }
    project["cases"][steady_case] = steady

    original = copy.deepcopy(
        project["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]
    )
    if time_grid == "hybrid":
        original["timesteps"] = hybrid_timesteps(original["timesteps"])
        original["output_intervals"] = [1] * len(original["timesteps"])
        original["inner_circuit_step_commit"] = True
        original["solver"] = {
            **original["solver"],
            "nonlinear_max_iterations": 120,
            "nonlinear_convergence_tolerance": 1e-8,
            "nonlinear_relaxation_factor": 1.0,
        }
    if smoke_steps is not None:
        original["timesteps"] = truncate_timesteps(
            original["timesteps"], smoke_steps
        )
        # Save only the terminal smoke step; VTU output remains disabled.
        original["output_intervals"] = [999999] * len(original["timesteps"])
        original["output_intervals"][-1] = 1
    original.update(
        {
            "mesh": mesh_name,
            "restart_from": steady_case,
            "restart_time": 0.02,
            "restart_file_path": f"../work/meshes/{mesh_dir}/{steady_case}.result",
            "state_file": f"work/meshes/{mesh_dir}/{steady_case}.state",
            "series_file": f"{pulse_case}_series.csv",
            "iteration_series_file": f"{pulse_case}_iterations.csv",
            "output_file_path": f"../work/meshes/{mesh_dir}/{pulse_case}.result",
            "vtu": False,
            "comparison_time_grid": {
                "source_case": "case_tes_mpi_comsol_grid_full_uniform_continuous",
                "source_project": SOURCE_PROJECT.name,
                "description": (
                    "Optimized SinglePixel mesh with a validated fine rise grid "
                    "and the original COMSOL-grid tail."
                ),
            },
        }
    )
    original["solver"] = dict(original["solver"])
    original["solver"]["linear_system"] = linear_system
    if reuse_steady:
        # The steady MPI restart already exists.  Remove the dependency edge
        # so run.py executes only the pulse case and does not rerun the
        # multi-hour steady solve.
        original["restart_from"] = None
        original["restart_file_base"] = steady_case
    elif reuse_known_serial_steady:
        # This result was produced successfully with the same production-v2
        # mesh and physics.  It is a serial restart, so this mode is only
        # valid when the pulse itself is also run with one MPI rank.
        original["restart_from"] = None
        original["restart_file_base"] = KNOWN_SERIAL_STEADY
        original["restart_file_path"] = (
            f"../work/meshes/{BASE_MESH}/{KNOWN_SERIAL_STEADY}.result"
        )
        # The historical validated result predates persisted circuit-state
        # files.  Reconstruct its matching electrical state from the final
        # converged iteration.  Use the run-specific path because transient
        # checkpointing overwrites it as timesteps commit.
        state_path = ROOT / str(original["state_file"])
        write_known_steady_state(state_path)
    project["cases"][pulse_case] = original
    OUTPUT_PROJECT.write_text(
        json.dumps(project, indent=2) + "\n", encoding="utf-8"
    )
    print(f"prepared {OUTPUT_PROJECT.relative_to(ROOT)}")
    print(f"mesh={mesh_dir}, steady={steady_case}, pulse={pulse_case}")
    print(f"time_grid={time_grid}, timestep stages={len(original['timesteps'])}")
    return pulse_case


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mpi-procs",
        type=int,
        default=4,
        help="MPI ranks; default 4, use 1 for the serial run",
    )
    parser.add_argument(
        "--elmer-solver",
        default=str(DEFAULT_SOLVER),
        help="ElmerSolver executable",
    )
    parser.add_argument(
        "--linear-system",
        choices=["iterative_hypre_boomeramg", "iterative", "mumps", "direct"],
        default="iterative_hypre_boomeramg",
        help=(
            "linear solver for steady and pulse; default is the MPI-safe "
            "HYPRE/BoomerAMG solver"
        ),
    )
    parser.add_argument(
        "--smoke-steps",
        type=int,
        help=(
            "run only the first N transient steps under a distinct case name; "
            "intended for GPU convergence and memory checks"
        ),
    )
    parser.add_argument(
        "--smoke-label",
        choices=["amgx", "cpu"],
        default="amgx",
        help="label for a truncated benchmark case name (default: amgx)",
    )
    parser.add_argument(
        "--time-grid",
        choices=["original", "hybrid"],
        default="original",
        help=(
            "original COMSOL grid or hybrid grid with the validated fine "
            "rise-time steps; default original"
        ),
    )
    parser.add_argument(
        "--reuse-steady",
        action="store_true",
        help="reuse the existing successful MPI steady restart and run pulse only",
    )
    parser.add_argument(
        "--reuse-known-serial-steady",
        action="store_true",
        help=(
            "reuse the previously validated serial production-v2 steady "
            "restart; requires --mpi-procs 1"
        ),
    )
    parser.add_argument(
        "--runtime-bin",
        default=str(DEFAULT_RUNTIME_BIN),
        help="runtime DLL directory for the selected Elmer build",
    )
    parser.add_argument(
        "--amgx-config",
        help=(
            "AMGX JSON config injected into the execution-only SIF; useful "
            "when launching the GPU build under WSL"
        ),
    )
    parser.add_argument(
        "--amgx-constraint-mode",
        choices=[
            "default", "no-scaling", "slave", "master", "slave-transpose", "master-transpose",
            "dual-lagrange", "penalty", "schur", "stabilized",
        ],
        default="default",
    )
    parser.add_argument(
        "--amgx-constraint-penalty", type=float, default=1.0e4,
    )
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument(
        "--force-deps",
        action="store_true",
        help="recompute restart dependencies before running the target case",
    )
    args = parser.parse_args()
    if args.mpi_procs < 1:
        parser.error("--mpi-procs must be >= 1")
    if args.smoke_steps is not None and args.smoke_steps < 1:
        parser.error("--smoke-steps must be >= 1")
    if args.reuse_steady and args.reuse_known_serial_steady:
        parser.error("--reuse-steady and --reuse-known-serial-steady are mutually exclusive")
    if args.reuse_known_serial_steady and args.mpi_procs != 1:
        parser.error("--reuse-known-serial-steady requires --mpi-procs 1")
    if not SOURCE_PROJECT.is_file():
        raise FileNotFoundError(SOURCE_PROJECT)

    mesh_dir = ensure_partition(args.mpi_procs)
    pulse_case = make_project(
        mesh_dir,
        args.linear_system,
        args.reuse_steady,
        args.reuse_known_serial_steady,
        args.time_grid,
        args.smoke_steps,
        args.smoke_label,
    )

    command = [
        sys.executable,
        str(ROOT / "run.py"),
        pulse_case,
        "--project",
        str(OUTPUT_PROJECT),
        "--mpi-procs",
        str(args.mpi_procs),
        "--elmer-solver",
        args.elmer_solver,
        "--runtime-bin",
        args.runtime_bin,
    ]
    if args.amgx_config:
        command += ["--amgx-config", args.amgx_config]
        command += ["--amgx-constraint-mode", args.amgx_constraint_mode]
        command += ["--amgx-constraint-penalty", str(args.amgx_constraint_penalty)]
    if args.dry_run:
        command.append("--dry-run")
    if args.force_deps:
        command.append("--force-deps")
    print("starting:", " ".join(command))
    return subprocess.run(command, cwd=ROOT).returncode


if __name__ == "__main__":
    raise SystemExit(main())
