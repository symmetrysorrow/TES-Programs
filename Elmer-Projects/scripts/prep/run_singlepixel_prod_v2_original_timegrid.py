"""Run the optimized SinglePixel mesh with the original COMSOL-grid tail.

This keeps the production-v2 spatial mesh and physics, but replaces its
provisional 10 us / 100 us / 1 ms tail schedule with the exact timestep
schedule used by ``case_tes_mpi_comsol_grid_full_uniform_continuous``.

One-command usage from the repository root::

    python scripts/prep/run_singlepixel_prod_v2_original_timegrid.py

The default linear solver is the MPI-safe HYPRE/BoomerAMG iterative solver.
Use ``--linear-system mumps`` only for a deliberate MUMPS retry.
"""

from __future__ import annotations

import argparse
import copy
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
KNOWN_SERIAL_STEADY = "case_tes_steady_singlepixel_prod_v2"

DEFAULT_SOLVER = Path(
    r"D:\Github\TES-Programs\tools\elmer-hypre\install-phase13-step-commit\bin\ElmerSolver.exe"
)
DEFAULT_RUNTIME_BIN = Path(r"C:\msys64\ucrt64\bin")


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
) -> None:
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
        "state_file": f"{mesh_dir}/{STEADY}.state",
        "series_file": f"{STEADY}_series.csv",
        "iteration_series_file": f"{STEADY}_iterations.csv",
        "output_file_path": f"../work/meshes/{mesh_dir}/{STEADY}.result",
    }
    project["cases"][STEADY] = steady

    original = copy.deepcopy(
        project["cases"]["case_tes_mpi_comsol_grid_full_uniform_continuous"]
    )
    original.update(
        {
            "mesh": mesh_name,
            "restart_from": STEADY,
            "restart_time": 0.02,
            # MPI writes restart files as <base>.result.0, .1, ... .  The
            # path must point into the active mesh directory; otherwise
            # Elmer looks for the files in the repository root.
            "restart_file_path": f"../work/meshes/{mesh_dir}/{STEADY}.result",
            "state_file": f"{mesh_dir}/{STEADY}.state",
            "series_file": f"{PULSE}_series.csv",
            "iteration_series_file": f"{PULSE}_iterations.csv",
            "output_file_path": f"../work/meshes/{mesh_dir}/{PULSE}.result",
            "vtu": False,
            "comparison_time_grid": {
                "source_case": "case_tes_mpi_comsol_grid_full_uniform_continuous",
                "source_project": SOURCE_PROJECT.name,
                "description": (
                    "Optimized SinglePixel mesh with the original COMSOL-grid "
                    "timesteps, including the original tail."
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
        original["restart_file_base"] = STEADY
    elif reuse_known_serial_steady:
        # This result was produced successfully with the same production-v2
        # mesh and physics.  It is a serial restart, so this mode is only
        # valid when the pulse itself is also run with one MPI rank.
        original["restart_from"] = None
        original["restart_file_base"] = KNOWN_SERIAL_STEADY
        original["restart_file_path"] = (
            f"../work/meshes/{BASE_MESH}/{KNOWN_SERIAL_STEADY}.result"
        )
    project["cases"][PULSE] = original
    OUTPUT_PROJECT.write_text(
        json.dumps(project, indent=2) + "\n", encoding="utf-8"
    )
    print(f"prepared {OUTPUT_PROJECT.relative_to(ROOT)}")
    print(f"mesh={mesh_dir}, steady={STEADY}, pulse={PULSE}")
    print(f"timestep stages={len(original['timesteps'])}")


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
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.mpi_procs < 1:
        parser.error("--mpi-procs must be >= 1")
    if args.reuse_steady and args.reuse_known_serial_steady:
        parser.error("--reuse-steady and --reuse-known-serial-steady are mutually exclusive")
    if args.reuse_known_serial_steady and args.mpi_procs != 1:
        parser.error("--reuse-known-serial-steady requires --mpi-procs 1")
    if not SOURCE_PROJECT.is_file():
        raise FileNotFoundError(SOURCE_PROJECT)

    mesh_dir = ensure_partition(args.mpi_procs)
    make_project(
        mesh_dir,
        args.linear_system,
        args.reuse_steady,
        args.reuse_known_serial_steady,
    )

    command = [
        sys.executable,
        str(ROOT / "run.py"),
        PULSE,
        "--project",
        str(OUTPUT_PROJECT),
        "--mpi-procs",
        str(args.mpi_procs),
        "--elmer-solver",
        args.elmer_solver,
        "--runtime-bin",
        args.runtime_bin,
    ]
    if args.dry_run:
        command.append("--dry-run")
    print("starting:", " ".join(command))
    return subprocess.run(command, cwd=ROOT).returncode


if __name__ == "__main__":
    raise SystemExit(main())
