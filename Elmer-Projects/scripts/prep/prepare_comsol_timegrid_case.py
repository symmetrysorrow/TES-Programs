"""Create a non-destructive Elmer case that uses COMSOL's exact time grid.

The COMSOL table is an adaptive-time export.  Elmer's BDF solver advances one
step per ``timesteps`` item, so each consecutive COMSOL time difference is
emitted as a one-step stage.  This makes the current CSV samples coincide with
the COMSOL timestamps while preserving the existing pulse case unchanged.

Usage:
    python scripts/prep/prepare_comsol_timegrid_case.py
    python run.py case_tes_pulse_20ms_3x_comsol_grid --project elmer_project_comsol_timegrid.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PROJECT = ROOT / "elmer_project.json"
COMSOL_TABLE = ROOT / "docs" / "Single-Pixel.txt"
OUTPUT_PROJECT = ROOT / "elmer_project_comsol_timegrid.json"
SOURCE_CASE = "case_tes_pulse_20ms_3x_refined"
SOURCE_STEADY_CASE = "case_tes_steady_3x_refined"
STEADY_CASE = "case_tes_steady_3x_comsol_grid"
TARGET_CASE = "case_tes_pulse_20ms_3x_comsol_grid_from_steady"
FAST_TARGET_CASE = "case_tes_pulse_20ms_3x_fast_compare"
HYPRE_STEADY_CASE = "case_tes_steady_3x_hypre_boomeramg"
HYPRE_FAST_CASE = "case_tes_pulse_20ms_3x_fast_hypre_boomeramg"
MUMPS_STEADY_CASE = "case_tes_steady_3x_mumps_mpi"
MUMPS_FAST_CASE = "case_tes_pulse_20ms_3x_fast_mumps_mpi"
MUMPS_PARALLEL_CIRCUIT_PILOT = "case_tes_steady_3x_mumps_parallel_circuit_pilot"
MUMPS_PARALLEL_CIRCUIT_PROBE = "case_tes_pulse_20ms_3x_parallel_circuit_probe"
MUMPS_INNER_CIRCUIT_STEADY = "case_tes_steady_3x_mumps_inner_circuit"
MUMPS_INNER_CIRCUIT_FAST = "case_tes_pulse_20ms_3x_mumps_inner_circuit_fast"
UMFPACK_STEADY_CASE = "case_tes_steady_3x_umfpack_recheck"


def main() -> None:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    times_ms = np.loadtxt(COMSOL_TABLE, comments="%", encoding="utf-8", usecols=0)
    if times_ms.ndim != 1 or len(times_ms) < 2 or not np.all(np.diff(times_ms) > 0):
        raise ValueError("COMSOL time column must be a strictly increasing series with at least two rows")

    # Run a dedicated steady case first.  It writes both the thermal restart
    # field and the converged UDF circuit state used by the pulse case.
    state_file = f"mesh_refined_3x/{STEADY_CASE}.state"
    steady = dict(project["cases"][SOURCE_STEADY_CASE])
    steady["output_result"] = True
    steady["state_file"] = state_file
    steady["series_file"] = "tes_steady_3x_comsol_grid_series.csv"
    project["cases"][STEADY_CASE] = steady

    case = dict(project["cases"][SOURCE_CASE])
    steps_s = np.diff(times_ms) * 1e-3
    # A separate one-step stage keeps the sequence of adaptive COMSOL times;
    # setting Output Intervals=1 also makes Elmer's standard time output exact.
    case["timesteps"] = [[f"{dt:.17g}[s]", 1] for dt in steps_s]
    case["output_intervals"] = [1] * len(steps_s)
    case["series_file"] = "tes_pulse_20ms_3x_comsol_grid_series.csv"
    case["output_result"] = True
    case["state_file"] = state_file
    case["restart_from"] = STEADY_CASE
    # The current CSV is the comparison deliverable.  Prevent 1,614 full-field
    # VTU files (tens of GB) while retaining a final field snapshot for QA.
    case["vtu"] = "after_simulation"
    case["comparison_time_grid"] = {
        "source": "docs/Single-Pixel.txt",
        "samples": int(len(times_ms)),
        "start_ms": float(times_ms[0]),
        "end_ms": float(times_ms[-1]),
        "description": "Each Elmer BDF step ends at the corresponding COMSOL probe-table time.",
    }
    project["cases"][TARGET_CASE] = case

    # Solver grid for a fast comparison.  It resolves the 1 ns deposition
    # window exactly, uses 5 us through the current-rise region, then expands
    # to 50 us and 0.5 ms.  COMSOL values are compared by interpolation in
    # post-processing rather than forcing Elmer to take COMSOL's tiny adaptive
    # steps (down to ~1e-13 s).
    fast_case = dict(case)
    fast_case["series_file"] = "tes_pulse_20ms_3x_fast_compare_series.csv"
    fast_case["timesteps"] = [
        # The 1 ns deposition step must immediately follow t=20 ms.  Putting
        # a coarse step first makes the overlap-based pulse UDF skip the
        # entire 1 ns energy deposition.
        ["1[ms]", 20], ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9],
        ["5[us]", 394], ["50[us]", 160], ["0.5[ms]", 340],
    ]
    fast_case["output_intervals"] = [20, 1, 10, 9, 9, 50, 20, 20]
    fast_case["vtu"] = "after_simulation"
    # The completed steady field/state remains the restart point.  Avoid a
    # full-field checkpoint at every transient step in this speed experiment.
    fast_case.pop("output_result", None)
    fast_case["comparison_time_grid"] = {
        "mode": "fast_solver_grid; interpolate Elmer current to COMSOL timestamps",
        "steps": 946,
        "smallest_step": "1 ns",
        "pulse_rise_step": "5 us",
    }
    project["cases"][FAST_TARGET_CASE] = fast_case

    # HYPRE/BoomerAMG speed-test pair.  The separate steady result/state
    # avoids mixing a direct-solver restart with the iterative experiment.
    hypre_state_file = f"mesh_refined_3x/{HYPRE_STEADY_CASE}.state"
    hypre_steady = dict(steady)
    hypre_steady["series_file"] = "tes_steady_3x_hypre_boomeramg_series.csv"
    hypre_steady["state_file"] = hypre_state_file
    hypre_steady["solver"] = dict(hypre_steady["solver"])
    hypre_steady["solver"]["linear_system"] = "iterative_hypre_boomeramg"
    hypre_steady["solver"]["linear_solver_backend"] = "HYPRE BoomerAMG"
    project["cases"][HYPRE_STEADY_CASE] = hypre_steady

    hypre_fast = dict(fast_case)
    hypre_fast["series_file"] = "tes_pulse_20ms_3x_fast_hypre_boomeramg_series.csv"
    hypre_fast["state_file"] = hypre_state_file
    hypre_fast["restart_from"] = HYPRE_STEADY_CASE
    hypre_fast["solver"] = dict(hypre_fast["solver"])
    hypre_fast["solver"]["linear_system"] = "iterative_hypre_boomeramg"
    hypre_fast["solver"]["linear_solver_backend"] = "HYPRE BoomerAMG"
    hypre_fast["comparison_time_grid"] = {
        **fast_case["comparison_time_grid"],
        "linear_solver": "BiCGStab + HYPRE BoomerAMG",
    }
    project["cases"][HYPRE_FAST_CASE] = hypre_fast

    # Parallel MUMPS retains a direct solve while distributing the matrix
    # across MPI ranks.  It therefore provides the closest result-preserving
    # speed comparison against the existing serial UMFPACK setup.
    mumps_state_file = f"mesh_refined_3x/{MUMPS_STEADY_CASE}.state"
    mumps_steady = dict(steady)
    mumps_steady["series_file"] = "tes_steady_3x_mumps_mpi_series.csv"
    mumps_steady["state_file"] = mumps_state_file
    mumps_steady["solver"] = dict(mumps_steady["solver"])
    mumps_steady["solver"]["linear_system"] = "mumps"
    # The parallel partitioned temperature update oscillates at ~1e-6 while
    # the TES current agrees with serial UMFPACK to 0.005%.  A 1e-5 nonlinear
    # criterion is therefore sufficient for the current-comparison workflow.
    mumps_steady["solver"]["nonlinear_convergence_tolerance"] = 1.0e-5
    project["cases"][MUMPS_STEADY_CASE] = mumps_steady

    # Prototype: the circuit is evaluated by a synchronized MPI solver once
    # per outer steady iteration, while the body-force UDF merely reads its
    # common power.  This avoids collective calls from element callbacks.
    parallel_circuit_pilot = dict(mumps_steady)
    parallel_circuit_pilot.pop("series_file", None)
    parallel_circuit_pilot.pop("state_file", None)
    parallel_circuit_pilot["heat_source"] = "circuit_parallel"
    parallel_circuit_pilot["steady_state_max_iterations"] = 30
    parallel_circuit_pilot["output_intervals"] = 30
    parallel_circuit_pilot["solver"] = dict(mumps_steady["solver"])
    parallel_circuit_pilot["solver"]["nonlinear_max_iterations"] = 8
    parallel_circuit_pilot["solver"]["nonlinear_convergence_tolerance"] = 1.0e-4
    project["cases"][MUMPS_PARALLEL_CIRCUIT_PILOT] = parallel_circuit_pilot

    # Production-oriented prototype: a single HeatSolve invocation with the
    # circuit update embedded in each of its nonlinear iterations.
    inner_circuit_steady = dict(mumps_steady)
    inner_circuit_steady.pop("series_file", None)
    inner_circuit_steady.pop("state_file", None)
    inner_circuit_steady["heat_source"] = "circuit_inner"
    inner_circuit_steady["solver"] = dict(mumps_steady["solver"])
    inner_circuit_steady["solver"]["nonlinear_convergence_tolerance"] = 1.0e-6
    project["cases"][MUMPS_INNER_CIRCUIT_STEADY] = inner_circuit_steady

    inner_circuit_fast = dict(fast_case)
    inner_circuit_fast["heat_source"] = "circuit_inner"
    inner_circuit_fast["restart_from"] = MUMPS_INNER_CIRCUIT_STEADY
    inner_circuit_fast["series_file"] = "tes_pulse_20ms_3x_mumps_inner_circuit_fast_series.csv"
    inner_circuit_fast.pop("state_file", None)
    inner_circuit_fast["solver"] = dict(inner_circuit_steady["solver"])
    project["cases"][MUMPS_INNER_CIRCUIT_FAST] = inner_circuit_fast

    # Short pulse-window validation of the synchronized circuit.  It reaches
    # 20.51 ms, covering the 20 ms deposition and the COMSOL current-rise
    # region, without committing to the full 200 ms comparison run.
    parallel_circuit_probe = dict(fast_case)
    parallel_circuit_probe["heat_source"] = "circuit_parallel"
    parallel_circuit_probe["restart_from"] = MUMPS_PARALLEL_CIRCUIT_PILOT
    parallel_circuit_probe["restart_time"] = 0.020
    parallel_circuit_probe["series_file"] = "tes_pulse_20ms_3x_parallel_circuit_probe_series.csv"
    parallel_circuit_probe["timesteps"] = [
        ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9], ["5[us]", 96],
    ]
    parallel_circuit_probe["output_intervals"] = [1, 10, 9, 9, 24]
    parallel_circuit_probe["vtu"] = False
    parallel_circuit_probe.pop("state_file", None)
    parallel_circuit_probe["solver"] = dict(parallel_circuit_pilot["solver"])
    # Four circuit/heat sweeps per time step provide an implicit-like
    # electrothermal coupling while keeping this diagnostic run tractable.
    parallel_circuit_probe["parallel_circuit_iterations"] = 4
    parallel_circuit_probe["solver"]["nonlinear_max_iterations"] = 1
    project["cases"][MUMPS_PARALLEL_CIRCUIT_PROBE] = parallel_circuit_probe

    mumps_fast = dict(fast_case)
    mumps_fast["series_file"] = "tes_pulse_20ms_3x_fast_mumps_mpi_series.csv"
    mumps_fast["state_file"] = mumps_state_file
    mumps_fast["restart_from"] = MUMPS_STEADY_CASE
    mumps_fast["solver"] = dict(mumps_fast["solver"])
    mumps_fast["solver"]["linear_system"] = "mumps"
    mumps_fast["solver"]["nonlinear_convergence_tolerance"] = 1.0e-5
    mumps_fast["comparison_time_grid"] = {
        **fast_case["comparison_time_grid"],
        "linear_solver": "parallel MUMPS direct",
    }
    project["cases"][MUMPS_FAST_CASE] = mumps_fast

    # Fresh serial direct-solver reference, isolated from prior restart and
    # UDF-state experiments.  This is used to verify the apparent 148 uA
    # baseline in the interrupted pulse CSV.
    umfpack_state_file = f"mesh_refined_3x/{UMFPACK_STEADY_CASE}.state"
    umfpack_steady = dict(steady)
    umfpack_steady["series_file"] = "tes_steady_3x_umfpack_recheck_series.csv"
    umfpack_steady["state_file"] = umfpack_state_file
    project["cases"][UMFPACK_STEADY_CASE] = umfpack_steady
    OUTPUT_PROJECT.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT_PROJECT}")
    print(f"{TARGET_CASE}: {len(steps_s)} BDF steps, {times_ms[0]:.9g}-{times_ms[-1]:.9g} ms")


if __name__ == "__main__":
    main()
