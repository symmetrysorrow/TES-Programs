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
MUMPS_INNER_CIRCUIT_STEADY_REPART_X = "case_tes_steady_3x_mumps_inner_circuit_repart_x"
MUMPS_INNER_CIRCUIT_FAST = "case_tes_pulse_20ms_3x_mumps_inner_circuit_fast"
MUMPS_INNER_CIRCUIT_FAST_ALIGNED = "case_tes_pulse_20ms_3x_mumps_inner_circuit_pulse_aligned"
UMFPACK_STEADY_CASE = "case_tes_steady_3x_umfpack_recheck"


def geometric_ramp(
    t_start: float, t_end: float, entry_dt: float, max_growth: float
) -> np.ndarray:
    """Timestamps from *t_start* to *t_end*, each step at most *max_growth*
    times the previous, growing from *entry_dt*. Ignores any original
    intermediate timestamps in between and lands exactly on *t_end*, so both
    the entry and exit step-size ratios stay bounded regardless of how
    abruptly the original grid changed cadence inside the window.
    """
    total = t_end - t_start
    dts = []
    dt = entry_dt * max_growth
    remaining = total
    while remaining > dt * max_growth:
        dts.append(dt)
        remaining -= dt
        dt *= max_growth
    dts.append(remaining)
    times = [t_start]
    for d in dts:
        times.append(times[-1] + d)
    times[-1] = t_end
    return np.asarray(times)


def main() -> None:
    project = json.loads(SOURCE_PROJECT.read_text(encoding="utf-8"))
    project["meshes"]["mesh_refined_3x_repart_x"] = {
        "geometry": "single_pixel",
        "dir": "mesh_refined_3x_repart_x",
        "recipe": {
            "generator": "validated geometric partition of mesh_refined_3x",
            "commands": [
                "ElmerGrid 2 2 mesh_refined_3x -partition 4 1 1 -out mesh_refined_3x_repart_x"
            ],
        },
        "notes": (
            "Validated 4-rank x partition; keeps nearly all non-conforming "
            "mortar pairs on the same rank."
        ),
    }
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

    # Solver grid for a fast comparison.  The pulse begins at 20.02 ms, so
    # first advance 20 us from the 20 ms steady restart, then make the 1 ns
    # deposition interval an explicit BDF step.  This guarantees non-zero
    # overlap in AbsorberWindowPulseHeatSource before switching to the faster
    # 10 ns--0.5 ms stages.
    fast_case = dict(case)
    fast_case["series_file"] = "tes_pulse_20ms_3x_fast_compare_series.csv"
    fast_case["timesteps"] = [
        ["1[ms]", 20], ["10[us]", 2], ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9],
        ["5[us]", 394], ["50[us]", 160], ["0.5[ms]", 340],
    ]
    fast_case["output_intervals"] = [20, 2, 1, 10, 9, 9, 50, 20, 20]
    fast_case["vtu"] = "after_simulation"
    # The completed steady field/state remains the restart point.  Avoid a
    # full-field checkpoint at every transient step in this speed experiment.
    fast_case.pop("output_result", None)
    fast_case["comparison_time_grid"] = {
        "mode": "fast_solver_grid; interpolate Elmer current to COMSOL timestamps",
        "steps": 945,
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

    # Production-oriented prototype: HeatSolver invokes the synchronized
    # circuit as a nonlinear pre-solver before every assembly sweep.
    inner_circuit_state_file = f"mesh_refined_3x/{MUMPS_INNER_CIRCUIT_STEADY}.state"
    inner_circuit_steady = dict(mumps_steady)
    inner_circuit_steady.pop("series_file", None)
    inner_circuit_steady["state_file"] = inner_circuit_state_file
    inner_circuit_steady["heat_source"] = "circuit_inner"
    inner_circuit_steady["solver"] = dict(mumps_steady["solver"])
    inner_circuit_steady["solver"]["nonlinear_convergence_tolerance"] = 1.0e-6
    project["cases"][MUMPS_INNER_CIRCUIT_STEADY] = inner_circuit_steady

    inner_circuit_steady_repart_x = dict(inner_circuit_steady)
    inner_circuit_steady_repart_x["mesh"] = "mesh_refined_3x_repart_x"
    inner_circuit_steady_repart_x["state_file"] = (
        f"mesh_refined_3x_repart_x/{MUMPS_INNER_CIRCUIT_STEADY_REPART_X}.state"
    )
    project["cases"][MUMPS_INNER_CIRCUIT_STEADY_REPART_X] = inner_circuit_steady_repart_x

    inner_circuit_fast = dict(fast_case)
    inner_circuit_fast["heat_source"] = "circuit_inner"
    inner_circuit_fast["restart_from"] = MUMPS_INNER_CIRCUIT_STEADY
    inner_circuit_fast["series_file"] = "tes_pulse_20ms_3x_mumps_inner_circuit_fast_series.csv"
    inner_circuit_fast["state_file"] = inner_circuit_state_file
    inner_circuit_fast["solver"] = dict(inner_circuit_steady["solver"])
    project["cases"][MUMPS_INNER_CIRCUIT_FAST] = inner_circuit_fast

    # Keep the completed, pulse-missed run intact and create an explicitly
    # named replacement for the quantitative COMSOL comparison.
    inner_circuit_fast_aligned = dict(inner_circuit_fast)
    inner_circuit_fast_aligned["series_file"] = "tes_pulse_20ms_3x_mumps_inner_circuit_pulse_aligned_series.csv"
    project["cases"][MUMPS_INNER_CIRCUIT_FAST_ALIGNED] = inner_circuit_fast_aligned

    # MPI/serial regression window: retain the 20.020 ms physical pulse time
    # and a pre-pulse baseline, but stop shortly after the expected peak.
    inner_circuit_regression = dict(inner_circuit_fast_aligned)
    inner_circuit_regression["mesh"] = "mesh_refined_3x_repart_x"
    # Use the same converged thermal field as the frozen direct reference.
    # The circuit state intentionally remains on the legacy T0 initialization
    # path because that reference did not persist a circuit state file.
    inner_circuit_regression["restart_from"] = SOURCE_STEADY_CASE
    inner_circuit_regression["restart_time"] = 0.0
    inner_circuit_regression.pop("state_file", None)
    inner_circuit_regression["series_file"] = "tes_mpi_legacy_regression_series.csv"
    inner_circuit_regression["iteration_series_file"] = "tes_mpi_legacy_regression_iterations.csv"
    inner_circuit_regression["solver"] = dict(project["cases"][SOURCE_CASE]["solver"])
    inner_circuit_regression["solver"]["linear_system"] = "mumps"
    # Match the frozen direct series exactly through, and slightly beyond,
    # its 20.520001 ms current minimum.  A different pre-pulse ramp or
    # post-pulse dt changes the nonlinear fixed point and is not a regression.
    inner_circuit_regression["timesteps"] = [
        ["1[ms]", 20], ["18[us]", 1], ["1[us]", 2], ["1[ns]", 1],
        ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9], ["10[us]", 9],
        ["100[us]", 6],
    ]
    inner_circuit_regression["output_intervals"] = [1] * len(inner_circuit_regression["timesteps"])
    inner_circuit_regression["vtu"] = False
    project["cases"]["case_tes_mpi_legacy_regression"] = inner_circuit_regression

    # MPI run on the exact COMSOL adaptive timestamps through the pulse peak.
    # Keeping the same repartitioned mesh and restart as the accepted legacy
    # regression isolates time-grid effects from mesh and circuit effects.
    comsol_grid_mpi = dict(inner_circuit_regression)
    comsol_grid_times_ms = [float(times_ms[0])]
    for _time_ms in times_ms[1:]:
        # COMSOL writes a dense cluster of round-off-level timestamps exactly at
        # the pulse switch. Keep the physical transition, but avoid asking
        # Elmer to take sub-microsecond steps that only reproduce duplicates.
        if _time_ms - comsol_grid_times_ms[-1] >= 1.0e-3:
            comsol_grid_times_ms.append(float(_time_ms))
    comsol_grid_mpi["series_file"] = "tes_mpi_comsol_grid_series.csv"
    comsol_grid_mpi["iteration_series_file"] = "tes_mpi_comsol_grid_iterations.csv"
    comsol_grid_mpi["timesteps"] = [
        [f"{dt:.17g}[s]", 1]
        for dt in np.diff(np.asarray(comsol_grid_times_ms)[np.asarray(comsol_grid_times_ms) <= 20.620001]) * 1.0e-3
    ]
    comsol_grid_mpi["output_intervals"] = [1] * len(comsol_grid_mpi["timesteps"])
    comsol_grid_mpi["comparison_time_grid"] = {
        "source": "docs/Single-Pixel.txt",
        "samples": int(np.count_nonzero(np.asarray(comsol_grid_times_ms) <= 20.620001)),
        "start_ms": float(comsol_grid_times_ms[0]),
        "end_ms": float(comsol_grid_times_ms[-1]),
        "minimum_step_ms": 1.0e-3,
        "description": "MPI BDF steps follow COMSOL timestamps through 20.62 ms, with round-off-level pulse duplicates filtered at 1 us.",
    }
    project["cases"]["case_tes_mpi_comsol_grid"] = comsol_grid_mpi

    # Same COMSOL-exact grid through pulse+500 us, then every 5th COMSOL
    # sample (~5x larger steps) out to the COMSOL end time.  The mostly
    # featureless settling tail does not need the full adaptive resolution,
    # so this shortens the run without touching the frozen fine-grid
    # comparison above.  A single full-field checkpoint (Output File) is
    # written at the fine/coarse boundary -- this transient case normally
    # persists no restart-capable state past the pre-pulse steady solve, so a
    # crash or interruption during the (short) coarse tail would otherwise
    # force redoing the whole ~500 us fine phase.
    CHECKPOINT_REL_US = 500.0
    PULSE_MS = 20.02
    checkpoint_ms = PULSE_MS + CHECKPOINT_REL_US * 1.0e-3
    grid_ms = np.asarray(comsol_grid_times_ms)
    grid_ms = grid_ms[grid_ms <= 20.620001]
    fine_ms = grid_ms[grid_ms <= checkpoint_ms]
    tail_ms = grid_ms[grid_ms > checkpoint_ms]
    tail_coarse_ms = tail_ms[4::5]
    if tail_coarse_ms.size == 0 or tail_coarse_ms[-1] != tail_ms[-1]:
        tail_coarse_ms = np.append(tail_coarse_ms, tail_ms[-1])
    coarse_tail_ms = np.concatenate([fine_ms, tail_coarse_ms])
    coarse_tail_steps_s = np.diff(coarse_tail_ms) * 1.0e-3
    # 0-based stage index of the step that lands exactly on fine_ms[-1].
    checkpoint_stage_index = len(fine_ms) - 2

    COARSE_TAIL_MPI = "case_tes_mpi_comsol_grid_coarse_tail"
    coarse_tail = dict(comsol_grid_mpi)
    coarse_tail["series_file"] = "tes_mpi_comsol_grid_coarse_tail_series.csv"
    coarse_tail["iteration_series_file"] = "tes_mpi_comsol_grid_coarse_tail_iterations.csv"
    coarse_tail["timesteps"] = [[f"{dt:.17g}[s]", 1] for dt in coarse_tail_steps_s]
    # Every stage has exactly one physical step (Timestep Intervals = 1), so
    # an Output Intervals value > 1 never fires within that stage; only the
    # checkpoint stage uses 1, giving exactly one saved snapshot.
    coarse_tail["output_intervals"] = [
        1 if i == checkpoint_stage_index else 999999
        for i in range(len(coarse_tail_steps_s))
    ]
    coarse_tail["output_result"] = True
    coarse_tail["comparison_time_grid"] = {
        "source": "docs/Single-Pixel.txt",
        "samples": int(len(coarse_tail_ms)),
        "start_ms": float(coarse_tail_ms[0]),
        "end_ms": float(coarse_tail_ms[-1]),
        "description": (
            "Matches COMSOL timestamps exactly through pulse+500 us, then "
            "keeps every 5th COMSOL sample to the COMSOL end time. A "
            f"full-field checkpoint is written at stage {checkpoint_stage_index} "
            f"(pulse+{CHECKPOINT_REL_US:g} us)."
        ),
    }
    project["cases"][COARSE_TAIL_MPI] = coarse_tail

    # Diagnostic: does a genuinely small, uniform absolute step size suppress
    # the pulse+10-25 ms ringing seen in both case_tes_mpi_comsol_grid_full
    # and its smoothed variant? Neither run's schedule changed abruptly
    # there (COMSOL's own cadence, or a max-1.2x geometric ramp, both grow
    # smoothly), yet the ringing persists with similar amplitude in both, so
    # it looks like a BDF1 stability/accuracy threshold tied to absolute step
    # size rather than to how smoothly the step size changes. Reuses the
    # pulse+500 us checkpoint already on disk (no need to redo the fine
    # phase) and covers pulse+500 us to pulse+25 ms at a uniform 50 us --
    # short and cheap on purpose, to confirm the fix before committing to a
    # full pulse+500 us-180 ms rerun at this resolution.
    RINGING_DIAG_MPI = "case_tes_mpi_comsol_grid_ringing_diag"
    ringing_diag = dict(coarse_tail)
    ringing_diag.pop("output_result", None)
    ringing_diag.pop("restart_time", None)
    ringing_diag["series_file"] = "tes_mpi_comsol_grid_ringing_diag_series.csv"
    ringing_diag["iteration_series_file"] = "tes_mpi_comsol_grid_ringing_diag_iterations.csv"
    ringing_diag["timesteps"] = [["50e-6[s]", 490]]
    ringing_diag["output_intervals"] = [1]
    ringing_diag["restart_from"] = COARSE_TAIL_MPI
    ringing_diag["comparison_time_grid"] = {
        "description": (
            "Restarts from the pulse+500 us checkpoint; 490 uniform 50 us "
            "steps to pulse+25 ms (ignores COMSOL's own sample times in "
            "this stretch), to check whether a small fixed step size removes "
            "the ringing seen with COMSOL-cadence or geometric-ramp steps."
        ),
    }
    project["cases"][RINGING_DIAG_MPI] = ringing_diag

    # Same diagnostic at 100 us, to see if a coarser (cheaper) uniform step
    # still suppresses the ringing well enough to use for the full
    # pulse+500 us-180 ms rerun instead of paying for 50 us throughout.
    RINGING_DIAG_100US_MPI = "case_tes_mpi_comsol_grid_ringing_diag_100us"
    ringing_diag_100us = dict(ringing_diag)
    ringing_diag_100us["series_file"] = "tes_mpi_comsol_grid_ringing_diag_100us_series.csv"
    ringing_diag_100us["iteration_series_file"] = (
        "tes_mpi_comsol_grid_ringing_diag_100us_iterations.csv"
    )
    ringing_diag_100us["timesteps"] = [["100e-6[s]", 245]]
    ringing_diag_100us["comparison_time_grid"] = dict(ringing_diag["comparison_time_grid"])
    project["cases"][RINGING_DIAG_100US_MPI] = ringing_diag_100us

    # Full pulse+500 us-180 ms reruns at each candidate resolution, both
    # restarting from the same pulse+500 us checkpoint (fine phase reused,
    # not recomputed). Ignores COMSOL's own sample times throughout this
    # stretch in favor of a uniform step -- pick whichever the ringing
    # diagnostics justify. A single checkpoint is written at pulse+90 ms
    # (roughly the midpoint) so a crash in either multi-hour run only costs
    # the second half.
    for label, dt_us, n_steps in (("100us", 100, 1795), ("50us", 50, 3590)):
        case_name = f"case_tes_mpi_comsol_grid_full_uniform_{label}"
        case = dict(coarse_tail)
        case.pop("restart_time", None)
        case["series_file"] = f"tes_mpi_comsol_grid_full_uniform_{label}_series.csv"
        case["iteration_series_file"] = f"tes_mpi_comsol_grid_full_uniform_{label}_iterations.csv"
        case["timesteps"] = [[f"{dt_us}e-6[s]", n_steps]]
        checkpoint_step = n_steps // 2
        case["output_intervals"] = [checkpoint_step]
        case["output_result"] = True
        case["restart_from"] = COARSE_TAIL_MPI
        case["comparison_time_grid"] = {
            "description": (
                f"Restarts from the pulse+500 us checkpoint; {n_steps} uniform "
                f"{dt_us} us steps to pulse+{500 + n_steps * dt_us:g} us "
                "(ignores COMSOL's own sample times throughout), to remove the "
                f"BDF1 ringing confirmed by case_tes_mpi_comsol_grid_ringing_diag"
                f"{'_100us' if label == '100us' else ''}. Checkpoint at step "
                f"{checkpoint_step}."
            ),
        }
        project["cases"][case_name] = case

    # Single continuous run combining both fixes with no restart boundary at
    # all: the case_tes_mpi_comsol_grid_full_uniform_100us cases above still
    # restart from the pulse+500 us checkpoint, and merging their series with
    # coarse_tail's own pulse+0-500 us series showed a small seam artifact at
    # that boundary (same root cause as the earlier checkpoint-boundary
    # humps: the TES circuit UDF's Aitken relaxation / previous-current
    # history lives in Fortran module SAVE variables that a fresh restarted
    # process does not inherit). This case instead restarts once from the
    # pre-pulse steady state (matching case_tes_mpi_comsol_grid_full) and
    # runs pulse+0-500 us on COMSOL's exact grid immediately followed by
    # pulse+500 us-180 ms at a uniform 100 us -- one process, no seam. A
    # checkpoint is still written mid-tail for crash resilience; unlike a
    # restart_from chain, an unused checkpoint does not affect this run's
    # own output.
    FULL_UNIFORM_CONTINUOUS_MPI = "case_tes_mpi_comsol_grid_full_uniform_continuous"
    continuous_tail_n_steps = 1795
    continuous_case = dict(inner_circuit_regression)
    continuous_case["series_file"] = "tes_mpi_comsol_grid_full_uniform_continuous_series.csv"
    continuous_case["iteration_series_file"] = (
        "tes_mpi_comsol_grid_full_uniform_continuous_iterations.csv"
    )
    continuous_case["timesteps"] = [
        [f"{dt:.17g}[s]", 1] for dt in np.diff(fine_ms) * 1.0e-3
    ] + [["100e-6[s]", continuous_tail_n_steps]]
    fine_stage_count = len(fine_ms) - 1
    checkpoint_step = continuous_tail_n_steps // 2
    continuous_case["output_intervals"] = [999999] * fine_stage_count + [checkpoint_step]
    continuous_case["output_result"] = True
    continuous_case["comparison_time_grid"] = {
        "description": (
            "Single continuous run (no restart mid-way): COMSOL's exact grid "
            "from t=0 through pulse+500 us, then "
            f"{continuous_tail_n_steps} uniform 100 us steps to pulse+"
            f"{500 + continuous_tail_n_steps * 100:g} us. Combines the "
            "ringing fix (case_tes_mpi_comsol_grid_ringing_diag_100us) with "
            "the checkpoint-seam fix in one run."
        ),
    }
    project["cases"][FULL_UNIFORM_CONTINUOUS_MPI] = continuous_case

    # Resume-from-checkpoint case: only the remaining coarse tail, restarting
    # from the single snapshot saved above. Use this if the coarse tail alone
    # needs to be redone after an interruption.
    COARSE_TAIL_RESUME_MPI = "case_tes_mpi_comsol_grid_coarse_tail_resume"
    resume_case = dict(coarse_tail)
    resume_case.pop("output_result", None)
    # Continue from the checkpoint's own embedded physical time (~pulse+500 us)
    # instead of resetting the clock to 0 like the steady->transient restarts
    # above.
    resume_case.pop("restart_time", None)
    resume_case["series_file"] = "tes_mpi_comsol_grid_coarse_tail_resume_series.csv"
    resume_case["iteration_series_file"] = "tes_mpi_comsol_grid_coarse_tail_resume_iterations.csv"
    resume_case["timesteps"] = [
        [f"{dt:.17g}[s]", 1] for dt in coarse_tail_steps_s[checkpoint_stage_index + 1 :]
    ]
    resume_case["output_intervals"] = [999999] * len(resume_case["timesteps"])
    resume_case["restart_from"] = COARSE_TAIL_MPI
    project["cases"][COARSE_TAIL_RESUME_MPI] = resume_case

    # Extend the resume case all the way to COMSOL's own end time (pulse+180 ms,
    # 200 ms absolute), restarting from the same pulse+500 us checkpoint so the
    # expensive fine-grid phase already computed above is not repeated. COMSOL
    # samples every ~5-50 us out to a few ms and only then relaxes to ~5 ms
    # steps; following every sample here would add ~1000 more MUMPS solves
    # (many hours). The post-transient current is essentially settled by
    # pulse+600 us, so every 10th COMSOL sample keeps the far tail in ~100
    # steps while still landing on real COMSOL timestamps. A second checkpoint
    # is written where the far tail begins (pulse+600 us) so a crash during
    # the long tail does not require repeating the pulse+500-600 us segment.
    FULL_TAIL_MPI = "case_tes_mpi_comsol_grid_full_tail"
    far_tail_ms = grid_ms_full = np.asarray(comsol_grid_times_ms)
    far_tail_ms = far_tail_ms[far_tail_ms > 20.620001]
    far_tail_coarse_ms = far_tail_ms[9::10]
    if far_tail_coarse_ms.size == 0 or far_tail_coarse_ms[-1] != far_tail_ms[-1]:
        far_tail_coarse_ms = np.append(far_tail_coarse_ms, far_tail_ms[-1])
    full_tail_ms = np.concatenate([coarse_tail_ms[checkpoint_stage_index + 1 :], far_tail_coarse_ms])
    full_tail_steps_s = np.diff(full_tail_ms) * 1.0e-3
    # index of the step landing on pulse+600 us (end of the already-validated
    # coarse_tail schedule, start of the newly added far tail). The reused
    # point array coarse_tail_ms[checkpoint_stage_index + 1 :] has
    # (len(coarse_tail_ms) - checkpoint_stage_index - 1) points, i.e. one
    # fewer step than that, and the boundary step is the last one in it.
    far_checkpoint_stage_index = len(coarse_tail_ms) - checkpoint_stage_index - 3

    full_tail = dict(resume_case)
    full_tail["series_file"] = "tes_mpi_comsol_grid_full_tail_series.csv"
    full_tail["iteration_series_file"] = "tes_mpi_comsol_grid_full_tail_iterations.csv"
    full_tail["timesteps"] = [[f"{dt:.17g}[s]", 1] for dt in full_tail_steps_s]
    full_tail["output_intervals"] = [
        1 if i == far_checkpoint_stage_index else 999999
        for i in range(len(full_tail_steps_s))
    ]
    full_tail["output_result"] = True
    full_tail["comparison_time_grid"] = {
        "source": "docs/Single-Pixel.txt",
        "samples": int(len(full_tail_ms)) + 1,
        "start_ms": float(checkpoint_ms),
        "end_ms": float(full_tail_ms[-1]),
        "description": (
            "Restarts from the pulse+500 us checkpoint, replays the pulse+"
            "500-600 us coarse tail, then keeps every 10th COMSOL sample out "
            "to COMSOL's own end time (pulse+180 ms). A full-field checkpoint "
            f"is written at stage {far_checkpoint_stage_index} (pulse+600 us)."
        ),
    }
    project["cases"][FULL_TAIL_MPI] = full_tail

    # Single continuous, no-restart run following every COMSOL timestamp from
    # t=0 to COMSOL's own end time (200 ms absolute). The checkpointed
    # coarse_tail/full_tail chain above showed small spurious humps exactly
    # at its restart boundary and in its most aggressively downsampled
    # region: the TES circuit UDF keeps history (last committed current,
    # Aitken relaxation state) in Fortran module SAVE variables that a freshly
    # restarted process does not inherit, and downsampling also loses some of
    # the decay's curvature. One uninterrupted process avoids both causes, at
    # the cost of no crash-resilience -- roughly 1262 steps, an estimated 5-7
    # hours based on the coarse_tail/full_tail per-step timings.
    FULL_MPI = "case_tes_mpi_comsol_grid_full"
    full_mpi = dict(inner_circuit_regression)
    full_mpi["series_file"] = "tes_mpi_comsol_grid_full_series.csv"
    full_mpi["iteration_series_file"] = "tes_mpi_comsol_grid_full_iterations.csv"
    full_grid_ms = np.asarray(comsol_grid_times_ms)
    full_mpi["timesteps"] = [[f"{dt:.17g}[s]", 1] for dt in np.diff(full_grid_ms) * 1.0e-3]
    full_mpi["output_intervals"] = [1] * len(full_mpi["timesteps"])
    full_mpi["comparison_time_grid"] = {
        "source": "docs/Single-Pixel.txt",
        "samples": int(len(full_grid_ms)),
        "start_ms": float(full_grid_ms[0]),
        "end_ms": float(full_grid_ms[-1]),
        "minimum_step_ms": 1.0e-3,
        "description": (
            "Single continuous run (no restart checkpoints) following every "
            "COMSOL timestamp from t=0 to COMSOL's own end time (200 ms "
            "absolute), with round-off-level pulse duplicates filtered at "
            "1 us. Avoids the restart-boundary and coarse-tail artifacts "
            "seen in the checkpointed coarse_tail/full_tail chain."
        ),
    }
    project["cases"][FULL_MPI] = full_mpi

    # Same full run, but with pulse+9.98-19.98 ms smoothed. COMSOL's own grid
    # repeats several consecutive ~2x step-size doublings there (10 us ->
    # ... -> 1 ms over ~17 samples) while the current is still actively
    # decaying; BDF1 is comparatively sensitive to that compounding growth
    # and case_tes_mpi_comsol_grid_full shows a small wiggle there that
    # COMSOL's own curve does not have. Replace just that window with a
    # geometric ramp (each step <=1.2x the previous) between the same
    # boundary timestamps, ignoring the original intermediate COMSOL samples
    # in between -- an earlier version kept subdividing those intermediate
    # gaps and, because growth was capped per *original* interval rather than
    # per synthesized step, needed ~150 steps and still only reached ~125 us
    # by the window's end, leaving a fresh 8x jump into the next untouched
    # 1000 us step. A clean ramp reaches ~1.8 ms by the end (28 steps) and
    # exits into that same 1000 us cadence as a step-size *decrease*.
    FULL_SMOOTH_MPI = "case_tes_mpi_comsol_grid_full_smooth"
    smooth_window_start_ms = PULSE_MS + 9.980
    smooth_window_end_ms = PULSE_MS + 19.980
    full_smooth_ms = np.asarray(comsol_grid_times_ms)
    win_idx = np.where(
        (full_smooth_ms >= smooth_window_start_ms - 1.0e-9)
        & (full_smooth_ms <= smooth_window_end_ms + 1.0e-9)
    )[0]
    entry_dt_ms = full_smooth_ms[win_idx[0]] - full_smooth_ms[win_idx[0] - 1]
    ramp = geometric_ramp(
        full_smooth_ms[win_idx[0]], full_smooth_ms[win_idx[-1]], entry_dt_ms, max_growth=1.2
    )
    full_smooth_ms = np.concatenate(
        [full_smooth_ms[: win_idx[0]], ramp, full_smooth_ms[win_idx[-1] + 1 :]]
    )
    full_smooth = dict(inner_circuit_regression)
    full_smooth["series_file"] = "tes_mpi_comsol_grid_full_smooth_series.csv"
    full_smooth["iteration_series_file"] = "tes_mpi_comsol_grid_full_smooth_iterations.csv"
    full_smooth["timesteps"] = [
        [f"{dt:.17g}[s]", 1] for dt in np.diff(full_smooth_ms) * 1.0e-3
    ]
    full_smooth["output_intervals"] = [1] * len(full_smooth["timesteps"])
    full_smooth["comparison_time_grid"] = {
        "source": "docs/Single-Pixel.txt",
        "samples": int(len(full_smooth_ms)),
        "start_ms": float(full_smooth_ms[0]),
        "end_ms": float(full_smooth_ms[-1]),
        "description": (
            "Same as case_tes_mpi_comsol_grid_full, with pulse+9.98-19.98 ms "
            "replaced by a max-1.2x-growth geometric ramp between the same "
            "boundary timestamps (ignoring COMSOL's intermediate samples in "
            "that window), to remove the BDF1 step-doubling wiggle seen "
            "there without introducing a new jump at the window exit."
        ),
    }
    project["cases"][FULL_SMOOTH_MPI] = full_smooth

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
    # The regression probe deliberately uses one circuit solver and the same
    # nonlinear iterations as the serial reference.  Duplicating circuit
    # solvers changes the update order and is therefore not a valid MPI/serial
    # equivalence test.
    parallel_circuit_probe["parallel_circuit_iterations"] = 1
    parallel_circuit_probe["solver"]["nonlinear_max_iterations"] = 25
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
