"""Create the small Phase24 production-path before/after smoke project."""
from __future__ import annotations

import copy
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_hypre_gpu_phase19.json"
OUTPUT = ROOT / "elmer_project_phase24_production.json"
BASE_CASE = "case_p19_hypre_flexgmres_boomeramg_cpu_time5us_smoke_7step"
CASE = "case_phase24_hypre_cpu_smoke_7step"
GPU_CASE = "case_phase24_hypre_gpu_smoke_7step"
ADAPTIVE_CASE = "case_phase24_adaptive_output_smoke"
ADAPTIVE_DEBUG_CASE = "case_phase24_adaptive_debug_strict_4us"
BDF2_SMOKE_CASE = "case_phase24_adaptive_bdf2_smoke_4us"
ADAPTIVE_POST_EVENT_CASE = "case_phase24_adaptive_post_event_cooldown_2ns"
ADAPTIVE_REJECTION_CASE = "case_phase24_adaptive_event_rejection_progression_2ns"


def main() -> None:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    base = project["cases"][BASE_CASE]
    candidate = copy.deepcopy(base)
    candidate["restart_from"] = None
    candidate["preexisting_restart"] = True
    candidate["series_file"] = f"{CASE}_series.csv"
    candidate["iteration_series_file"] = f"{CASE}_iterations.csv"
    candidate["output_file_path"] = (
        f"../work/meshes/{candidate['mesh']}/{CASE}.result"
    )
    candidate["phase24_smoke"] = {
        "purpose": "same production HYPRE CPU smoke with Phase24 opt-in",
        "baseline_case": BASE_CASE,
        "path": "scalar 3D P1 tetra, static/dynamic split, cached CSR insertion",
    }
    candidate["phase24_vector_assembly"] = True
    candidate["phase24_wall_profiling"] = True
    candidate["phase24_static_matrix_reuse"] = True
    # Stage 8 keeps A current and moves reuse to the AMG preconditioner.
    candidate["phase24_operator_lagging"] = "disabled"
    candidate["phase24_operator_lag_threshold"] = 1.0e-4
    candidate["phase24_operator_lag_max_reuse"] = 3
    candidate["phase24_operator_lag_residual_growth"] = 1.25
    candidate["phase24_preconditioner_lagging"] = "adaptive"
    candidate["phase24_preconditioner_lag_threshold"] = 1.0e-4
    candidate["phase24_preconditioner_lag_max_age"] = 2
    candidate["phase24_preconditioner_lag_krylov_relative_growth"] = 0.50
    candidate["phase24_preconditioner_lag_nonlinear_threshold"] = 1.0e-3
    candidate["phase24_hypre_reuse"] = True
    candidate["solver"] = dict(candidate["solver"])
    candidate["solver"]["matrix_dump_prefix"] = CASE
    candidate["solver"]["boomer_amg_strong_threshold"] = 0.5
    # Stage 4 uses the physical linear-solve policy accepted in Stage 3:
    # avoid spending the run budget on an artificial 1e-11 residual target.
    candidate["solver"]["linear_system_convergence_tolerance"] = 5e-7
    project["cases"][CASE] = candidate

    # Keep GPU selection explicit while sharing the exact Phase24 lifecycle
    # policy and production inputs.  The WSL runner chooses CUDA, HIP, or CPU;
    # this case only selects HYPRE's device-capable solver entry.
    gpu_candidate = copy.deepcopy(candidate)
    gpu_candidate["series_file"] = f"{GPU_CASE}_series.csv"
    gpu_candidate["iteration_series_file"] = f"{GPU_CASE}_iterations.csv"
    gpu_candidate["output_file_path"] = (
        f"../work/meshes/{gpu_candidate['mesh']}/{GPU_CASE}.result"
    )
    gpu_candidate["phase24_smoke"] = {
        "purpose": "same Phase24 production HYPRE lifecycle with device-resident linear algebra",
        "baseline_case": BASE_CASE,
        "path": "CPU FEM assembly, persistent IJ/ParCSR GPU linear algebra, explicit transfer counters",
        "backend_selection": "runner-selected CUDA/HIP/CPU",
    }
    gpu_candidate["phase24_hypre_backend"] = "device"
    gpu_candidate["solver"] = dict(gpu_candidate["solver"])
    gpu_candidate["solver"]["linear_system"] = "iterative_hypre_flexgmres_boomeramg_gpu"
    gpu_candidate["solver"]["matrix_dump_prefix"] = GPU_CASE
    project["cases"][GPU_CASE] = gpu_candidate

    # Stage 11 uses one outer interval for the physical window.  The native
    # driver chooses accepted internal steps inside it, while this schedule
    # is only sampled through dense output.  The count is intentionally a
    # case default, not a solver invariant; callers may replace the schedule
    # with 50, 500, 5000, or an explicit/nonuniform array.
    adaptive_candidate = copy.deepcopy(candidate)
    adaptive_candidate["series_file"] = f"{ADAPTIVE_CASE}_series.csv"
    adaptive_candidate["iteration_series_file"] = f"{ADAPTIVE_CASE}_iterations.csv"
    adaptive_candidate["output_file_path"] = (
        f"../work/meshes/{adaptive_candidate['mesh']}/{ADAPTIVE_CASE}.result"
    )
    adaptive_candidate["timesteps"] = [["31[us]", 1]]
    adaptive_candidate["output_intervals"] = [0]
    adaptive_candidate["bdf_order"] = 2
    adaptive_candidate["adaptive_time"] = {
        "start": "20[ms]",
        "end": "20.031[ms]",
        "requested_output_times": {
            "mode": "uniform",
            "start": "20[ms]",
            "end": "20.031[ms]",
            "count": 64,
        },
        "dt_initial": "1[us]",
        "dt_min": "0.5[ns]",
        "dt_max": "100[us]",
        # Native production calibration: the strict debug tolerance is useful
        # for exercising rollback, but the fixed-step reference shows that
        # 0.2% relative temporal control is the first viable CPU candidate.
        "relative_tolerance": 2.0e-3,
        "absolute_tolerance": 1.0e-8,
        "r_min": 0.5,
        "r_max": 2.0,
        "max_growth": 1.5,
        "max_shrink": 0.5,
        # The production window crosses two physical events.  Keep the
        # rejection budget finite, but large enough to cover both event
        # neighborhoods without changing the temporal error target.
        "max_rejected": 64,
        "physical_event_times": ["20.02[ms]", "20.020001[ms]"],
    }
    adaptive_candidate["phase24_smoke"] = {
        "purpose": "Stage 11 output/internal-step decoupling smoke",
        "reference_case": BASE_CASE,
        "path": "persistent Phase24 HYPRE lifecycle under adaptive BDF1/BDF2 stepping with dense output",
    }
    project["cases"][ADAPTIVE_CASE] = adaptive_candidate

    # Same production mesh/material/UDF path, shortened to two adaptive
    # intervals so the native debug gate reaches the first BDF2 trial.
    adaptive_debug = copy.deepcopy(adaptive_candidate)
    adaptive_debug["series_file"] = f"{ADAPTIVE_DEBUG_CASE}_series.csv"
    adaptive_debug["iteration_series_file"] = f"{ADAPTIVE_DEBUG_CASE}_iterations.csv"
    adaptive_debug["output_file_path"] = (
        f"../work/meshes/{adaptive_debug['mesh']}/{ADAPTIVE_DEBUG_CASE}.result"
    )
    adaptive_debug["timesteps"] = [["2[us]", 2]]
    adaptive_debug["adaptive_time"] = {
        "start": "20[ms]",
        "end": "20.004[ms]",
        "requested_output_times": {
            "mode": "uniform",
            "start": "20[ms]",
            "end": "20.004[ms]",
            "count": 5,
        },
        "dt_initial": "2[us]",
        "dt_min": "1[ns]",
        "dt_max": "2[us]",
        "relative_tolerance": 1.0e-4,
        "absolute_tolerance": 1.0e-8,
        "r_min": 0.5,
        "r_max": 2.0,
        "max_growth": 1.5,
        "max_shrink": 0.5,
        "max_rejected": 2,
        "debug": True,
    }
    adaptive_debug["phase24_smoke"] = {
        "purpose": "short same-mesh Stage 11 trial sequencing diagnostic",
        "reference_case": ADAPTIVE_CASE,
        "path": "same Phase24 HeatSolve/HYPRE/TES path, two outer microsecond intervals",
    }
    project["cases"][ADAPTIVE_DEBUG_CASE] = adaptive_debug

    # Reproducible production-tolerance BDF2 gate.  Keep this separate from
    # the strict temporal-tolerance diagnostic so the successful policy does
    # not depend on an unrecorded SIF override.
    adaptive_bdf2_smoke = copy.deepcopy(adaptive_debug)
    adaptive_bdf2_smoke["series_file"] = f"{BDF2_SMOKE_CASE}_series.csv"
    adaptive_bdf2_smoke["iteration_series_file"] = f"{BDF2_SMOKE_CASE}_iterations.csv"
    adaptive_bdf2_smoke["output_file_path"] = (
        f"../work/meshes/{adaptive_bdf2_smoke['mesh']}/{BDF2_SMOKE_CASE}.result"
    )
    adaptive_bdf2_smoke["adaptive_time"]["relative_tolerance"] = 2.0e-3
    adaptive_bdf2_smoke["phase24_smoke"] = {
        "purpose": "committed production-tolerance BDF2 cache-invalidation gate",
        "reference_case": ADAPTIVE_CASE,
        "path": "same production mesh/material/TES/HYPRE path, 4-us adaptive interval",
    }
    project["cases"][BDF2_SMOKE_CASE] = adaptive_bdf2_smoke

    # Bounded reproducible diagnostic that reaches the production event
    # neighborhood and retains a 2 ns tail after the second event.
    # The restart state is the same 20 ms steady state as production; the
    # short tail makes the post-event controller behavior inspectable without
    # changing the production tolerances, mesh, TES coupling, or HYPRE path.
    adaptive_post_event = copy.deepcopy(adaptive_candidate)
    adaptive_post_event["series_file"] = f"{ADAPTIVE_POST_EVENT_CASE}_series.csv"
    adaptive_post_event["iteration_series_file"] = f"{ADAPTIVE_POST_EVENT_CASE}_iterations.csv"
    adaptive_post_event["output_file_path"] = (
        f"../work/meshes/{adaptive_post_event['mesh']}/{ADAPTIVE_POST_EVENT_CASE}.result"
    )
    adaptive_post_event["output_result"] = True
    adaptive_post_event["timesteps"] = [["20.002[us]", 1]]
    adaptive_post_event["adaptive_time"] = {
        "start": "20[ms]",
        "end": "20.020002[ms]",
        "requested_output_times": {
            "mode": "explicit",
            "times": [
                "20[ms]",
                "20.02[ms]",
                "20.020001[ms]",
                "20.020002[ms]",
            ],
        },
        "dt_initial": "20[us]",
        "dt_min": "1[ns]",
        "dt_max": "20[us]",
        "relative_tolerance": 2.0e-3,
        "absolute_tolerance": 1.0e-8,
        "r_min": 0.5,
        "r_max": 2.0,
        "max_growth": 1.5,
        "max_shrink": 0.5,
        "max_rejected": 8,
        "physical_event_times": ["20.02[ms]", "20.020001[ms]"],
    }
    adaptive_post_event["phase24_smoke"] = {
        "purpose": "bounded post-event adaptive rejection/cooldown diagnostic",
        "reference_case": ADAPTIVE_CASE,
        "path": "same production restart, mesh/material/TES/HYPRE path; 2-ns tail after the second event",
    }
    project["cases"][ADAPTIVE_POST_EVENT_CASE] = adaptive_post_event

    # Bounded lifecycle regression for the production 0.5-ns floor.  This
    # reaches the same two physical events as production, enables the native
    # trial trace, and stops immediately after the 2-ns local tail.
    adaptive_rejection = copy.deepcopy(adaptive_candidate)
    adaptive_rejection["series_file"] = f"{ADAPTIVE_REJECTION_CASE}_series.csv"
    adaptive_rejection["iteration_series_file"] = f"{ADAPTIVE_REJECTION_CASE}_iterations.csv"
    adaptive_rejection["output_file_path"] = (
        f"../work/meshes/{adaptive_rejection['mesh']}/{ADAPTIVE_REJECTION_CASE}.result"
    )
    adaptive_rejection["output_result"] = True
    adaptive_rejection["timesteps"] = [["20.002[us]", 1]]
    adaptive_rejection["adaptive_time"] = {
        "start": "20[ms]",
        "end": "20.020002[ms]",
        "requested_output_times": {
            "mode": "explicit",
            "times": [
                "20[ms]",
                "20.02[ms]",
                "20.020001[ms]",
                "20.020002[ms]",
            ],
        },
        "dt_initial": "20[us]",
        "dt_min": "0.5[ns]",
        "dt_max": "20[us]",
        "relative_tolerance": 2.0e-3,
        "absolute_tolerance": 1.0e-8,
        "r_min": 0.5,
        "r_max": 2.0,
        "max_growth": 1.5,
        "max_shrink": 0.5,
        "max_rejected": 8,
        "physical_event_times": ["20.02[ms]", "20.020001[ms]"],
        "debug": True,
    }
    adaptive_rejection["phase24_smoke"] = {
        "purpose": "bounded post-event adaptive rejection lifecycle regression",
        "reference_case": ADAPTIVE_CASE,
        "path": "same production restart, mesh/material/TES/HYPRE path; 2-ns tail after the second event",
    }
    project["cases"][ADAPTIVE_REJECTION_CASE] = adaptive_rejection

    # Fixed-dt local scaling probes.  All probes restart from output position
    # 3 of the bounded case (the accepted state at 20.020001 ms), so their
    # embedded estimators are directly comparable across dt.
    local_dt_values_ns = (1.5, 1.25, 1.0, 0.875, 0.75, 0.625, 0.5)
    local_state_time = 0.020020001
    for dt_ns in local_dt_values_ns:
        tag = str(dt_ns).replace(".", "p")
        local_case_name = f"case_phase24_adaptive_local_dt_{tag}ns"
        local_dt = dt_ns * 1.0e-9
        local_end = local_state_time + local_dt
        local_case = copy.deepcopy(adaptive_post_event)
        local_case["restart_from"] = ADAPTIVE_POST_EVENT_CASE
        local_case["restart_file_path"] = (
            f"../work/meshes/{local_case['mesh']}/{ADAPTIVE_POST_EVENT_CASE}.result"
        )
        local_case["restart_position"] = 3
        local_case["restart_time"] = local_state_time
        local_case["series_file"] = f"{local_case_name}_series.csv"
        local_case["iteration_series_file"] = f"{local_case_name}_iterations.csv"
        local_case["output_file_path"] = (
            f"../work/meshes/{local_case['mesh']}/{local_case_name}.result"
        )
        local_case["timesteps"] = [[f"{dt_ns:g}[ns]", 1]]
        local_case["adaptive_time"] = {
            "start": f"{local_state_time * 1000:.12f}[ms]",
            "end": f"{local_end * 1000:.12f}[ms]",
            "requested_output_times": {
                "mode": "explicit",
                "times": [
                    f"{local_state_time * 1000:.12f}[ms]",
                    f"{local_end * 1000:.12f}[ms]",
                ],
            },
            "dt_initial": f"{dt_ns:g}[ns]",
            "dt_min": f"{dt_ns:g}[ns]",
            "dt_max": f"{dt_ns:g}[ns]",
            "relative_tolerance": 2.0e-3,
            "absolute_tolerance": 1.0e-8,
            "r_min": 0.5,
            "r_max": 2.0,
            "max_growth": 1.0,
            "max_shrink": 0.5,
            "max_rejected": 2,
            "physical_event_times": [],
            "debug": True,
        }
        local_case["phase24_smoke"] = {
            "purpose": "same-state adaptive estimator local dt scaling probe",
            "reference_case": ADAPTIVE_POST_EVENT_CASE,
            "restart_position": 3,
            "dt_ns": dt_ns,
        }
        project["cases"][local_case_name] = local_case

    # dt_min sensitivity probes use the same accepted restart state and the
    # same 1.5 ns proposed step.  Only the floor is changed, allowing the
    # forced-accept path to be compared at 1, 0.5, and 0.25 ns.
    for dt_min_ns in (1.0, 0.5, 0.25):
        tag = str(dt_min_ns).replace(".", "p")
        diag_case_name = f"case_phase24_adaptive_dtmin_{tag}ns"
        diag_case = copy.deepcopy(adaptive_post_event)
        diag_case["restart_from"] = ADAPTIVE_POST_EVENT_CASE
        diag_case["restart_file_path"] = (
            f"../work/meshes/{diag_case['mesh']}/{ADAPTIVE_POST_EVENT_CASE}.result"
        )
        diag_case["restart_position"] = 3
        diag_case["restart_time"] = local_state_time
        diag_case["series_file"] = f"{diag_case_name}_series.csv"
        diag_case["iteration_series_file"] = f"{diag_case_name}_iterations.csv"
        diag_case["output_file_path"] = (
            f"../work/meshes/{diag_case['mesh']}/{diag_case_name}.result"
        )
        diag_case["output_result"] = False
        diag_case["timesteps"] = [["1.5[ns]", 1]]
        diag_end = local_state_time + 1.5e-9
        diag_case["adaptive_time"] = {
            "start": f"{local_state_time * 1000:.12f}[ms]",
            "end": f"{diag_end * 1000:.12f}[ms]",
            "requested_output_times": {
                "mode": "explicit",
                "times": [
                    f"{local_state_time * 1000:.12f}[ms]",
                    f"{diag_end * 1000:.12f}[ms]",
                ],
            },
            "dt_initial": "1.5[ns]",
            "dt_min": f"{dt_min_ns:g}[ns]",
            "dt_max": "1.5[ns]",
            "relative_tolerance": 2.0e-3,
            "absolute_tolerance": 1.0e-8,
            "r_min": 0.5,
            "r_max": 2.0,
            "max_growth": 1.5,
            "max_shrink": 0.5,
            "max_rejected": 8,
            "physical_event_times": [],
            "debug": True,
        }
        diag_case["phase24_smoke"] = {
            "purpose": "same-state dt_min forced-accept sensitivity diagnostic",
            "reference_case": ADAPTIVE_POST_EVENT_CASE,
            "restart_position": 3,
            "proposed_dt_ns": 1.5,
            "dt_min_ns": dt_min_ns,
        }
        project["cases"][diag_case_name] = diag_case
    OUTPUT.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
