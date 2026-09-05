"""Assemble Phase21 host-overhead and replacement benchmark artifacts.

The report keeps measurement clocks explicit.  HYPRE and UDF markers are
solver/UDF CPU clocks, while ``MAIN: Elapsed time`` and ``WALL_SECONDS`` are
wall clocks.  HeatSolve's assembly counter is never silently presented as
wall time.  Uninstrumented FEM subterms therefore remain visible as named
``not_independently_timed`` fields instead of being hidden in a misleading
large ``other`` bucket.
"""
from __future__ import annotations

import csv
import json
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.assemble_phase20_control_reports import body_temperature_average
from scripts.analysis.evaluate_physical_parity import result_values


ROOT = Path(__file__).resolve().parents[2]
RESULT_ROOT = ROOT / "results"
ARTIFACT_ROOT = ROOT / "artifacts" / "phase21_host_optimization"
PHASE20_ARTIFACT_ROOT = ROOT / "artifacts" / "phase20_conformal"
FLOAT = r"([0-9.Ee+\-]+)"


def floats(pattern: str, text: str) -> list[float]:
    return [float(value) for value in re.findall(pattern, text)]


def parse_log(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    solver_total = re.search(
        rf"SOLVER TOTAL TIME\(CPU,REAL\):\s*{FLOAT}\s+{FLOAT}", text
    )
    wall = re.search(rf"WALL_SECONDS\s+{FLOAT}", text)
    markers = list(re.finditer(r"MAIN:\s+Time:\s*(\d+)/(\d+):", text))
    elapsed = [float(value) for value in re.findall(rf"MAIN:\s+Elapsed time:\s*{FLOAT}\s+seconds", text)]
    assembly = [
        (float(step), float(total))
        for step, total in re.findall(
            rf"HeatSolve: iter:\s*\d+ Assembly: \(s\)\s*{FLOAT}\s+{FLOAT}", text
        )
    ]
    heat_wall_breakdown = []
    heat_wall_re = re.compile(
        rf"HeatSolveWallBreakdown:\s*step=(?P<step>\d+)\s+iter=(?P<iter>\d+)"
        rf"\s+total_wall_s=\s*(?P<total>{FLOAT})"
        rf"\s+element_traversal_wall_s=\s*(?P<traversal>{FLOAT})"
        rf"\s+local_fem_wall_s=\s*(?P<local>{FLOAT})"
        rf"\s+local_stiffness_wall_s=\s*(?P<stiffness>{FLOAT})"
        rf"\s+local_mass_wall_s=\s*(?P<mass>{FLOAT})"
        rf"\s+nonlinear_material_wall_s=\s*(?P<material>{FLOAT})"
        rf"\s+global_insertion_wall_s=\s*(?P<global>{FLOAT})"
        rf"\s+rhs_assembly_wall_s=\s*(?P<rhs>{FLOAT})"
        rf"\s+matrix_conversion_wall_s=\s*(?P<conversion>{FLOAT})"
        rf"\s+boundary_assembly_wall_s=\s*(?P<boundary>{FLOAT})"
    )
    for match in heat_wall_re.finditer(text):
        heat_wall_breakdown.append(
            {
                "time_step": int(match.group("step")),
                "nonlinear_iter": int(match.group("iter")),
                **{
                    key: float(match.group(name))
                    for key, name in (
                        ("total_wall_s", "total"),
                        ("element_traversal_wall_s", "traversal"),
                        ("local_fem_wall_s", "local"),
                        ("local_stiffness_wall_s", "stiffness"),
                        ("local_mass_wall_s", "mass"),
                        ("nonlinear_material_wall_s", "material"),
                        ("global_insertion_wall_s", "global"),
                        ("rhs_assembly_wall_s", "rhs"),
                        ("matrix_conversion_wall_s", "conversion"),
                        ("boundary_assembly_wall_s", "boundary"),
                    )
                },
            }
        )
    solve = [
        (float(step), float(total))
        for step, total in re.findall(
            rf"HeatSolve: iter:\s*\d+ Solve:\s*\(s\)\s*{FLOAT}\s+{FLOAT}", text
        )
    ]
    setup = floats(rf"SolveHypre: setup time \(method \d+\):\s*{FLOAT}", text)
    hypre_solve = floats(rf"SolveHypre: Solution time \(method \d+\):\s*{FLOAT}", text)
    profile = []
    profile_re = re.compile(
        rf"TESParallelCircuitProfile:\s*step=(\d+)\s+iter=(\d+)"
        rf"\s+integration_cpu_s=\s*{FLOAT}\s+circuit_output_cpu_s=\s*{FLOAT}"
        rf"\s+total_cpu_s=\s*{FLOAT}"
        rf"(?:\s+integration_wall_s=\s*{FLOAT}\s+circuit_output_wall_s=\s*{FLOAT}"
        rf"\s+total_wall_s=\s*{FLOAT})?"
        rf"\s+cached_elements=\s*(\d+)\s+cached_nodes=\s*(\d+)"
    )
    for match in profile_re.finditer(text):
        profile.append(
            {
                "time_step": int(match.group(1)),
                "nonlinear_iter": int(match.group(2)),
                "integration_cpu_s": float(match.group(3)),
                "circuit_and_output_cpu_s": float(match.group(4)),
                "total_cpu_s": float(match.group(5)),
                "integration_wall_s": float(match.group(6)) if match.group(6) else None,
                "circuit_and_output_wall_s": float(match.group(7)) if match.group(7) else None,
                "total_wall_s": float(match.group(8)) if match.group(8) else None,
                "cached_elements": int(match.group(9)),
                "cached_nodes": int(match.group(10)),
            }
        )
    step_rows = []
    previous_elapsed = 0.0
    for index, marker in enumerate(markers):
        end = markers[index + 1].start() if index + 1 < len(markers) else len(text)
        segment = text[marker.start():end]
        step = int(marker.group(1))
        step_assembly = [value for _, value in re.findall(
            rf"HeatSolve: iter:\s*\d+ Assembly: \(s\)\s*{FLOAT}\s+{FLOAT}", segment
        )]
        step_setup = floats(rf"SolveHypre: setup time \(method \d+\):\s*{FLOAT}", segment)
        step_solve = floats(rf"SolveHypre: Solution time \(method \d+\):\s*{FLOAT}", segment)
        step_elapsed = [float(value) for value in re.findall(
            rf"MAIN:\s+Elapsed time:\s*{FLOAT}\s+seconds", segment
        )]
        cumulative = step_elapsed[-1] if step_elapsed else None
        step_rows.append(
            {
                "time_step": step,
                "time_step_count": int(marker.group(2)),
                "wall_seconds_marker_cumulative": cumulative,
                "wall_seconds_marker_delta": (
                    max(cumulative - previous_elapsed, 0.0)
                    if cumulative is not None
                    else None
                ),
                "hypre_setup_s": sum(step_setup),
                "hypre_solve_s": sum(step_solve),
                "assembly_reported_cpu_s": float(step_assembly[-1]) if step_assembly else None,
                "circuit_profile_cpu_s": sum(
                    item["total_cpu_s"] for item in profile if item["time_step"] == step
                ),
            }
        )
        if cumulative is not None:
            previous_elapsed = cumulative
    manifest = path.with_name("manifest.json")
    manifest_data = {}
    if manifest.exists():
        try:
            manifest_data = json.loads(manifest.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            manifest_data = {}
    return {
        "case": path.parent.name,
        "log": str(path.resolve()),
        "all_done": "MAIN: *** Elmer Solver: ALL DONE ***" in text,
        "wall_seconds": float(wall.group(1)) if wall else manifest_data.get("wall_seconds"),
        "solver_wall_seconds": manifest_data.get("solver_wall_seconds"),
        "output_io_wall_seconds": manifest_data.get("output_io_wall_seconds"),
        "process_cpu_seconds": manifest_data.get("process_cpu_seconds"),
        "thread_cpu_seconds": manifest_data.get("thread_cpu_seconds"),
        "clock_policy": manifest_data.get("clock_policy", {}),
        "solver_cpu_seconds": float(solver_total.group(1)) if solver_total else None,
        "solver_real_seconds": float(solver_total.group(2)) if solver_total else None,
        "hypre_setup_total_s": sum(setup),
        "hypre_solution_total_s": sum(hypre_solve),
        "hypre_setup_count": len(setup),
        "hypre_solution_count": len(hypre_solve),
        "assembly_reported_cpu_s": assembly[-1][1] if assembly else None,
        "heat_wall_breakdown": heat_wall_breakdown,
        "solve_reported_cpu_s": solve[-1][1] if solve else None,
        "circuit_profile": profile,
        "time_steps": step_rows,
        "time_step_markers": len(markers),
        "gpu_requested": "HYPRE GPU requested" in text,
        "gpu_matrix_migration_logged": "migrated HYPRE IJ matrices" in text,
        "gpu_vector_migration_logged": "migrated HYPRE IJ vectors" in text,
        "collected_outputs": manifest_data.get("collected_outputs", []),
        "manifest": manifest_data,
    }


def case_log(name: str) -> Path:
    return RESULT_ROOT / name / "solver.log"


def available(name: str) -> dict | None:
    path = case_log(name)
    return parse_log(path) if path.exists() else None


def stats(values: list[float]) -> dict:
    if not values:
        return {"count": 0, "median": None, "minimum": None, "maximum": None}
    return {
        "count": len(values),
        "median": statistics.median(values),
        "minimum": min(values),
        "maximum": max(values),
    }


def transient_wall_breakdown() -> dict:
    rows = []
    for backend in ("cpu", "gpu"):
        for phase, prefix in (("phase20_baseline", "50step"), ("phase21_cached", "50step")):
            name = (
                f"case_phase20_perf_transient_{backend}_{prefix}"
                if phase == "phase20_baseline"
                else f"case_phase21_host_transient_{backend}_{prefix}"
            )
            parsed = available(name)
            if not parsed:
                continue
            step_wall = [
                row["wall_seconds_marker_delta"]
                for row in parsed["time_steps"]
                if row["wall_seconds_marker_delta"] is not None
            ]
            wall = parsed["wall_seconds"]
            udf_wall_values = [
                item["total_wall_s"] for item in parsed["circuit_profile"]
                if item.get("total_wall_s") is not None
            ]
            measured_categories = {
                "assembly_wall": sum(row["total_wall_s"] for row in parsed["heat_wall_breakdown"]),
                "hypre_setup_wall": parsed["hypre_setup_total_s"],
                "hypre_krylov_wall": parsed["hypre_solution_total_s"],
                "circuit_udf_wall": sum(udf_wall_values) if udf_wall_values else None,
                "output_io_wall": parsed.get("output_io_wall_seconds"),
            }
            known_wall = sum(value for value in measured_categories.values() if value is not None)
            unclassified_wall = max(wall - known_wall, 0.0) if wall is not None else None

            def wall_category(seconds: float | None, clock: str, source: str = "") -> dict:
                return {
                    "seconds": seconds,
                    "percentage": seconds / wall * 100.0 if seconds is not None and wall else None,
                    "clock": clock,
                    **({"source": source} if source else {}),
                }

            rows.append(
                {
                    "case": name,
                    "backend": backend,
                    "phase": phase,
                    "steps": parsed["time_step_markers"],
                    "wall_seconds": wall,
                    "solver_wall_seconds": parsed.get("solver_wall_seconds"),
                    "output_io_wall_seconds": parsed.get("output_io_wall_seconds"),
                    "cpu_profile": {
                        "solver_reported_cpu_seconds": parsed.get("solver_cpu_seconds"),
                        "launcher_child_process_cpu_seconds": parsed.get("process_cpu_seconds"),
                        "launcher_thread_cpu_seconds": parsed.get("thread_cpu_seconds"),
                        "assembly_reported_cpu_seconds": parsed.get("assembly_reported_cpu_s"),
                    },
                    "wall_breakdown_sum_seconds": known_wall + (unclassified_wall or 0.0),
                    "wall_breakdown_residual_seconds": (
                        wall - (known_wall + (unclassified_wall or 0.0))
                        if wall is not None else None
                    ),
                    "wall_breakdown_non_overlapping": True,
                    "categories": {
                        "assembly_wall": wall_category(
                            measured_categories["assembly_wall"] or None,
                            "HeatSolve_RealTime" if parsed["heat_wall_breakdown"] else "not_instrumented",
                            "native HeatSolve wall marker; subcomponents are reported in assembly_profile",
                        ),
                        "hypre_setup_wall": wall_category(parsed["hypre_setup_total_s"], "native_wall_timer"),
                        "hypre_krylov_wall": wall_category(parsed["hypre_solution_total_s"], "native_wall_timer"),
                        "circuit_udf_wall": wall_category(
                            measured_categories["circuit_udf_wall"],
                            "UDF_SYSTEM_CLOCK" if udf_wall_values else "not_instrumented",
                            "UDF profile wall timer; CPU profile remains separate",
                        ),
                        "output_io_wall": wall_category(
                            measured_categories["output_io_wall"],
                            "launcher_perf_counter" if measured_categories["output_io_wall"] is not None else "not_instrumented",
                            "post-solver output collection only",
                        ),
                        "unclassified_wall": wall_category(
                            unclassified_wall,
                            "wall_residual",
                            "includes assembly until native HeatSolve marker is available, MPI synchronization, and other uninstrumented work",
                        ),
                    },
                    "step_wall_seconds": step_wall,
                    "step_rows": parsed["time_steps"],
                    "reported_heat_assembly_cpu_s": parsed["assembly_reported_cpu_s"],
                }
            )
    return {
        "status": "PASS" if rows else "NOT_RUN",
        "schema_version": "phase22.1",
        "clock_policy": {
            "total_wall": "run.py time.perf_counter from solver launch through output collection",
            "solver_wall": "run.py time.perf_counter around the solver process",
            "process_cpu": "child process user+system CPU where resource.RUSAGE_CHILDREN exists",
            "thread_cpu": "launcher thread time.thread_time; never mixed into wall breakdown",
            "native_solver": "HYPRE setup and Krylov timers",
            "udf_wall": "Fortran SYSTEM_CLOCK markers; UDF CPU_TIME markers remain a separate profile",
            "assembly_wall": "reserved for HeatSolveWallBreakdown native marker",
        },
        "runs": rows,
        "limitations": [
            "When the rebuilt HeatSolve marker is present, assembly_wall is a top-level non-overlapping bucket; its local/global subterms are a nested diagnostic and must not be added again to the top-level sum.",
            "HYPRE setup currently includes IJ construction and GPU migration in the native integration timer.",
            "The category sum is non-overlapping by construction; unavailable categories are not imputed from CPU counters.",
        ],
    }


def assembly_profile() -> dict:
    runs = []
    for backend in ("cpu", "gpu"):
        for name in (
            f"case_phase21_host_transient_{backend}_50step",
            f"case_phase20_perf_transient_{backend}_50step",
        ):
            parsed = available(name)
            if parsed:
                runs.append(
                    {
                        "case": name,
                        "backend": backend,
                        "phase": "phase21_cached" if "phase21" in name else "phase20_baseline",
                        "reported_heat_assembly_cpu_s": parsed["assembly_reported_cpu_s"],
                        "reported_heat_solve_cpu_s": parsed["solve_reported_cpu_s"],
                        "wall_clock_breakdown": parsed["heat_wall_breakdown"],
                        "per_step": parsed["time_steps"],
                    }
                )
    return {
        "status": "PASS" if runs else "NOT_RUN",
        "classification": {
            "geometry": "static",
            "tetra_shape_function_gradients": "static",
            "element_connectivity": "static",
            "density": "static_if_material_constant",
            "nominal_heat_capacity": "static_if_material_constant",
            "constant_k_material_contribution": "static_for_constant_k_bodies",
            "mass_matrix": "static_for lumped/constant-capacity terms; timestep coefficient remains dynamic",
            "sparsity_pattern": "static",
            "temperature_dependent_membrane_conductivity": "dynamic",
            "timestep_coefficient": "dynamic",
            "nonlinear_source_terms": "dynamic",
            "TES_Joule_power": "dynamic but UDF circuit metadata is cached",
            "pulse_source": "dynamic in time; spatial geometry/normalization is cached",
        },
        "runs": runs,
        "prototype": {
            "implemented": [
                "TES body element/node connectivity cache",
                "TES and pulse immutable constant cache",
                "pulse temporal factor cache per timestep",
            ],
            "not_implemented": [
                "native HeatSolver A_static + A_dynamic matrix split",
                "AMG hierarchy reuse",
            ],
            "reason": "matrix split needs native HeatSolver evidence and must preserve temperature-dependent coefficients",
        },
        "wall_subcomponent_policy": {
            "local_fem_calculation": "DiffuseConvective(Gen)Compose wall timer",
            "global_sparse_insertion": "DefaultUpdateEquations wall timer",
            "mass_contribution": "DefaultUpdateMass/Default1stOrderTime wall timer",
            "element_traversal": "bulk residual after measured local/mass/global calls",
            "local_stiffness_formation": "not independently separable in Compose; nested diagnostic only",
            "nonlinear_material_evaluation": "not independently separable in stock HeatSolve; remains in traversal residual",
            "rhs_assembly": "included in global insertion until a lower-level RHS hook is added",
            "matrix_format_conversion": "not independently timed; remains in assembly residual",
            "boundary_assembly": "boundary loop wall timer",
        },
    }


def circuit_profile() -> dict:
    runs = []
    for backend in ("cpu", "gpu"):
        for phase, name in (
            ("phase20_baseline", f"case_phase20_perf_transient_{backend}_50step"),
            ("phase21_cached", f"case_phase21_host_transient_{backend}_50step"),
        ):
            parsed = available(name)
            if parsed:
                profiles = parsed["circuit_profile"]
                runs.append(
                    {
                        "case": name,
                        "backend": backend,
                        "phase": phase,
                        "profile_samples": len(profiles),
                        "total_cpu_seconds": sum(row["total_cpu_s"] for row in profiles),
                        "integration_cpu_seconds": sum(row["integration_cpu_s"] for row in profiles),
                        "circuit_and_output_cpu_seconds": sum(row["circuit_and_output_cpu_s"] for row in profiles),
                        "cache_contract": {
                            "cached_elements": sorted({row["cached_elements"] for row in profiles}),
                            "cached_nodes": sorted({row["cached_nodes"] for row in profiles}),
                        },
                        "samples": profiles,
                    }
                )
    return {
        "status": "PASS" if runs else "NOT_RUN",
        "runs": runs,
        "cache_scope": [
            "body element discovery",
            "node connectivity traversal",
            "TES material/electrical constants",
            "pulse spatial geometry and discrete normalization",
            "pulse temporal interval factor per timestep",
        ],
    }


def udf_cache_benchmark() -> dict:
    rows = []
    for backend in ("cpu", "gpu"):
        baseline = available(f"case_phase20_perf_transient_{backend}_50step")
        optimized = available(f"case_phase21_host_transient_{backend}_50step")
        if baseline and optimized and baseline["wall_seconds"] and optimized["wall_seconds"]:
            rows.append(
                {
                    "backend": backend,
                    "baseline_case": baseline["case"],
                    "optimized_case": optimized["case"],
                    "baseline_wall_s": baseline["wall_seconds"],
                    "optimized_wall_s": optimized["wall_seconds"],
                    "wall_speedup": baseline["wall_seconds"] / optimized["wall_seconds"],
                    "baseline_circuit_cpu_s": sum(row["total_cpu_s"] for row in baseline["circuit_profile"]),
                    "optimized_circuit_cpu_s": sum(row["total_cpu_s"] for row in optimized["circuit_profile"]),
                }
            )
    return {
        "status": "PASS" if rows else "WAITING_FOR_PHASE21_RUNS",
        "comparison": rows,
        "physics_gate": "must compare TES T/current/resistance/power and absorber temperature before accepting speedup",
        "baseline_provenance": "Phase20 50-step logs were captured before the cache revision",
    }


def _last_series_row(case: str) -> dict[str, str] | None:
    candidates = sorted((RESULT_ROOT / case).glob("*_series.csv"))
    if not candidates:
        return None
    with candidates[0].open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, skipinitialspace=True))
    return rows[-1] if rows else None


def regression_gate() -> dict:
    baseline_case = "case_phase20_perf_transient_cpu_50step"
    optimized_case = "case_phase21_host_transient_cpu_50step"
    baseline = _last_series_row(baseline_case)
    optimized = _last_series_row(optimized_case)
    mesh = ROOT / "work/meshes/mesh_singlepixel_conformal_gpu_fine"
    baseline_result = mesh / f"{baseline_case}.result"
    optimized_result = mesh / f"{optimized_case}.result"
    fields = ("tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W")
    electrical = {}
    if baseline and optimized:
        for field in fields:
            electrical[field] = {
                "baseline": float(baseline[field]),
                "optimized": float(optimized[field]),
                "absolute_difference": abs(float(baseline[field]) - float(optimized[field])),
            }
    temperature = {}
    if baseline_result.exists() and optimized_result.exists():
        base_values = result_values(baseline_result, field_index=0)
        opt_values = result_values(optimized_result, field_index=0)
        common = sorted(set(base_values) & set(opt_values))
        deltas = [abs(base_values[node] - opt_values[node]) for node in common]
        temperature = {
            "common_nodes": len(common),
            "max_abs_difference_K": max(deltas) if deltas else None,
            "rms_difference_K": (sum(value * value for value in deltas) / len(deltas)) ** 0.5 if deltas else None,
            "baseline_absorber_temperature_K": body_temperature_average(mesh, baseline_result, "abs"),
            "optimized_absorber_temperature_K": body_temperature_average(mesh, optimized_result, "abs"),
        }
        temperature["absorber_temperature_absolute_difference_K"] = abs(
            temperature["baseline_absorber_temperature_K"]
            - temperature["optimized_absorber_temperature_K"]
        )
    tolerances = {
        "tes_temperature_K": 1.0e-5,
        "tes_current_A": 1.0e-7,
        "tes_resistance_ohm": 1.0e-5,
        "tes_power_W": 1.0e-12,
        "absorber_temperature_K": 1.0e-5,
        "temperature_field_max_abs_K": 1.0e-5,
    }
    checks = {
        field: item["absolute_difference"] <= tolerances[field]
        for field, item in electrical.items()
    }
    if temperature:
        checks["absorber_temperature_K"] = temperature["absorber_temperature_absolute_difference_K"] <= tolerances["absorber_temperature_K"]
        checks["temperature_field_max_abs_K"] = temperature["max_abs_difference_K"] <= tolerances["temperature_field_max_abs_K"]
    return {
        "status": "PASS" if checks and all(checks.values()) else "NOT_RUN_OR_FAIL",
        "baseline_case": baseline_case,
        "optimized_case": optimized_case,
        "electrical_observables": electrical,
        "temperature_observables": temperature,
        "tolerances": tolerances,
        "checks": checks,
        "scope": "CPU 50-step cache regression; GPU rerun unavailable because the current WSL session has no CUDA-capable device",
    }


def fine_transient() -> dict:
    rows = []
    for prefix in ("7step", "20step"):
        pair = {}
        for backend in ("cpu", "gpu"):
            parsed = available(f"case_phase21_fine_transient_{backend}_{prefix}")
            if parsed:
                pair[backend] = parsed
        if pair:
            entry = {"prefix": prefix}
            for backend, parsed in pair.items():
                entry[backend] = {
                    "case": parsed["case"],
                    "wall_seconds": parsed["wall_seconds"],
                    "hypre_setup_s": parsed["hypre_setup_total_s"],
                    "hypre_solve_s": parsed["hypre_solution_total_s"],
                    "steps": parsed["time_step_markers"],
                }
            if "cpu" in pair and "gpu" in pair:
                entry["gpu_wall_speedup"] = pair["cpu"]["wall_seconds"] / pair["gpu"]["wall_seconds"]
            rows.append(entry)
    gpu_probe = available("case_phase21_host_transient_gpu_7step")
    gpu_status = None
    if gpu_probe:
        gpu_status = {
            "status": "PASS" if gpu_probe["all_done"] else "UNAVAILABLE",
            "all_done": gpu_probe["all_done"],
            "log": gpu_probe["log"],
            "reason": "no CUDA-capable device is detected" if not gpu_probe["all_done"] else None,
        }
    return {
        "status": "PASS" if rows else "NOT_RUN",
        "runs": rows,
        "gpu_probe": gpu_status,
    }


def io_benchmark() -> dict:
    rows = []
    modes = ("full_io", "no_vtu", "no_result", "no_iteration_csv", "no_series_csv")
    for backend in ("cpu", "gpu"):
        parsed = {mode: available(f"case_phase21_io_{backend}_{mode}") for mode in modes}
        if any(parsed.values()):
            entry = {"backend": backend, "modes": {}}
            for mode, row in parsed.items():
                if row:
                    entry["modes"][mode] = {
                        "case": row["case"],
                        "wall_seconds": row["wall_seconds"],
                        "outputs": row["collected_outputs"],
                    }
            if "full_io" in entry["modes"]:
                full = entry["modes"]["full_io"]["wall_seconds"]
                for mode, item in entry["modes"].items():
                    item["speedup_vs_full_io"] = full / item["wall_seconds"]
            rows.append(entry)
    return {
        "status": "PASS" if rows else "NOT_RUN",
        "definition": "production 7-step, one output class removed at a time; full_io enables VTU",
        "runs": rows,
        "policy": "only full_io/no_vtu with preserved validation outputs may be used for production speedup claims",
    }


def replacement_benchmark() -> dict:
    # The strict common-mesh Mortar runs are prepared separately because the
    # global interface override differs between the two project JSONs.  Keep a
    # state row for every planned case: a partial benchmark is evidence, not
    # an empty PREPARED_NOT_RUN placeholder.
    historical = PHASE20_ARTIFACT_ROOT / "phase20_performance_acceptance.json"
    historical_data = json.loads(historical.read_text(encoding="utf-8")) if historical.exists() else {}
    windows = {"steady": "steady", "1step": "1step", "7step": "7step"}
    case_states = []

    def state_for(case: str, backend: str, blocked_by: str | None = None) -> dict:
        parsed = available(case)
        if parsed is not None:
            log_text = Path(parsed["log"]).read_text(encoding="utf-8", errors="replace")
            if parsed["all_done"]:
                status = "DONE"
            elif "failed to converge" in log_text.lower() or "not converged" in log_text.lower():
                status = "FAILED_NONCONVERGENCE"
            else:
                status = "NOT_RUN"
            details = {
                "all_done": parsed["all_done"],
                "wall_seconds": parsed["wall_seconds"],
                "log": parsed["log"],
            }
            if status == "FAILED_NONCONVERGENCE":
                failures = [line.strip() for line in log_text.splitlines() if "converge" in line.lower() or "residual" in line.lower()]
                details["failure_tail"] = failures[-5:]
        else:
            status = "BLOCKED_NO_CUDA" if backend == "gpu" and gpu_runtime_blocked() else "NOT_RUN"
            details = {"all_done": False, "wall_seconds": None}
            if blocked_by:
                details["blocked_by"] = blocked_by
        return {"case": case, "backend": backend, "status": status, **details}

    def gpu_runtime_blocked() -> bool:
        for candidate in RESULT_ROOT.glob("*/solver.log"):
            try:
                text = candidate.read_text(encoding="utf-8", errors="replace").lower()
            except OSError:
                continue
            if "no cuda-capable device" in text or "cuda error code=100" in text:
                return True
        return False

    for suffix in windows.values():
        case_states.append(state_for(f"case_phase21_mortar_{suffix}", "mortar"))
    cpu_steady = "case_phase21_conformal_cpu_steady"
    cpu_steady_state = state_for(cpu_steady, "cpu")
    case_states.append(cpu_steady_state)
    for suffix in ("1step", "7step"):
        case_states.append(state_for(f"case_phase21_conformal_cpu_{suffix}", "cpu", cpu_steady))
    for suffix in windows.values():
        case_states.append(state_for(f"case_phase21_conformal_gpu_{suffix}", "gpu", "CUDA runtime unavailable"))

    measured = []
    for window, suffix in windows.items():
        mortar = available(f"case_phase21_mortar_{suffix}")
        cpu = available(f"case_phase21_conformal_cpu_{suffix}")
        gpu = available(f"case_phase21_conformal_gpu_{suffix}")
        if not (mortar and cpu):
            continue
        mortar_t = mortar["wall_seconds"] or mortar["solver_real_seconds"]
        cpu_t = cpu["wall_seconds"] or cpu["solver_real_seconds"]
        gpu_t = (gpu["wall_seconds"] or gpu["solver_real_seconds"]) if gpu else None
        measured.append(
            {
                "window": window,
                "status": "DONE" if gpu_t else "BLOCKED_NO_CUDA",
                "mortar_cpu_s": mortar_t,
                "conformal_cpu_s": cpu_t,
                "conformal_gpu_s": gpu_t,
                "S_algorithmic": mortar_t / cpu_t if cpu_t else None,
                "S_gpu": cpu_t / gpu_t if gpu_t else None,
                "S_total": mortar_t / gpu_t if gpu_t else None,
                "timing_source": "WALL_SECONDS when available, otherwise solver real clock",
                "field_gate": "physical-reference comparison is recorded separately; replacement speedup requires validated observables",
            }
        )
    done_count = sum(row["status"] == "DONE" for row in case_states)
    failed_count = sum(row["status"].startswith("FAILED") for row in case_states)
    blocked_count = sum(row["status"].startswith("BLOCKED") for row in case_states)
    return {
        "status": "DONE" if all(row["status"] == "DONE" for row in case_states) else "PARTIAL",
        "case_state_counts": {
            "DONE": done_count,
            "FAILED_NONCONVERGENCE": failed_count,
            "BLOCKED_NO_CUDA": blocked_count,
            "NOT_RUN": len(case_states) - done_count - failed_count - blocked_count,
        },
        "definition": {
            "window": ["steady", "one-step", "7-step"],
            "same_geometry": "physical-parity single-pixel geometry",
            "same_time_grid": "shared pulse grid from the physical-parity project",
            "mortar": "validated CPU Mortar control",
            "conformal_cpu": "conformal shared-node CPU HYPRE",
            "conformal_gpu": "conformal shared-node GPU HYPRE",
        },
        "speedup_formulas": {
            "S_algorithmic": "t_MortarCPU / t_ConformalCPU",
            "S_gpu": "t_ConformalCPU / t_ConformalGPU",
            "S_total": "t_MortarCPU / t_conformalGPU",
        },
        "tolerance_policy": {
            "production_candidate": 1.0e-7,
            "strict_reference": 1.0e-8,
            "selection_artifact": "artifacts/phase22_physical_tolerance/hypre_physical_tolerance_study.json",
            "status": "PROVISIONAL",
        },
        "physical_reference_gate": {
            "status": "INCOMPLETE",
            "reason": "retained Mortar run has no machine-readable TES electrical/pulse series and shows a route/mesh temperature offset; CPU timing is retained but not declared a validated physical replacement",
        },
        "case_states": case_states,
        "runs": measured,
        "historical_non_strict_reference": historical_data.get("cpu_mortar_to_conformal_gpu"),
        "reason": "GPU cases are BLOCKED_NO_CUDA; Mortar and conformal CPU timing results are retained, while the physical replacement gate remains incomplete.",
    }


def backend_recommendation(fine: dict, cache: dict) -> dict:
    phase20_acceptance = PHASE20_ARTIFACT_ROOT / "phase20_performance_acceptance.json"
    previous = json.loads(phase20_acceptance.read_text(encoding="utf-8")) if phase20_acceptance.exists() else {}
    fine_rows = [row for row in fine.get("runs", []) if row.get("gpu_wall_speedup") is not None]
    cache_rows = cache.get("comparison", [])
    if cache_rows:
        cpu = next((row for row in cache_rows if row["backend"] == "cpu"), None)
        gpu = next((row for row in cache_rows if row["backend"] == "gpu"), None)
    else:
        cpu = gpu = None
    if cpu and gpu:
        production_backend = "gpu" if gpu["optimized_wall_s"] < cpu["optimized_wall_s"] else "cpu"
        production_reason = "Phase21 cached-UDF 50-step measured wall time"
    elif cpu:
        production_backend = "cpu"
        production_reason = "CPU 50-step measured; GPU optimized rerun is unavailable because WSL has no CUDA-capable device"
    else:
        production_backend = "undetermined"
        production_reason = "Phase21 cached-UDF 50-step pair not run"
    if fine_rows and all(row["gpu_wall_speedup"] > 1.0 for row in fine_rows):
        classification = "GPU_EFFECTIVE_AT_LARGE_SCALE"
    elif production_backend == "cpu":
        # A CPU run without a GPU runtime comparison is a measured CPU
        # datapoint, not a final backend verdict.
        classification = "GPU_VERDICT_PENDING_RUNTIME"
    else:
        classification = "GPU_VERDICT_PENDING_RUNTIME"
    return {
        "status": "PROVISIONAL",
        "classification": classification,
        "backend_policy": {
            "small_medium": "conformal CPU HYPRE unless a measured case-specific result overrides this",
            "large_fine": "conformal GPU HYPRE only when the transient wall benchmark is >1.0x",
            "50_step_priority": True,
        },
        "phase21_50step_choice": {"backend": production_backend, "reason": production_reason},
        "previous_phase20_classification": previous.get("classification"),
        "strict_mortar_replacement_gate": "pending",
        "optimization_classification": "SAFE_OPTIMIZATION_WITH_NO_MEASURED_SPEEDUP",
        "note": "The UDF immutable-data cache remains enabled and physics-exact, but its wall speedup is not claimed. GPU verdict remains pending because the optimized GPU rerun and common replacement GPU case are blocked by CUDA availability.",
    }


def runtime_blockers() -> dict:
    probes = {}
    for name in (
        "case_phase21_host_transient_gpu_7step",
        "case_phase21_conformal_cpu_steady",
    ):
        parsed = available(name)
        if parsed:
            text = Path(parsed["log"]).read_text(encoding="utf-8", errors="replace")
            probes[name] = {
                "all_done": parsed["all_done"],
                "wall_seconds": parsed["wall_seconds"],
                "failure_lines": [
                    line.strip()
                    for line in text.splitlines()
                    if "failed" in line.lower() or "no cuda-capable" in line.lower()
                ][-5:],
            }
    planned_gpu_runs = [
        "case_phase21_host_transient_gpu_50step",
        "case_phase21_fine_transient_gpu_7step",
        "case_phase21_conformal_gpu_7step",
    ]
    cuda_blocked = any(
        "no cuda-capable device" in line.lower()
        for probe in probes.values()
        for line in probe.get("failure_lines", [])
    )
    for case in planned_gpu_runs:
        parsed = available(case)
        probes.setdefault(
            case,
            {
                "all_done": parsed["all_done"] if parsed else False,
                "wall_seconds": parsed["wall_seconds"] if parsed else None,
                "status": "DONE" if parsed and parsed["all_done"] else ("BLOCKED_NO_CUDA" if cuda_blocked else "NOT_RUN"),
                "failure_lines": [],
            },
        )
    return {
        "status": "RECORDED",
        "probes": probes,
        "gpu_recovery_first_runs": planned_gpu_runs,
        "policy": "failed GPU/device or convergence probes remain blockers; they are not replaced by relaxed tolerances or synthetic timings",
    }


def write(name: str, data: dict) -> None:
    ARTIFACT_ROOT.mkdir(parents=True, exist_ok=True)
    (ARTIFACT_ROOT / name).write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    breakdown = transient_wall_breakdown()
    assembly = assembly_profile()
    circuit = circuit_profile()
    cache = udf_cache_benchmark()
    fine = fine_transient()
    io = io_benchmark()
    replacement = replacement_benchmark()
    gate = regression_gate()
    write("transient_wall_breakdown.json", breakdown)
    write("assembly_profile.json", assembly)
    write("circuit_profile.json", circuit)
    write("static_dynamic_matrix_analysis.json", assembly)
    write("udf_cache_benchmark.json", cache)
    write("phase21_regression_gate.json", gate)
    write("fine_transient_cpu_gpu.json", fine)
    write("io_overhead_benchmark.json", io)
    write("mortar_conformal_replacement_benchmark.json", replacement)
    write("production_backend_recommendation.json", backend_recommendation(fine, cache))
    write("runtime_blockers.json", runtime_blockers())
    print(ARTIFACT_ROOT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
