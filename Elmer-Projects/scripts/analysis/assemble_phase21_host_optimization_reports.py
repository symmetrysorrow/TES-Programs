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
        rf"\s+total_cpu_s=\s*{FLOAT}\s+cached_elements=\s*(\d+)\s+cached_nodes=\s*(\d+)"
    )
    for match in profile_re.finditer(text):
        profile.append(
            {
                "time_step": int(match.group(1)),
                "nonlinear_iter": int(match.group(2)),
                "integration_cpu_s": float(match.group(3)),
                "circuit_and_output_cpu_s": float(match.group(4)),
                "total_cpu_s": float(match.group(5)),
                "cached_elements": int(match.group(6)),
                "cached_nodes": int(match.group(7)),
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
        "solver_cpu_seconds": float(solver_total.group(1)) if solver_total else None,
        "solver_real_seconds": float(solver_total.group(2)) if solver_total else None,
        "hypre_setup_total_s": sum(setup),
        "hypre_solution_total_s": sum(hypre_solve),
        "hypre_setup_count": len(setup),
        "hypre_solution_count": len(hypre_solve),
        "assembly_reported_cpu_s": assembly[-1][1] if assembly else None,
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
            measured_wall = parsed["hypre_setup_total_s"] + parsed["hypre_solution_total_s"]
            wall = parsed["wall_seconds"]
            rows.append(
                {
                    "case": name,
                    "backend": backend,
                    "phase": phase,
                    "steps": parsed["time_step_markers"],
                    "wall_seconds": wall,
                    "categories": {
                        "element_matrix_assembly": {
                            "seconds": None,
                            "percentage": None,
                            "clock": "not_independently_timed",
                            "source": "HeatSolve assembly CPU counter is reported separately",
                        },
                        "mass_matrix_contribution": {"seconds": None, "percentage": None, "clock": "not_instrumented"},
                        "conductivity_matrix_contribution": {"seconds": None, "percentage": None, "clock": "not_instrumented"},
                        "rhs_source_assembly": {"seconds": None, "percentage": None, "clock": "not_instrumented"},
                        "tes_circuit_evaluation": {
                            "seconds": sum(item["total_cpu_s"] for item in parsed["circuit_profile"]),
                            "percentage": None,
                            "clock": "UDF_CPU_TIME",
                        },
                        "tes_temperature_integration": {
                            "seconds": sum(item["integration_cpu_s"] for item in parsed["circuit_profile"]),
                            "percentage": None,
                            "clock": "UDF_CPU_TIME",
                        },
                        "absorber_pulse_integration": {
                            "seconds": None,
                            "percentage": None,
                            "clock": "cache_marker_only",
                            "source": "pulse temporal factor is cached per timestep; callback wall timer is not yet in HeatSolver",
                        },
                        "hypre_ij_fill_update": {
                            "seconds": None,
                            "percentage": None,
                            "clock": "included_in_hypre_setup",
                        },
                        "host_device_transfer": {
                            "seconds": None,
                            "percentage": None,
                            "clock": "included_in_hypre_setup_for_gpu",
                        },
                        "boomeramg_setup": {
                            "seconds": parsed["hypre_setup_total_s"],
                            "percentage": parsed["hypre_setup_total_s"] / wall * 100.0 if wall else None,
                            "clock": "native_wall_timer",
                        },
                        "krylov_solve": {
                            "seconds": parsed["hypre_solution_total_s"],
                            "percentage": parsed["hypre_solution_total_s"] / wall * 100.0 if wall else None,
                            "clock": "native_wall_timer",
                        },
                        "result_series_output": {
                            "seconds": None,
                            "percentage": None,
                            "clock": "not_independently_timed",
                            "source": "run.py collects output after solver exit; solver-side file time is not exposed",
                        },
                        "miscellaneous_host_residual": {
                            "seconds": max(wall - measured_wall, 0.0) if wall else None,
                            "percentage": max(wall - measured_wall, 0.0) / wall * 100.0 if wall else None,
                            "clock": "wall_residual",
                            "note": "contains FEM assembly, circuit, output, MPI synchronization, and other host work; not a hidden zero",
                        },
                    },
                    "step_wall_seconds": step_wall,
                    "step_rows": parsed["time_steps"],
                    "reported_heat_assembly_cpu_s": parsed["assembly_reported_cpu_s"],
                }
            )
    return {
        "status": "PASS" if rows else "NOT_RUN",
        "schema_version": "phase21.1",
        "clock_policy": {
            "wall": "WALL_SECONDS / MAIN Elapsed time",
            "native_solver": "HYPRE setup and solution timers",
            "udf": "Fortran CPU_TIME markers",
            "assembly": "HeatSolve reported CPU counter; never added to wall categories",
        },
        "runs": rows,
        "limitations": [
            "Element/mass/conductivity/RHS subterms require a HeatSolver instrumentation build; current report preserves them as not_independently_timed.",
            "HYPRE setup currently includes IJ construction and GPU migration in the native integration timer.",
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
    # global interface override differs between the two project JSONs.
    historical = PHASE20_ARTIFACT_ROOT / "phase20_performance_acceptance.json"
    historical_data = json.loads(historical.read_text(encoding="utf-8")) if historical.exists() else {}
    windows = {"steady": "steady", "1step": "1step", "7step": "7step"}
    measured = []
    for window, suffix in windows.items():
        mortar = available(f"case_phase21_mortar_{suffix}")
        cpu = available(f"case_phase21_conformal_cpu_{suffix}")
        gpu = available(f"case_phase21_conformal_gpu_{suffix}")
        if not (mortar and cpu and gpu):
            continue
        mortar_t = mortar["wall_seconds"] or mortar["solver_real_seconds"]
        cpu_t = cpu["wall_seconds"] or cpu["solver_real_seconds"]
        gpu_t = gpu["wall_seconds"] or gpu["solver_real_seconds"]
        measured.append(
            {
                "window": window,
                "mortar_cpu_s": mortar_t,
                "conformal_cpu_s": cpu_t,
                "conformal_gpu_s": gpu_t,
                "S_conformalCPU": mortar_t / cpu_t if cpu_t else None,
                "S_GPU": cpu_t / gpu_t if gpu_t else None,
                "S_total": mortar_t / gpu_t if gpu_t else None,
                "timing_source": "WALL_SECONDS when available, otherwise solver real clock",
                "field_gate": "pending parser integration for TES/absorber observables",
            }
        )
    return {
        "status": "PASS" if measured else "PREPARED_NOT_RUN",
        "definition": {
            "window": ["steady", "one-step", "7-step"],
            "same_geometry": "physical-parity single-pixel geometry",
            "same_time_grid": "shared pulse grid from the physical-parity project",
            "mortar": "validated CPU Mortar control",
            "conformal_cpu": "conformal shared-node CPU HYPRE",
            "conformal_gpu": "conformal shared-node GPU HYPRE",
        },
        "speedup_formulas": {
            "S_conformalCPU": "t_MortarCPU / t_conformalCPU",
            "S_GPU": "t_conformalCPU / t_conformalGPU",
            "S_total": "t_MortarCPU / t_conformalGPU",
        },
        "runs": measured,
        "historical_non_strict_reference": historical_data.get("cpu_mortar_to_conformal_gpu"),
        "reason": None if measured else "Phase21 common-window Mortar and conformal cases are prepared but require execution on the paired physical-parity meshes.",
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
        classification = "CPU_CONFORMAL_BEST"
    else:
        classification = "ARCHITECTURE_LIMITED"
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
        "note": "Production transient GPU acceleration is not claimed until the optimized 50-step and common Mortar window both pass.",
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
    return {
        "status": "RECORDED",
        "probes": probes,
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
