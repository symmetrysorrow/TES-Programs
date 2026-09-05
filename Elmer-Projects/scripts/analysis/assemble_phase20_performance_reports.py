"""Assemble Phase20 timing, crossover, transient, and baseline reports.

The Elmer log exposes two clocks: the solver's reported CPU/REAL totals and
the external wall marker written by ``run.py``.  This report preserves both
and never treats the CPU-clock assembly counters as wall time.  HYPRE setup
time is emitted by the native integration at residual verbosity >= 6 and is
kept separate from Krylov solution time.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.assemble_phase20_control_reports import body_temperature_average
from scripts.analysis.evaluate_physical_parity import result_values


ROOT = Path(__file__).resolve().parents[2]
MESH_ROOT = ROOT / "work/meshes"
RESULT_ROOT = ROOT / "results"

FLOAT = r"([0-9.Ee+\-]+)"


def _floats(pattern: str, text: str) -> list[float]:
    return [float(item) for item in re.findall(pattern, text)]


def parse_log(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    solver_total = re.search(
        rf"SOLVER TOTAL TIME\(CPU,REAL\):\s*{FLOAT}\s+{FLOAT}", text
    )
    wall = re.search(rf"WALL_SECONDS\s+{FLOAT}", text)
    manifest_wall = None
    manifest = path.with_name("manifest.json")
    if not wall and manifest.exists():
        try:
            manifest_wall = json.loads(manifest.read_text(encoding="utf-8")).get("wall_seconds")
        except (OSError, json.JSONDecodeError):
            manifest_wall = None
    solve = _floats(rf"SolveHypre: Solution time \(method \d+\):\s*{FLOAT}", text)
    setup = _floats(rf"SolveHypre: setup time \(method \d+\):\s*{FLOAT}", text)
    iterations = [
        int(value)
        for value in re.findall(r"SolveHypre: Required iterations (\d+) \(method", text)
    ]
    assembly = [
        (float(step), float(total))
        for step, total in re.findall(
            rf"HeatSolve: iter:\s*\d+ Assembly: \(s\)\s*{FLOAT}\s+{FLOAT}",
            text,
        )
    ]
    heat_solve = [
        (float(step), float(total))
        for step, total in re.findall(
            rf"HeatSolve: iter:\s*\d+ Solve:\s*\(s\)\s*{FLOAT}\s+{FLOAT}",
            text,
        )
    ]
    return {
        "log": str(path.resolve()),
        "all_done": "MAIN: *** Elmer Solver: ALL DONE ***" in text,
        "wall_seconds": float(wall.group(1)) if wall else manifest_wall,
        "solver_cpu_seconds": float(solver_total.group(1)) if solver_total else None,
        "solver_real_seconds": float(solver_total.group(2)) if solver_total else None,
        "hypre_solution_times_s": solve,
        "hypre_solution_total_s": sum(solve),
        "hypre_setup_times_s": setup,
        "hypre_setup_total_s": sum(setup),
        "hypre_required_iterations": iterations,
        "reported_heat_assembly_cpu_s": assembly[-1][1] if assembly else None,
        "reported_heat_solve_cpu_s": heat_solve[-1][1] if heat_solve else None,
        "gpu_requested": "HYPRE GPU requested" in text,
        "gpu_matrix_migration_logged": "migrated HYPRE IJ matrices" in text,
        "gpu_vector_migration_logged": "migrated HYPRE IJ vectors" in text,
        "setup_timer_logged": bool(setup),
    }


def mesh_counts(mesh: Path) -> tuple[int, int]:
    fields = (mesh / "mesh.header").read_text(encoding="utf-8").split()
    return int(fields[0]), int(fields[1])


def enrich(case: str, size: str, backend: str, repeat: int | None) -> dict:
    log = RESULT_ROOT / case / "solver.log"
    mesh_name = {
        "small": "mesh_physical_parity_conformal",
        "production": "mesh_singlepixel_conformal_gpu_fine",
        "medium": "mesh_stycast_convergence_medium",
        "fine": "mesh_stycast_convergence_fine",
    }[size]
    mesh = MESH_ROOT / mesh_name
    parsed = parse_log(log)
    nodes, tets = mesh_counts(mesh)
    parsed.update(
        {
            "case": case,
            "size": size,
            "backend": backend,
            "repeat": repeat,
            "mesh": str(mesh.resolve()),
            "node_count": nodes,
            "tetrahedron_count": tets,
        }
    )
    if parsed["wall_seconds"] is not None and parsed["solver_real_seconds"] is not None:
        parsed["outside_solver_wall_s"] = max(
            parsed["wall_seconds"] - parsed["solver_real_seconds"], 0.0
        )
    return parsed


def stats(rows: list[dict], key: str) -> dict:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    if not values:
        return {"count": 0, "first": None, "median": None, "minimum": None, "stdev": None}
    return {
        "count": len(values),
        "first": values[0],
        "median": statistics.median(values),
        "minimum": min(values),
        "stdev": statistics.stdev(values) if len(values) > 1 else 0.0,
    }


def crossover(rows: list[dict]) -> list[dict]:
    by_size: dict[str, dict[str, list[dict]]] = {}
    for row in rows:
        by_size.setdefault(row["size"], {}).setdefault(row["backend"], []).append(row)
    result = []
    for size in ("small", "production", "medium", "fine"):
        groups = by_size.get(size, {})
        cpu = groups.get("cpu", [])
        gpu = groups.get("gpu", [])
        if not cpu or not gpu:
            continue
        cpu_wall = stats(cpu, "wall_seconds")
        gpu_wall = stats(gpu, "wall_seconds")
        cpu_linear = stats(cpu, "hypre_solution_total_s")
        gpu_linear = stats(gpu, "hypre_solution_total_s")
        cpu_result = MESH_ROOT / cpu[0]["mesh"].split("work\\meshes\\")[-1] / f"{cpu[0]['case']}.result"
        gpu_result = MESH_ROOT / gpu[0]["mesh"].split("work\\meshes\\")[-1] / f"{gpu[0]['case']}.result"
        parity = compare_results(cpu_result, gpu_result)
        result.append(
            {
                "size": size,
                "node_count": cpu[0]["node_count"],
                "tetrahedron_count": cpu[0]["tetrahedron_count"],
                "cpu": {"wall_s": cpu_wall, "linear_solve_s": cpu_linear},
                "gpu": {"wall_s": gpu_wall, "linear_solve_s": gpu_linear},
                "gpu_wall_speedup_median": (
                    cpu_wall["median"] / gpu_wall["median"]
                    if cpu_wall["median"] is not None and gpu_wall["median"] is not None
                    else None
                ),
                "gpu_linear_speedup_median": (
                    cpu_linear["median"] / gpu_linear["median"]
                    if cpu_linear["median"] is not None and gpu_linear["median"] is not None
                    else None
                ),
                "result_parity": parity,
            }
        )
    return result


def parse_transient_log(path: Path, backend: str, prefix: str) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    # One HeatSolve assembly/solve pair can represent one nonlinear iteration;
    # the iteration CSV is the authoritative timestep count when available.
    parsed = parse_log(path)
    markers = list(re.finditer(r"MAIN:\s+Time:\s*(\d+)/(\d+):", text))
    steps = []
    for index, marker in enumerate(markers):
        segment = text[marker.start(): markers[index + 1].start() if index + 1 < len(markers) else len(text)]
        assembly = [
            (float(step), float(total))
            for step, total in re.findall(
                rf"HeatSolve: iter:\s*\d+ Assembly: \(s\)\s*{FLOAT}\s+{FLOAT}",
                segment,
            )
        ]
        solve = _floats(rf"SolveHypre: Solution time \(method \d+\):\s*{FLOAT}", segment)
        setup = _floats(rf"SolveHypre: setup time \(method \d+\):\s*{FLOAT}", segment)
        iterations = [int(value) for value in re.findall(r"SolveHypre: Required iterations (\d+) \(method", segment)]
        elapsed = re.search(rf"MAIN: Elapsed time:\s*{FLOAT}\s+seconds", segment)
        steps.append(
            {
                "time_step": int(marker.group(1)),
                "time_step_count": int(marker.group(2)),
                "setup_s": sum(setup),
                "solve_s": sum(solve),
                "required_iterations": iterations,
                "assembly_reported_cpu_s": assembly[-1][1] if assembly else None,
                "elapsed_seconds_marker": float(elapsed.group(1)) if elapsed else None,
            }
        )
    parsed.update({"backend": backend, "prefix": prefix, "time_step_markers": len(markers), "time_steps": steps})
    return parsed


def compare_results(cpu_result: Path, gpu_result: Path) -> dict:
    if not cpu_result.exists() or not gpu_result.exists():
        return {"status": "NOT_AVAILABLE", "cpu_result": str(cpu_result), "gpu_result": str(gpu_result)}
    cpu = result_values(cpu_result, field_index=0)
    gpu = result_values(gpu_result, field_index=0)
    common = sorted(set(cpu) & set(gpu))
    if not common:
        return {"status": "NOT_AVAILABLE", "reason": "no common result nodes"}
    deltas = [abs(cpu[node] - gpu[node]) for node in common]
    return {
        "status": "PASS" if max(deltas) <= 1.0e-5 else "FAIL",
        "common_nodes": len(common),
        "temperature_max_abs_difference_K": max(deltas),
        "temperature_rms_difference_K": (sum(delta * delta for delta in deltas) / len(deltas)) ** 0.5,
        "cpu_result": str(cpu_result.resolve()),
        "gpu_result": str(gpu_result.resolve()),
    }


def compare_series(cpu_series: Path, gpu_series: Path) -> dict:
    if not cpu_series.exists() or not gpu_series.exists():
        return {"status": "NOT_AVAILABLE"}
    with cpu_series.open(encoding="utf-8", newline="") as handle:
        cpu_rows = list(csv.DictReader(handle))
    with gpu_series.open(encoding="utf-8", newline="") as handle:
        gpu_rows = list(csv.DictReader(handle))
    if not cpu_rows or not gpu_rows:
        return {"status": "NOT_AVAILABLE"}
    cpu = cpu_rows[-1]
    gpu = gpu_rows[-1]
    fields = ("tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W")
    deltas = {field: abs(float(cpu[field]) - float(gpu[field])) for field in fields}
    tolerances = {
        "tes_temperature_K": 1.0e-6,
        "tes_current_A": 1.0e-7,
        "tes_resistance_ohm": 1.0e-5,
        "tes_power_W": 1.0e-12,
    }
    return {
        "status": "PASS" if all(deltas[field] <= tolerances[field] for field in fields) else "FAIL",
        "final_time_s_cpu": cpu.get("time_s"),
        "final_time_s_gpu": gpu.get("time_s"),
        "absolute_differences": deltas,
        "absolute_tolerances": tolerances,
    }


def baseline() -> dict:
    # This is the previously validated CPU Mortar steady reference.  Its
    # solver log is retained in the sibling production worktree; keep the
    # provenance explicit instead of silently mixing runs.
    external = Path(r"D:/github/TES-Programs/Elmer-Projects/results/case_tes_steady_singlepixel_conformal_mortar_reference/solver.log")
    if not external.exists():
        return {"status": "NOT_AVAILABLE", "reason": str(external)}
    return {
        "status": "AVAILABLE",
        "source": str(external),
        "timing": parse_log(external),
        "comparison_note": "Historical validated CPU Mortar run; mesh/time-grid provenance must be checked before using replacement speedup as a strict apples-to-apples claim.",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "artifacts/phase20_conformal")
    args = parser.parse_args()

    rows = []
    for size in ("small", "production", "medium", "fine"):
        for backend in ("cpu", "gpu"):
            for repeat in range(1, 4):
                case = f"case_phase20_perf_{size}_{backend}_r{repeat}"
                log = RESULT_ROOT / case / "solver.log"
                if log.exists():
                    rows.append(enrich(case, size, backend, repeat))
    transient_rows = []
    for backend in ("cpu", "gpu"):
        for prefix in ("7step", "50step"):
            case = f"case_phase20_perf_transient_{backend}_{prefix}"
            log = RESULT_ROOT / case / "solver.log"
            if log.exists():
                transient_rows.append(parse_transient_log(log, backend, prefix))
    reuse_rows = []
    for backend in ("cpu", "gpu"):
        case = f"case_phase20_perf_transient_{backend}_7step_reuse_probe"
        log = RESULT_ROOT / case / "solver.log"
        if log.exists():
            row = parse_transient_log(log, backend, "7step_reuse_probe")
            row["case"] = case
            row["convergence_failure"] = "failed to converge" in log.read_text(encoding="utf-8", errors="replace").lower()
            baseline_case = RESULT_ROOT / f"case_phase20_perf_transient_{backend}_7step" / "solver.log"
            row["rebuild_baseline_log"] = str(baseline_case.resolve())
            result_dir = MESH_ROOT / "mesh_singlepixel_conformal_gpu_fine"
            row["result_parity_vs_rebuild"] = compare_results(
                result_dir / f"case_phase20_perf_transient_{backend}_7step.result",
                result_dir / f"{case}.result",
            )
            reuse_rows.append(row)

    breakdown = []
    for row in rows:
        breakdown.append(
            {
                "case": row["case"],
                "size": row["size"],
                "backend": row["backend"],
                "repeat": row["repeat"],
                "wall_seconds": row["wall_seconds"],
                "solver_real_seconds": row["solver_real_seconds"],
                "outside_solver_wall_s": row.get("outside_solver_wall_s"),
                "hypre_setup_total_s": row["hypre_setup_total_s"],
                "hypre_solution_total_s": row["hypre_solution_total_s"],
                "host_solver_remainder_s": (
                    row["solver_real_seconds"]
                    - row["hypre_setup_total_s"]
                    - row["hypre_solution_total_s"]
                    if row["solver_real_seconds"] is not None
                    else None
                ),
                "reported_heat_assembly_cpu_s": row["reported_heat_assembly_cpu_s"],
                "reported_heat_solve_cpu_s": row["reported_heat_solve_cpu_s"],
                "unmeasured_components": [
                    "HYPRE IJ construction versus device migration are included in native setup timer",
                    "TES circuit/UDF and result I/O are not independently timestamped by current runtime",
                ],
            }
        )
    report = {
        "status": "PASS" if rows else "NOT_RUN",
        "measurement_contract": {
            "wall": "external WALL_SECONDS marker",
            "solver_real": "Elmer SOLVER TOTAL TIME(CPU,REAL) real clock",
            "assembly": "HeatSolve cumulative reported CPU clock; not added to wall clock",
            "hypre_setup": "SolveHypre setup timer, includes IJ construction/migration/AMG setup in current binary",
            "hypre_solve": "SolveHypre solution timer",
        },
        "wall_time_breakdown": breakdown,
        "crossover": crossover(rows),
        "raw_runs": rows,
        "cpu_mortar_reference": baseline(),
        "limitations": [
            "Separate H2D migration and BoomerAMG setup timers require a rebuilt SolveHypre integration; current setup timer is the conservative combined bucket.",
            "GPU utilization requires an Nsight Systems summary when the profiler run is available.",
        ],
    }
    transient_speedups = {}
    for prefix in ("7step", "50step"):
        cpu_wall = next(
            (r["wall_seconds"] for r in transient_rows
             if r["backend"] == "cpu" and r["prefix"] == prefix),
            None,
        )
        gpu_wall = next(
            (r["wall_seconds"] for r in transient_rows
             if r["backend"] == "gpu" and r["prefix"] == prefix),
            None,
        )
        transient_speedups[prefix] = (
            cpu_wall / gpu_wall if cpu_wall is not None and gpu_wall else None
        )
    transient_report = {
        "status": "PASS" if transient_rows else "NOT_RUN",
        "runs": transient_rows,
        "reuse_probes": reuse_rows,
        "speedups": transient_speedups,
        "parity": {
            prefix: {
                "field": compare_results(
                    MESH_ROOT / "mesh_singlepixel_conformal_gpu_fine" / f"case_phase20_perf_transient_cpu_{prefix}.result",
                    MESH_ROOT / "mesh_singlepixel_conformal_gpu_fine" / f"case_phase20_perf_transient_gpu_{prefix}.result",
                ),
                "electrical_series": compare_series(
                    RESULT_ROOT / f"case_phase20_perf_transient_cpu_{prefix}" / f"case_phase20_perf_transient_cpu_{prefix}_series.csv",
                    RESULT_ROOT / f"case_phase20_perf_transient_gpu_{prefix}" / f"case_phase20_perf_transient_gpu_{prefix}_series.csv",
                ),
            }
            for prefix in ("7step", "50step")
        },
        "setup_reuse_observation": "No Precondition Recompute was not enabled in the accepted benchmark cases; explicit refactorize=False probes failed to converge at timestep 2 and are rejected.",
    }
    production = next((item for item in report["crossover"] if item["size"] == "production"), None)
    medium = next((item for item in report["crossover"] if item["size"] == "medium"), None)
    fine = next((item for item in report["crossover"] if item["size"] == "fine"), None)
    acceptance = {
        "status": "PASS" if rows and transient_rows else "NOT_RUN",
        "classification": "GPU_SOLVER_EFFECTIVE_BUT_HOST_BOUND",
        "speedups": {
            "production_wall_median": production["gpu_wall_speedup_median"] if production else None,
            "production_linear_median": production["gpu_linear_speedup_median"] if production else None,
            "medium_wall": medium["gpu_wall_speedup_median"] if medium else None,
            "fine_wall": fine["gpu_wall_speedup_median"] if fine else None,
            "transient_7step_wall": transient_speedups["7step"],
            "transient_50step_wall": transient_speedups["50step"],
        },
        "cpu_mortar_reference": report["cpu_mortar_reference"],
        "cpu_mortar_to_conformal_gpu": {
            "status": "NOT_STRICTLY_COMPARABLE",
            "historical_cpu_mortar_solver_real_s": report["cpu_mortar_reference"].get("timing", {}).get("solver_real_seconds") if report["cpu_mortar_reference"].get("status") == "AVAILABLE" else None,
            "reason": "The retained Mortar reference uses a different mesh/runtime and is not an apples-to-apples production transient. A common-mesh Mortar transient is required before claiming replacement speedup.",
        },
        "limiting_component": "host-side assembly/circuit/integration remainder plus repeated HYPRE setup; GPU HYPRE solve is effective but is not the dominant wall bucket at 90k transient scale.",
        "recommended_next_optimization": "Implement a correctness-guarded matrix-update/preconditioner reuse path (SolveHypre3-equivalent) or parallelize host FEM assembly; do not enable unconditional No Precondition Recompute for the temperature-dependent matrix.",
        "gpu_utilization_artifact": "gpu_utilization_summary.json",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "wall_time_breakdown.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "phase20_performance_crossover.json").write_text(json.dumps({"status": report["status"], "crossover": report["crossover"]}, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "phase20_transient_performance.json").write_text(json.dumps(transient_report, indent=2) + "\n", encoding="utf-8")
    (args.output_dir / "phase20_performance_acceptance.json").write_text(json.dumps(acceptance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"steady_rows": len(rows), "transient_rows": len(transient_rows), "status": report["status"]}, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
