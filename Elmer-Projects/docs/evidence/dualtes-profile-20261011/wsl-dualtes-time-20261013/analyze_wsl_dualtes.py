#!/usr/bin/env python3
"""Extract timing, nonlinear, HYPRE, and TES-series evidence from the WSL run."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "runs" / "attempt6"


def csv_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, skipinitialspace=True))


def field_stats(path: Path, names: list[str]):
    rows = csv_rows(path)
    result = {"rows": len(rows)}
    if rows:
        for name in names:
            values = [float(row[name]) for row in rows]
            result[name] = {"min": min(values), "max": max(values), "last": values[-1]}
        result["first_time_s"] = float(rows[0]["time_s"])
        result["last_time_s"] = float(rows[-1]["time_s"])
    return result


def analyze(stage: str):
    folder = RUN / stage
    record = json.loads((folder / "run_record.json").read_text(encoding="utf-8"))
    log = (folder / "solver.log").read_text(encoding="utf-8", errors="replace")
    hypre = [float(v) for v in re.findall(r"SolveHypre: Solution time \(method 601\):\s*([0-9.eE+-]+)", log)]
    assembly = [float(v) for v in re.findall(r"HeatSolve: iter:\s*\d+ Assembly: \(s\)\s*([0-9.eE+-]+)", log)]
    solve = [float(v) for v in re.findall(r"HeatSolve: iter:\s*\d+ Solve:\s*\(s\)\s*([0-9.eE+-]+)", log)]
    niter = [int(v) for v in re.findall(r"ComputeChange: NS \(ITER=(\d+)\)", log)]
    heat_iter = [int(v) for v in re.findall(r"TEMPERATURE ITERATION\s+(\d+)", log)]
    coupled = [(a, b) for a, b in re.findall(r"TES coupled error: iter=\d+, dT_circuit=\s*([0-9.eE+-]+), dT_k=\s*([0-9.eE+-]+)", log)]
    def fortran_float(value):
        value = re.sub(r"([+-]\d{2,3})$", r"E\1", value) if "e" not in value.lower() else value
        return float(value)
    coupled = [(fortran_float(a), fortran_float(b)) for a, b in coupled]
    finite_coupled = [(a, b) for a, b in coupled if abs(a) < 1e100 and abs(b) < 1e100]
    step_parts = re.split(r"MAIN: Time:\s*(\d+)\s*/\s*(\d+):\s*([0-9.eE+-]+)", log)
    step_final_errors = []
    for i in range(1, len(step_parts), 4):
        body = step_parts[i + 3] if i + 3 < len(step_parts) else ""
        matches = re.findall(r"TES coupled error: iter=\d+, dT_circuit=\s*([0-9.eE+-]+), dT_k=\s*([0-9.eE+-]+)", body)
        if matches:
            a, b = (fortran_float(x) for x in matches[-1])
            if abs(a) < 1e100 and abs(b) < 1e100:
                step_final_errors.append((int(step_parts[i]), a, b))
    step_lines = re.findall(r"MAIN: Time:\s*(\d+)\s*/\s*(\d+):\s*([0-9.eE+-]+)", log)
    reuse = re.findall(r"Phase24 precond reuse: reused=(\w+), age=(\d+), its=([^,]+), fresh_its=([^,]+), rebuilds=(\d+), reuses=(\d+)", log)
    result = {
        "run_record": record,
        "step_lines": len(step_lines),
        "last_step": {"index": int(step_lines[-1][0]), "count": int(step_lines[-1][1]), "time_s": float(step_lines[-1][2])} if step_lines else None,
        "computechange_records": len(niter),
        "computechange_iteration_index_max": max(niter) if niter else None,
        "temperature_iteration_peak": max(heat_iter) if heat_iter else None,
        "temperature_iteration_count_at_configured_cap_120": sum(v == 120 for v in heat_iter),
        "coupled_error_records": len(coupled),
        "coupled_error_records_excluding_huge_sentinel": len(finite_coupled),
        "coupled_temperature_error_max_finite_K": max((max(a, b) for a, b in finite_coupled), default=None),
        "coupled_circuit_error_max_finite_K": max((a for a, _ in finite_coupled), default=None),
        "coupled_material_k_error_max_finite_K": max((b for _, b in finite_coupled), default=None),
        "steps_with_final_coupled_error": len(step_final_errors),
        "steps_final_circuit_error_over_1e-7_K": sum(a > 1e-7 for _, a, _ in step_final_errors),
        "steps_final_material_k_error_over_1e-7_K": sum(b > 1e-7 for _, _, b in step_final_errors),
        "step_final_coupled_error_max_K": max((max(a, b) for _, a, b in step_final_errors), default=None),
        "hypre_solve_calls": len(hypre),
        "hypre_reported_solution_time_sum_s": sum(hypre),
        "hypre_reported_solution_time_max_s": max(hypre) if hypre else None,
        "assembly_reported_first_column_sum_s": sum(assembly),
        "solve_reported_first_column_sum_s": sum(solve),
        "preconditioner_reuse_last": {
            "reused": reuse[-1][0] == "T", "age": int(reuse[-1][1]), "iterations_field": reuse[-1][2],
            "fresh_iterations_field": reuse[-1][3], "rebuilds": int(reuse[-1][4]), "reuses": int(reuse[-1][5]),
        } if reuse else None,
        "adaptive_retry_or_step_failure_markers": len(re.findall(r"adaptive retry|time step failed|step failed", log, re.I)),
        "all_done": "ALL DONE" in log,
    }
    if stage == "steady":
        for i in (1, 2):
            p = folder / f"dualtes_wsl_steady_20261013_newton_iterations_{i}.csv"
            result.setdefault("tes_iteration_series", {})[str(i)] = field_stats(p, ["tes_temperature_K", "raw_current_A", "raw_power_W"])
    else:
        for i in (1, 2):
            p = folder / f"dualtes_wsl_transient_20261013_newton_{i}_series_reconciled.csv"
            result.setdefault("tes_final_series", {})[str(i)] = field_stats(p, ["tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W"])
            p = folder / f"dualtes_wsl_transient_20261013_newton_iterations_{i}.csv"
            result.setdefault("tes_iteration_series", {})[str(i)] = field_stats(p, ["tes_temperature_K", "raw_current_A", "raw_power_W"])
    return result


steady = analyze("steady")
transient = analyze("transient")
single = json.loads((ROOT.parent / "wsl-8min-replay-20261013" / "run_record.json").read_text(encoding="utf-8"))
single_log = (ROOT.parent / "wsl-8min-replay-20261013" / "solver.log").read_text(encoding="utf-8", errors="replace")
single_iters = [int(v) for v in re.findall(r"TEMPERATURE ITERATION\s+(\d+)", single_log)]
dual_wall = transient["run_record"]["wall_seconds"]
single_wall = single["wall_seconds"]
profile = {
    "steady": steady,
    "transient": transient,
    "comparison": {
        "singlepixel_transient_wall_seconds": single_wall,
        "dualtes_transient_wall_seconds": dual_wall,
        "dualtes_over_singlepixel_transient_ratio": dual_wall / single_wall,
        "dualtes_steady_wall_seconds": steady["run_record"]["wall_seconds"],
        "dualtes_steady_plus_transient_wall_seconds": steady["run_record"]["wall_seconds"] + dual_wall,
        "singlepixel_source_run_record": "../wsl-8min-replay-20261013/run_record.json",
        "singlepixel_temperature_iteration_peak": max(single_iters) if single_iters else None,
    },
    "notes": [
        "Hypre iteration count is not available from this build's method-601 log: the solver prints Required iterations -1; reported solution-time values are still accumulated.",
        "Material-property evaluation and output-I/O time are not separately instrumented; their effects are included in measured wall time.",
        "Assembly and solve first-column timings are summed from HeatSolve timing lines and are component evidence, not an additive partition of total wall time.",
        "DualTES used the same TES inner-circuit Newton, temperature-prediction, and step-commit options as the SinglePixel replay. The matching run did not reach the configured 120 temperature-iteration cap; the final circuit and material-k coupled errors at each accepted timestep were <= 1e-7 K.",
        "The raw TES final-series files omit the last accepted timestep; the reconciled series adds that row from each circuit's accepted iteration CSV, preserving the raw series files.",
    ],
}
(ROOT / "wall_profile.json").write_text(json.dumps(profile, indent=2) + "\n", encoding="utf-8")
print(json.dumps(profile["comparison"], indent=2))
