"""Validate and summarize Phase20 block-Schur probe CSVs.

The native probe records measurements, not inferred solver state.  In
particular an unavailable outer residual stays unavailable; it is never
converted to zero.  Cumulative timers are sampled at their last row while
per-call timers are summed, so the two meanings cannot be double counted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable

TRUE = {"true", "t", "1", "yes"}
FALSE = {"false", "f", "0", "no"}
OUTER_REQUIRED = {"outer_iteration", "preconditioner_application", "solver_reported_iteration"}
SCHUR_REQUIRED = {"schur_solve", "iterations", "initial_residual", "final_residual", "reached_tolerance", "hit_maxiter", "breakdown", "nonfinite"}
SCHUR_RE = re.compile(r"Matrix-free Schur GMRES iterations:\s*(\d+)\s+residual:\s*([0-9.Ee+-]+)")
WARN_RE = re.compile(r"Inner Schur GMRES reached its limit")


def _read_csv(path: Path | None, required: set[str]) -> tuple[list[dict[str, str]], dict[str, Any]]:
    if path is None:
        return [], {"state": "NOT PROVIDED", "path": None}
    if not path.exists():
        return [], {"state": "MISSING", "path": str(path)}
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            headers = set(reader.fieldnames or [])
            if not headers:
                return [], {"state": "EMPTY", "path": str(path)}
            missing = sorted(required - headers)
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        return [], {"state": "READ_ERROR", "path": str(path), "error": str(exc)}
    if missing:
        return rows, {"state": "SCHEMA_ERROR", "path": str(path), "missing_columns": missing, "headers": sorted(headers)}
    if not rows:
        return [], {"state": "HEADER_ONLY", "path": str(path), "headers": sorted(headers)}
    return rows, {"state": "OK", "path": str(path), "headers": sorted(headers)}


def validate_csv_schema(path: Path, kind: str) -> dict[str, Any]:
    """Return a machine-readable schema result without raising on bad output."""
    required = OUTER_REQUIRED if kind == "outer" else SCHUR_REQUIRED
    _, status = _read_csv(path, required)
    return status


def _float(row: dict[str, str], *names: str) -> float | None:
    for name in names:
        raw = row.get(name, "")
        if raw not in (None, ""):
            try:
                return float(raw)
            except (TypeError, ValueError):
                return None
    return None


def _bool(row: dict[str, str], name: str) -> bool | None:
    value = str(row.get(name, "")).strip().lower()
    if value in TRUE:
        return True
    if value in FALSE:
        return False
    return None


def _finite_nonnegative(value: float | None) -> bool:
    return value is not None and math.isfinite(value) and value >= 0.0


def _last(rows: Iterable[dict[str, str]], *names: str) -> float | None:
    values = [_float(row, *names) for row in rows]
    values = [value for value in values if value is not None]
    return values[-1] if values else None


def _sum(rows: Iterable[dict[str, str]], *names: str) -> float | None:
    values = [_float(row, *names) for row in rows]
    values = [value for value in values if value is not None]
    return sum(values) if values else None


def _k_stage_totals(rows: list[dict[str, str]]) -> dict[str, Any]:
    stages = {
        "primal_block_solve": ("k_actions_primal_block_solve", "k_actions_primal"),
        "matrix_free_schur_action": ("k_actions_matrix_free_schur", "k_actions_schur"),
        "full_factorization_upper_correction": ("k_actions_full_upper_correction", "k_actions_upper_correction"),
        "setup_rebuild": ("k_actions_setup_rebuild",),
    }
    result: dict[str, float | None] = {}
    for label, names in stages.items():
        result[label] = _sum(rows, *names)
    explicit = _sum(rows, "k_actions_total", "k_actions")
    stage_values = [result[key] for key in stages]
    mismatch = False
    if explicit is not None and all(value is not None for value in stage_values):
        mismatch = not math.isclose(explicit, sum(stage_values), rel_tol=1.0e-9, abs_tol=1.0e-12)
    if explicit is not None:
        result["total"] = explicit
    else:
        values = [value for key, value in result.items()
                  if key not in {"total", "total_mismatch", "stage_sum"} and value is not None]
        result["total"] = sum(values) if values else None
    result["total_mismatch"] = mismatch
    result["stage_sum"] = sum(stage_values) if all(value is not None for value in stage_values) else None
    return result


def _log_reduction(initial: float | None, final: float | None) -> tuple[float | None, str | None]:
    if not _finite_nonnegative(initial) or not _finite_nonnegative(final):
        return None, "missing_or_invalid_residual"
    if initial == 0.0:
        return None, "zero_initial_residual"
    if final == 0.0:
        return None, "zero_final_residual"
    return math.log10(initial / final), None


def _workload_signature(rows: list[dict[str, str]]) -> str | None:
    if not rows:
        return None
    # Strategy, backend, case name, and probe prefix are not workload identity.
    fields = ("mesh_id", "restart_id", "restart_file", "timestep", "nonlinear_iteration",
              "rhs_id", "matrix_fingerprint", "outer_limit", "outer_tolerance",
              "linear_system_tolerance", "schur_tolerance", "schur_max_iterations",
              "schur_restart")
    values = {field: rows[0].get(field, "") for field in fields if rows[0].get(field, "") != ""}
    if not values and rows[0].get("workload_id", ""):
        # Legacy fixtures only. Native output should provide canonical fields.
        values = {"legacy_workload_id": rows[0]["workload_id"]}
    if not values:
        return None
    encoded = json.dumps(values, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def validate_same_workload(*summaries: dict[str, Any]) -> dict[str, Any]:
    signatures = [summary.get("workload_signature") for summary in summaries]
    available = all(signature is not None for signature in signatures)
    return {"same_workload": available and len(set(signatures)) == 1, "status": "PASS" if available and len(set(signatures)) == 1 else "INCOMPLETE", "signatures": signatures}


def summarize(outer_path: Path | None = None, schur_path: Path | None = None, log_path: Path | None = None) -> dict[str, Any]:
    outer, outer_status = _read_csv(outer_path, OUTER_REQUIRED)
    schur, schur_status = _read_csv(schur_path, SCHUR_REQUIRED)
    log_solves: list[dict[str, str]] = []
    log_hit_maxiter = 0
    if log_path and log_path.exists():
        for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
            match = SCHUR_RE.search(line)
            if match:
                log_solves.append({"iterations": match.group(1), "initial_residual": "", "final_residual": match.group(2)})
            if WARN_RE.search(line):
                log_hit_maxiter += 1
    schur_data = schur or log_solves
    initial = _float(outer[0], "initial_residual") if outer else None
    final = _float(outer[-1], "current_residual", "outer_residual", "solver_reported_residual") if outer else None
    invalid_rows = []
    invalid_timing_rows = []
    for index, row in enumerate(outer + schur):
        for field in ("initial_residual", "current_residual", "outer_residual", "solver_reported_residual", "final_residual"):
            if row.get(field, "") != "":
                value = _float(row, field)
                if not _finite_nonnegative(value):
                    invalid_rows.append({"row": index, "field": field, "value": row.get(field)})
        for field in ("elapsed_wall_seconds_cumulative", "elapsed_wall_seconds_per_call",
                      "elapsed_wall_seconds", "k_apply_seconds_cumulative", "k_apply_seconds_per_call",
                      "schur_action_seconds_cumulative", "schur_action_seconds_per_call"):
            if row.get(field, "") != "":
                value = _float(row, field)
                if not _finite_nonnegative(value):
                    invalid_timing_rows.append({"row": index, "field": field, "value": row.get(field)})
    log_reduction, log_reason = _log_reduction(initial, final)
    outer_iterations = [int(value) for value in (_float(row, "solver_reported_iteration") for row in outer) if value is not None]
    unique_outer_iterations = list(dict.fromkeys(outer_iterations))
    duplicate_outer_iterations = sorted({value for value in outer_iterations if outer_iterations.count(value) > 1})
    nonmonotonic_outer_iterations = any(b <= a for a, b in zip(outer_iterations, outer_iterations[1:]))
    sorted_outer_iterations = sorted(unique_outer_iterations)
    skipped_outer_iterations = [value for a, b in zip(sorted_outer_iterations, sorted_outer_iterations[1:])
                                for value in range(a + 1, b)]
    actual_outer = len(unique_outer_iterations) if unique_outer_iterations else None
    counter_rows = outer if any(row.get("k_actions_total", "") != "" for row in outer) else schur
    k = _k_stage_totals(counter_rows)
    elapsed_cumulative = _last(outer, "elapsed_wall_seconds_cumulative", "elapsed_wall_seconds")
    elapsed_per_call = _sum(outer, "elapsed_wall_seconds_per_call")
    k_time_cumulative = _last(schur, "k_apply_seconds_cumulative", "k_apply_seconds")
    k_time_per_call = _sum(schur, "k_apply_seconds_per_call")
    schur_time_cumulative = _last(schur, "schur_action_seconds_cumulative", "schur_action_seconds")
    cumulative_violations = []
    for label, rows, field in (
        ("outer_elapsed", outer, "elapsed_wall_seconds_cumulative"),
        ("schur_elapsed", schur, "elapsed_wall_seconds_cumulative"),
        ("k_apply", schur, "k_apply_seconds_cumulative"),
        ("schur_action", schur, "schur_action_seconds_cumulative"),
    ):
        previous = None
        for index, row in enumerate(rows):
            value = _float(row, field)
            if value is None:
                continue
            if previous is not None and value < previous:
                cumulative_violations.append({"series": label, "row": index,
                                              "previous": previous, "value": value})
            previous = value
    reduction_ratio = final / initial if _finite_nonnegative(initial) and _finite_nonnegative(final) and initial > 0 else None
    incomplete_reasons = []
    if outer_status["state"] != "OK":
        incomplete_reasons.append(f"outer_csv:{outer_status['state']}")
    if schur_path is not None and schur_status["state"] != "OK" and not log_solves:
        incomplete_reasons.append(f"schur_csv:{schur_status['state']}")
    if outer and final is None:
        incomplete_reasons.append("outer_residual_missing")
    if outer and any(_float(row, "initial_residual") is None for row in outer):
        incomplete_reasons.append("outer_initial_residual_missing")
    if schur and any(_float(row, "initial_residual") is None for row in schur):
        incomplete_reasons.append("schur_initial_residual_missing")
    if outer and not outer_iterations:
        incomplete_reasons.append("outer_iteration_missing")
    if duplicate_outer_iterations or skipped_outer_iterations or nonmonotonic_outer_iterations:
        incomplete_reasons.append("outer_iteration_sequence_invalid")
    if invalid_rows or invalid_timing_rows or cumulative_violations or k["total_mismatch"]:
        status = "FAIL"
    elif incomplete_reasons or not outer:
        status = "INCOMPLETE"
    else:
        status = "PASS"
    per_outer = log_reduction / actual_outer if log_reduction is not None and actual_outer and actual_outer > 0 else None
    per_k = log_reduction / k["total"] if log_reduction is not None and k["total"] and k["total"] > 0 else None
    per_second = log_reduction / elapsed_cumulative if log_reduction is not None and elapsed_cumulative and elapsed_cumulative > 0 else None
    return {
        "status": status,
        "incomplete_reasons": incomplete_reasons,
        "invalid_residual_rows": invalid_rows,
        "invalid_timing_rows": invalid_timing_rows,
        "cumulative_timing_violations": cumulative_violations,
        "outer_rows": len(outer),
        "schur_solves": len(schur_data),
        "schur_iteration_counts": [int(float(row["iterations"])) for row in schur_data if row.get("iterations", "") != ""],
        "schur_reached_tolerance_csv": sum(_bool(row, "reached_tolerance") is True for row in schur),
        "schur_hit_maxiter_csv": sum(_bool(row, "hit_maxiter") is True for row in schur),
        "schur_hit_maxiter_log": log_hit_maxiter,
        "schur_hit_maxiter": (sum(_bool(row, "hit_maxiter") is True for row in schur)
                              if schur else log_hit_maxiter) > 0,
        "schur_breakdown": sum(_bool(row, "breakdown") is True for row in schur),
        "schur_nonfinite": sum(_bool(row, "nonfinite") is True for row in schur),
        "schur_final_residual_min": min((value for value in (_float(row, "final_residual") for row in schur_data) if value is not None), default=None),
        "schur_final_residual_max": max((value for value in (_float(row, "final_residual") for row in schur_data) if value is not None), default=None),
        "k_actions_by_stage": k,
        "k_actions_total": k["total"],
        "k_actions_per_actual_outer_iteration": k["total"] / actual_outer if k["total"] is not None and actual_outer else None,
        "outer_first_residual": initial,
        "outer_last_residual": final,
        "outer_residual_reduction": reduction_ratio,
        "log10_reduction": log_reduction,
        "log10_reduction_reason": log_reason,
        "actual_outer_iterations": actual_outer,
        "unique_outer_iterations": unique_outer_iterations,
        "duplicate_outer_iterations": duplicate_outer_iterations,
        "skipped_outer_iterations": skipped_outer_iterations,
        "nonmonotonic_outer_iterations": nonmonotonic_outer_iterations,
        "k_stage_total_mismatch": k["total_mismatch"],
        "log10_reduction_per_outer_iteration": per_outer,
        "log10_reduction_per_k_action": per_k,
        "log10_reduction_per_second": per_second,
        "elapsed_wall_seconds": elapsed_cumulative,
        "elapsed_wall_seconds_cumulative_last": elapsed_cumulative,
        "elapsed_wall_seconds_per_call_sum": elapsed_per_call,
        "k_apply_seconds_cumulative_last": k_time_cumulative,
        "k_apply_seconds_per_call_sum": k_time_per_call,
        "schur_action_seconds_cumulative_last": schur_time_cumulative,
        "timing_definition": {
            "elapsed_wall_seconds": "last cumulative outer wall time; never a sum of cumulative samples",
            "*_per_call_sum": "sum of rows explicitly named per_call",
            "k_apply_seconds": "K solve/application inclusive time, excluding Schur vector algebra",
            "schur_action_seconds": "matrix-free Schur action inclusive time, including its K actions",
        },
        "workload_signature": _workload_signature(outer or schur),
        "csv_schema": {"outer": outer_status, "schur": schur_status},
        "source": {"outer_csv": str(outer_path) if outer_path else None, "schur_csv": str(schur_path) if schur_path else None, "log": str(log_path) if log_path else None},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-csv", type=Path)
    parser.add_argument("--schur-csv", type=Path)
    parser.add_argument("--log", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not any((args.outer_csv, args.schur_csv, args.log)):
        parser.error("at least one of --outer-csv, --schur-csv, or --log is required")
    result = summarize(args.outer_csv, args.schur_csv, args.log)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["status"] != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
