"""Evaluate diagnostic and production acceptance profiles independently."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

DEFAULT_POLICY = {
    "absolute_residual_max": 1.0e-12,
    "relative_residual_max": 1.0e-10,
    "backward_error_max": 1.0e-12,
    "constraint_absolute_residual_max": 1.0e-12,
    "relative_primal_error_max": 1.0e-4,
    "temperature_difference_max": 1.0e-3,
    "current_difference_max": 1.0e-3,
    "relative_residual_warning": 1.0e-11,
    "numerical_floor": 1.0e-14,
}


def _first(data: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in data:
            return data[name]
    return None


def _finite(value: Any) -> bool:
    try:
        return value is not None and math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _status(ok: bool, missing: bool = False) -> str:
    return "PASS" if ok else "INCOMPLETE" if missing else "FAIL"


def evaluate(metrics: dict[str, Any], policy: dict[str, float] | None = None,
             *, profile: str = "diagnostic") -> dict[str, Any]:
    if profile not in {"diagnostic", "production"}:
        raise ValueError("profile must be 'diagnostic' or 'production'")
    limits = dict(DEFAULT_POLICY)
    if policy:
        limits.update(policy)
    absolute = _first(metrics, "original_system_absolute_residual", "full_absolute_residual")
    relative = _first(metrics, "original_system_relative_residual", "full_relative_residual")
    constraint = _first(metrics, "constraint_absolute_residual", "absolute_constraint_residual")
    backward = _first(metrics, "backward_error", "original_system_backward_error")
    primal = _first(metrics, "relative_primal_error_vs_MUMPS", "relative_primal_agreement_with_mumps")
    temperature = _first(metrics, "tes_temperature_difference", "temperature_difference")
    current = _first(metrics, "tes_current_difference", "current_difference")
    constrained = bool(metrics.get("constrained", metrics.get("mortar", True)))
    no_mortar = bool(metrics.get("no_mortar", False))
    constraint_required = bool(metrics.get("constraint_metric_required", constrained and not no_mortar))
    nonfinite = bool(metrics.get("nonfinite", False) or metrics.get("nonfinite_or_breakdown", False))
    breakdown = bool(metrics.get("breakdown", False))
    hard_fail_reasons: list[str] = []
    if nonfinite:
        hard_fail_reasons.append("nonfinite")
    if breakdown:
        hard_fail_reasons.append("breakdown")
    values = {"absolute_residual": absolute, "relative_residual": relative, "constraint_residual": constraint,
              "backward_error": backward, "primal_comparison": primal, "tes_temperature_difference": temperature,
              "tes_current_difference": current}
    for label, value in values.items():
        if value is not None and not _finite(value):
            hard_fail_reasons.append(f"nonfinite {label}")
        elif value is not None and float(value) < 0.0:
            hard_fail_reasons.append(f"negative {label}")

    timing_fields = ("elapsed_wall_seconds", "elapsed_wall_seconds_cumulative",
                     "elapsed_wall_seconds_per_call_sum", "k_apply_seconds",
                     "k_apply_seconds_cumulative_last", "k_apply_seconds_per_call_sum",
                     "schur_action_seconds_cumulative_last", "k_actions_total")
    invalid_timing = []
    for field in timing_fields:
        value = metrics.get(field)
        if value is not None and (not _finite(value) or float(value) < 0.0):
            invalid_timing.append(field)
            hard_fail_reasons.append(f"invalid performance metric {field}")

    numerical_checks = {
        "absolute_residual": _finite(absolute) and float(absolute) <= limits["absolute_residual_max"],
        "relative_residual": _finite(relative) and float(relative) <= limits["relative_residual_max"],
        "backward_error": _finite(backward) and float(backward) <= limits["backward_error_max"],
        "constraint_absolute_residual": (not constraint_required) or (_finite(constraint) and float(constraint) <= limits["constraint_absolute_residual_max"]),
    }
    numerical_missing = [name for name, value in (("absolute_residual", absolute), ("relative_residual", relative), ("backward_error", backward)) if value is None]
    if constraint_required and constraint is None:
        numerical_missing.append("constraint_absolute_residual")
    physical_checks = {
        "primal_comparison": _finite(primal) and float(primal) <= limits["relative_primal_error_max"],
        "tes_temperature_difference": _finite(temperature) and float(temperature) <= limits["temperature_difference_max"],
        "tes_current_difference": _finite(current) and float(current) <= limits["current_difference_max"],
    }
    physical_missing = [name for name, value in (("primal_comparison", primal), ("tes_temperature_difference", temperature), ("tes_current_difference", current)) if value is None]
    floor_warning = _finite(absolute) and _finite(relative) and float(absolute) <= limits["numerical_floor"] and float(relative) > limits["relative_residual_warning"]
    relative_warning = _finite(relative) and float(relative) > limits["relative_residual_warning"]
    implementation_ok = not hard_fail_reasons
    numerical_ok = all(numerical_checks.values()) and not numerical_missing
    physical_ok = all(physical_checks.values()) and not physical_missing
    numerical_status = _status(numerical_ok, bool(numerical_missing))
    physical_status = _status(physical_ok, bool(physical_missing))
    performance_values_present = [metrics.get(field) for field in timing_fields]
    performance_missing = not any(value is not None for value in performance_values_present)
    performance_status = "FAIL" if invalid_timing else "INCOMPLETE" if performance_missing else "PASS"
    if profile == "diagnostic":
        production_ready = False
        overall = ("FAIL" if hard_fail_reasons or (not numerical_ok and not numerical_missing)
                   else "INCOMPLETE" if numerical_missing else "PASS")
    else:
        production_ready = implementation_ok and numerical_ok and physical_ok
        overall = (
            "FAIL" if hard_fail_reasons or (not numerical_ok and not numerical_missing)
            or (not physical_ok and not physical_missing)
            else "INCOMPLETE" if not production_ready else "PASS"
        )
    return {
        "profile": profile,
        "status": overall,
        "production_ready": production_ready,
        "implementation_correctness": {"status": "FAIL" if hard_fail_reasons else "PASS", "hard_fail_reasons": hard_fail_reasons},
        "numerical_convergence": {"status": numerical_status, "checks": numerical_checks, "missing_metrics": numerical_missing},
        "physical_acceptance": {"status": physical_status, "checks": physical_checks, "missing_metrics": physical_missing},
        "performance_readiness": {"status": performance_status, "invalid_metrics": invalid_timing},
        "gpu_acceleration": {"status": str(metrics.get("gpu_acceleration_status", "NOT RUN")), "speedup": metrics.get("gpu_speedup")},
        "constraint_metric_required": constraint_required,
        "no_mortar_constraint_exemption": no_mortar and not constraint_required,
        "metrics": values,
        "backward_error_norm_definition": metrics.get("backward_error_norm_definition", "componentwise max: |b-Ax|/(|A||x|+|b|), denominator elementwise"),
        "relative_residual_warning": relative_warning,
        "relative_residual_is_numerical_floor_warning": floor_warning,
        "policy": limits,
        "interpretation": "relative residual is diagnostic; numerical-floor warning does not override absolute/backward gates" if floor_warning else "profiles keep numerical convergence, physical acceptance, and performance separate",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", choices=("diagnostic", "production"), default="diagnostic")
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    result = evaluate(json.loads(args.input.read_text(encoding="utf-8")), profile=args.profile)
    if args.require_complete and result["status"] == "INCOMPLETE":
        result["status"] = "FAIL"
        result["implementation_correctness"]["hard_fail_reasons"].append("missing metrics")
        result["implementation_correctness"]["status"] = "FAIL"
        for category in ("numerical_convergence", "physical_acceptance", "performance_readiness"):
            if result[category]["status"] == "INCOMPLETE":
                result[category]["status"] = "FAIL"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if result["status"] != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
