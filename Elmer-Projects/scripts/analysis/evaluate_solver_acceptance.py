"""Evaluate Phase19/20 solver metrics with separate numerical and physical gates.

The evaluator intentionally does not make a small right-hand-side relative
residual the sole production criterion.  It accepts either a metrics JSON
object on stdin/file or a JSON artifact produced by a solver comparison.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


DEFAULT_POLICY = {
    "absolute_residual_max": 1.0e-12,
    "backward_error_max": 1.0e-12,
    "constraint_absolute_residual_max": 1.0e-12,
    "relative_primal_error_max": 1.0e-4,
    "relative_residual_warning": 1.0e-11,
    "numerical_floor": 1.0e-14,
}


def _first(data: dict[str, Any], *names: str) -> Any:
    for name in names:
        if name in data:
            return data[name]
    return None


def _finite(value: Any) -> bool:
    return value is not None and math.isfinite(float(value))


def evaluate(metrics: dict[str, Any], policy: dict[str, float] | None = None) -> dict[str, Any]:
    limits = dict(DEFAULT_POLICY)
    if policy:
        limits.update(policy)

    absolute = _first(metrics, "original_system_absolute_residual", "full_absolute_residual")
    relative = _first(metrics, "original_system_relative_residual", "full_relative_residual")
    constraint = _first(metrics, "constraint_absolute_residual", "absolute_constraint_residual")
    backward = _first(metrics, "backward_error", "original_system_backward_error")
    primal = _first(metrics, "relative_primal_error_vs_MUMPS", "relative_primal_agreement_with_mumps")

    nonfinite = bool(metrics.get("nonfinite", False) or metrics.get("nonfinite_or_breakdown", False))
    breakdown = bool(metrics.get("breakdown", False))
    hard_fail_reasons: list[str] = []
    if nonfinite:
        hard_fail_reasons.append("nonfinite")
    if breakdown:
        hard_fail_reasons.append("breakdown")
    for label, value in (("absolute residual", absolute), ("constraint residual", constraint),
                         ("backward error", backward), ("primal comparison", primal)):
        if value is not None and not _finite(value):
            hard_fail_reasons.append(f"nonfinite {label}")

    numerical_checks = {
        "absolute_residual": absolute is not None and _finite(absolute)
        and float(absolute) <= limits["absolute_residual_max"],
        "constraint_absolute_residual": constraint is not None and _finite(constraint)
        and float(constraint) <= limits["constraint_absolute_residual_max"],
        "backward_error": backward is not None and _finite(backward)
        and float(backward) <= limits["backward_error_max"],
    }
    # A missing backward error is a diagnostic limitation, not a fabricated
    # zero.  Callers can require it with require_complete=True.
    missing_primary = [name for name, value in (
        ("absolute_residual", absolute), ("constraint_absolute_residual", constraint),
        ("backward_error", backward),
    ) if value is None]
    primal_check = primal is None or float(primal) <= limits["relative_primal_error_max"]
    floor_warning = (
        _finite(absolute) and _finite(relative)
        and float(absolute) <= limits["numerical_floor"]
        and float(relative) > limits["relative_residual_warning"]
    )
    relative_warning = _finite(relative) and float(relative) > limits["relative_residual_warning"]

    complete = not missing_primary
    numerical_pass = all(numerical_checks.values()) and primal_check
    production_ready = not hard_fail_reasons and complete and numerical_pass
    return {
        "status": "PASS" if production_ready else "FAIL" if hard_fail_reasons or (complete and not numerical_pass) else "INCOMPLETE",
        "production_ready": production_ready,
        "hard_fail_reasons": hard_fail_reasons,
        "metrics": {
            "original_system_absolute_residual": absolute,
            "original_system_relative_residual": relative,
            "constraint_absolute_residual": constraint,
            "backward_error": backward,
            "relative_primal_error_vs_MUMPS": primal,
        },
        "checks": numerical_checks | {"primal_comparison": primal_check},
        "missing_primary_metrics": missing_primary,
        "relative_residual_warning": relative_warning,
        "relative_residual_is_numerical_floor_warning": floor_warning,
        "policy": limits,
        "interpretation": (
            "relative residual exceeded its diagnostic threshold, but the absolute residual "
            "is at the configured numerical floor; this is not a hard failure"
            if floor_warning else
            "relative residual is diagnostic only; production readiness is determined by "
            "complete absolute/backward/constraint and primal checks"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    payload = json.loads(args.input.read_text(encoding="utf-8"))
    result = evaluate(payload)
    if args.require_complete and result["missing_primary_metrics"] and result["status"] == "INCOMPLETE":
        result["status"] = "FAIL"
        result["hard_fail_reasons"].append("missing primary metrics")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
