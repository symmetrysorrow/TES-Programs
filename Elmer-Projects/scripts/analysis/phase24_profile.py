"""Build and validate Phase24 wall-clock profiling artifacts.

The profiler consumes explicit wall-clock events exported by an instrumented
Elmer run.  It deliberately does not infer wall time from CPU time.  Events
may be nested; child intervals are removed from their parent before totals are
reported so the accounting is additive.

Example::

    python scripts/analysis/phase24_profile.py profile \
      --events phase24_events.json --total-wall 12.5 \
      --case production_smoke --output artifacts/assembly_wall_profile.json

    python scripts/analysis/phase24_profile.py validate \
      --baseline baseline_gate.json --candidate cached_gate.json \
      --output artifacts/phase24_regression_gate.json
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


SCHEMA_VERSION = 1
NOISE_FRACTION = 0.02

COMPONENTS = (
    "element_traversal",
    "local_element_calculation",
    "material_property_evaluation",
    "mass_contribution",
    "stiffness_contribution",
    "boundary_element_assembly",
    "global_sparse_insertion",
    "rhs_insertion",
    "matrix_zero_reset",
    "matrix_format_conversion",
    "hypre_matrix_transfer_fill",
    "hypre_setup",
    "hypre_solve",
    "udf_circuit",
    "output_io",
    "unclassified",
)

ALIASES = {
    "element traversal": "element_traversal",
    "element_traversal": "element_traversal",
    "local element matrix/vector calculation": "local_element_calculation",
    "local element calculation": "local_element_calculation",
    "local_element_calculation": "local_element_calculation",
    "material-property evaluation": "material_property_evaluation",
    "material property evaluation": "material_property_evaluation",
    "material_property_evaluation": "material_property_evaluation",
    "mass contribution": "mass_contribution",
    "mass_contribution": "mass_contribution",
    "stiffness/conductivity contribution": "stiffness_contribution",
    "stiffness contribution": "stiffness_contribution",
    "stiffness_contribution": "stiffness_contribution",
    "boundary element assembly": "boundary_element_assembly",
    "boundary_element_assembly": "boundary_element_assembly",
    "global sparse insertion": "global_sparse_insertion",
    "global insertion": "global_sparse_insertion",
    "global_sparse_insertion": "global_sparse_insertion",
    "rhs insertion": "rhs_insertion",
    "rhs_insertion": "rhs_insertion",
    "matrix zero/reset": "matrix_zero_reset",
    "matrix zero reset": "matrix_zero_reset",
    "matrix_zero_reset": "matrix_zero_reset",
    "matrix format/conversion": "matrix_format_conversion",
    "matrix format conversion": "matrix_format_conversion",
    "matrix_format_conversion": "matrix_format_conversion",
    "hypre matrix transfer/fill": "hypre_matrix_transfer_fill",
    "hypre matrix transfer fill": "hypre_matrix_transfer_fill",
    "hypre_matrix_transfer_fill": "hypre_matrix_transfer_fill",
    "hypre setup": "hypre_setup",
    "hypre_setup": "hypre_setup",
    "hypre solve": "hypre_solve",
    "hypre_solve": "hypre_solve",
    "udf/circuit": "udf_circuit",
    "udf circuit": "udf_circuit",
    "udf_circuit": "udf_circuit",
    "output/io": "output_io",
    "output io": "output_io",
    "output_io": "output_io",
    "unclassified": "unclassified",
}


class ProfileError(ValueError):
    """Raised when an artifact would be misleading or internally invalid."""


def _read_json(path: Path) -> Any:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ProfileError(f"cannot read JSON {path}: {exc}") from exc


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def canonical_component(name: str) -> str:
    key = " ".join(name.strip().lower().replace("-", " ").split())
    return ALIASES.get(key, key.replace(" ", "_"))


def _interval(event: dict[str, Any]) -> tuple[float, float] | None:
    if "start_s" not in event or "end_s" not in event:
        return None
    start = float(event["start_s"])
    end = float(event["end_s"])
    if not (math.isfinite(start) and math.isfinite(end)) or end < start:
        raise ProfileError(f"invalid event interval: {event!r}")
    return start, end


def _union_length(intervals: Iterable[tuple[float, float]]) -> float:
    ordered = sorted((a, b) for a, b in intervals if b > a)
    total = 0.0
    if not ordered:
        return total
    left, right = ordered[0]
    for start, end in ordered[1:]:
        if start <= right:
            right = max(right, end)
        else:
            total += right - left
            left, right = start, end
    return total + right - left


def exclusive_event_seconds(events: list[dict[str, Any]]) -> dict[str, float]:
    """Return exclusive component totals from wall-clock intervals.

    Each event's descendants are subtracted by interval union.  Events with a
    direct ``wall_seconds`` field are accepted for exporters that do not expose
    timestamps, but such events must not be marked as CPU-only.
    """
    totals: defaultdict[str, float] = defaultdict(float)
    for index, event in enumerate(events):
        if not isinstance(event, dict):
            raise ProfileError(f"event {index} is not an object")
        if "cpu_seconds" in event and "wall_seconds" not in event and _interval(event) is None:
            raise ProfileError("CPU-only timing is not valid for a wall profile")
        name = canonical_component(str(event.get("component", event.get("name", "unclassified"))))
        if name not in COMPONENTS:
            raise ProfileError(f"unknown component {event.get('component', event.get('name'))!r}")
        iv = _interval(event)
        if iv is not None:
            children = []
            for child in events:
                child_iv = _interval(child)
                if child is event or child_iv is None:
                    continue
                if child_iv[0] >= iv[0] and child_iv[1] <= iv[1]:
                    children.append(child_iv)
            seconds = (iv[1] - iv[0]) - _union_length(children)
        elif "wall_seconds" in event:
            seconds = float(event["wall_seconds"])
        else:
            raise ProfileError(f"event has neither interval nor wall_seconds: {event!r}")
        if not math.isfinite(seconds) or seconds < 0:
            raise ProfileError(f"invalid wall_seconds in event: {event!r}")
        totals[name] += seconds
    return dict(totals)


def make_profile(events_payload: Any, *, case: str, total_wall: float, metadata: dict[str, Any] | None = None) -> dict[str, Any]:
    if isinstance(events_payload, dict):
        events = events_payload.get("events")
    else:
        events = events_payload
    if not isinstance(events, list):
        raise ProfileError("events JSON must be a list or an object containing 'events'")
    total_wall = float(total_wall)
    if not math.isfinite(total_wall) or total_wall < 0:
        raise ProfileError("total wall time must be finite and non-negative")
    components = {name: 0.0 for name in COMPONENTS}
    components.update(exclusive_event_seconds(events))
    measured = sum(components[name] for name in COMPONENTS if name != "unclassified")
    unclassified = max(total_wall - measured, 0.0)
    components["unclassified"] = unclassified
    accounting_error = measured + unclassified - total_wall
    result = {
        "schema_version": SCHEMA_VERSION,
        "artifact": "assembly_wall_profile",
        "case": case,
        "clock": "wall",
        "cpu_time_used_for_accounting": False,
        "total_solver_wall_seconds": total_wall,
        "components_seconds": components,
        "measured_components_seconds": measured,
        "unclassified_seconds": unclassified,
        "accounting_error_seconds": accounting_error,
        "accounting_ok": abs(accounting_error) <= max(1.0e-9, total_wall * 1.0e-9),
        "metadata": metadata or {},
    }
    return result


PHASE24_LOG_EVENT = re.compile(
    r"component=(?P<component>[A-Za-z0-9_ /-]+?)\s+wall_seconds=\s*(?P<seconds>[-+0-9.eE]+)"
)


def parse_phase24_log(path: Path) -> list[dict[str, Any]]:
    """Extract additive PHASE24_EVENT records emitted by HeatSolve.F90."""
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise ProfileError(f"cannot read log {path}: {exc}") from exc
    events: list[dict[str, Any]] = []
    for line in lines:
        if "PHASE24_EVENT" not in line:
            continue
        match = PHASE24_LOG_EVENT.search(line)
        if not match:
            raise ProfileError(f"malformed PHASE24_EVENT line: {line!r}")
        events.append({"component": match.group("component"), "wall_seconds": float(match.group("seconds"))})
    if not events:
        raise ProfileError(f"no PHASE24_EVENT records found in {path}")
    return events


def _number(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ProfileError(f"{label} must be numeric") from exc
    if not math.isfinite(result):
        raise ProfileError(f"{label} must be finite")
    return result


def _extract_gate(payload: dict[str, Any]) -> dict[str, Any]:
    gate = payload.get("gate", payload)
    if not isinstance(gate, dict):
        raise ProfileError("gate artifact must be an object")
    return gate


def validate_gate(baseline_payload: dict[str, Any], candidate_payload: dict[str, Any], tolerance: float) -> dict[str, Any]:
    baseline = _extract_gate(baseline_payload)
    candidate = _extract_gate(candidate_payload)
    required = ("matrix_sparsity_identical", "matrix_values_identical", "rhs_identical", "temperature_field_identical", "tes_observables_identical", "solver_converged")
    checks = {name: bool(candidate.get(name, False)) for name in required}
    numeric = candidate.get("numeric_parity", {})
    if not isinstance(numeric, dict):
        raise ProfileError("numeric_parity must be an object")
    numeric_errors = {str(k): _number(v, f"numeric_parity.{k}") for k, v in numeric.items()}
    numeric_ok = all(abs(value) <= tolerance for value in numeric_errors.values())
    checks["numeric_parity_within_tolerance"] = numeric_ok
    passed = all(checks.values())
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact": "phase24_regression_gate",
        "baseline_case": baseline.get("case"),
        "candidate_case": candidate.get("case"),
        "tolerance": tolerance,
        "checks": checks,
        "numeric_parity_errors": numeric_errors,
        "passed": passed,
        "rejection_reason": None if passed else "correctness gate failed",
    }


def classify(baseline_profile: dict[str, Any], candidate_profile: dict[str, Any], gate: dict[str, Any]) -> dict[str, Any]:
    if not bool(gate.get("passed", False)):
        classification = "REJECTED_CORRECTNESS_GATE"
    else:
        base = baseline_profile.get("components_seconds", {})
        cand = candidate_profile.get("components_seconds", {})
        def speedup(name: str) -> float:
            old = float(base.get(name, 0.0))
            new = float(cand.get(name, 0.0))
            return (old - new) / old if old > 0 else 0.0
        total_old = float(baseline_profile["total_solver_wall_seconds"])
        total_new = float(candidate_profile["total_solver_wall_seconds"])
        total_gain = (total_old - total_new) / total_old if total_old > 0 else 0.0
        gains = {name: speedup(name) for name in COMPONENTS}
        dominant = max((name for name in COMPONENTS if name != "unclassified"), key=lambda name: float(base.get(name, 0.0)), default="unclassified")
        if total_gain <= NOISE_FRACTION:
            classification = "SPARSE_INSERTION_LIMITED" if dominant == "global_sparse_insertion" else "ELMER_ARCHITECTURE_LIMITED"
        elif gains.get("global_sparse_insertion", 0.0) > NOISE_FRACTION:
            classification = "ASSEMBLY_CACHE_EFFECTIVE"
        elif gains.get("local_element_calculation", 0.0) > NOISE_FRACTION:
            classification = "CPU_PARALLEL_ASSEMBLY_EFFECTIVE"
        elif gains.get("matrix_zero_reset", 0.0) > NOISE_FRACTION or gains.get("hypre_matrix_transfer_fill", 0.0) > NOISE_FRACTION:
            classification = "STATIC_REUSE_EFFECTIVE"
        else:
            classification = "ASSEMBLY_AND_REUSE_EFFECTIVE" if gains.get("global_sparse_insertion", 0.0) > 0 or gains.get("local_element_calculation", 0.0) > 0 else "MEMORY_BANDWIDTH_LIMITED"
    return {
        "schema_version": SCHEMA_VERSION,
        "artifact": "phase24_backend_recommendation",
        "classification": classification,
        "production_backend": "cpu_conformal_hypre",
        "gpu_path_preserved": True,
        "noise_fraction": NOISE_FRACTION,
    }


def command_profile(args: argparse.Namespace) -> None:
    profile = make_profile(_read_json(Path(args.events)), case=args.case, total_wall=args.total_wall, metadata={"source": str(args.events)})
    _write_json(Path(args.output), profile)


def command_log(args: argparse.Namespace) -> None:
    events = parse_phase24_log(Path(args.log))
    profile = make_profile(events, case=args.case, total_wall=args.total_wall, metadata={"source": str(args.log)})
    _write_json(Path(args.output), profile)


def command_validate(args: argparse.Namespace) -> None:
    result = validate_gate(_read_json(Path(args.baseline)), _read_json(Path(args.candidate)), args.tolerance)
    _write_json(Path(args.output), result)
    if not result["passed"]:
        raise SystemExit(2)


def command_classify(args: argparse.Namespace) -> None:
    result = classify(_read_json(Path(args.baseline_profile)), _read_json(Path(args.candidate_profile)), _read_json(Path(args.gate)))
    _write_json(Path(args.output), result)


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__)
    sub = root.add_subparsers(dest="command", required=True)
    p = sub.add_parser("profile", help="create an additive wall-clock profile")
    p.add_argument("--events", required=True, type=Path)
    p.add_argument("--total-wall", required=True, type=float)
    p.add_argument("--case", required=True)
    p.add_argument("--output", required=True, type=Path)
    p.set_defaults(func=command_profile)
    l = sub.add_parser("log", help="create a profile from HeatSolve PHASE24_EVENT log lines")
    l.add_argument("--log", required=True, type=Path)
    l.add_argument("--total-wall", required=True, type=float)
    l.add_argument("--case", required=True)
    l.add_argument("--output", required=True, type=Path)
    l.set_defaults(func=command_log)
    v = sub.add_parser("validate", help="run the matrix/RHS/field correctness gate")
    v.add_argument("--baseline", required=True, type=Path)
    v.add_argument("--candidate", required=True, type=Path)
    v.add_argument("--tolerance", type=float, default=1.0e-10)
    v.add_argument("--output", required=True, type=Path)
    v.set_defaults(func=command_validate)
    c = sub.add_parser("classify", help="classify a validated before/after result")
    c.add_argument("--baseline-profile", required=True, type=Path)
    c.add_argument("--candidate-profile", required=True, type=Path)
    c.add_argument("--gate", required=True, type=Path)
    c.add_argument("--output", required=True, type=Path)
    c.set_defaults(func=command_classify)
    return root


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    try:
        args.func(args)
    except ProfileError as exc:
        print(f"phase24_profile: error: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
