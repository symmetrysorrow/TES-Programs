from __future__ import annotations

import ast
from copy import deepcopy
from typing import Any

from scripts.support.vendored.dimensioned_expression import (
    _replace_unit_literals,
    evaluate_dimensioned_expression,
)
from scripts.support.vendored.expression import evaluate_expression

_DIMENSIONED_BUILTIN_NAMES = {"abs", "min", "max", "pow", "__unit__"}


def _round_noise(value: float) -> float:
    # Round away float64 arithmetic noise (e.g. 715*1e-6 != 0.000715) while
    # preserving far more precision than this simulation needs.
    return float(f"{value:.15g}")


def _dimensioned_dependencies(expression: str) -> set[str]:
    processed = _replace_unit_literals(expression)
    tree = ast.parse(processed, mode="eval")
    names = {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
    return names - _DIMENSIONED_BUILTIN_NAMES


def resolve_parameters(parameter_expressions: dict[str, str]) -> dict[str, float]:
    """Evaluate parameter_expressions into SI floats, resolving cross-references
    (e.g. I_0 = "I_bias*R_sh/(R_0+R_sh)") in dependency order via fixed-point
    iteration. Raises ValueError on a circular or missing dependency."""
    resolved: dict[str, float] = {}
    remaining = dict(parameter_expressions)

    progress = True
    while remaining and progress:
        progress = False
        for name, expr in list(remaining.items()):
            try:
                needed = _dimensioned_dependencies(expr)
            except SyntaxError as exc:
                raise ValueError(f"Could not parse expression for '{name}': {expr!r} ({exc})") from exc
            if needed <= resolved.keys():
                value = evaluate_dimensioned_expression(expr, variables=resolved).value
                resolved[name] = _round_noise(value)
                del remaining[name]
                progress = True

    if remaining:
        raise ValueError(
            "Could not resolve parameter expressions (circular or missing dependency): "
            f"{sorted(remaining)}"
        )
    return resolved


def _evaluate_material_scalar(expression: str, context: dict[str, float]) -> float | None:
    try:
        return _round_noise(float(evaluate_expression(expression, variables=context)))
    except Exception:
        return None


def _reconcile_materials(materials: dict[str, Any], context: dict[str, float]) -> None:
    for material in materials.values():
        for key in ("k", "rho", "cp"):
            prop = material.get(key)
            if not isinstance(prop, dict):
                continue
            expression = str(prop.get("expression", "") or "").strip()
            if not expression:
                continue
            value = _evaluate_material_scalar(expression, context)
            if value is not None:
                prop["nominal"] = value


def _reconcile_expr_fields(node: Any, context: dict[str, float]) -> None:
    if isinstance(node, dict):
        for key in list(node.keys()):
            if not key.endswith("_expr"):
                continue
            base_key = key[: -len("_expr")]
            expression = str(node.get(key) or "").strip()
            if not expression:
                continue
            try:
                value = evaluate_dimensioned_expression(expression, variables=context).value
            except Exception:
                continue
            node[base_key] = _round_noise(value)
        for value in node.values():
            _reconcile_expr_fields(value, context)
    elif isinstance(node, list):
        for item in node:
            _reconcile_expr_fields(item, context)


def reconcile_project(model: dict[str, Any]) -> dict[str, Any]:
    """Return a deep copy of *model* with every numeric field re-derived from
    its paired expression field (parameter_expressions -> parameters,
    materials.*.expression -> materials.*.nominal, geometries.*.*_expr -> the
    matching numeric sibling, for every tree in the `geometries` registry).
    expression/_expr fields are treated as the
    single source of truth; fields whose expression cannot be evaluated
    (e.g. temperature-dependent k(T) laws with a free variable T) are left
    untouched rather than erroring out.

    NOTE: since Phase 0 of the redesign this function no longer writes
    anything back to disk. Generation must use the returned copy; the stored
    numeric fields are informational and checked by `numeric_drift_report`."""
    reconciled = deepcopy(model)

    parameter_expressions = reconciled.get("parameter_expressions", {})
    if parameter_expressions:
        resolved_parameters = resolve_parameters(parameter_expressions)
        reconciled["parameters"] = resolved_parameters
    else:
        resolved_parameters = dict(reconciled.get("parameters", {}))

    _reconcile_materials(reconciled.get("materials", {}), resolved_parameters)
    for geometry_tree in reconciled.get("geometries", {}).values():
        _reconcile_expr_fields(geometry_tree, resolved_parameters)

    return reconciled


def _diff_paths(stored: Any, derived: Any, path: str, out: list[str]) -> None:
    if isinstance(stored, dict) and isinstance(derived, dict):
        for key in derived:
            # Keys absent from the file are derived-only by design
            # (expressions-only schema); only real mismatches are drift.
            if key in stored:
                _diff_paths(stored[key], derived[key], f"{path}.{key}" if path else str(key), out)
        return
    if isinstance(stored, list) and isinstance(derived, list):
        for i, (s, d) in enumerate(zip(stored, derived)):
            _diff_paths(s, d, f"{path}[{i}]", out)
        return
    if isinstance(stored, float) and isinstance(derived, float):
        if stored != derived:
            out.append(f"{path}: file has {stored!r}, expressions give {derived!r}")
        return
    if stored != derived:
        out.append(f"{path}: file has {stored!r}, expressions give {derived!r}")


def numeric_drift_report(model: dict[str, Any], reconciled: dict[str, Any]) -> list[str]:
    """List numeric fields whose stored value disagrees with the value derived
    from its expression. Empty list means the file is internally consistent."""
    out: list[str] = []
    for section in ("parameters", "materials"):
        if section in model and section in reconciled:
            _diff_paths(model[section], reconciled[section], section, out)
    for geometry_name in model.get("geometries", {}):
        if geometry_name in reconciled.get("geometries", {}):
            _diff_paths(
                model["geometries"][geometry_name],
                reconciled["geometries"][geometry_name],
                f"geometries.{geometry_name}",
                out,
            )
    return out
