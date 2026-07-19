"""Dimensioned expression evaluation.

Vendored 2026-07-14 from Thermal-and-Electoric-Sim
(core/project/dimensioned_expression.py); only this header and the import
path changed.
"""
from __future__ import annotations

import ast
import re
from dataclasses import dataclass
from typing import Any, Mapping

from scripts.support.vendored.unit_utils import (
    dimension_name_from_signature,
    dimension_signature_for_name,
    normalize_unit,
    unit_factor_and_dimension,
)

_UNIT_LITERAL_PATTERN = re.compile(
    r"(?<![\w.])"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"\s*\[\s*([a-zA-Z0-9Ωµμ/\^\*\(\)\-]+)\s*\]"
)

_BARE_UNIT_LITERAL_PATTERN = re.compile(
    r"(?<![\w.\[])"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?!\s*\[)"
    r"(?![eE][+-]?\d)"
    r"\s*([a-zA-ZΩµμ][a-zA-Z0-9Ωµμ]*)"
    r"(?![\w])"
)


def _format_dimensions(dimensions: tuple[tuple[str, int], ...]) -> str:
    if not dimensions:
        return "dimensionless"
    named = dimension_name_from_signature(dimensions)
    if named is not None:
        return named
    parts: list[str] = []
    for dim, power in dimensions:
        if power == 1:
            parts.append(dim)
        else:
            parts.append(f"{dim}^{power}")
    return "*".join(parts)


def _format_quantity(quantity: "DimensionedQuantity") -> str:
    return f"{quantity.value:.6g}[{_format_dimensions(quantity.dimensions)}]"


def _dimension_mismatch_message(operation: str, left: "DimensionedQuantity", right: "DimensionedQuantity") -> str:
    return (
        f"Cannot {operation} values with different dimensions: "
        f"left={_format_quantity(left)}, right={_format_quantity(right)}."
    )


@dataclass(frozen=True)
class DimensionedQuantity:
    value: float
    dimensions: tuple[tuple[str, int], ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "value", float(self.value))
        object.__setattr__(self, "dimensions", _normalize_dimensions(self.dimensions))

    @staticmethod
    def dimensionless(value: float) -> "DimensionedQuantity":
        return DimensionedQuantity(float(value), ())

    @staticmethod
    def from_value(value: Any) -> "DimensionedQuantity":
        if isinstance(value, DimensionedQuantity):
            return value
        if isinstance(value, tuple) and len(value) == 2:
            raw_value, raw_dimension = value
            dims = _normalize_dimensions(raw_dimension)
            return DimensionedQuantity(float(raw_value), dims)
        if isinstance(value, (int, float, bool)):
            return DimensionedQuantity.dimensionless(float(value))
        return DimensionedQuantity.dimensionless(float(value))

    def _combine(self, other: Any, value: float, *, scale_self: float = 1.0, scale_other: float = 1.0) -> "DimensionedQuantity":
        rhs = DimensionedQuantity.from_value(other)
        dims = _merge_dimensions(self.dimensions, rhs.dimensions, scale_self=scale_self, scale_other=scale_other)
        return DimensionedQuantity(value, dims)

    def __add__(self, other: Any) -> "DimensionedQuantity":
        rhs = DimensionedQuantity.from_value(other)
        if self.dimensions != rhs.dimensions:
            raise ValueError(_dimension_mismatch_message("add", self, rhs))
        return DimensionedQuantity(self.value + rhs.value, self.dimensions)

    def __radd__(self, other: Any) -> "DimensionedQuantity":
        return self.__add__(other)

    def __sub__(self, other: Any) -> "DimensionedQuantity":
        rhs = DimensionedQuantity.from_value(other)
        if self.dimensions != rhs.dimensions:
            raise ValueError(_dimension_mismatch_message("subtract", self, rhs))
        return DimensionedQuantity(self.value - rhs.value, self.dimensions)

    def __rsub__(self, other: Any) -> "DimensionedQuantity":
        lhs = DimensionedQuantity.from_value(other)
        return lhs.__sub__(self)

    def __mul__(self, other: Any) -> "DimensionedQuantity":
        rhs = DimensionedQuantity.from_value(other)
        return DimensionedQuantity(self.value * rhs.value, _merge_dimensions(self.dimensions, rhs.dimensions))

    def __rmul__(self, other: Any) -> "DimensionedQuantity":
        return self.__mul__(other)

    def __truediv__(self, other: Any) -> "DimensionedQuantity":
        rhs = DimensionedQuantity.from_value(other)
        return DimensionedQuantity(self.value / rhs.value, _merge_dimensions(self.dimensions, rhs.dimensions, scale_other=-1))

    def __rtruediv__(self, other: Any) -> "DimensionedQuantity":
        lhs = DimensionedQuantity.from_value(other)
        return lhs.__truediv__(self)

    def __pow__(self, other: Any) -> "DimensionedQuantity":
        rhs = DimensionedQuantity.from_value(other)
        if rhs.dimensions:
            raise ValueError("Exponent must be dimensionless.")
        exponent = float(rhs.value)
        if self.dimensions and not exponent.is_integer():
            raise ValueError("Only integer powers are supported for dimensioned values.")
        if self.dimensions:
            scaled = tuple((dim, int(power * exponent)) for dim, power in self.dimensions)
            return DimensionedQuantity(self.value ** exponent, _normalize_dimensions(scaled))
        return DimensionedQuantity(self.value ** exponent, ())

    def __neg__(self) -> "DimensionedQuantity":
        return DimensionedQuantity(-self.value, self.dimensions)

    def __pos__(self) -> "DimensionedQuantity":
        return self

    def __abs__(self) -> "DimensionedQuantity":
        return DimensionedQuantity(abs(self.value), self.dimensions)

    @property
    def is_dimensionless(self) -> bool:
        return not self.dimensions


def _normalize_dimensions(raw: Any) -> tuple[tuple[str, int], ...]:
    if raw is None:
        return ()
    if isinstance(raw, str):
        text = raw.strip().lower()
        if not text or text in {"dimensionless", "none", "-", "[]"}:
            return ()
        signature = dimension_signature_for_name(text)
        if signature is not None:
            return signature
        return ((text, 1),)
    if isinstance(raw, Mapping):
        merged: dict[str, int] = {}
        for key, value in raw.items():
            power = int(value)
            if power == 0:
                continue
            key_text = str(key).strip().lower()
            signature = dimension_signature_for_name(key_text)
            if signature is not None and signature and len(signature) == 1 and signature[0][1] == 1:
                base_dim = signature[0][0]
                merged[base_dim] = merged.get(base_dim, 0) + power
            elif signature is not None:
                for dim_name, dim_power in signature:
                    merged[dim_name] = merged.get(dim_name, 0) + dim_power * power
            else:
                merged[key_text] = merged.get(key_text, 0) + power
        items = [(dim, power) for dim, power in merged.items() if power != 0]
        return tuple(sorted(items))
    if isinstance(raw, tuple):
        merged: dict[str, int] = {}
        for item in raw:
            if isinstance(item, tuple) and len(item) == 2:
                dim, power = item
                power_i = int(power)
                if power_i != 0:
                    dim_text = str(dim).strip().lower()
                    signature = dimension_signature_for_name(dim_text)
                    if signature is not None and signature and len(signature) == 1 and signature[0][1] == 1:
                        base_dim = signature[0][0]
                        merged[base_dim] = merged.get(base_dim, 0) + power_i
                    elif signature is not None:
                        for dim_name, dim_power in signature:
                            merged[dim_name] = merged.get(dim_name, 0) + dim_power * power_i
                    else:
                        merged[dim_text] = merged.get(dim_text, 0) + power_i
        items = [(dim, power) for dim, power in merged.items() if power != 0]
        return tuple(sorted(items))
    return ()


def _merge_dimensions(
    left: tuple[tuple[str, int], ...],
    right: tuple[tuple[str, int], ...],
    *,
    scale_self: float = 1.0,
    scale_other: float = 1.0,
) -> tuple[tuple[str, int], ...]:
    merged: dict[str, int] = {dim: power for dim, power in left}
    for dim, power in right:
        merged[dim] = merged.get(dim, 0) + int(power * scale_other if scale_other != 1.0 else power)
    if scale_self != 1.0:
        merged = {dim: int(power * scale_self) for dim, power in merged.items()}
    return tuple(sorted((dim, power) for dim, power in merged.items() if power != 0))


def _dimension_name(dimensions: tuple[tuple[str, int], ...]) -> str | None:
    return dimension_name_from_signature(dimensions)


def _check_expected_dimension(result: DimensionedQuantity, expected_dimension: str | None) -> None:
    if expected_dimension is None or not str(expected_dimension).strip():
        return
    expected = str(expected_dimension).strip().lower()
    if expected in {"dimensionless", "none", "-", "[]"}:
        if result.dimensions:
            raise ValueError("Expected a dimensionless expression.")
        return
    dim_name = _dimension_name(result.dimensions)
    if dim_name is None:
        raise ValueError(f"Expression has composite dimensions; expected {expected}.")
    if dim_name != expected:
        raise ValueError(f"Expression has {dim_name} dimensions; expected {expected}.")


def _replace_unit_literals(expression: str) -> str:
    raw = str(expression or "")
    if not raw.strip():
        return raw

    raw_without_bracket_units = _UNIT_LITERAL_PATTERN.sub("", raw)
    bare_match = _BARE_UNIT_LITERAL_PATTERN.search(raw_without_bracket_units)
    if bare_match:
        unit = bare_match.group(2).strip().replace("µ", "u").replace("μ", "u").lower()
        raise ValueError(f"Units must use [unit] notation; found bare unit '{unit}'.")

    def repl(match: re.Match[str]) -> str:
        number = match.group(1)
        unit = normalize_unit(match.group(2))
        return f"__unit__({number!s}, {unit!r})"

    return _UNIT_LITERAL_PATTERN.sub(repl, raw)


class _DimensionedEvaluator(ast.NodeVisitor):
    def __init__(self, variables: Mapping[str, Any]) -> None:
        self.variables = dict(variables)

    def visit_Expression(self, node: ast.Expression) -> DimensionedQuantity:  # pragma: no cover - entrypoint
        return self.visit(node.body)

    def visit_Constant(self, node: ast.Constant) -> DimensionedQuantity:
        if isinstance(node.value, (int, float, bool)):
            return DimensionedQuantity.dimensionless(float(node.value))
        raise ValueError(f"Unsupported constant type: {type(node.value).__name__}")

    def visit_Name(self, node: ast.Name) -> DimensionedQuantity:
        if node.id in self.variables:
            value = self.variables[node.id]
            if isinstance(value, Mapping):
                return value  # type: ignore[return-value]
            return DimensionedQuantity.from_value(value)
        raise ValueError(f"Unknown symbol: {node.id}")

    def visit_Attribute(self, node: ast.Attribute) -> Any:
        value = self.visit(node.value)
        if isinstance(value, DimensionedQuantity):
            raise ValueError(f"Unknown attribute: {node.attr}")
        if isinstance(value, Mapping) and node.attr in value:
            resolved = value[node.attr]
            if isinstance(resolved, Mapping):
                return resolved
            return DimensionedQuantity.from_value(resolved)
        if hasattr(value, node.attr):
            return DimensionedQuantity.from_value(getattr(value, node.attr))
        raise ValueError(f"Unknown attribute: {node.attr}")

    def visit_BinOp(self, node: ast.BinOp) -> DimensionedQuantity:
        left = self.visit(node.left)
        right = self.visit(node.right)
        if isinstance(node.op, ast.Add):
            return left + right
        if isinstance(node.op, ast.Sub):
            return left - right
        if isinstance(node.op, ast.Mult):
            return left * right
        if isinstance(node.op, ast.Div):
            return left / right
        if isinstance(node.op, ast.Pow):
            return left**right
        raise ValueError(f"Unsupported binary operator: {type(node.op).__name__}")

    def visit_UnaryOp(self, node: ast.UnaryOp) -> DimensionedQuantity:
        operand = self.visit(node.operand)
        if isinstance(node.op, ast.UAdd):
            return +operand
        if isinstance(node.op, ast.USub):
            return -operand
        raise ValueError(f"Unsupported unary operator: {type(node.op).__name__}")

    def visit_Call(self, node: ast.Call) -> DimensionedQuantity:
        if not isinstance(node.func, ast.Name):
            raise ValueError("Only direct function calls are supported.")
        func_name = node.func.id
        if func_name == "__unit__":
            if len(node.args) != 2:
                raise ValueError("Invalid unit literal.")
            value = self.visit(node.args[0])
            unit_name = node.args[1]
            if not isinstance(unit_name, ast.Constant) or not isinstance(unit_name.value, str):
                raise ValueError("Invalid unit literal.")
            unit = unit_name.value
            factor_and_dim = unit_factor_and_dimension(unit)
            if factor_and_dim is None:
                raise ValueError(f"Unknown unit: {unit}")
            dimensions, factor = factor_and_dim
            if value.dimensions:
                raise ValueError("Unit literals must be attached to plain numbers only.")
            return DimensionedQuantity(float(value.value) * factor, dimensions)
        if func_name == "abs":
            args = [self.visit(arg) for arg in node.args]
            if len(args) != 1:
                raise ValueError("abs() takes exactly one argument.")
            return abs(args[0])
        if func_name in {"min", "max"}:
            args = [self.visit(arg) for arg in node.args]
            if not args:
                raise ValueError(f"{func_name}() requires at least one argument.")
            first = args[0]
            for arg in args[1:]:
                if arg.dimensions != first.dimensions:
                    formatted = ", ".join(
                        f"arg{i + 1}={_format_quantity(arg_value)}"
                        for i, arg_value in enumerate(args)
                    )
                    raise ValueError(
                        f"{func_name}() arguments must have the same dimensions: {formatted}."
                    )
            chosen = min(args, key=lambda q: q.value) if func_name == "min" else max(args, key=lambda q: q.value)
            return DimensionedQuantity(chosen.value, chosen.dimensions)
        if func_name == "pow":
            args = [self.visit(arg) for arg in node.args]
            if len(args) != 2:
                raise ValueError("pow() takes exactly two arguments.")
            return args[0] ** args[1]
        raise ValueError(f"Unsupported function: {func_name}")

    def generic_visit(self, node: ast.AST) -> DimensionedQuantity:  # pragma: no cover - guardrail
        raise ValueError(f"Unsupported expression node: {type(node).__name__}")


def evaluate_dimensioned_expression(expression: str | None, *, variables: Mapping[str, Any] | None = None) -> DimensionedQuantity:
    expr = str(expression or "").strip()
    if not expr:
        raise ValueError("Expression is empty.")
    tree = ast.parse(_replace_unit_literals(expr), mode="eval")
    result = _DimensionedEvaluator(variables or {}).visit(tree)
    if not isinstance(result, DimensionedQuantity):
        result = DimensionedQuantity.from_value(result)
    return result


def dimension_name_of(quantity: DimensionedQuantity) -> str | None:
    return _dimension_name(quantity.dimensions)


def has_only_numeric_terms(expression: str | None) -> bool:
    expr = str(expression or "").strip()
    if not expr:
        return True
    processed = _replace_unit_literals(expr)
    tree = ast.parse(processed, mode="eval")
    for node in ast.walk(tree):
        if isinstance(node, (ast.Name, ast.Call)):
            return False
    return True
