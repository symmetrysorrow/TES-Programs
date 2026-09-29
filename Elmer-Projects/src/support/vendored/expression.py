"""Plain (dimensionless) expression evaluation.

Vendored 2026-07-14 from Thermal-and-Electoric-Sim (core/circuit/expression.py),
trimmed to `evaluate_expression` and its helpers (autodiff/explain utilities
removed). Only length units are accepted as [unit] literals here, matching
upstream behaviour.
"""
from __future__ import annotations

import ast
import math
import re
from typing import Any, Mapping

_ALLOWED_FUNCTIONS = {
    "abs": abs,
    "min": min,
    "max": max,
    "pow": pow,
    "sqrt": math.sqrt,
    "sin": math.sin,
    "cos": math.cos,
    "tan": math.tan,
    "asin": math.asin,
    "acos": math.acos,
    "atan": math.atan,
    "exp": math.exp,
    "log": math.log,
    "log10": math.log10,
    "tanh": math.tanh,
    "sinh": math.sinh,
    "cosh": math.cosh,
    "floor": math.floor,
    "ceil": math.ceil,
    "pi": math.pi,
    "e": math.e,
}

_LENGTH_UNIT_FACTORS = {
    "km": 1000.0,
    "m": 1.0,
    "mm": 1.0e-3,
    "um": 1.0e-6,
    "nm": 1.0e-9,
    "pm": 1.0e-12,
}

_UNIT_LITERAL_PATTERN = re.compile(
    r"(?<![\w.])"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"\s*\[\s*([a-zA-Zµμ]+)\s*\]"
    r"(?![\w])"
)

_BARE_UNIT_LITERAL_PATTERN = re.compile(
    r"(?<![\w.])"
    r"([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)"
    r"(?![eE][+-]?\d)"
    r"\s*([a-zA-Zµμ]+)"
    r"(?![\w])"
)


def _normalize_unit_literals(expression: str | None) -> str:
    raw = str(expression or "")
    if not raw.strip():
        return raw

    bare_unit_match = _BARE_UNIT_LITERAL_PATTERN.search(raw)
    if bare_unit_match:
        unit = bare_unit_match.group(2).strip().replace("µ", "u").replace("μ", "u").lower()
        raise ValueError(f"Units must use [unit] notation; found bare unit '{unit}'.")

    def repl(match: re.Match[str]) -> str:
        value = float(match.group(1))
        unit = (match.group(2) or "").strip().replace("µ", "u").replace("μ", "u").lower()
        factor = _LENGTH_UNIT_FACTORS.get(unit)
        if factor is None:
            raise ValueError(f"Unknown unit: {unit}")
        return format(value * factor, ".15g")

    return _UNIT_LITERAL_PATTERN.sub(repl, raw)


def _resolve_attribute_value(value: Any, attr: str) -> Any:
    if isinstance(value, Mapping) and attr in value:
        return value[attr]
    if hasattr(value, attr):
        return getattr(value, attr)
    raise ValueError(f"Unknown attribute: {attr}")


class _ExpressionEvaluator(ast.NodeVisitor):
    def __init__(self, variables: Mapping[str, Any]) -> None:
        self.variables = dict(variables)

    def visit_Expression(self, node: ast.Expression) -> Any:  # pragma: no cover - ast entrypoint
        return self.visit(node.body)

    def visit_Constant(self, node: ast.Constant) -> Any:
        if isinstance(node.value, (int, float, bool)):
            return node.value
        raise ValueError(f"Unsupported constant type: {type(node.value).__name__}")

    def visit_Name(self, node: ast.Name) -> Any:
        if node.id in self.variables:
            return self.variables[node.id]
        if node.id in _ALLOWED_FUNCTIONS:
            return _ALLOWED_FUNCTIONS[node.id]
        raise ValueError(f"Unknown symbol: {node.id}")

    def visit_Attribute(self, node: ast.Attribute) -> Any:
        value = self.visit(node.value)
        return _resolve_attribute_value(value, node.attr)

    def visit_BinOp(self, node: ast.BinOp) -> Any:
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
        if isinstance(node.op, ast.FloorDiv):
            return left // right
        if isinstance(node.op, ast.Mod):
            return left % right
        if isinstance(node.op, ast.Pow):
            return left**right
        raise ValueError(f"Unsupported binary operator: {type(node.op).__name__}")

    def visit_UnaryOp(self, node: ast.UnaryOp) -> Any:
        operand = self.visit(node.operand)
        if isinstance(node.op, ast.UAdd):
            return +operand
        if isinstance(node.op, ast.USub):
            return -operand
        if isinstance(node.op, ast.Not):
            return not operand
        raise ValueError(f"Unsupported unary operator: {type(node.op).__name__}")

    def visit_BoolOp(self, node: ast.BoolOp) -> Any:
        values = [self.visit(v) for v in node.values]
        if isinstance(node.op, ast.And):
            result = values[0]
            for value in values[1:]:
                result = result and value
            return result
        if isinstance(node.op, ast.Or):
            result = values[0]
            for value in values[1:]:
                result = result or value
            return result
        raise ValueError(f"Unsupported boolean operator: {type(node.op).__name__}")

    def visit_Compare(self, node: ast.Compare) -> Any:
        left = self.visit(node.left)
        for op, comparator in zip(node.ops, node.comparators, strict=True):
            right = self.visit(comparator)
            if isinstance(op, ast.Eq):
                ok = left == right
            elif isinstance(op, ast.NotEq):
                ok = left != right
            elif isinstance(op, ast.Lt):
                ok = left < right
            elif isinstance(op, ast.LtE):
                ok = left <= right
            elif isinstance(op, ast.Gt):
                ok = left > right
            elif isinstance(op, ast.GtE):
                ok = left >= right
            else:
                raise ValueError(f"Unsupported comparison operator: {type(op).__name__}")
            if not ok:
                return False
            left = right
        return True

    def visit_IfExp(self, node: ast.IfExp) -> Any:
        return self.visit(node.body) if self.visit(node.test) else self.visit(node.orelse)

    def visit_Call(self, node: ast.Call) -> Any:
        if not isinstance(node.func, ast.Name):
            raise ValueError("Only direct function calls are supported.")
        func = self.visit_Name(node.func)
        if func not in _ALLOWED_FUNCTIONS.values():
            raise ValueError(f"Unsupported function: {node.func.id}")
        args = [self.visit(arg) for arg in node.args]
        kwargs = {
            kw.arg: self.visit(kw.value)
            for kw in node.keywords
            if kw.arg is not None
        }
        return func(*args, **kwargs)

    def generic_visit(self, node: ast.AST) -> Any:  # pragma: no cover - guardrail
        raise ValueError(f"Unsupported expression node: {type(node).__name__}")


def evaluate_expression(expression: str | None, *, variables: Mapping[str, Any]) -> float:
    if expression is None or not str(expression).strip():
        raise ValueError("Expression is empty.")
    tree = ast.parse(_normalize_unit_literals(expression), mode="eval")
    result = _ExpressionEvaluator(variables).visit(tree)
    if isinstance(result, bool):
        return float(result)
    return float(result)
