"""Numerically audit native outer HeatSolve captures without rerunning physics.

The TES average deliberately mirrors ``TESInnerCircuitUpdate``: every owned TES
element has equal weight, and each of its nodal Temperature DOFs has weight
``1 / element_node_count``.  It is not a nodal or volume average.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def indexed_vector(path: Path, count: int) -> np.ndarray:
    values = np.zeros(count, dtype=np.float64)
    with path.open(encoding="utf-8") as source:
        for line in source:
            fields = line.split()
            if len(fields) < 2:
                continue
            index = int(fields[0])
            if 1 <= index <= count:
                values[index - 1] = float(fields[1])
    return values


def matrix_and_mask(path: Path) -> tuple[object, np.ndarray, int, int]:
    rows: list[int] = []
    columns: list[int] = []
    values: list[float] = []
    max_index = 0
    with path.open(encoding="utf-8") as source:
        for line in source:
            fields = line.split()
            if len(fields) < 3:
                continue
            row, column = int(fields[0]), int(fields[1])
            rows.append(row - 1)
            columns.append(column - 1)
            values.append(float(fields[2]))
            max_index = max(max_index, row, column)
    diagonal = np.zeros(max_index, dtype=bool)
    for row, column, value in zip(rows, columns, values):
        if row == column and value != 0.0:
            diagonal[row] = True
    matrix = coo_matrix((values, (rows, columns)), shape=(max_index, max_index)).tocsc()
    return matrix, diagonal, max_index, len(values)


def residual_stats(residual: np.ndarray, primal: np.ndarray, rhs: np.ndarray) -> dict[str, object]:
    def block(mask: np.ndarray) -> dict[str, float | int]:
        values = residual[mask]
        rhs_values = rhs[mask]
        norm = float(np.linalg.norm(values))
        rhs_norm = float(np.linalg.norm(rhs_values))
        return {
            "rows": int(np.count_nonzero(mask)),
            "l2": norm,
            "max_abs": float(np.max(np.abs(values))) if values.size else 0.0,
            "relative_to_block_rhs_l2": norm / max(rhs_norm, 1.0e-300),
        }

    norm = float(np.linalg.norm(residual))
    return {
        "full": {
            "rows": int(residual.size),
            "l2": norm,
            "max_abs": float(np.max(np.abs(residual))),
            "relative_to_rhs_l2": norm / max(float(np.linalg.norm(rhs)), 1.0e-300),
        },
        "primal": block(primal),
        "constraint": block(~primal),
    }


def load_temp_permutation(result: Path) -> np.ndarray:
    """Read the first Temperature ``Perm`` table: node id -> saved system DOF."""
    with result.open(encoding="utf-8") as source:
        lines = iter(source)
        for line in lines:
            if line.strip().lower() != "temperature":
                continue
            perm_header = next(lines).split()
            if not perm_header or perm_header[0].lower().rstrip(":") != "perm":
                continue
            # Elmer writes ``Perm: node_count dof_count``.  Nonconforming
            # meshes can contain hanging nodes with no temperature DOF, so
            # these two counts need not be equal.
            node_count = int(perm_header[1])
            dof_count = int(perm_header[2]) if len(perm_header) > 2 else node_count
            permutation = np.zeros(node_count, dtype=np.int64)
            for _ in range(dof_count):
                node, dof = next(lines).split()[:2]
                permutation[int(node) - 1] = int(dof)
            return permutation
    raise ValueError(f"Temperature Perm table not found in {result}")


def tes_element_weights(elements: Path, permutation: np.ndarray, body_id: int) -> np.ndarray:
    """Return the exact serial/MPI-reduced TESInnerCircuitUpdate DOF weights."""
    weights = np.zeros(permutation.size, dtype=np.float64)
    element_count = 0
    with elements.open(encoding="utf-8") as source:
        for line in source:
            fields = line.split()
            if len(fields) < 4 or int(fields[1]) != body_id:
                continue
            nodes = [int(value) for value in fields[3:]]
            if not nodes:
                continue
            element_count += 1
            for node in nodes:
                dof = int(permutation[node - 1])
                if dof > 0:
                    weights[dof - 1] += 1.0 / len(nodes)
    if element_count == 0:
        raise ValueError(f"no body {body_id} elements in {elements}")
    return weights / element_count


def capture_paths(directory: Path) -> dict[str, Path]:
    return {
        "before_a": directory / "outer_before_a.dat",
        "before_b": directory / "outer_before_b.dat",
        "before_x": directory / "outer_before_sol.dat",
        "after_a": directory / "outer_after_a.dat",
        "after_b": directory / "outer_after_b.dat",
        "after_x": directory / "outer_after_sol.dat",
    }


def analyze_one(directory: Path, weights: np.ndarray) -> dict[str, object]:
    paths = capture_paths(directory)
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(", ".join(missing))
    before_a, before_primal, count, before_nnz = matrix_and_mask(paths["before_a"])
    after_a, after_primal, after_count, after_nnz = matrix_and_mask(paths["after_a"])
    if count != after_count or weights.size != count:
        raise ValueError("saved-system / Temperature Perm dimensions differ")
    before_b = indexed_vector(paths["before_b"], count)
    after_b = indexed_vector(paths["after_b"], count)
    x_before = indexed_vector(paths["before_x"], count)
    x_after = indexed_vector(paths["after_x"], count)
    x_direct = np.asarray(spsolve(before_a, before_b), dtype=np.float64)
    if not np.all(np.isfinite(x_direct)):
        raise ValueError("direct solution is non-finite")

    native_after_residual = after_a.dot(x_after) - after_b
    native_after_on_before = before_a.dot(x_after) - before_b
    null_before = before_a.dot(x_after - x_direct)
    null_after = after_a.dot(x_after - x_direct)
    averages = {
        "x_before_K": float(weights.dot(x_before)),
        "x_direct_K": float(weights.dot(x_direct)),
        "x_after_K": float(weights.dot(x_after)),
    }
    averages["direct_minus_before_mK"] = (averages["x_direct_K"] - averages["x_before_K"]) * 1.0e3
    averages["after_minus_before_mK"] = (averages["x_after_K"] - averages["x_before_K"]) * 1.0e3
    averages["after_minus_direct_mK"] = (averages["x_after_K"] - averages["x_direct_K"]) * 1.0e3
    return {
        "directory": str(directory),
        "files_sha256": {name: sha256(path) for name, path in paths.items()},
        "systems": {
            "before": {"rows": count, "matrix_records": before_nnz},
            "after": {"rows": after_count, "matrix_records": after_nnz},
            "A_raw_bytes_identical": sha256(paths["before_a"]) == sha256(paths["after_a"]),
            "b_raw_bytes_identical": sha256(paths["before_b"]) == sha256(paths["after_b"]),
            "primal_mask_identical": bool(np.array_equal(before_primal, after_primal)),
        },
        "native_outer_after_residual": residual_stats(native_after_residual, after_primal, after_b),
        "native_outer_after_on_before_system_residual": residual_stats(native_after_on_before, before_primal, before_b),
        "direct_vs_native_null_residual": {
            "before_system": residual_stats(null_before, before_primal, before_b),
            "after_system": residual_stats(null_after, after_primal, after_b),
        },
        "tes_inner_circuit_element_equal_average": averages,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--mesh-elements", type=Path, required=True)
    parser.add_argument("--restart-result", type=Path, required=True)
    parser.add_argument("--tes-body-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=3)
    args = parser.parse_args()
    permutation = load_temp_permutation(args.restart_result)
    weights = tes_element_weights(args.mesh_elements, permutation, args.tes_body_id)
    report = {
        "method": {
            "direct_solver": "scipy.sparse.linalg.spsolve (SuperLU)",
            "tes_average": "TESInnerCircuitUpdate element-equal nodal average",
            "primal_constraint_partition": "nonzero diagonal rows are primal; others constraint",
        },
        "tes_weight_audit": {
            "body_id": args.tes_body_id,
            "system_dofs": int(weights.size),
            "nonzero_weight_dofs": int(np.count_nonzero(weights)),
            "sum": float(weights.sum()),
        },
        "iterations": {},
    }
    for iteration in range(1, args.iterations + 1):
        directory = args.root / f"ts0001_nl{iteration:04d}"
        report["iterations"][f"nl{iteration}"] = analyze_one(directory, weights)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
