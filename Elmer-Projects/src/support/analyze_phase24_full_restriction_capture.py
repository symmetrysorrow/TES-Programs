"""Independently audit the full mortar system captured at MUMPS' call site."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


def indexed_vector(path: Path, count: int) -> np.ndarray:
    result = np.zeros(count, dtype=np.float64)
    with path.open(encoding="utf-8", errors="replace") as source:
        for line in source:
            fields = line.split()
            if len(fields) >= 2:
                result[int(fields[0]) - 1] = float(fields[1])
    return result


def load_matrix(path: Path, shape: tuple[int, int]):
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    with path.open(encoding="utf-8", errors="replace") as source:
        for line in source:
            fields = line.split()
            if len(fields) < 3:
                continue
            rows.append(int(fields[0]) - 1)
            cols.append(int(fields[1]) - 1)
            values.append(float(fields[2]))
    return coo_matrix((values, (rows, cols)), shape=shape).tocsc()


def norm_stats(values: np.ndarray, rhs: np.ndarray | None = None) -> dict[str, float]:
    l2 = float(np.linalg.norm(values))
    rhs_norm = float(np.linalg.norm(rhs)) if rhs is not None else 0.0
    return {
        "absolute_l2": l2,
        "relative_l2": l2 / max(rhs_norm, 1.0e-300),
        "max_abs": float(np.max(np.abs(values))) if values.size else 0.0,
    }


def block_residuals(matrix, rhs: np.ndarray, solution: np.ndarray, primal_rows: int) -> dict[str, object]:
    n = primal_rows
    k = matrix[:n, :n]
    bt = matrix[:n, n:]
    b = matrix[n:, :n]
    d = matrix[n:, n:]
    xp, lam = solution[:n], solution[n:]
    r_primal = np.asarray(k @ xp + bt @ lam - rhs[:n]).ravel()
    r_constraint = np.asarray(b @ xp + d @ lam - rhs[n:]).ravel()
    constraint_rhs_norm = float(np.linalg.norm(rhs[n:]))
    constraint_scale = float(max(
        np.linalg.norm(b.data) * np.linalg.norm(xp),
        np.linalg.norm(d.data) * np.linalg.norm(lam),
        constraint_rhs_norm,
    ))
    constraint = norm_stats(r_constraint, rhs[n:])
    if constraint_rhs_norm <= 1.0e-14:
        constraint["relative_l2"] = None
    constraint["relative_l2_to_constraint_action_scale"] = (
        constraint["absolute_l2"] / max(constraint_scale, 1.0e-300)
    )
    return {
        "primal": norm_stats(r_primal, rhs[:n]),
        "constraint": constraint,
        "constraint_rhs_l2": constraint_rhs_norm,
        "constraint_action_scale": constraint_scale,
        "constraint_matrix_solution_scale": constraint_scale,
    }


def one(directory: Path, weights: np.ndarray) -> tuple[dict[str, object], object, np.ndarray]:
    metadata = json.loads((directory / "metadata.json").read_text(encoding="utf-8"))
    runtime = metadata["runtime"]
    total = int(runtime["total_rows"])
    primal = int(runtime["primal_rows"])
    constraints = total - primal
    before = load_matrix(directory / "full_A_before.dat", (total, total))
    after = load_matrix(directory / "full_A_after.dat", (total, total))
    b_before = indexed_vector(directory / "full_b_before.dat", total)
    b_after = indexed_vector(directory / "full_b_after.dat", total)
    x_before = indexed_vector(directory / "full_x_before.dat", total)
    x_native = indexed_vector(directory / "full_x_after.dat", total)
    lu = splu(before)
    x_direct = np.asarray(lu.solve(b_before), dtype=np.float64)
    direct_residual = before @ x_direct - b_before
    native_residual = before @ x_native - b_before
    native_after_residual = after @ x_native - b_after
    xp_direct = x_direct[:primal]
    xp_native = x_native[:primal]
    k = before[:primal, :primal]
    bt = before[:primal, primal:]
    b = before[primal:, :primal]
    d = before[primal:, primal:]
    k_direct = np.asarray(splu(k).solve(b_before[:primal]), dtype=np.float64)
    delta_mortar = xp_direct - k_direct
    reaction_direct = np.asarray(bt @ x_direct[primal:]).ravel()
    outer_path = directory / "outer_after_sol.dat"
    outer = indexed_vector(outer_path, primal) if outer_path.is_file() else None
    diff = x_native - x_direct
    reaction_abs = np.abs(reaction_direct)
    top = np.argsort(reaction_abs)[-10:][::-1]
    averages: dict[str, float | None] = {
        "x_before_TES": float(weights.dot(x_before[:primal])),
        "full_direct_primal_TES": float(weights.dot(xp_direct)),
        "full_native_primal_TES": float(weights.dot(xp_native)),
        "heat_solve_outer_after_TES": None if outer is None else float(weights.dot(outer)),
        "primal_only_direct_TES": float(weights.dot(k_direct)),
        "mortar_delta_direct_TES": float(weights.dot(delta_mortar)),
    }
    report: dict[str, object] = {
        "directory": str(directory),
        "dimension": {"total_rows": total, "primal_rows": primal, "constraint_rows": constraints},
        "runtime": runtime,
        "sha256": metadata.get("sha256", {}),
        "raw_consistency": {
            "matrix_before_after_identical": metadata.get("sha256", {}).get("matrix_before") == metadata.get("sha256", {}).get("matrix_after"),
            "rhs_before_after_identical": metadata.get("sha256", {}).get("rhs_before") == metadata.get("sha256", {}).get("rhs_after"),
            "matrix_before_after_max_abs": float(np.max(np.abs((before - after).data))) if (before - after).nnz else 0.0,
            "rhs_before_after_max_abs": float(np.max(np.abs(b_before - b_after))),
        },
        "direct_full_residual": norm_stats(direct_residual, b_before),
        "native_full_residual": norm_stats(native_residual, b_before),
        "native_after_system_residual": norm_stats(native_after_residual, b_after),
        "native_vs_direct": {
            "full_vector_max_abs": float(np.max(np.abs(diff))),
            "primal_vector_max_abs": float(np.max(np.abs(diff[:primal]))),
            "multiplier_vector_max_abs": float(np.max(np.abs(diff[primal:]))) if constraints else 0.0,
            "full_l2_relative": float(np.linalg.norm(diff) / max(np.linalg.norm(x_direct), 1.0e-300)),
        },
        "block_residuals": {
            "direct": block_residuals(before, b_before, x_direct, primal),
            "native": block_residuals(before, b_before, x_native, primal),
        },
        "tes_average": averages,
        "heat_solve_outer_after_vs_native": None if outer is None else {
            "max_abs": float(np.max(np.abs(xp_native - outer))),
            "l2_relative": float(np.linalg.norm(xp_native - outer) / max(np.linalg.norm(outer), 1.0e-300)),
            "tes_mK": (averages["full_native_primal_TES"] - averages["heat_solve_outer_after_TES"]) * 1.0e3,
        },
        "mortar_reaction": {
            "direct_Bt_lambda_l2": float(np.linalg.norm(reaction_direct)),
            "direct_Bt_lambda_max_abs": float(np.max(reaction_abs)) if reaction_abs.size else 0.0,
            "direct_Bt_lambda_TES_weighted": float(weights.dot(reaction_direct[:primal])),
            "direct_Bt_lambda_top_dofs_1based": [int(i + 1) for i in top],
            "delta_full_primal_minus_primal_only_l2": float(np.linalg.norm(delta_mortar)),
            "delta_full_primal_minus_primal_only_max_abs": float(np.max(np.abs(delta_mortar))),
            "delta_full_primal_minus_primal_only_TES": float(weights.dot(delta_mortar)),
            "K_nnz": int(k.nnz),
            "Bt_nnz": int(bt.nnz),
            "B_nnz": int(b.nnz),
            "D_nnz": int(d.nnz),
            "D_frobenius": float(np.linalg.norm(d.data)),
        },
    }
    return report, lu, b_before


def load_weights(restart_result: Path, elements: Path, body_id: int) -> np.ndarray:
    permutation = np.zeros(1, dtype=np.int64)
    lines = iter(restart_result.read_text(encoding="utf-8", errors="replace").splitlines())
    for line in lines:
        if line.strip().lower() != "temperature":
            continue
        header = next(lines).split()
        if header[0].lower().rstrip(":") != "perm":
            continue
        permutation = np.zeros(int(header[1]), dtype=np.int64)
        for _ in range(permutation.size):
            node, dof = next(lines).split()[:2]
            permutation[int(node) - 1] = int(dof)
        break
    if permutation.size == 1:
        raise ValueError(f"Temperature Perm table not found in {restart_result}")
    weights = np.zeros(int(permutation.max()), dtype=np.float64)
    count = 0
    for line in elements.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 4 or int(fields[1]) != body_id:
            continue
        nodes = [int(value) for value in fields[3:]]
        count += 1
        for node in nodes:
            dof = int(permutation[node - 1])
            if dof > 0:
                weights[dof - 1] += 1.0 / len(nodes)
    if not count:
        raise ValueError(f"no body {body_id} elements in {elements}")
    return weights / count


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--restart-result", type=Path, required=True)
    parser.add_argument("--mesh-elements", type=Path, required=True)
    parser.add_argument("--tes-body-id", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iterations", type=int, default=3)
    args = parser.parse_args()
    weights = load_weights(args.restart_result, args.mesh_elements, args.tes_body_id)
    reports: dict[str, object] = {}
    lus: list[object] = []
    rhs: list[np.ndarray] = []
    for iteration in range(1, args.iterations + 1):
        report, factor, vector = one(args.root / f"ts0001_nl{iteration:04d}", weights)
        reports[f"nl{iteration}"] = report
        lus.append(factor)
        rhs.append(vector)
    crossed: dict[str, object] = {}
    if len(lus) >= 2:
        base = reports["nl1"]
        base_tes = float(base["tes_average"]["full_direct_primal_TES"])
        for index in range(1, len(lus)):
            current = reports[f"nl{index + 1}"]
            n = int(current["dimension"]["primal_rows"])
            w = weights[:n]
            a11 = base_tes
            a12 = float(w.dot(lus[index - 1].solve(rhs[0])[:n]))
            a21 = float(w.dot(lus[index].solve(rhs[0])[:n]))
            a22 = float(current["tes_average"]["full_direct_primal_TES"])
            crossed[f"nl1_to_nl{index + 1}"] = {
                "rhs_only_mK": (a12 - a11) * 1.0e3,
                "operator_only_mK": (a21 - a11) * 1.0e3,
                "interaction_mK": (a22 - a12 - a21 + a11) * 1.0e3,
                "total_mK": (a22 - a11) * 1.0e3,
            }
    result = {
        "method": {
            "direct_solver": "scipy.sparse.linalg.splu (SuperLU)",
            "tes_average": "TESInnerCircuitUpdate element-equal nodal average",
            "block_partition": "runtime primal_rows; tail rows are explicit restriction equations",
        },
        "iterations": reports,
        "crossed_sensitivity": crossed,
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
