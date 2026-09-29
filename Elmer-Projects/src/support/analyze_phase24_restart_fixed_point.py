"""Compare the Gate3 steady field with the first transient saddle system."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import lsqr


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def indexed_vector(path: Path, count: int) -> np.ndarray:
    result = np.zeros(count, dtype=np.float64)
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                result[int(fields[0]) - 1] = float(fields[1])
    return result


def load_matrix(path: Path, shape: tuple[int, int]) -> csr_matrix:
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 3:
                rows.append(int(fields[0]) - 1)
                cols.append(int(fields[1]) - 1)
                values.append(float(fields[2]))
    matrix = coo_matrix((values, (rows, cols)), shape=shape).tocsr()
    matrix.eliminate_zeros()
    return matrix


def matrix_shape(path: Path) -> tuple[int, int, int]:
    records = 0
    extent = 0
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 3:
                continue
            records += 1
            extent = max(extent, int(fields[0]), int(fields[1]))
    return extent, extent, records


def norm_stats(values: np.ndarray, reference: np.ndarray | None = None) -> dict[str, float]:
    l2 = float(np.linalg.norm(values))
    ref = float(np.linalg.norm(reference)) if reference is not None else 0.0
    return {
        "l2": l2,
        "relative_l2": l2 / max(ref, 1.0e-300),
        "max_abs": float(np.max(np.abs(values))) if values.size else 0.0,
    }


def matrix_diff(left: csr_matrix, right: csr_matrix) -> dict[str, float | int | str]:
    delta = (left - right).tocsr()
    data = delta.data
    return {
        "shape": [int(left.shape[0]), int(left.shape[1])],
        "nnz_left": int(left.nnz),
        "nnz_right": int(right.nnz),
        "max_abs_difference": float(np.max(np.abs(data))) if data.size else 0.0,
        "l2_difference": float(np.linalg.norm(data)),
        "relative_l2_difference": float(np.linalg.norm(data) / max(np.linalg.norm(right.data), 1.0e-300)),
        "differing_records": int(np.count_nonzero(data)),
    }


def vector_diff(left: np.ndarray, right: np.ndarray) -> dict[str, float | int]:
    delta = left - right
    return {
        "rows": int(delta.size),
        "max_abs_difference": float(np.max(np.abs(delta))) if delta.size else 0.0,
        "l2_difference": float(np.linalg.norm(delta)),
        "relative_l2_difference": float(np.linalg.norm(delta) / max(np.linalg.norm(right), 1.0e-300)),
    }


def load_result_field(result: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = result.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
    header = lines[marker + 1].split()
    count = int(header[1])
    permutation = np.zeros(count, dtype=np.int64)
    for offset in range(count):
        node, dof = lines[marker + 2 + offset].split()[:2]
        permutation[int(node) - 1] = int(dof)
    values = np.asarray(
        [float(line.split()[0]) for line in lines[marker + 2 + count: marker + 2 + 2 * count]],
        dtype=np.float64,
    )
    if values.size != count:
        raise ValueError(f"expected {count} result values, got {values.size}")
    field = np.zeros(count, dtype=np.float64)
    field[permutation - 1] = values
    return field, permutation


def tes_weights(elements: Path, permutation: np.ndarray, body_id: int) -> np.ndarray:
    weights = np.zeros(permutation.size, dtype=np.float64)
    elements_seen = 0
    with elements.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 4 or int(fields[1]) != body_id:
                continue
            nodes = [int(value) for value in fields[3:]]
            elements_seen += 1
            for node in nodes:
                dof = int(permutation[node - 1])
                if dof > 0:
                    weights[dof - 1] += 1.0 / len(nodes)
    if not elements_seen:
        raise ValueError(f"no body {body_id} elements in {elements}")
    return weights / elements_seen


def node_coordinates(nodes_path: Path, node_ids: set[int]) -> dict[int, list[float]]:
    found: dict[int, list[float]] = {}
    if not node_ids:
        return found
    with nodes_path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 5 and int(fields[0]) in node_ids:
                found[int(fields[0])] = [float(fields[2]), float(fields[3]), float(fields[4])]
                if len(found) == len(node_ids):
                    break
    return found


def interface_rows(residual: np.ndarray, b: csr_matrix, permutation: np.ndarray, nodes: Path) -> list[dict[str, object]]:
    order = np.argsort(np.abs(residual))[::-1][:10]
    dof_to_node = np.zeros(permutation.size, dtype=np.int64)
    dof_to_node[permutation - 1] = np.arange(1, permutation.size + 1)
    chosen: list[tuple[int, np.ndarray]] = []
    node_ids: set[int] = set()
    for row in order:
        cols = b.getrow(int(row)).indices
        chosen.append((int(row) + 1, cols))
        node_ids.update(int(dof_to_node[col]) for col in cols if dof_to_node[col] > 0)
    coords = node_coordinates(nodes, node_ids)
    return [
        {
            "constraint_row_1based": row,
            "residual": float(residual[row - 1]),
            "columns_1based": [int(col + 1) for col in cols[:20]],
            "interface_nodes": [int(dof_to_node[col]) for col in cols[:20] if dof_to_node[col] > 0],
            "coordinates_m": [coords.get(int(dof_to_node[col]), []) for col in cols[:20] if dof_to_node[col] > 0],
        }
        for row, cols in chosen
    ]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--transient", type=Path, required=True)
    parser.add_argument("--steady", type=Path, required=True)
    parser.add_argument("--gate3-result", type=Path, required=True)
    parser.add_argument("--mesh-elements", type=Path, required=True)
    parser.add_argument("--mesh-nodes", type=Path, required=True)
    parser.add_argument("--tes-body-id", type=int, default=101)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    transient_n, _, transient_records = matrix_shape(args.transient / "full_A_before.dat")
    steady_a_path = args.steady / "full_A.dat"
    steady_b_path = args.steady / "full_b.dat"
    steady_n, _, steady_records = matrix_shape(steady_a_path)
    transient_meta = json.loads((args.transient / "metadata.json").read_text(encoding="utf-8"))
    transient_p = int(transient_meta["runtime"]["primal_rows"])
    steady_p = int(transient_meta["runtime"]["primal_rows"])
    at = load_matrix(args.transient / "full_A_before.dat", (transient_n, transient_n))
    bt_rhs = indexed_vector(args.transient / "full_b_before.dat", transient_n)
    x_restart = indexed_vector(args.transient / "full_x_before.dat", transient_n)[:transient_p]
    steady_a = load_matrix(steady_a_path, (steady_n, steady_n))
    steady_rhs = indexed_vector(steady_b_path, steady_n)
    x_steady, permutation = load_result_field(args.gate3_result)
    weights = tes_weights(args.mesh_elements, permutation, args.tes_body_id)

    k_t = at[:transient_p, :transient_p].tocsr()
    bt = at[:transient_p, transient_p:].tocsr()
    b = at[transient_p:, :transient_p].tocsr()
    g = bt_rhs[transient_p:]
    k_s = steady_a[:steady_p, :steady_p].tocsr()
    bt_s = steady_a[:steady_p, steady_p:].tocsr()
    b_s = steady_a[steady_p:, :steady_p].tocsr()
    g_s = steady_rhs[steady_p:]
    r_constraint = np.asarray(b @ x_restart - g).ravel()
    r0 = np.asarray(k_t @ x_restart - bt_rhs[:transient_p]).ravel()
    lsqr_result = lsqr(bt, -r0, atol=1.0e-14, btol=1.0e-14, iter_lim=2000)
    lam = np.asarray(lsqr_result[0], dtype=np.float64)
    r_best = r0 + np.asarray(bt @ lam).ravel()
    residual_mask = weights != 0.0

    xs = x_steady[:steady_p]
    steady_constraint_residual = np.asarray(b_s @ xs - g_s).ravel()
    steady_r0 = np.asarray(k_s @ xs - steady_rhs[:steady_p]).ravel()
    steady_lsqr = lsqr(bt_s, -steady_r0, atol=1.0e-14, btol=1.0e-14, iter_lim=2000)
    steady_lam = np.asarray(steady_lsqr[0], dtype=np.float64)
    steady_best = steady_r0 + np.asarray(bt_s @ steady_lam).ravel()
    steady_field_residual = np.asarray(steady_a @ np.r_[xs, steady_lam] - steady_rhs).ravel()
    fixed_constraint = np.asarray(b @ xs - g).ravel()
    fixed_r0 = np.asarray(k_t @ xs - bt_rhs[:transient_p]).ravel()
    fixed_lsqr = lsqr(bt, -fixed_r0, atol=1.0e-14, btol=1.0e-14, iter_lim=2000)
    fixed_lam = np.asarray(fixed_lsqr[0], dtype=np.float64)
    fixed_best = fixed_r0 + np.asarray(bt @ fixed_lam).ravel()
    delta_k = (k_t - k_s).tocsr()
    delta_rhs = bt_rhs[:transient_p] - steady_rhs[:steady_p]
    history_proxy = np.asarray(delta_k @ xs).ravel()
    source_boundary_unexplained = delta_rhs - history_proxy

    report: dict[str, object] = {
        "method": {
            "transient_system": "SolveWithLinearRestriction full_A_before at timestep 1 nonlinear iteration 1",
            "steady_system": "SolveWithLinearRestriction full_A_after at Gate3 steady final nonlinear iteration 4",
            "steady_capture_caveat": "The steady capture selector uses the solver's internal timestep=1 even though the simulation is steady-state.",
            "least_squares": "scipy.sparse.linalg.lsqr on Bt with primal field fixed",
        },
        "dimensions": {
            "transient": {"total_rows": transient_n, "primal_rows": transient_p, "constraint_rows": transient_n - transient_p, "matrix_records": transient_records},
            "steady": {"total_rows": steady_n, "primal_rows": steady_p, "constraint_rows": steady_n - steady_p, "matrix_records": steady_records},
        },
        "transient_constraint_residual_restart": {
            "row_count": int(r_constraint.size),
            "g_norm": float(np.linalg.norm(g)),
            "Bx_restart_norm": float(np.linalg.norm(np.asarray(b @ x_restart).ravel())),
            "stats": norm_stats(r_constraint, g),
            "relative_to_Bx_scale": float(np.linalg.norm(r_constraint) / max(np.linalg.norm(np.asarray(b @ x_restart).ravel()), 1.0e-300)),
            "worst_rows": interface_rows(r_constraint, b, permutation, args.mesh_nodes),
        },
        "best_fit_multiplier_restart": {
            "lambda_norm": float(np.linalg.norm(lam)),
            "lsqr_istop": int(lsqr_result[1]),
            "lsqr_iterations": int(lsqr_result[2]),
            "primal_residual": {
                **norm_stats(r_best, bt_rhs[:transient_p]),
                "tes_region": norm_stats(r_best[residual_mask], bt_rhs[:transient_p][residual_mask]),
                "tes_weighted_signed": float(weights.dot(r_best)),
            },
            "Bt_lambda": norm_stats(np.asarray(bt @ lam).ravel()),
            "Kx_minus_b": norm_stats(r0, bt_rhs[:transient_p]),
        },
        "block_comparison_steady_vs_transient": {
            "K": {
                **matrix_diff(k_t, k_s),
                "transient_raw_sha256": sha256(args.transient / "full_A_before.dat"),
                "steady_raw_sha256": sha256(steady_a_path),
            },
            "B": {
                "transient_shape": [int(b.shape[0]), int(b.shape[1])],
                "transient_nnz": int(b.nnz),
                "steady_shape": [int(b_s.shape[0]), int(b_s.shape[1])],
                "steady_nnz": int(b_s.nnz),
                "steady_has_constraint_block": bool(b_s.shape[0] > 0),
                "changed_rows": int((b - b_s).getnnz()),
                "changed_interface_dofs": int(np.count_nonzero(np.asarray((b - b_s).getnnz(axis=0)).ravel())),
                "max_abs_transient": float(np.max(np.abs(b.data))) if b.nnz else 0.0,
                "l2_transient": float(np.linalg.norm(b.data)),
                "relative_l2_transient_to_K": float(np.linalg.norm(b.data) / max(np.linalg.norm(k_t.data), 1.0e-300)),
                "max_abs_difference": float(np.max(np.abs((b - b_s).data))) if (b - b_s).nnz else 0.0,
                "l2_difference": float(np.linalg.norm((b - b_s).data)),
            },
            "Bt": {
                "transient_shape": [int(bt.shape[0]), int(bt.shape[1])],
                "transient_nnz": int(bt.nnz),
                "transpose_consistent": bool((bt - b.transpose()).nnz == 0),
                "steady_shape": [int(bt_s.shape[0]), int(bt_s.shape[1])],
                "steady_nnz": int(bt_s.nnz),
                "max_abs_difference": float(np.max(np.abs((bt - bt_s).data))) if (bt - bt_s).nnz else 0.0,
                "l2_difference": float(np.linalg.norm((bt - bt_s).data)),
            },
            "g": {
                "transient_norm": float(np.linalg.norm(g)),
                "steady_norm": float(np.linalg.norm(g_s)),
                "steady_shape": [int(g_s.size)],
                "steady_nonzero_rows": int(np.count_nonzero(g_s)),
                "transient_nonzero_rows": int(np.count_nonzero(g)),
                "difference_l2": float(np.linalg.norm(g - g_s)),
            },
            "primal_rhs": vector_diff(bt_rhs[:transient_p], steady_rhs[:steady_p]),
            "steady_field_residual": norm_stats(steady_field_residual, steady_rhs),
            "steady_best_fit_primal": {
                "lambda_norm": float(np.linalg.norm(steady_lam)),
                "residual": norm_stats(steady_best, steady_rhs[:steady_p]),
            },
            "transient_minus_steady_expected_bdf1": {
                "definition": "(b_t-b_s) - (K_t-K_s)*x_s; proxy for non-mass/source/history mismatch",
                "stats": norm_stats(source_boundary_unexplained, bt_rhs[:transient_p]),
                "tes_weighted_signed": float(weights.dot(source_boundary_unexplained)),
            },
        },
        "fixed_point_gate3_field_in_transient": {
            "steady_field_tes_average_K": float(weights.dot(xs)),
            "restart_field_tes_average_K": float(weights.dot(x_restart)),
            "restart_minus_steady_field": norm_stats(x_restart - xs, xs),
            "steady_field_residual": norm_stats(steady_field_residual, steady_rhs),
            "constraint": {
                "stats": norm_stats(fixed_constraint, g),
                "g_norm": float(np.linalg.norm(g)),
                "Bx_s_norm": float(np.linalg.norm(np.asarray(b @ xs).ravel())),
                "worst_rows": interface_rows(fixed_constraint, b, permutation, args.mesh_nodes),
            },
            "steady_constraint": {
                "stats": norm_stats(steady_constraint_residual, g_s),
                "g_norm": float(np.linalg.norm(g_s)),
                "Bx_s_norm": float(np.linalg.norm(np.asarray(b_s @ xs).ravel())),
            },
            "best_fit_primal": {
                "lambda_norm": float(np.linalg.norm(fixed_lam)),
                "lsqr_istop": int(fixed_lsqr[1]),
                "lsqr_iterations": int(fixed_lsqr[2]),
                "residual": {
                    **norm_stats(fixed_best, bt_rhs[:transient_p]),
                    "tes_region": norm_stats(fixed_best[residual_mask], bt_rhs[:transient_p][residual_mask]),
                    "tes_weighted_signed": float(weights.dot(fixed_best)),
                },
            },
            "expected_rhs_difference": {
                "definition": "actual transient primal RHS minus (K_t*x_s + Bt*lambda_best)",
                "stats": norm_stats(-fixed_best, bt_rhs[:transient_p]),
                "tes_weighted_signed": float(-weights.dot(fixed_best)),
            },
        },
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
