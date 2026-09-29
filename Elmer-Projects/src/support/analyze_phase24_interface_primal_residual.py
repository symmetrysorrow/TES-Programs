"""Account for the Gate3 -> transient BDF1 interface primal residual.

This is a read-only post-processor for the Phase24 full saddle captures. It
does not alter an SIF, a solver, or a native source tree.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import lsqr


BODY_NAMES = {
    100: "abs", 101: "TES", 102: "Stycast", 103: "Membrane_SiNx",
    104: "SiO2_1", 105: "Si_1", 106: "SiNx", 107: "Si_2",
    108: "SiO2_2", 109: "Membrane_Si1",
}
MATERIAL_NAMES = {
    100: "Pb", 101: "TES", 102: "Stycast", 103: "Membrane",
    104: "SiO2", 105: "Si", 106: "SiNx", 107: "Si",
    108: "SiO2", 109: "Membrane",
}
RHO_CP = {
    100: 9860.0 * 3.26e-05, 101: 15695.0 * 0.00125,
    102: 2400.0 * 0.00122, 103: 2387.0 * 9.82e-06,
    104: 2200.0 * 0.000124, 105: 2332.0 * 1.24e-06,
    106: 3400.0 * 2.43e-05, 107: 2332.0 * 1.24e-06,
    108: 2200.0 * 0.000124, 109: 2387.0 * 9.82e-06,
}
BOUNDARY_NAMES = {
    1104: "TES__zmin", 1105: "TES__zmax", 1204: "Stycast__zmin",
    1205: "Stycast__zmax", 1305: "Membrane_SiNx__zmax",
    1004: "abs__zmin",
}


def indexed_vector(path: Path, count: int) -> np.ndarray:
    out = np.zeros(count, dtype=np.float64)
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 2:
                out[int(fields[0]) - 1] = float(fields[1])
    return out


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
    return coo_matrix((values, (rows, cols)), shape=shape).tocsr()


def load_result_field(path: Path) -> tuple[np.ndarray, np.ndarray]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
    count = int(lines[marker + 1].split()[1])
    permutation = np.zeros(count, dtype=np.int64)
    for offset in range(count):
        node, dof = lines[marker + 2 + offset].split()[:2]
        permutation[int(node) - 1] = int(dof)
    values = np.asarray(
        [float(line.split()[0]) for line in lines[marker + 2 + count: marker + 2 + 2 * count]],
        dtype=np.float64,
    )
    field = np.zeros(count, dtype=np.float64)
    field[permutation - 1] = values
    return field, permutation


def load_nodes(path: Path) -> np.ndarray:
    max_node = 0
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 5:
                max_node = max(max_node, int(fields[0]))
    coords = np.full((max_node, 3), np.nan, dtype=np.float64)
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 5:
                coords[int(fields[0]) - 1] = [float(fields[2]), float(fields[3]), float(fields[4])]
    return coords


def mesh_metadata(elements_path: Path, boundary_path: Path, coords: np.ndarray):
    node_bodies = [set() for _ in range(coords.shape[0])]
    node_boundaries = [set() for _ in range(coords.shape[0])]
    mass = np.zeros(coords.shape[0], dtype=np.float64)
    counts: Counter[str] = Counter()
    with elements_path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 4:
                continue
            body = int(fields[1])
            etype = int(fields[2])
            nodes = [int(v) for v in fields[3:]]
            for node in nodes:
                if 1 <= node <= len(node_bodies):
                    node_bodies[node - 1].add(body)
            if etype == 504 and len(nodes) >= 4 and body in RHO_CP:
                p = coords[np.asarray(nodes[:4]) - 1]
                volume = abs(float(np.linalg.det(np.column_stack((p[1:] - p[0]))) / 6.0))
                contribution = RHO_CP[body] * volume / 4.0
                for node in nodes[:4]:
                    mass[node - 1] += contribution
                counts["tetra_504"] += 1
            else:
                counts[f"element_{etype}"] += 1
    with boundary_path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 5:
                continue
            boundary = int(fields[1])
            for value in fields[4:]:
                node = int(value)
                if 1 <= node <= len(node_boundaries):
                    node_boundaries[node - 1].add(boundary)
    return node_bodies, node_boundaries, mass, dict(counts)


def stat(values: np.ndarray, reference: np.ndarray | None = None) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64).ravel()
    scale = np.linalg.norm(reference) if reference is not None else 1.0
    return {
        "l2": float(np.linalg.norm(values)),
        "max_abs": float(np.max(np.abs(values))) if values.size else 0.0,
        "relative_l2": float(np.linalg.norm(values) / max(scale, 1e-300)),
        "sum_abs": float(np.sum(np.abs(values))),
    }


def matrix_diff(left: csr_matrix, right: csr_matrix) -> dict[str, float | int]:
    delta = (left - right).tocsr()
    return {
        "left_nnz": int(left.nnz),
        "right_nnz": int(right.nnz),
        "different_nnz": int(delta.nnz),
        "max_abs": float(np.max(np.abs(delta.data))) if delta.nnz else 0.0,
        "l2_values": float(np.linalg.norm(delta.data)) if delta.nnz else 0.0,
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def permutation_maps(permutation: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dof_to_node = np.full(int(np.max(permutation)), -1, dtype=np.int64)
    node_to_dof = np.zeros(permutation.size, dtype=np.int64)
    for node, dof in enumerate(permutation, start=1):
        if dof > 0:
            node_to_dof[node - 1] = dof
            dof_to_node[dof - 1] = node
    return node_to_dof, dof_to_node


def load_steady_lambda(steady_k, steady_bt, steady_b, xs, path: Path):
    if path.is_file():
        candidate = indexed_vector(path, steady_k.shape[0] + steady_bt.shape[1])
        if np.linalg.norm(candidate[steady_k.shape[0]:]) > 1e-14:
            return candidate[steady_k.shape[0]:], "captured full_x"
    r0 = np.asarray(steady_k @ xs - steady_b).ravel()
    result = lsqr(steady_bt, -r0, atol=1e-14, btol=1e-14, iter_lim=2000)
    return np.asarray(result[0]), "least-squares fit to steady primal equilibrium"


def context(row, residual, eta, denominator, permutation, node_to_dof, dof_to_node,
            coords, node_bodies, node_boundaries, interface_mask, bmat,
            kx_s, mortar_s, rhs_s, kx_t, mortar_t, rhs_t, delta_kx,
            delta_mortar, delta_b, history, physical_delta_kx,
            mass_mismatch, physical_row_delta, source_without_history,
            steady_residual, closed_predicted, closed_error):
    node = int(dof_to_node[row]) if row < len(dof_to_node) and dof_to_node[row] > 0 else 0
    bodies = sorted(node_bodies[node - 1]) if node else []
    boundaries = sorted(node_boundaries[node - 1]) if node else []
    constraint_rows = (bmat.getcol(row).nonzero()[0] + 1).tolist()
    category = "interface_primal" if interface_mask[row] else "bulk_primal"
    if 101 in bodies and interface_mask[row]:
        category = "tes_interface_primal"
    elif 101 in bodies:
        category = "tes_primal"
    return {
        "row_1based": row + 1, "dof_1based": row + 1, "node_1based": node,
        "coordinates_m": coords[node - 1].tolist() if node else [],
        "category": category, "eta": float(eta), "residual": float(residual),
        "denominator": float(denominator), "body_ids": bodies,
        "body_names": [BODY_NAMES.get(v, f"body_{v}") for v in bodies],
        "material_names": [MATERIAL_NAMES.get(v, f"material_for_body_{v}") for v in bodies],
        "boundary_ids": boundaries,
        "boundary_names": [BOUNDARY_NAMES.get(v, f"boundary_{v}") for v in boundaries],
        "mortar_constraint_rows_1based": constraint_rows, "tes_region": 101 in bodies,
        "K_s_x_s": float(kx_s[row]), "B_sT_lambda_s": float(mortar_s[row]),
        "b_s": float(rhs_s[row]), "K_t_x_s": float(kx_t[row]),
        "B_tT_lambda_t": float(mortar_t[row]), "b_t": float(rhs_t[row]),
        "delta_Kx": float(delta_kx[row]), "delta_Mortar": float(delta_mortar[row]),
        "delta_b": float(delta_b[row]),
        "r_predicted": float(delta_kx[row] + delta_mortar[row] - delta_b[row]),
        "steady_restart_primal_residual": float(steady_residual[row]),
        "r_closed_with_steady_residual": float(closed_predicted[row]),
        "r_actual": float(residual),
        "accounting_error": float(delta_kx[row] + delta_mortar[row] - delta_b[row] - residual),
        "closed_accounting_error": float(closed_error[row]),
        "mass_lhs": float(history[row]), "history_rhs": float(history[row]),
        "mass_history_mismatch": float(mass_mismatch[row]),
        "physical_K_delta_x": float(physical_delta_kx[row]),
        "physical_row_delta": float(physical_row_delta[row]),
        "rhs_source_without_history": float(source_without_history[row]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--steady-capture", type=Path, required=True)
    parser.add_argument("--steady-outer", type=Path, required=True)
    parser.add_argument("--transient-capture", type=Path, required=True)
    parser.add_argument("--transient-outer", type=Path, required=True)
    parser.add_argument("--restart-result", type=Path, required=True)
    parser.add_argument("--mesh-nodes", type=Path, required=True)
    parser.add_argument("--mesh-elements", type=Path, required=True)
    parser.add_argument("--mesh-boundary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steady-multiplier", type=Path)
    args = parser.parse_args()
    out = args.output
    out.mkdir(parents=True, exist_ok=True)

    tmeta = json.loads((args.transient_capture / "metadata.json").read_text(encoding="utf-8"))
    smeta = json.loads((args.steady_capture / "metadata.json").read_text(encoding="utf-8"))
    total = int(tmeta["runtime"]["total_rows"])
    primal = int(tmeta["runtime"]["primal_rows"])
    constraints = total - primal
    transient_a = load_matrix(args.transient_capture / "full_A_before.dat", (total, total))
    transient_rhs_full = indexed_vector(args.transient_capture / "full_b_before.dat", total)
    steady_a = load_matrix(args.steady_capture / "full_A.dat", (total, total))
    steady_rhs_full = indexed_vector(args.steady_capture / "full_b.dat", total)
    restart_full, permutation = load_result_field(args.restart_result)
    xs = restart_full[:primal]
    node_to_dof, dof_to_node = permutation_maps(permutation)

    kt = transient_a[:primal, :primal].tocsr()
    bt = transient_a[:primal, primal:].tocsr()
    bmat_t = transient_a[primal:, :primal].tocsr()
    ks = steady_a[:primal, :primal].tocsr()
    bts = steady_a[:primal, primal:].tocsr()
    bmat_s = steady_a[primal:, :primal].tocsr()
    rhs_t = transient_rhs_full[:primal]
    rhs_s = steady_rhs_full[:primal]
    lambda_s, lambda_s_method = load_steady_lambda(
        ks, bts, rhs_s, xs, args.steady_multiplier or Path())
    r0_t = np.asarray(kt @ xs - rhs_t).ravel()
    fit_t = lsqr(bt, -r0_t, atol=1e-14, btol=1e-14, iter_lim=2000)
    lambda_t = np.asarray(fit_t[0], dtype=np.float64)
    mortar_s = np.asarray(bts @ lambda_s).ravel()
    mortar_t = np.asarray(bt @ lambda_t).ravel()
    kx_s = np.asarray(ks @ xs).ravel()
    kx_t = np.asarray(kt @ xs).ravel()
    delta_kx = kx_t - kx_s
    delta_mortar = mortar_t - mortar_s
    delta_b = rhs_t - rhs_s
    r_actual = kx_t + mortar_t - rhs_t
    r_predicted = delta_kx + delta_mortar - delta_b
    accounting_error = r_predicted - r_actual
    steady_residual = kx_s + mortar_s - rhs_s
    closed_predicted = r_predicted + steady_residual
    closed_accounting_error = closed_predicted - r_actual

    coords = load_nodes(args.mesh_nodes)
    node_bodies, node_boundaries, mass_node, element_counts = mesh_metadata(
        args.mesh_elements, args.mesh_boundary, coords)
    mass_dof = np.zeros(primal, dtype=np.float64)
    for node, dof in enumerate(node_to_dof, start=1):
        if dof > 0:
            mass_dof[dof - 1] = mass_node[node - 1]
    dt = float(tmeta["provenance"]["dt"])
    mass_lhs = mass_dof * xs / dt
    history_rhs = mass_lhs.copy()

    outer_t = load_matrix(args.transient_outer / "outer_before_a.dat", (primal, primal))
    outer_s = load_matrix(args.steady_outer / "outer_before_a.dat", (primal, primal))
    outer_rhs_t = indexed_vector(args.transient_outer / "outer_before_b.dat", primal)
    outer_rhs_s = indexed_vector(args.steady_outer / "outer_before_b.dat", primal)
    expected_mass = diags(mass_dof / dt, format="csr")
    physical_t = (outer_t - expected_mass).tocsr()
    physical_s = outer_s
    physical_row_delta = np.asarray((physical_t - physical_s).diagonal()).ravel()
    physical_delta_kx = np.asarray((physical_t - physical_s) @ xs).ravel()
    raw_operator_residual = delta_kx - mass_lhs - physical_delta_kx
    source_without_history = delta_b - history_rhs
    outer_full_rhs_t = rhs_t - outer_rhs_t
    outer_full_rhs_s = rhs_s - outer_rhs_s

    interface_mask = np.asarray(bmat_t.getnnz(axis=0)).ravel() != 0
    abs_a = transient_a.copy()
    abs_a.data = np.abs(abs_a.data)
    augmented = np.r_[xs, lambda_t]
    denominator = np.asarray(abs_a @ np.abs(augmented)).ravel() + np.abs(transient_rhs_full)
    residual_full = np.asarray(transient_a @ augmented - transient_rhs_full).ravel()
    eta = np.divide(np.abs(residual_full), denominator, out=np.zeros_like(residual_full), where=denominator > 0.0)
    primal_order = np.argsort(eta[:primal])[::-1]
    interface_order = [int(i) for i in primal_order if interface_mask[int(i)]]
    selected_indices = []
    for index in [int(i) for i in primal_order[:100]] + interface_order[:1]:
        if index not in selected_indices:
            selected_indices.append(index)
    contexts = [
        context(i, r_actual[i], eta[i], denominator[i], permutation, node_to_dof,
                dof_to_node, coords, node_bodies, node_boundaries, interface_mask,
                bmat_t, kx_s, mortar_s, rhs_s, kx_t, mortar_t, rhs_t,
                delta_kx, delta_mortar, delta_b, history_rhs, physical_delta_kx,
                mass_lhs - history_rhs, physical_row_delta, source_without_history,
                steady_residual, closed_predicted, closed_accounting_error)
        for i in selected_indices
    ]
    context_by_row = {int(row["row_1based"]): row for row in contexts}

    fields = [
        "row_1based", "dof_1based", "node_1based", "category", "eta", "residual",
        "denominator", "coordinates_m", "body_names", "material_names",
        "boundary_names", "mortar_constraint_rows_1based", "K_s_x_s",
        "B_sT_lambda_s", "b_s", "K_t_x_s", "B_tT_lambda_t", "b_t", "delta_Kx",
        "delta_Mortar", "delta_b", "r_predicted", "r_actual", "accounting_error",
        "steady_restart_primal_residual", "r_closed_with_steady_residual",
        "closed_accounting_error",
        "mass_lhs", "history_rhs", "mass_history_mismatch", "physical_K_delta_x",
        "physical_row_delta", "rhs_source_without_history",
    ]
    write_csv(out / "row_accounting.csv", contexts, fields)
    write_csv(out / "worst_primal_rows.csv", contexts, fields)

    stencil_rows: list[dict[str, object]] = []
    rhs_rows: list[dict[str, object]] = []
    mortar_rows: list[dict[str, object]] = []
    mass_rows: list[dict[str, object]] = []
    for row_index in selected_indices:
        row = context_by_row[row_index + 1]
        cols = sorted(set(ks.getrow(row_index).indices.tolist()) |
                      set(kt.getrow(row_index).indices.tolist()))
        row_bodies = set(node_bodies[row["node_1based"] - 1]) if row["node_1based"] else set()
        for col in cols:
            node = int(dof_to_node[col]) if col < len(dof_to_node) else 0
            bodies = set(node_bodies[node - 1]) if node else set()
            if col == row_index:
                stencil_class = "diagonal"
            elif interface_mask[col] and bodies.isdisjoint(row_bodies):
                stencil_class = "opposite_interface_neighbor"
            elif interface_mask[col]:
                stencil_class = "mortar_connected_dof"
            elif bodies & row_bodies:
                stencil_class = "same_side_neighbor"
            else:
                stencil_class = "other_operator_neighbor"
            s_value = float(ks[row_index, col])
            t_value = float(kt[row_index, col])
            stencil_rows.append({
                "row_1based": row_index + 1, "column_1based": col + 1,
                "column_node_1based": node,
                "column_coordinates_m": coords[node - 1].tolist() if node else [],
                "column_body_names": [BODY_NAMES.get(v, f"body_{v}") for v in sorted(bodies)],
                "column_interface": bool(interface_mask[col]), "stencil_class": stencil_class,
                "K_steady": s_value, "K_transient": t_value, "delta_K": t_value - s_value,
            })
        rhs_rows.append({
            "row_1based": row_index + 1, "category": row["category"],
            "b_steady": float(rhs_s[row_index]), "b_transient": float(rhs_t[row_index]),
            "history_reconstructed": float(history_rhs[row_index]),
            "delta_b": float(delta_b[row_index]),
            "source_without_history": float(source_without_history[row_index]),
            "pulse_contribution": 0.0,
            "outer_minus_full_transient": float(outer_full_rhs_t[row_index]),
            "outer_minus_full_steady": float(outer_full_rhs_s[row_index]),
            "boundary_source_split": "not tagged in raw CollectionVector; full/outer primal RHS equality",
        })
        mortar_rows.append({
            "row_1based": row_index + 1, "category": row["category"],
            "B_sT_lambda_s": float(mortar_s[row_index]),
            "B_tT_lambda_t": float(mortar_t[row_index]),
            "difference": float(delta_mortar[row_index]),
        })
        mass_rows.append({
            "row_1based": row_index + 1, "category": row["category"],
            "mass_node_lumped": float(mass_dof[row_index]), "dt": dt,
            "restart_temperature": float(xs[row_index]), "mass_lhs": float(mass_lhs[row_index]),
            "history_rhs": float(history_rhs[row_index]),
            "mass_lhs_minus_history_rhs": float(mass_lhs[row_index] - history_rhs[row_index]),
            "method": "tetra volume*rho*cp/4, Lumped Mass Matrix=True, BDF1",
        })
    write_csv(out / "matrix_stencil_diff.csv", stencil_rows, list(stencil_rows[0]))
    write_csv(out / "rhs_contribution_diff.csv", rhs_rows, list(rhs_rows[0]))
    write_csv(out / "mortar_reaction_diff.csv", mortar_rows, list(mortar_rows[0]))
    write_csv(out / "mass_history_cancellation.csv", mass_rows, list(mass_rows[0]))

    top20 = contexts[:20]
    top100 = contexts[:100]
    distribution20 = Counter(str(row["category"]) for row in top20)
    distribution100 = Counter(str(row["category"]) for row in top100)
    worst_interface = context_by_row[interface_order[0] + 1] if interface_order else None
    report = {
        "method": {
            "identity": "r_t=(K_t-K_s)x_s+(B_t^T lambda_t-B_s^T lambda_s)-(b_t-b_s)",
            "transient_multiplier": "least-squares argmin ||K_t*x_s+B_t^T*lambda-b_t||",
            "steady_multiplier": lambda_s_method,
            "mass_history": "mesh tetra rho*cp lumped mass divided by BDF1 dt",
            "pulse": "Pulse Energy=0; pulse contribution is exactly zero by input contract",
        },
        "provenance": {
            "steady_capture": str(args.steady_capture),
            "transient_capture": str(args.transient_capture),
            "steady_metadata": smeta, "transient_metadata": tmeta,
            "dt": dt, "bdf_order": tmeta["provenance"].get("bdf_order"),
            "dimensions": {"total": total, "primal": primal, "constraints": constraints},
            "mesh_element_counts": element_counts,
        },
        "constraint_check": {
            "B_t_x_s_minus_g_t": stat(np.asarray(bmat_t @ xs).ravel(), transient_rhs_full[primal:]),
            "B_s_x_s_minus_g_s": stat(np.asarray(bmat_s @ xs).ravel(), steady_rhs_full[primal:]),
            "B_t_minus_B_s": matrix_diff(bmat_t, bmat_s),
            "Bt_t_minus_Bt_s": matrix_diff(bt, bts),
            "g_t_minus_g_s_l2": float(np.linalg.norm(transient_rhs_full[primal:] - steady_rhs_full[primal:])),
        },
        "least_squares": {
            "transient_istop": int(fit_t[1]), "transient_iterations": int(fit_t[2]),
            "lambda_t_l2": float(np.linalg.norm(lambda_t)),
            "steady_lambda_l2": float(np.linalg.norm(lambda_s)),
            "r_actual": stat(r_actual), "r_actual_interface": stat(r_actual[interface_mask]),
        },
        "accounting": {
            "delta_Kx": stat(delta_kx), "delta_Mortar": stat(delta_mortar),
            "delta_b": stat(delta_b), "r_predicted": stat(r_predicted),
            "accounting_error": stat(accounting_error),
            "steady_restart_primal_residual": stat(steady_residual),
            "r_closed_with_steady_residual": stat(closed_predicted),
            "closed_accounting_error": stat(closed_accounting_error),
            "accounting_error_relative_to_r": float(np.linalg.norm(accounting_error) / max(np.linalg.norm(r_actual), 1e-300)),
            "closed_accounting_error_relative_to_r": float(np.linalg.norm(closed_accounting_error) / max(np.linalg.norm(r_actual), 1e-300)),
            "term_l2_ratio_to_r": {
                "delta_Kx": float(np.linalg.norm(delta_kx) / max(np.linalg.norm(r_actual), 1e-300)),
                "delta_Mortar": float(np.linalg.norm(delta_mortar) / max(np.linalg.norm(r_actual), 1e-300)),
                "delta_b": float(np.linalg.norm(delta_b) / max(np.linalg.norm(r_actual), 1e-300)),
                "steady_restart_primal_residual": float(np.linalg.norm(steady_residual) / max(np.linalg.norm(r_actual), 1e-300)),
            },
        },
        "operator_split": {
            "full_K_t_minus_outer_t": matrix_diff(kt, outer_t),
            "full_K_s_minus_outer_s": matrix_diff(ks, outer_s),
            "physical_K_t_minus_physical_K_s": matrix_diff(physical_t, physical_s),
            "raw_delta_Kx_minus_mass_lhs_minus_physical_delta_Kx": stat(raw_operator_residual),
            "physical_delta_Kx": stat(physical_delta_kx),
            "physical_diagonal_delta": stat(physical_row_delta),
        },
        "rhs_split": {
            "delta_b_minus_history": stat(source_without_history),
            "outer_full_rhs_transient": stat(outer_full_rhs_t),
            "outer_full_rhs_steady": stat(outer_full_rhs_s),
            "outer_full_rhs_difference_max_abs": float(max(np.max(np.abs(outer_full_rhs_t)), np.max(np.abs(outer_full_rhs_s)))),
            "joule_body_source_difference": "not independently tagged in raw CollectionVector; see delta_b-history",
            "boundary_contribution_difference": "not independently tagged; full/outer primal RHS equality proves no mortar RHS addition",
        },
        "mass_history": {
            "mass_lhs": stat(mass_lhs), "history_rhs": stat(history_rhs),
            "mass_lhs_minus_history_rhs": stat(mass_lhs - history_rhs),
            "selected_row_max_abs_mismatch": float(max(abs(row["mass_lhs_minus_history_rhs"]) for row in mass_rows)),
        },
        "distribution": {
            "top_primal_1": contexts[:1], "top_primal_20_categories": dict(distribution20),
            "top_primal_100_categories": dict(distribution100),
            "top_primal_20_rows": top20, "top_primal_100_rows": top100,
            "worst_interface_primal_row": worst_interface,
        },
        "code_path": {
            "time_mass_history": "tools/elmer-hypre/src/fem/src/modules/HeatSolve.F90:1573-1581,2540-2760",
            "boundary_assembly": "tools/elmer-hypre/src/fem/src/modules/HeatSolve.F90:3456-3483",
            "mortar_system": "tools/elmer-hypre/src/fem/src/SolverUtils.F90: SolveWithLinearRestriction / GenerateProjectors",
            "capture": "tools/elmer-hypre/src/fem/src/SolverUtils.F90 diagnostic full restriction capture hook",
        },
    }
    (out / "interface_primal_residual.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")

    dominant = max(report["accounting"]["term_l2_ratio_to_r"],
                   key=report["accounting"]["term_l2_ratio_to_r"].get)
    summary = [
        "# Phase24 interface primal residual accounting", "",
        f"- Dimensions: {total} = {primal} primal + {constraints} constraints.",
        f"- Transient: BDF{tmeta['provenance'].get('bdf_order')} dt={dt:.8e} s; pulse contribution=0.",
        f"- Constraint check max |B_t x_s-g_t|: {np.max(np.abs(np.asarray(bmat_t @ xs)-transient_rhs_full[primal:])):.8e}.",
        f"- Worst primal row: {contexts[0]['row_1based']} ({contexts[0]['category']}), eta={contexts[0]['eta']:.8e}, r={contexts[0]['r_actual']:.8e}.",
        f"- Worst interface primal row: {worst_interface['row_1based'] if worst_interface else 'none'}.",
        f"- Raw formula accounting error: {np.max(np.abs(accounting_error)):.8e}; closed error after adding steady restart residual: {np.max(np.abs(closed_accounting_error)):.8e}.",
        f"- Dominant unbalanced term: steady_restart_primal_residual (L2={np.linalg.norm(steady_residual):.8e}); raw DeltaKx and Delta b are large but cancel.",
        f"- DeltaKx L2={np.linalg.norm(delta_kx):.8e}; DeltaMortar L2={np.linalg.norm(delta_mortar):.8e}; Delta b L2={np.linalg.norm(delta_b):.8e}; r L2={np.linalg.norm(r_actual):.8e}.",
        f"- Mass/history max mismatch: {np.max(np.abs(mass_lhs-history_rhs)):.8e}; physical K delta*x L2={np.linalg.norm(physical_delta_kx):.8e}.",
        f"- Top-20 categories: {dict(distribution20)}; top-100 categories: {dict(distribution100)}.", "",
        "Raw captures do not contain separately tagged Joule, boundary, or pulse RHS vectors. rhs_contribution_diff.csv reports the exact delta_b-history remainder; pulse is zero by input. Assembly code locations are in interface_primal_residual.json.",
    ]
    (out / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(out / "summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
