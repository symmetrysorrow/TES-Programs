"""Compare Elmer's emitted K solves with the independent SciPy/SuperLU solve."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from validate_matrix_free_schur import deterministic_vectors, read_vector, triplets


def relative_error(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right) / max(np.linalg.norm(right), np.finfo(float).tiny))


def relative_norm(value: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(value) / max(np.linalg.norm(reference), np.finfo(float).tiny))


def componentwise_backward_error(matrix, rhs: np.ndarray, solution: np.ndarray) -> float:
    residual = rhs - matrix @ solution
    denominator = np.asarray(np.abs(matrix).dot(np.abs(solution))).reshape(-1) + np.abs(rhs)
    denominator = np.maximum(denominator,
                            np.finfo(float).eps * max(float(np.max(denominator)),
                                                      float(np.linalg.norm(rhs)),
                                                      np.finfo(float).tiny))
    return float(np.max(np.abs(residual) /
                        np.maximum(denominator, np.finfo(float).tiny)))


def read_matrix_triplets(path: Path) -> tuple[tuple[int, int, int], np.ndarray,
                                             np.ndarray, np.ndarray]:
    with path.open(encoding="utf-8") as stream:
        header = stream.readline().strip().split()
    if len(header) != 4 or header[0] != "#":
        raise ValueError(f"expected Elmer matrix header in {path}")
    shape = (int(header[1]), int(header[2]))
    stored_nnz = int(header[3])
    if stored_nnz == 0:
        data = np.empty((0, 3), dtype=np.float64)
    else:
        data = np.loadtxt(path, dtype=np.float64, comments="#")
    if data.size == 0:
        data = np.empty((0, 3), dtype=np.float64)
    elif data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 3:
        raise ValueError(f"expected row, column, value triplets in {path}")
    if len(data) != stored_nnz:
        raise ValueError(f"header nnz does not match data in {path}")
    return shape, data[:, 0].astype(np.int64) - 1, data[:, 1].astype(np.int64) - 1, data[:, 2]


def fingerprint(shape: tuple[int, int], rows: np.ndarray, cols: np.ndarray,
                values: np.ndarray) -> str:
    order = np.lexsort((cols, rows))
    digest = hashlib.sha256()
    digest.update(f"{shape[0]} {shape[1]} {len(values)}\n".encode("ascii"))
    for index in order:
        digest.update(f"{int(rows[index])} {int(cols[index])} ".encode("ascii"))
        digest.update(np.float64(values[index]).tobytes())
    return digest.hexdigest()


def compare_block(name: str, expected_shape: tuple[int, int], expected_rows: np.ndarray,
                 expected_cols: np.ndarray, expected_values: np.ndarray,
                 actual_path: Path) -> dict[str, object]:
    if not actual_path.exists():
        return {"name": name, "match": False, "failure": f"missing {actual_path}"}
    actual_shape, actual_rows, actual_cols, actual_values = read_matrix_triplets(actual_path)
    expected_nonzero = expected_values != 0.0
    actual_nonzero = actual_values != 0.0
    expected_map = {(int(row), int(col)): float(value)
                    for row, col, value in zip(expected_rows[expected_nonzero],
                                                expected_cols[expected_nonzero],
                                                expected_values[expected_nonzero])}
    actual_map = {(int(row), int(col)): float(value)
                  for row, col, value in zip(actual_rows[actual_nonzero],
                                              actual_cols[actual_nonzero],
                                              actual_values[actual_nonzero])}
    common = set(expected_map) & set(actual_map)
    missing = set(expected_map) - set(actual_map)
    extra = set(actual_map) - set(expected_map)
    differences = np.array([actual_map[key] - expected_map[key] for key in common])
    expected_common = np.array([expected_map[key] for key in common])
    extra_values = np.array([actual_map[key] for key in extra])
    max_abs = float(np.max(np.abs(differences))) if len(differences) else 0.0
    max_rel = float(np.max(np.abs(differences) /
                           np.maximum(np.abs(expected_common), np.finfo(float).tiny))) \
        if len(differences) else 0.0
    canonical_expected_rows = np.array(list(expected_map.keys()), dtype=np.int64)[:, 0] \
        if expected_map else np.empty(0, dtype=np.int64)
    canonical_expected_cols = np.array(list(expected_map.keys()), dtype=np.int64)[:, 1] \
        if expected_map else np.empty(0, dtype=np.int64)
    canonical_expected_values = np.array(list(expected_map.values()), dtype=np.float64)
    canonical_actual_rows = np.array(list(actual_map.keys()), dtype=np.int64)[:, 0] \
        if actual_map else np.empty(0, dtype=np.int64)
    canonical_actual_cols = np.array(list(actual_map.keys()), dtype=np.int64)[:, 1] \
        if actual_map else np.empty(0, dtype=np.int64)
    canonical_actual_values = np.array(list(actual_map.values()), dtype=np.float64)
    actual_fingerprint = fingerprint(actual_shape, actual_rows, actual_cols, actual_values)
    expected_fingerprint = fingerprint(expected_shape, expected_rows, expected_cols, expected_values)
    canonical_expected_fingerprint = fingerprint(
        expected_shape, canonical_expected_rows, canonical_expected_cols,
        canonical_expected_values)
    canonical_actual_fingerprint = fingerprint(
        expected_shape, canonical_actual_rows, canonical_actual_cols,
        canonical_actual_values)
    return {
        "name": name,
        "expected_shape": list(expected_shape),
        "actual_shape": list(actual_shape),
        "expected_nnz": int(len(expected_values)),
        "actual_nnz": int(len(actual_values)),
        "expected_sha256": expected_fingerprint,
        "actual_sha256": actual_fingerprint,
        "expected_canonical_sha256": canonical_expected_fingerprint,
        "actual_canonical_sha256_padded": canonical_actual_fingerprint,
        "missing_nonzero_entries": len(missing),
        "extra_nonzero_entries": len(extra),
        "extra_nonzero_max_abs": float(np.max(np.abs(extra_values))) if len(extra_values) else 0.0,
        "extra_nonzero_l2": float(np.linalg.norm(extra_values)) if len(extra_values) else 0.0,
        "max_abs_value_difference": max_abs,
        "max_relative_value_difference": max_rel,
        "numerical_match_after_zero_pruning": bool(not missing and not extra and
                                                     max_abs == 0.0 and
                                                     actual_shape[1] == expected_shape[1]),
        "match": bool(actual_shape == expected_shape and
                      expected_fingerprint == actual_fingerprint),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--c-start", type=int, required=True)
    parser.add_argument("--elmer-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows, cols, values = triplets(args.matrix)
    n_f = args.c_start - 1
    n_c = args.rows - n_f
    matrix = coo_matrix((values, (rows, cols)), shape=(args.rows, args.rows)).tocsr()
    k = matrix[:n_f, :n_f].tocsc()
    bt = matrix[:n_f, n_f:].tocsr()
    b = matrix[n_f:, :n_f].tocsr()
    d = matrix[n_f:, n_f:].tocsr()
    lu = splu(k, permc_spec="COLAMD")

    raw_rows = rows
    raw_cols = cols
    raw_values = values
    block_specs = {
        "K": ((n_f, n_f), raw_rows < n_f, raw_cols < n_f, raw_rows, raw_cols),
        "B": ((n_c, n_f), raw_rows >= n_f, raw_cols < n_f, raw_rows - n_f, raw_cols),
        "Bt": ((n_f, n_c), raw_rows < n_f, raw_cols >= n_f, raw_rows, raw_cols - n_f),
        "D": ((n_c, n_c), raw_rows >= n_f, raw_cols >= n_f,
              raw_rows - n_f, raw_cols - n_f),
    }
    block_fingerprints = []
    for name, (shape, row_mask, col_mask, local_rows, local_cols) in block_specs.items():
        mask = row_mask & col_mask
        block_fingerprints.append(compare_block(
            name, shape, local_rows[mask], local_cols[mask], raw_values[mask],
            Path(f"{args.elmer_prefix}_{name}.triplets")))

    records: list[dict[str, object]] = []
    for number, (name, vector) in enumerate(deterministic_vectors(n_c), start=1):
        rhs = bt @ vector
        scipy_ku = lu.solve(rhs)
        ku_path = Path(f"{args.elmer_prefix}_ku{number}.dat")
        y_path = Path(f"{args.elmer_prefix}_y{number}.dat")
        if not ku_path.exists() or not y_path.exists():
            records.append({"name": name, "pass": False,
                            "failure": f"missing parity output: {ku_path} or {y_path}"})
            continue
        elmer_ku = read_vector(ku_path, n_f)
        emitted_y = read_vector(y_path, n_c)
        scipy_action = d @ vector - b @ scipy_ku
        action_from_elmer_ku = d @ vector - b @ elmer_ku
        elmer_residual = rhs - k @ elmer_ku
        scipy_residual = rhs - k @ scipy_ku
        bt_path = Path(f"{args.elmer_prefix}_bt{number}.dat")
        bku_path = Path(f"{args.elmer_prefix}_bku{number}.dat")
        dv_path = Path(f"{args.elmer_prefix}_dv{number}.dat")
        stage = {}
        if bt_path.exists() and bku_path.exists() and dv_path.exists():
            emitted_bt = read_vector(bt_path, n_f)
            emitted_bku = read_vector(bku_path, n_c)
            emitted_dv = read_vector(dv_path, n_c)
            stage = {
                "bt_v_relative_error": relative_error(emitted_bt, rhs),
                "b_ku_relative_error": relative_error(emitted_bku, b @ elmer_ku),
                "d_v_relative_error": relative_error(emitted_dv, d @ vector),
                "y_reconstruction_relative_error": relative_error(
                    emitted_y, emitted_dv - emitted_bku),
                "y_vs_explicit_relative_error": relative_error(emitted_y, action_from_elmer_ku),
                "stage_outputs_present": True,
            }
        else:
            stage = {"stage_outputs_present": False,
                     "stage_failure": f"missing {bt_path}, {bku_path}, or {dv_path}"}
        records.append({
            "name": name,
            "ku_relative_error": relative_error(elmer_ku, scipy_ku),
            "action_relative_error_vs_scipy": relative_error(emitted_y, scipy_action),
            "action_relative_error_vs_elmer_ku": relative_error(emitted_y, action_from_elmer_ku),
            "elmer_backward_residual_norm": float(np.linalg.norm(elmer_residual)),
            "elmer_backward_residual_relative_to_rhs": relative_norm(elmer_residual, rhs),
            "elmer_componentwise_backward_error": componentwise_backward_error(k, rhs, elmer_ku),
            "scipy_backward_residual_norm": float(np.linalg.norm(scipy_residual)),
            "scipy_backward_residual_relative_to_rhs": relative_norm(scipy_residual, rhs),
            **stage,
            "scipy_ku_norm": float(np.linalg.norm(scipy_ku)),
            "elmer_ku_norm": float(np.linalg.norm(elmer_ku)),
            "pass": bool(relative_error(elmer_ku, scipy_ku) <= 1.0e-10 and
                          relative_error(emitted_y, action_from_elmer_ku) <= 1.0e-12 and
                          stage.get("stage_outputs_present", False) and
                          stage.get("bt_v_relative_error", np.inf) <= 1.0e-12 and
                          stage.get("b_ku_relative_error", np.inf) <= 1.0e-12 and
                          stage.get("d_v_relative_error", np.inf) <= 1.0e-12 and
                          stage.get("y_reconstruction_relative_error", np.inf) <= 1.0e-12),
        })

    result = {
        "n_f": n_f,
        "n_c": n_c,
        "scipy_solver": "SuperLU via scipy.sparse.linalg.splu, COLAMD",
        "elmer_solver": "system SuperLU emitted by the diagnostic",
        "block_fingerprints": block_fingerprints,
        "block_fingerprints_match": all(item["match"] for item in block_fingerprints),
        "records": records,
        "passed": bool(all(record["pass"] for record in records) and
                       all(item["match"] for item in block_fingerprints)),
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
