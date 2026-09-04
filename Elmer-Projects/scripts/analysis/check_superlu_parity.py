"""Separate Elmer-block self consistency from monolithic block equivalence."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu

from validate_matrix_free_schur import deterministic_vectors, read_vector, triplets


def relative_error(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right) /
                 max(np.linalg.norm(right), np.finfo(float).tiny))


def relative_norm(value: np.ndarray, reference: np.ndarray) -> float:
    return float(np.linalg.norm(value) /
                 max(np.linalg.norm(reference), np.finfo(float).tiny))


def componentwise_backward_error(matrix: csr_matrix, rhs: np.ndarray,
                                 solution: np.ndarray) -> float:
    residual = rhs - matrix @ solution
    denominator = np.asarray(np.abs(matrix).dot(np.abs(solution))).reshape(-1) + np.abs(rhs)
    floor = np.finfo(float).eps * max(float(np.max(denominator)),
                                      float(np.linalg.norm(rhs)),
                                      np.finfo(float).tiny)
    return float(np.max(np.abs(residual) / np.maximum(denominator, floor)))


def read_matrix_triplets(path: Path) -> tuple[tuple[int, int], np.ndarray,
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
        if data.ndim == 1:
            data = data.reshape(1, -1)
    if data.shape[1] != 3 or len(data) != stored_nnz:
        raise ValueError(f"invalid triplet count in {path}")
    return shape, data[:, 0].astype(np.int64) - 1, data[:, 1].astype(np.int64) - 1, data[:, 2]


def canonical_triplets(shape: tuple[int, int], rows: np.ndarray, cols: np.ndarray,
                       values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    matrix = coo_matrix((values, (rows, cols)), shape=shape).tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    coo = matrix.tocoo()
    order = np.lexsort((coo.col, coo.row))
    return coo.row[order], coo.col[order], coo.data[order]


def fingerprint(shape: tuple[int, int], rows: np.ndarray, cols: np.ndarray,
                values: np.ndarray) -> str:
    order = np.lexsort((cols, rows))
    digest = hashlib.sha256()
    digest.update(f"{shape[0]} {shape[1]} {len(values)}\n".encode("ascii"))
    for index in order:
        digest.update(f"{int(rows[index])} {int(cols[index])} ".encode("ascii"))
        digest.update(np.float64(values[index]).tobytes())
    return digest.hexdigest()


def sparse_from_triplets(shape: tuple[int, int], rows: np.ndarray, cols: np.ndarray,
                         values: np.ndarray) -> csr_matrix:
    matrix = coo_matrix((values, (rows, cols)), shape=shape).tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return matrix


def load_elmer_block(path: Path, logical_shape: tuple[int, int]) -> tuple[csr_matrix, dict[str, object]]:
    storage_shape, rows, cols, values = read_matrix_triplets(path)
    if len(rows) and (rows.min() < 0 or cols.min() < 0 or
                      rows.max() >= logical_shape[0] or cols.max() >= logical_shape[1]):
        raise ValueError(f"out-of-range entry in {path}")
    return sparse_from_triplets(logical_shape, rows, cols, values), {
        "path": str(path),
        "storage_shape": list(storage_shape),
        "logical_shape": list(logical_shape),
        "stored_nnz": int(len(values)),
        "stored_zero_entries": int(np.count_nonzero(values == 0.0)),
    }


def compare_block(name: str, expected_shape: tuple[int, int], expected_rows: np.ndarray,
                 expected_cols: np.ndarray, expected_values: np.ndarray,
                 actual_path: Path) -> dict[str, object]:
    if not actual_path.exists():
        return {"name": name, "raw_match": False, "canonical_match": False,
                "failure": f"missing {actual_path}"}
    actual_shape, actual_rows, actual_cols, actual_values = read_matrix_triplets(actual_path)
    exp_raw_sha = fingerprint(expected_shape, expected_rows, expected_cols, expected_values)
    act_raw_sha = fingerprint(actual_shape, actual_rows, actual_cols, actual_values)
    exp_r, exp_c, exp_v = canonical_triplets(expected_shape, expected_rows, expected_cols,
                                              expected_values)
    act_r, act_c, act_v = canonical_triplets(expected_shape, actual_rows, actual_cols,
                                              actual_values)
    expected_map = {(int(r), int(c)): float(v) for r, c, v in zip(exp_r, exp_c, exp_v)}
    actual_map = {(int(r), int(c)): float(v) for r, c, v in zip(act_r, act_c, act_v)}
    common = set(expected_map) & set(actual_map)
    missing = set(expected_map) - set(actual_map)
    extra = set(actual_map) - set(expected_map)
    differences = np.array([actual_map[key] - expected_map[key] for key in common])
    common_expected = np.array([expected_map[key] for key in common])
    extra_values = np.array([actual_map[key] for key in extra])
    actual_matrix = sparse_from_triplets(expected_shape, actual_rows, actual_cols, actual_values)
    expected_matrix = sparse_from_triplets(expected_shape, expected_rows, expected_cols, expected_values)
    delta = actual_matrix - expected_matrix
    expected_frob = float(np.sqrt(np.dot(expected_matrix.data, expected_matrix.data)))
    delta_frob = float(np.sqrt(np.dot(delta.data, delta.data)))
    max_abs = float(np.max(np.abs(differences))) if len(differences) else 0.0
    max_rel = float(np.max(np.abs(differences) /
                           np.maximum(np.abs(common_expected), np.finfo(float).tiny))) \
        if len(differences) else 0.0
    canonical_match = bool(actual_shape == expected_shape and len(missing) == 0 and
                           len(extra) == 0 and max_abs == 0.0)
    return {
        "name": name,
        "expected_shape": list(expected_shape),
        "actual_storage_shape": list(actual_shape),
        "expected_stored_nnz": int(len(expected_values)),
        "actual_stored_nnz": int(len(actual_values)),
        "expected_raw_sha256": exp_raw_sha,
        "actual_raw_sha256": act_raw_sha,
        "raw_match": bool(actual_shape == expected_shape and exp_raw_sha == act_raw_sha),
        "expected_canonical_sha256": fingerprint(expected_shape, exp_r, exp_c, exp_v),
        "actual_canonical_sha256_padded": fingerprint(expected_shape, act_r, act_c, act_v),
        "canonical_match": canonical_match,
        "missing_nonzeros": int(len(missing)),
        "extra_nonzeros": int(len(extra)),
        "extra_max_abs": float(np.max(np.abs(extra_values))) if len(extra_values) else 0.0,
        "extra_l2_norm": float(np.linalg.norm(extra_values)) if len(extra_values) else 0.0,
        "common_entry_max_abs_diff": max_abs,
        "common_entry_max_relative_diff": max_rel,
        "relative_frobenius_difference": delta_frob / max(expected_frob, np.finfo(float).tiny),
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
    full_shape = (args.rows, args.rows)
    raw_matrix = coo_matrix((values, (rows, cols)), shape=full_shape).tocsr()
    k_mono = raw_matrix[:n_f, :n_f].tocsr()
    bt_mono = raw_matrix[:n_f, n_f:].tocsr()
    b_mono = raw_matrix[n_f:, :n_f].tocsr()
    d_mono = raw_matrix[n_f:, n_f:].tocsr()
    lu_mono = splu(k_mono.tocsc(), permc_spec="COLAMD")

    raw_specs = {
        "K": ((n_f, n_f), rows < n_f, cols < n_f, rows, cols),
        "B": ((n_c, n_f), rows >= n_f, cols < n_f, rows - n_f, cols),
        "Bt": ((n_f, n_c), rows < n_f, cols >= n_f, rows, cols - n_f),
        "D": ((n_c, n_c), rows >= n_f, cols >= n_f, rows - n_f, cols - n_f),
    }
    block_equivalence = []
    for name, (shape, row_mask, col_mask, local_rows, local_cols) in raw_specs.items():
        mask = row_mask & col_mask
        block_equivalence.append(compare_block(
            name, shape, local_rows[mask], local_cols[mask], values[mask],
            Path(f"{args.elmer_prefix}_{name}.triplets")))

    k_elmer, k_meta = load_elmer_block(Path(f"{args.elmer_prefix}_K.triplets"), (n_f, n_f))
    b_elmer, b_meta = load_elmer_block(Path(f"{args.elmer_prefix}_B.triplets"), (n_c, n_f))
    bt_elmer, bt_meta = load_elmer_block(Path(f"{args.elmer_prefix}_Bt.triplets"), (n_f, n_c))
    d_elmer, d_meta = load_elmer_block(Path(f"{args.elmer_prefix}_D.triplets"), (n_c, n_c))
    lu_elmer = splu(k_elmer.tocsc(), permc_spec="COLAMD")

    b_vs_bt = b_elmer - bt_elmer.transpose().tocsr()
    b_norm = float(np.sqrt(np.dot(b_elmer.data, b_elmer.data)))
    b_vs_bt_norm = float(np.sqrt(np.dot(b_vs_bt.data, b_vs_bt.data)))
    b_vs_bt_report = {
        "relative_frobenius_difference": b_vs_bt_norm / max(b_norm, np.finfo(float).tiny),
        "max_abs_difference": float(np.max(np.abs(b_vs_bt.data))) if b_vs_bt.nnz else 0.0,
        "nnz_difference": int(b_vs_bt.nnz),
        "pass": bool(b_vs_bt.nnz == 0),
    }
    tail_rows = int(np.count_nonzero((raw_specs["Bt"][3][raw_specs["Bt"][1] & raw_specs["Bt"][2]] >= bt_meta["storage_shape"][0])))
    short_storage = {
        "stored_bt_shape": bt_meta["storage_shape"],
        "logical_bt_shape": [n_f, n_c],
        "monolithic_nonzeros_in_rows_after_stored_bt": tail_rows,
        "zero_padding_required": bool(bt_meta["storage_shape"][0] < n_f),
        "zero_padding_is_safe_for_emitted_output": tail_rows == 0,
    }

    self_records = []
    mono_records = []
    for number, (name, vector) in enumerate(deterministic_vectors(n_c), start=1):
        v = read_vector(Path(f"{args.elmer_prefix}_v{number}.dat"), n_c)
        emitted_bt = read_vector(Path(f"{args.elmer_prefix}_bt{number}.dat"), n_f)
        emitted_ku = read_vector(Path(f"{args.elmer_prefix}_ku{number}.dat"), n_f)
        emitted_bku = read_vector(Path(f"{args.elmer_prefix}_bku{number}.dat"), n_c)
        emitted_dv = read_vector(Path(f"{args.elmer_prefix}_dv{number}.dat"), n_c)
        emitted_y = read_vector(Path(f"{args.elmer_prefix}_y{number}.dat"), n_c)

        bt_ref = np.asarray(bt_elmer @ v).reshape(-1)
        ku_ref = lu_elmer.solve(bt_ref)
        bku_ref = np.asarray(b_elmer @ ku_ref).reshape(-1)
        dv_ref = np.asarray(d_elmer @ v).reshape(-1)
        y_ref = dv_ref - bku_ref
        bku_from_emitted = np.asarray(b_elmer @ emitted_ku).reshape(-1)
        y_from_emitted_stages = emitted_dv - emitted_bku
        elmer_residual = bt_ref - k_elmer @ emitted_ku
        scipy_residual = bt_ref - k_elmer @ ku_ref
        self_record = {
            "name": name,
            "vector_relative_error": relative_error(v, vector),
            "bt_self_relative_error": relative_error(emitted_bt, bt_ref),
            "ku_solution_relative_difference": relative_error(emitted_ku, ku_ref),
            "elmer_k_backward_residual_norm": float(np.linalg.norm(elmer_residual)),
            "elmer_k_backward_residual_relative": relative_norm(elmer_residual, bt_ref),
            "elmer_k_componentwise_backward_error": componentwise_backward_error(k_elmer, bt_ref, emitted_ku),
            "scipy_k_backward_residual_norm": float(np.linalg.norm(scipy_residual)),
            "scipy_k_backward_residual_relative": relative_norm(scipy_residual, bt_ref),
            "bku_oracle_relative_error": relative_error(emitted_bku, bku_ref),
            "bku_emitted_ku_matvec_relative_error": relative_error(emitted_bku, bku_from_emitted),
            "dv_oracle_relative_error": relative_error(emitted_dv, dv_ref),
            "schur_oracle_relative_error": relative_error(emitted_y, y_ref),
            "schur_oracle_absolute_error": float(np.linalg.norm(emitted_y - y_ref)),
            "schur_emitted_stage_reconstruction_relative_error": relative_error(
                emitted_y, y_from_emitted_stages),
        }
        self_record["emitted_stage_pass"] = bool(
            self_record["vector_relative_error"] <= 1.0e-14 and
            self_record["bt_self_relative_error"] <= 1.0e-12 and
            self_record["bku_emitted_ku_matvec_relative_error"] <= 1.0e-12 and
            self_record["dv_oracle_relative_error"] <= 1.0e-12 and
            self_record["schur_emitted_stage_reconstruction_relative_error"] <= 1.0e-10)
        self_record["oracle_pass"] = bool(
            self_record["vector_relative_error"] <= 1.0e-14 and
            self_record["bt_self_relative_error"] <= 1.0e-12 and
            self_record["bku_oracle_relative_error"] <= 1.0e-12 and
            self_record["dv_oracle_relative_error"] <= 1.0e-12 and
            self_record["schur_oracle_relative_error"] <= 1.0e-10)
        self_record["pass"] = self_record["oracle_pass"]
        self_records.append(self_record)

        elmer_action = np.asarray(d_elmer @ v).reshape(-1) - b_elmer @ lu_elmer.solve(bt_elmer @ v)
        mono_action = np.asarray(d_mono @ v).reshape(-1) - b_mono @ lu_mono.solve(bt_mono @ v)
        mono_records.append({
            "name": name,
            "elmer_block_action_norm": float(np.linalg.norm(elmer_action)),
            "monolithic_action_norm": float(np.linalg.norm(mono_action)),
            "absolute_action_difference": float(np.linalg.norm(elmer_action - mono_action)),
            "relative_action_difference": relative_error(elmer_action, mono_action),
        })

    raw_match = all(item.get("raw_match", False) for item in block_equivalence)
    canonical_match = all(item.get("canonical_match", False) for item in block_equivalence)
    max_frob = max(item.get("relative_frobenius_difference", 0.0) for item in block_equivalence)
    max_extra = max(item.get("extra_max_abs", 0.0) for item in block_equivalence)
    if raw_match:
        classification = "RAW MATCH"
    elif all(item.get("missing_nonzeros", 1) == 0 and item.get("common_entry_max_abs_diff", 1.0) == 0.0
             and item.get("extra_max_abs", 1.0) <= 1.0e-12 for item in block_equivalence):
        classification = "NUMERICALLY CLOSE"
    else:
        classification = "MATERIAL MISMATCH"

    self_pass = all(item["oracle_pass"] for item in self_records)
    emitted_stage_pass = all(item["emitted_stage_pass"] for item in self_records)
    result = {
        "metadata": {"n_f": n_f, "n_c": n_c, "operator": "D - B K^-1 B^T",
                     "elmer_solver": "system SuperLU wrapper; temporary CSR->CSC",
                     "monolithic_solver": "SciPy SuperLU, COLAMD",
                     "strict_gates": {"vector": 1.0e-14, "bt": 1.0e-12,
                                      "bku": 1.0e-12, "dv": 1.0e-12, "schur": 1.0e-10}},
        "self_consistency": {"records": self_records, "passed": self_pass,
                              "oracle_passed": self_pass,
                              "emitted_stage_passed": emitted_stage_pass},
        "b_vs_bt_transpose": b_vs_bt_report,
        "short_bt_storage": short_storage,
        "block_equivalence": {"blocks": block_equivalence,
                               "raw_match": raw_match,
                               "canonical_match": canonical_match,
                               "classification": classification,
                               "max_relative_frobenius_difference": max_frob,
                               "max_extra_abs": max_extra,
                               "elmer_block_metadata": [k_meta, b_meta, bt_meta, d_meta]},
        "monolithic_comparison": {"records": mono_records,
                                   "max_relative_action_difference": max(
                                       item["relative_action_difference"] for item in mono_records)},
        "matrix_free_matvec_self_consistency": emitted_stage_pass,
        "matrix_free_block_oracle_consistency": self_pass,
        "matrix_free_implementation_valid": self_pass,
        "block_extraction_equivalent": canonical_match,
        "production_ready": False,
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
