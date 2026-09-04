"""Validate Elmer's matrix-free Schur action against an independent oracle.

The Elmer diagnostic writes one ``<prefix>_vN.dat`` and ``<prefix>_yN.dat``
pair for each deterministic constraint vector.  This script reconstructs
``S = D - B K^{-1} B^T`` with an independent SuperLU factorization of K and
reports action errors without solving the full saddle system.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu


def triplets(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 3:
        raise ValueError(f"expected row, column, value triplets in {path}")
    return data[:, 0].astype(np.int64) - 1, data[:, 1].astype(np.int64) - 1, data[:, 2]


def read_vector(path: Path, size: int) -> np.ndarray:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError(f"expected index and value columns in {path}")
    result = np.zeros(size, dtype=np.float64)
    indexes = data[:, 0].astype(np.int64) - 1
    result[indexes] = data[:, -1]
    return result


def deterministic_vectors(n: int) -> list[tuple[str, np.ndarray]]:
    vectors = [
        ("all_ones", np.ones(n)),
        ("alternating_sign", np.where(np.arange(n) % 2 == 0, 1.0, -1.0)),
        ("deterministic_sine", np.sin(0.017 * (np.arange(n) + 1))),
    ]
    basis = np.zeros(n)
    basis[min(16, n - 1)] = 1.0
    vectors.append(("basis_like_17", basis))
    return vectors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--c-start", type=int, required=True,
                        help="one-based first constraint row")
    parser.add_argument("--elmer-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows, cols, values = triplets(args.matrix)
    n_f = args.c_start - 1
    n_c = args.rows - n_f
    if not (0 < n_f < args.rows):
        raise ValueError("invalid block split")
    matrix = coo_matrix((values, (rows, cols)), shape=(args.rows, args.rows)).tocsr()
    k = matrix[:n_f, :n_f].tocsr()
    bt = matrix[:n_f, n_f:].tocsr()
    b = matrix[n_f:, :n_f].tocsr()
    d = matrix[n_f:, n_f:].tocsr()
    lu = splu(k.tocsc(), permc_spec="COLAMD")

    records: list[dict[str, object]] = []
    for number, (name, vector) in enumerate(deterministic_vectors(n_c), start=1):
        exact = d @ vector - b @ lu.solve(bt @ vector)
        v_path = Path(f"{args.elmer_prefix}_v{number}.dat")
        y_path = Path(f"{args.elmer_prefix}_y{number}.dat")
        if not v_path.exists() or not y_path.exists():
            raise FileNotFoundError(f"missing Elmer diagnostic pair: {v_path}, {y_path}")
        emitted_v = read_vector(v_path, n_c)
        emitted_y = read_vector(y_path, n_c)
        vector_error = np.linalg.norm(emitted_v - vector) / max(np.linalg.norm(vector), np.finfo(float).tiny)
        action_error = np.linalg.norm(emitted_y - exact) / max(np.linalg.norm(exact), np.finfo(float).tiny)
        records.append({
            "name": name,
            "vector_relative_error": float(vector_error),
            "action_relative_error": float(action_error),
            "exact_action_norm": float(np.linalg.norm(exact)),
            "elmer_action_norm": float(np.linalg.norm(emitted_y)),
            "pass": bool(action_error <= 1.0e-10 and vector_error <= 1.0e-14),
        })

    result = {
        "n_f": n_f,
        "n_c": n_c,
        "operator": "D - B K^-1 B^T",
        "k_solver": "independent SciPy SuperLU",
        "vectors": records,
        "max_action_relative_error": max(item["action_relative_error"] for item in records),
        "passed": all(item["pass"] for item in records),
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
