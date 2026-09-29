"""Analyze the saved one-step explicit mortar matrix by F/C row blocks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix


def load_triplets(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 3:
        raise ValueError(f"expected 3 columns in {path}, got {data.shape[1]}")
    rows = data[:, 0].astype(np.int64) - 1
    cols = data[:, 1].astype(np.int64) - 1
    vals = data[:, 2]
    return rows, cols, vals


def frob(a: csr_matrix) -> float:
    return float(np.sqrt(np.dot(a.data, a.data)))


def row_norm_range(a: csr_matrix) -> tuple[float, float]:
    norms = np.sqrt(np.asarray(a.multiply(a).sum(axis=1)).ravel())
    return float(norms.min()), float(norms.max())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--c-start", type=int, required=True,
                        help="one-based first C row in the assembled matrix")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows, cols, vals = load_triplets(args.matrix)
    n = args.rows
    if rows.min() < 0 or cols.min() < 0 or rows.max() >= n or cols.max() >= n:
        raise ValueError("matrix index outside requested size")
    n_f = args.c_start - 1
    n_c = n - n_f
    if not (0 < n_f < n):
        raise ValueError("invalid C split")

    a = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    k = a[:n_f, :n_f].tocsr()
    bt = a[:n_f, n_f:].tocsr()
    b = a[n_f:, :n_f].tocsr()
    d = a[n_f:, n_f:].tocsr()

    b_diff = b - bt.transpose()
    k_diff = k - k.transpose()
    diag = k.diagonal()
    result = {
        "n_rows": n,
        "n_f": n_f,
        "n_c": n_c,
        "nnz": {"A": int(a.nnz), "K": int(k.nnz), "B": int(b.nnz),
                "Bt": int(bt.nnz), "D": int(d.nnz)},
        "norm": {"A": frob(a), "K": frob(k), "B": frob(b), "Bt": frob(bt),
                 "D": frob(d)},
        "relative": {
            "B_minus_BtT_over_B": frob(b_diff) / max(frob(b), np.finfo(float).tiny),
            "K_symmetry_error": frob(k_diff) / max(frob(k), np.finfo(float).tiny),
        },
        "K_diagonal_min_max": [float(diag.min()), float(diag.max())],
        "K_row_norm_min_max": list(row_norm_range(k)),
        "B_row_norm_min_max": list(row_norm_range(b)),
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    if args.output:
        args.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
