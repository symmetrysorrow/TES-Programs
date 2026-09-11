"""Trusted sparse direct reference for an exact Phase24 A/b dump."""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("a", type=Path)
    parser.add_argument("b", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--solution-npy", type=Path)
    args = parser.parse_args()
    a = np.loadtxt(args.a, dtype=np.float64)
    b_dump = np.loadtxt(args.b, dtype=np.float64)
    n = 87534
    rows = a[:, 0].astype(np.int64) - 1
    cols = a[:, 1].astype(np.int64) - 1
    values = a[:, 2]
    rhs = np.zeros(n, dtype=np.float64)
    rhs[b_dump[:, 0].astype(np.int64) - 1] = b_dump[:, 1]
    matrix = coo_matrix((values, (rows, cols)), shape=(n, n)).tocsr()
    started = time.perf_counter()
    solution = spsolve(matrix, rhs)
    wall = time.perf_counter() - started
    residual = matrix @ solution - rhs
    rel_residual = float(np.linalg.norm(residual) / np.linalg.norm(rhs))
    result = {
        "solver": "scipy.sparse.linalg.spsolve (SuperLU backend)",
        "rows": n,
        "nnz": int(matrix.nnz),
        "rhs_l2": float(np.linalg.norm(rhs)),
        "relative_residual": rel_residual,
        "max_abs_residual": float(np.max(np.abs(residual))),
        "wall_seconds": wall,
        "solution_sha256": hashlib.sha256(np.ascontiguousarray(solution).tobytes()).hexdigest(),
        "solution_l2": float(np.linalg.norm(solution)),
        "solution_max_abs": float(np.max(np.abs(solution))),
    }
    if args.solution_npy:
        np.save(args.solution_npy, np.ascontiguousarray(solution))
        result["solution_npy"] = str(args.solution_npy)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
