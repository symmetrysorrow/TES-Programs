"""Compare Elmer's emitted K solves with the independent SciPy/SuperLU solve."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

from validate_matrix_free_schur import deterministic_vectors, read_vector, triplets


def relative_error(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right) / max(np.linalg.norm(right), np.finfo(float).tiny))


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
        records.append({
            "name": name,
            "ku_relative_error": relative_error(elmer_ku, scipy_ku),
            "action_relative_error_vs_scipy": relative_error(emitted_y, scipy_action),
            "action_relative_error_vs_elmer_ku": relative_error(emitted_y, action_from_elmer_ku),
            "scipy_ku_norm": float(np.linalg.norm(scipy_ku)),
            "elmer_ku_norm": float(np.linalg.norm(elmer_ku)),
            "pass": bool(relative_error(elmer_ku, scipy_ku) <= 1.0e-10 and
                          relative_error(emitted_y, action_from_elmer_ku) <= 1.0e-12),
        })

    result = {
        "n_f": n_f,
        "n_c": n_c,
        "scipy_solver": "SuperLU via scipy.sparse.linalg.splu, COLAMD",
        "elmer_solver": "system SuperLU emitted by the diagnostic",
        "records": records,
        "passed": all(record["pass"] for record in records),
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
