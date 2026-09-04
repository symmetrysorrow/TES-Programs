"""Compare explicit mortar Schur strategies on a saved one-rank system.

The script is a correctness/performance diagnostic, not the production
solver.  It evaluates the exact SuperLU Schur oracle, the cheap
``diag(K)^-1`` approximation, and a sparse-ILU K inverse action.  The latter
is a portable proxy for the matrix-free AMG action used by the native block
candidate; the actual HYPRE BoomerAMG path is validated by the Elmer smoke
case because SciPy cannot execute HYPRE cycles.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spilu, splu


def triplets(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] != 3:
        raise ValueError(f"expected three columns in {path}")
    return (data[:, 0].astype(np.int64) - 1,
            data[:, 1].astype(np.int64) - 1,
            data[:, 2])


def solve_with_action(name: str, K: csr_matrix, B: csr_matrix,
                      Bt: csr_matrix, D: csr_matrix, f: np.ndarray,
                      g: np.ndarray, A: csr_matrix, rhs: np.ndarray,
                      mumps: np.ndarray, ilu_drop: float,
                      ilu_fill: float) -> dict[str, object]:
    setup_start = time.perf_counter()
    if name == "exact":
        action = splu(K.tocsc(), permc_spec="COLAMD")
        z = action.solve(Bt.toarray())
        S = D.toarray() - B @ z
        ksolve = lambda q: action.solve(q)
    elif name == "diag":
        diag = K.diagonal()
        if np.any(diag == 0.0):
            raise ValueError("diag(K) contains zero entries")
        invdiag = 1.0 / diag
        z = invdiag[:, None] * Bt.toarray()
        S = D.toarray() - B @ z
        ksolve = lambda q: invdiag * q
    elif name == "ilu":
        action = spilu(K.tocsc(), drop_tol=ilu_drop, fill_factor=ilu_fill,
                       permc_spec="COLAMD")
        z = action.solve(Bt.toarray())
        S = D.toarray() - B @ z
        ksolve = action.solve
    else:
        raise ValueError(f"unknown strategy: {name}")
    setup_seconds = time.perf_counter() - setup_start

    solve_start = time.perf_counter()
    rhs_s = g - B @ ksolve(f)
    lam = np.linalg.solve(S, rhs_s)
    u = ksolve(f - Bt @ lam)
    solve_seconds = time.perf_counter() - solve_start

    x = np.concatenate((u, lam))
    residual = np.asarray(A @ x - rhs).reshape(-1)
    full_abs = float(np.linalg.norm(residual))
    full_rel = full_abs / max(float(np.linalg.norm(rhs)), np.finfo(float).tiny)
    backward = full_abs / max(
        float(np.linalg.norm(A.data)) * float(np.linalg.norm(x))
        + float(np.linalg.norm(rhs)), np.finfo(float).tiny
    )
    constraint = B @ u + D @ lam - g
    constraint_abs = float(np.linalg.norm(constraint))
    constraint_scale = max(float(np.linalg.norm(B @ u)),
                           float(np.linalg.norm(D @ lam)),
                           float(np.linalg.norm(g)))
    s_bytes = int(S.nbytes)
    return {
        "strategy": name,
        "interface_contract": {
            "apply_K_inverse": "q -> K^{-1}q",
            "apply_schur": "x -> D*x - B*(K^{-1}*(Bt*x))",
            "factorization": "setup once, then apply repeatedly",
            "explicit_schur": "optional materialized D - B*K^{-1}*Bt for this diagnostic only",
        },
        "setup_seconds": setup_seconds,
        "solve_seconds": solve_seconds,
        # The factorization is intentionally built once, then reused for all
        # columns of B^T and the final primal action.  This is the prototype
        # reuse boundary for a production explicit/frozen Schur candidate.
        "reuse_scope": "same matrix / same timestep",
        "reuse_boundaries": {
            "reuse_within_timestep": True,
            "rebuild_on_matrix_or_constraint_pattern_change": True,
            "rebuild_on_nonlinear_or_timestep_assembly": "required unless matrix fingerprint matches",
        },
        "cache_key_fields": ["matrix_fingerprint", "constraint_pattern", "backend", "preconditioner_settings"],
        "factorization_reused_for_schur_rhs": True,
        "schur_shape": list(S.shape),
        "schur_nnz": int(np.count_nonzero(S)),
        "schur_dense_bytes": s_bytes,
        "full_absolute_residual": full_abs,
        "full_relative_residual": float(full_rel),
        "backward_error": float(backward),
        "backward_error_norm_definition": "normwise ||A*x-b||_2 / (||A||_F*||x||_2 + ||b||_2)",
        "absolute_constraint_residual": constraint_abs,
        "relative_constraint_residual": (
            None if constraint_scale <= 1.0e-14 else float(constraint_abs / constraint_scale)
        ),
        "constraint_rhs_norm": float(np.linalg.norm(g)),
        "constraint_residual_scale": constraint_scale,
        "relative_primal_agreement_with_mumps": float(
            np.linalg.norm(u - mumps) / max(float(np.linalg.norm(mumps)),
                                             np.finfo(float).tiny)),
        "primal_norm": float(np.linalg.norm(u)),
        "multiplier_norm": float(np.linalg.norm(lam)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--rhs", type=Path, required=True)
    parser.add_argument("--mumps-sol", type=Path, required=True)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--c-start", type=int, required=True,
                        help="one-based first constraint row")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--strategies", default="exact,diag,ilu")
    parser.add_argument("--ilu-drop", type=float, default=1.0e-4)
    parser.add_argument("--ilu-fill", type=float, default=10.0)
    args = parser.parse_args()

    rows, cols, vals = triplets(args.matrix)
    n_f = args.c_start - 1
    if not (0 < n_f < args.rows):
        raise ValueError("invalid constraint split")
    A = coo_matrix((vals, (rows, cols)), shape=(args.rows, args.rows)).tocsr()
    K = A[:n_f, :n_f].tocsr()
    Bt = A[:n_f, n_f:].tocsr()
    B = A[n_f:, :n_f].tocsr()
    D = A[n_f:, n_f:].tocsr()
    rhs = np.loadtxt(args.rhs, dtype=np.float64)
    if rhs.ndim == 2:
        rhs = rhs[:, -1]
    rhs = rhs[:args.rows]
    f, g = rhs[:n_f], rhs[n_f:]
    mumps = np.loadtxt(args.mumps_sol, dtype=np.float64)
    if mumps.ndim == 2:
        mumps = mumps[:, -1]
    mumps = mumps[:n_f]

    result = {
        "n_f": n_f,
        "n_c": args.rows - n_f,
        "nnz": {"A": int(A.nnz), "K": int(K.nnz), "B": int(B.nnz),
                "Bt": int(Bt.nnz), "D": int(D.nnz)},
        "strategies": [
            solve_with_action(name.strip(), K, B, Bt, D, f, g, A, rhs,
                              mumps, args.ilu_drop, args.ilu_fill)
            for name in args.strategies.split(",") if name.strip()
        ],
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    args.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
