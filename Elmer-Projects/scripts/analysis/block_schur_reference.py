"""One-rank dense-Schur reference for an explicit mortar matrix.

This diagnostic is intentionally independent of HYPRE/MGR.  It factors the
primal K block with SuperLU, forms the 2,898-by-2,898 Schur matrix, and compares
the reconstructed primal field with the saved MUMPS temperature vector.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import splu


def triplets(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.loadtxt(path, dtype=np.float64)
    if d.ndim == 1:
        d = d.reshape(1, -1)
    return d[:, 0].astype(np.int64) - 1, d[:, 1].astype(np.int64) - 1, d[:, 2]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--rhs", type=Path, required=True)
    p.add_argument("--mumps-sol", type=Path, required=True)
    p.add_argument("--rows", type=int, required=True)
    p.add_argument("--c-start", type=int, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()

    rows, cols, vals = triplets(a.matrix)
    n_f = a.c_start - 1
    n_c = a.rows - n_f
    A = coo_matrix((vals, (rows, cols)), shape=(a.rows, a.rows)).tocsr()
    K = A[:n_f, :n_f].tocsr()
    Bt = A[:n_f, n_f:].tocsr()
    B = A[n_f:, :n_f].tocsr()
    D = A[n_f:, n_f:].tocsr()
    rhs = np.loadtxt(a.rhs, dtype=np.float64)
    if rhs.ndim == 2:
        rhs = rhs[:, -1]
    rhs = rhs[: a.rows]
    f, g = rhs[:n_f], rhs[n_f:]
    m = np.loadtxt(a.mumps_sol, dtype=np.float64)
    if m.ndim == 2:
        m = m[:, -1]
    m = m[:n_f]

    # SuperLU is used only to establish an independent one-rank reference;
    # no HYPRE/MGR object participates in this calculation.
    lu = splu(K.tocsc(), permc_spec="COLAMD")
    y = lu.solve(f)
    Z = lu.solve(Bt.toarray())
    S = D.toarray() - B @ Z
    rhs_s = g - B @ y
    lam = np.linalg.solve(S, rhs_s)
    u = y - Z @ lam
    x = np.concatenate((u, lam))
    ax_b = A @ x - rhs
    saddle_abs = float(np.linalg.norm(ax_b))
    saddle_res = saddle_abs / np.linalg.norm(rhs)
    constraint_abs = float(np.linalg.norm(B @ u + D @ lam - g))
    g_norm = float(np.linalg.norm(g))
    constraint_res = None if g_norm == 0.0 else constraint_abs / g_norm
    agreement = np.linalg.norm(u - m) / np.linalg.norm(m)
    result = {
        "n_f": n_f,
        "n_c": n_c,
        "nnz_K": int(K.nnz),
        "nnz_S": int(np.count_nonzero(np.abs(S) > 0.0)),
        "norm_S": float(np.linalg.norm(S)),
        "absolute_residual_Ax_minus_b": saddle_abs,
        "relative_residual_Ax_minus_b": float(saddle_res),
        "absolute_constraint_residual": constraint_abs,
        "relative_constraint_residual": (None if constraint_res is None
                                          else float(constraint_res)),
        "relative_primal_agreement_with_mumps": float(agreement),
        "mumps_primal_norm": float(np.linalg.norm(m)),
        "block_primal_norm": float(np.linalg.norm(u)),
    }
    encoded = json.dumps(result, indent=2)
    print(encoded)
    a.output.write_text(encoded + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
