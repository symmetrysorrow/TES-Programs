"""Element-level conductance reassembly, verified against the native matrix.

The native capture stores only the assembled operator.  To split the exact
discrete dissipation ``T^T K T = P (T_src - T_bath)`` into non-negative
element contributions we recompute the element conductance matrices
``K_e = int k grad N grad N^T`` for the linear tetra (504) and wedge (706)
elements with the SIF conductivities, then *require* that their assembly
reproduces the native ``K`` on every non-Dirichlet primal row.  Only when that
verification passes are the element dissipations used.
"""
from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp

from scripts.support.phase24_native_operator_lib import Mesh

T_INITIAL = 0.16857  # SIF Initial Condition 1; one nonlinear iteration => k(T_initial)


def membrane_k(t: float = T_INITIAL) -> float:
    tx = 0.5 * ((t + 1.0e-12) + abs(t - 1.0e-12))
    return 7.854e-06 * tx ** (4.252 - 1) * 0.4 * (0.0007 - 0.0005) / (8 * (1e-06 + 1.5e-05 + 1e-06) * 0.0005)


CONDUCTIVITY = {
    "abs": 0.0168, "TES": 68.0, "Stycast": 2.69094e-06, "SiO2_1": 0.000494, "SiO2_2": 0.000494,
    "Si_1": 0.37, "Si_2": 0.37, "SiNx": 0.000102, "Membrane_SiNx": None, "Membrane_Si1": None,
}


def conductivity(name: str) -> float:
    value = CONDUCTIVITY[name]
    return membrane_k() if value is None else value


def _tet(p: np.ndarray) -> np.ndarray:
    m = np.ones((4, 4))
    m[:, 1:] = p
    inv = np.linalg.inv(m)
    grad = inv[1:, :].T              # (4, 3)
    vol = abs(np.linalg.det(m)) / 6.0
    return vol * grad @ grad.T


_G3 = [(1 / 6, 1 / 6), (2 / 3, 1 / 6), (1 / 6, 2 / 3)]
_G2 = [-1 / math.sqrt(3), 1 / math.sqrt(3)]


def _wedge(p: np.ndarray) -> np.ndarray:
    k = np.zeros((6, 6))
    for (u, v) in _G3:
        for z in _G2:
            w = 1 / 6 * 1.0
            L = np.array([1 - u - v, u, v])
            dL = np.array([[-1, -1], [1, 0], [0, 1]], dtype=float)
            lo, hi = 0.5 * (1 - z), 0.5 * (1 + z)
            dN = np.zeros((6, 3))
            dN[:3, :2] = dL * lo
            dN[3:, :2] = dL * hi
            dN[:3, 2] = -0.5 * L
            dN[3:, 2] = 0.5 * L
            J = dN.T @ p
            detJ = np.linalg.det(J)
            g = np.linalg.solve(J, dN.T).T
            k += w * abs(detJ) * g @ g.T
    return k


def element_matrices(mesh: Mesh) -> list[np.ndarray]:
    out = []
    for e, nodes in enumerate(mesh.elem_nodes):
        kval = conductivity(mesh.body_names[int(mesh.elem_body[e])])
        p = mesh.coords[nodes]
        code = int(mesh.elem_type[e])
        if code == 504:
            out.append(kval * _tet(p))
        elif code == 706:
            out.append(kval * _wedge(p))
        else:
            raise NotImplementedError(code)
    return out


def assemble(mesh: Mesh, kes: list[np.ndarray], n: int) -> sp.csr_matrix:
    rows, cols, vals = [], [], []
    for nodes, ke in zip(mesh.elem_nodes, kes):
        d = mesh.node_to_dof[nodes]
        rows.append(np.repeat(d, len(d)))
        cols.append(np.tile(d, len(d)))
        vals.append(ke.ravel())
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


def verify(mesh: Mesh, kes: list[np.ndarray], K_native: sp.csr_matrix, dirichlet: np.ndarray) -> dict:
    n = K_native.shape[0]
    Kp = assemble(mesh, kes, n)
    mask = np.ones(n, dtype=bool)
    mask[dirichlet] = False
    D = sp.diags(mask.astype(float))
    diff = (D @ (Kp - K_native)).tocoo()
    ref = (D @ K_native).tocoo()
    scale = np.abs(K_native.diagonal())
    rel = np.abs(diff.data) / np.maximum(scale[diff.row], 1e-300) if diff.nnz else np.zeros(1)
    return {
        "max_abs_entry_diff": float(np.abs(diff.data).max()) if diff.nnz else 0.0,
        "max_entry_diff_relative_to_row_diagonal": float(rel.max()),
        "native_nnz_checked": int(ref.nnz),
        "passed": bool(rel.max() < 1e-8),
    }


def dissipation(mesh: Mesh, kes: list[np.ndarray], T: np.ndarray, t_ref: float) -> np.ndarray:
    """Per-element T_e^T K_e T_e (W*K).  Shift-invariant since K_e 1 = 0."""
    out = np.empty(len(kes))
    for e, (nodes, ke) in enumerate(zip(mesh.elem_nodes, kes)):
        t = T[mesh.node_to_dof[nodes]] - t_ref
        out[e] = float(t @ ke @ t)
    return out
