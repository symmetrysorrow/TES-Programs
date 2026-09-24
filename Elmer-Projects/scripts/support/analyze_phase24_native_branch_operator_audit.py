"""Phase24 native branch heat-flow and trace-to-flux operator audit.

Diagnostic only.  Meshes, SIFs, materials, TES law, circuit and solver are
untouched; no ElmerSolver run is needed.  The audit re-reads the existing
opt-in native captures (the exact MUMPS saddle-point system
``[[K, C^T], [C, 0]] [T; lambda] = [f; 0]`` and its solution) of

* the historical reference         (``mesh_singlepixel_prod_v2``),
* the best Phase24 control         (``mesh_phase24_trace_membrane_substrate_historical``),
* the refined mortar parent        (``mesh_phase24_stycast_density_10um``),

at 0.95/1.00/1.05 P0 (frozen power, CPU native HeatSolve MUMPS).

Heat-flow sources (no area splitting, no gradient proxy, no series summation):

* cut flow  = assembled native edge flow ``q_ij = K_ij (T_j - T_i)`` across a
  closed node cut, plus ``C_pair^T lambda`` for mortar-coupled pairs;
* bath flow = native edge flow into the Dirichlet rows;
* branch conductance of the parallel suspended-membrane sheets = exact element
  dissipation ``T_e^T K_e T_e / dT^2`` from element matrices whose assembly is
  verified entry-by-entry against the native ``K``.

Operators are condensed natively on sparse local subdomains (no dense global
matrix): ``K_cond = K_ff - K_fi K_ii^-1 K_if`` (with mortar rows kept as a
local saddle system where present).
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as sla
import scipy.sparse as sp
import scipy.sparse.linalg as spla

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.support.phase24_element_dissipation import (  # noqa: E402
    assemble, conductivity, dissipation, element_matrices, verify,
)
from scripts.support.phase24_native_operator_lib import (  # noqa: E402
    CASES, FRACTIONS, P0, TBATH, Capture, Mesh, body_mean_temperature_like_harness,
    dirichlet_rows, load_capture, load_mesh, sha256, tag,
)

OUT = ROOT / "artifacts/phase24_native_branch_operator_audit"
UM = 1e-6
ZTOL = 1e-3 * UM
STOL = 0.5 * UM
YC = 1000 * UM
A_TES = 250 * UM          # TES half width
A_WIN = 350 * UM          # suspended membrane window half width
GHIST_REPORTED = 1.8782944616314638e-8
GBEST_REPORTED = 1.959713948e-8
TOP = ("TES", "Stycast", "abs")
SHEETS = ("Membrane_SiNx", "Membrane_Si1", "SiO2_1")
ZONES = ((0, 230), (230, 250), (250, 270), (270, 330), (330, 350), (350, 370), (370, 1e9))

FACE_TABLE = {
    504: ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
    706: ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
}


# --------------------------------------------------------------------------
# small utilities


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        fields.extend(k for k in row if k not in fields)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (np.floating, float)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    return value


def square_radius(xy: np.ndarray) -> np.ndarray:
    """Chebyshev distance from the TES/window centre (the geometry is square)."""
    return np.maximum(np.abs(xy[:, 0]), np.abs(xy[:, 1] - YC))


def body_faces(mesh: Mesh, body: int) -> list[np.ndarray]:
    faces = []
    for e in np.where(mesh.elem_body == body)[0]:
        nodes = mesh.elem_nodes[e]
        for local in FACE_TABLE[int(mesh.elem_type[e])]:
            faces.append(nodes[list(local)])
    return faces


def face_area(mesh: Mesh, face: np.ndarray) -> float:
    p = mesh.coords[face]
    area = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
    if len(face) == 4:
        area += 0.5 * np.linalg.norm(np.cross(p[3] - p[0], p[2] - p[0]))
    return float(area)


def plane_faces(mesh: Mesh, body: int, z0: float, smax: float | None = None) -> list[np.ndarray]:
    out = {}
    for face in body_faces(mesh, body):
        p = mesh.coords[face]
        if np.all(np.abs(p[:, 2] - z0) < ZTOL):
            if smax is not None and square_radius(p.mean(axis=0, keepdims=True))[0] > smax + STOL:
                continue
            out[tuple(sorted(face.tolist()))] = face
    return list(out.values())


def lumped(mesh: Mesh, faces: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    weights: dict[int, float] = {}
    for face in faces:
        area = face_area(mesh, face)
        for node in face:
            weights[int(node)] = weights.get(int(node), 0.0) + area / len(face)
    keys = np.array(sorted(weights), dtype=np.int64)
    return keys, np.array([weights[k] for k in keys])


def aw(mesh: Mesh, T: np.ndarray, faces: list[np.ndarray]) -> float:
    if not faces:
        return float("nan")
    nodes, w = lumped(mesh, faces)
    return float(np.dot(w, T[mesh.dofs(nodes)]) / w.sum())


def node_xyz_by_dof(mesh: Mesh, n: int) -> np.ndarray:
    xyz = np.full((n, 3), np.nan)
    has = mesh.node_to_dof >= 0
    xyz[mesh.node_to_dof[has]] = mesh.coords[has]
    return xyz


def body_dofs(mesh: Mesh, name: str) -> np.ndarray:
    return mesh.dofs(mesh.body_node_sets[mesh.name_to_body(name)])


def mortar_pairs(mesh: Mesh, cap: Capture) -> dict[tuple[str, ...], np.ndarray]:
    """Constraint rows grouped by the set of bodies their columns touch."""
    n = cap.n_primal
    label = np.zeros(n, dtype=np.int64)
    names = sorted(mesh.body_names)
    for bit, body in enumerate(names):
        label[mesh.dofs(mesh.body_node_sets[body])] |= 1 << bit
    C = cap.C.tocsr()
    rows: dict[int, list[int]] = {}
    for r in range(C.shape[0]):
        cols = C.indices[C.indptr[r]:C.indptr[r + 1]]
        mask = int(np.bitwise_or.reduce(label[cols])) if len(cols) else 0
        rows.setdefault(mask, []).append(r)
    out = {}
    for mask, rr in rows.items():
        key = tuple(sorted(mesh.body_names[b] for bit, b in enumerate(names) if mask >> bit & 1))
        out[key] = np.array(rr)
    return out


# --------------------------------------------------------------------------
# native cut flows


def cut_flow(cap: Capture, dirichlet_mask: np.ndarray, inside: np.ndarray, mortar_m: np.ndarray | None = None) -> tuple[float, np.ndarray]:
    """Native flow out of node set ``inside`` (bool over primal DOFs).

    Returns total and the per-edge (row-node attributed) outflow vector.  Only
    non-Dirichlet rows are used; mortar reactions ``m = C^T lambda`` on inside
    nodes are added when given.
    """
    K = cap.K.tocoo()
    T = cap.T
    sel = inside[K.row] & ~inside[K.col] & ~dirichlet_mask[K.row]
    q = K.data[sel] * (T[K.col[sel]] - T[K.row[sel]])
    per = np.zeros(cap.n_primal)
    np.add.at(per, K.row[sel], q)
    total = float(q.sum())
    if mortar_m is not None:
        total += float(mortar_m[inside].sum())
    return total, per


def layer_attributed_cut(cap: Capture, dmask: np.ndarray, inside: np.ndarray, xyz: np.ndarray) -> dict[str, float]:
    """Split a closed lateral cut by the vertical position of each crossing edge."""
    K = cap.K.tocoo()
    T = cap.T
    sel = inside[K.row] & ~inside[K.col] & ~dmask[K.row]
    r, c = K.row[sel], K.col[sel]
    q = K.data[sel] * (T[c] - T[r])
    zr, zc = xyz[r, 2] / UM, xyz[c, 2] / UM
    zmid = 0.5 * (zr + zc)
    lo, hi = np.minimum(zr, zc), np.maximum(zr, zc)
    out = {
        "Membrane_SiNx_layer(z191-192)": 0.0, "Membrane_Si1_layer(z176-191)": 0.0,
        "SiO2_1_layer(z175-176)": 0.0, "plane_z191_junction": 0.0, "plane_z176_junction": 0.0,
        "plane_z175_junction": 0.0, "plane_z192_surface": 0.0, "other": 0.0,
    }
    for value, a, b, m in zip(q, lo, hi, zmid):
        if abs(a - b) < 1e-3:
            key = {191: "plane_z191_junction", 176: "plane_z176_junction", 175: "plane_z175_junction", 192: "plane_z192_surface"}.get(int(round(a)), "other")
        elif 191 - 1e-3 <= m <= 192 + 1e-3:
            key = "Membrane_SiNx_layer(z191-192)"
        elif 176 - 1e-3 <= m <= 191 + 1e-3:
            key = "Membrane_Si1_layer(z176-191)"
        elif 175 - 1e-3 <= m <= 176 + 1e-3:
            key = "SiO2_1_layer(z175-176)"
        else:
            key = "other"
        out[key] += float(value)
    out["total"] = float(q.sum())
    out["negative_edge_flow_sum"] = float(q[q < 0].sum())
    return out


# --------------------------------------------------------------------------
# condensed operators


class LocalOperator:
    """Condensed trace-to-flux operator of a sparse local subdomain.

    ``K`` is the subdomain conductance matrix (global DOF numbering), ``f`` the
    trace DOFs (prescribed), ``g`` grounded DOFs (T=0), everything else in the
    subdomain free.  Optional mortar rows ``C`` (rows x global DOFs) are kept as
    a saddle block so hanging/mortar traces are condensed exactly.
    """

    def __init__(self, K: sp.csr_matrix, f: np.ndarray, g: np.ndarray, domain: np.ndarray, C: sp.csr_matrix | None = None):
        n = K.shape[0]
        f = np.asarray(f)
        self.dropped_constraint_rows = 0
        if C is not None and C.shape[0]:
            # Constraint rows whose columns are all prescribed trace DOFs carry
            # no information for the condensed operator (their multiplier is
            # undetermined); hanging trace DOFs coupled only through such rows
            # have no path into the subdomain and are removed from the trace.
            fset = np.zeros(n, dtype=bool); fset[f] = True
            Cc = C.tocsr()
            informative = np.array([(~fset[Cc.indices[Cc.indptr[r]:Cc.indptr[r + 1]]]).any() for r in range(Cc.shape[0])])
            self.dropped_constraint_rows = int((~informative).sum())
            C = Cc[informative]
        Kr = K.tocsr()
        touched = np.asarray(Kr[f].getnnz(axis=1)) > 0
        if C is not None and C.shape[0]:
            touched |= np.asarray(C.tocsc()[:, f].getnnz(axis=0)).ravel() > 0
        self.keep = touched
        self.f = f[touched]
        in_dom = np.zeros(n, dtype=bool)
        in_dom[domain] = True
        in_dom[self.f] = True
        fixed = np.zeros(n, dtype=bool)
        fixed[self.f] = True
        fixed[g] = True
        free = np.where(in_dom & ~fixed)[0]
        self.free = free
        Kc = K.tocsc()
        self.K_ii = Kc[free][:, free].tocsc()
        self.K_if = Kc[free][:, self.f].tocsc()
        self.K_ff = Kc[self.f][:, self.f].tocsc()
        self.K_fi = Kc[self.f][:, free].tocsc()
        self.nc = 0
        if C is not None and C.shape[0]:
            C = C.tocsc()
            self.C_i = C[:, free]
            self.C_f = C[:, self.f]
            self.nc = C.shape[0]
            S = sp.bmat([[self.K_ii, self.C_i.T], [self.C_i, None]], format="csc")
        else:
            S = self.K_ii
        self.lu = spla.splu(S, permc_spec="COLAMD")

    def apply(self, t: np.ndarray) -> np.ndarray:
        t = np.atleast_2d(t.T).T if t.ndim == 1 else t
        rhs = -(self.K_if @ t)
        if self.nc:
            rhs = np.vstack([rhs, -(self.C_f @ t)])
        sol = self.lu.solve(np.asarray(rhs))
        u = sol[: len(self.free)]
        r = self.K_ff @ t + self.K_fi @ u
        if self.nc:
            r = r + self.C_f.T @ sol[len(self.free):]
        return np.asarray(r)

    def dense(self, block: int = 256) -> np.ndarray:
        m = len(self.f)
        out = np.zeros((m, m))
        for s in range(0, m, block):
            e = min(m, s + block)
            E = np.zeros((m, e - s))
            E[np.arange(s, e), np.arange(e - s)] = 1.0
            out[:, s:e] = self.apply(E)
        return 0.5 * (out + out.T)


def modes(xyz_f: np.ndarray) -> dict[str, np.ndarray]:
    xi = xyz_f[:, 0] / A_TES
    eta = (xyz_f[:, 1] - YC) / A_TES
    return {"constant": np.ones(len(xi)), "linear_x": xi, "linear_y": eta, "radial_quadratic": xi ** 2 + eta ** 2}


def operator_report(name: str, op: LocalOperator, xyz_f: np.ndarray, area_f: np.ndarray, dense: bool = True) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    xyz_f, area_f = xyz_f[op.keep], area_f[op.keep]
    rep: dict[str, Any] = {"operator": name, "trace_dof_count": int(len(op.f)), "interior_dof_count": int(len(op.free)),
                           "constraint_rows": int(op.nc), "trace_area_m2": float(area_f.sum()),
                           "dropped_hanging_trace_dofs": int((~op.keep).sum()), "dropped_uninformative_constraint_rows": op.dropped_constraint_rows}
    rows = []
    s = square_radius(xyz_f) / UM
    for mode, t in modes(xyz_f).items():
        r = op.apply(t).ravel()
        energy = float(t @ r)
        l2 = float(np.dot(area_f, t * t))
        rep[f"mode_{mode}"] = {
            "total_reaction_W_per_K": float(r.sum()),
            "energy_W_per_K": energy,
            "normalized_response_W_per_K_m2": energy / l2 if l2 else float("nan"),
        }
        edges = [0, 50, 100, 150, 200, 225, 240, 250.0001]
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (s >= lo) & (s < hi)
            rows.append({"operator": name, "mode": mode, "s_bin_um": f"{lo:g}-{min(hi, 250):g}", "trace_nodes": int(m.sum()),
                         "input_mean": float(np.average(t[m], weights=area_f[m])) if m.any() and area_f[m].sum() > 0 else float("nan"),
                         "reaction_W_per_K": float(r[m].sum()),
                         "reaction_density_W_per_K_m2": float(r[m].sum() / area_f[m].sum()) if area_f[m].sum() > 0 else float("nan"),
                         "area_m2": float(area_f[m].sum())})
    if dense:
        Kc = op.dense()
        rs = Kc.sum(axis=1)
        pos = area_f > 0
        rep["row_sum_total_W_per_K"] = float(rs.sum())
        rep["row_sum_per_area_percentiles_W_per_K_m2"] = np.percentile(rs[pos] / area_f[pos], [1, 10, 50, 90, 99]).tolist()
        rep["diag_per_area_percentiles_W_per_K_m2"] = np.percentile(np.diag(Kc)[pos] / area_f[pos], [1, 10, 50, 90, 99]).tolist()
        rep["row_sum_zero_area_trace_dofs"] = {"count": int((~pos).sum()), "row_sum_total": float(rs[~pos].sum())}
        # generalized spectrum K v = mu M v with lumped trace area (hanging DOFs
        # without own area get a tiny mass so the pencil stays definite)
        M = np.where(pos, area_f, area_f[pos].min() * 1e-6)
        Mi = 1.0 / np.sqrt(M)
        mu = sla.eigvalsh(Mi[:, None] * Kc * Mi[None, :])
        rep["generalized_eigenvalues_W_per_K_m2"] = {"smallest_6": mu[:6].tolist(), "median": float(np.median(mu)), "largest_3": mu[-3:].tolist()}
        rep["effective_total_conductance_W_per_K"] = float(rs.sum())
    return rep, rows


# --------------------------------------------------------------------------
# 2-D reference for the lateral sheet shape factor


def annulus_ntd_shape_factor(h_um: float) -> float:
    """Same annulus, uniform heat density on |s|<=250um, outer square at 0.

    Returns S_ntd = 1 / (k t R) with R = (q . u) for unit total heat."""
    n = int(round(2 * 350 / h_um)) + 1
    xs = np.linspace(-350, 350, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    S = np.maximum(np.abs(X), np.abs(Y))
    outer = S >= 350 - 1e-9
    free = ~outer
    idx = -np.ones((n, n), dtype=np.int64)
    idx[free] = np.arange(free.sum())
    fi, fj = np.where(free)
    k = idx[fi, fj]
    rows, cols, vals = [k], [k], [np.full(len(k), 4.0)]
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nb = idx[fi + di, fj + dj]
        m = nb >= 0
        rows.append(k[m]); cols.append(nb[m]); vals.append(-np.ones(m.sum()))
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(free.sum(),) * 2)
    # cell-area weight of each node inside the closed footprint
    wx = np.clip(np.minimum(250 - np.abs(xs) + h_um / 2, h_um), 0, None) / h_um
    W = np.outer(wx, wx)[free]
    q = W / W.sum()
    u = spla.spsolve(A.tocsc(), q)
    return 1.0 / float(q @ u)


def annulus_shape_factor(h_um: float) -> float:
    """5-point FD shape factor S of the square annulus a_in=250, a_out=350 um.

    G_sheet(exact) = k * t * S for a sheet with isothermal inner square and
    isothermal outer square.  Quarter symmetry is not used (robustness)."""
    n = int(round(2 * 350 / h_um)) + 1
    xs = np.linspace(-350, 350, n)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    S = np.maximum(np.abs(X), np.abs(Y))
    inner = S <= 250 + 1e-9
    outer = S >= 350 - 1e-9
    free = ~(inner | outer)
    idx = -np.ones((n, n), dtype=np.int64)
    idx[free] = np.arange(free.sum())
    rows, cols, vals = [], [], []
    rhs = np.zeros(free.sum())
    fi, fj = np.where(free)
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ni, nj = fi + di, fj + dj
        k = idx[fi, fj]
        rows.append(k); cols.append(k); vals.append(np.ones_like(k, dtype=float))
        nb = idx[ni, nj]
        m = nb >= 0
        rows.append(k[m]); cols.append(nb[m]); vals.append(-np.ones(m.sum()))
        rhs[k[inner[ni, nj]]] += 1.0
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(free.sum(),) * 2)
    u = spla.spsolve(A.tocsc(), rhs)
    U = np.zeros((n, n)); U[inner] = 1.0; U[free] = u
    # flux out of the inner square through the first ring of links
    flux = 0.0
    ii, jj = np.where(inner)
    for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ni, nj = ii + di, jj + dj
        ok = (ni >= 0) & (ni < n) & (nj >= 0) & (nj < n)
        ni, nj, i0, j0 = ni[ok], nj[ok], ii[ok], jj[ok]
        m = ~inner[ni, nj]
        flux += float((U[i0[m], j0[m]] - U[ni[m], nj[m]]).sum())
    return flux


def native_ntd(cap: Capture, mesh: Mesh, dmask: np.ndarray, xyz: np.ndarray) -> dict[str, Any]:
    """Low-order TES source modes through the complete native saddle matrix."""
    n = cap.n_primal
    src = cap.b[:n].copy()
    src[dmask] = 0.0
    P = float(src.sum())
    lu = spla.splu(cap.A.tocsc(), permc_spec="COLAMD")
    nz = src != 0
    out: dict[str, Any] = {"source_dofs": int(nz.sum())}
    tes = body_dofs(mesh, "TES")
    for mode, w in modes(xyz[nz]).items():
        q = np.zeros(cap.A.shape[0])
        q[np.where(nz)[0]] = src[nz] * w / P
        y = lu.solve(q)[:n]
        out[f"mode_{mode}"] = {
            "energy_K_per_W": float(q[:n] @ y),
            "TES_mean_temperature_response_K_per_W": float(y[tes].mean()),
            "net_source_W": float(q[:n].sum()),
        }
    # consistency: constant mode reproduces the native solution (shifted)
    q = np.zeros(cap.A.shape[0]); q[:n] = src
    y = lu.solve(q)[:n]
    out["constant_mode_reproduces_native_relative_error"] = float(np.abs(y - (cap.T - TBATH)).max() / (cap.T - TBATH).max())
    return out


def stack_ntd(mesh: Mesh, K_net: sp.csr_matrix, dir_rows: np.ndarray, net_dofs: np.ndarray, membrane_body: int) -> dict[str, Any]:
    n = K_net.shape[0]
    faces = plane_faces(mesh, membrane_body, 192 * UM)
    rng = np.random.default_rng(1)
    bary = rng.dirichlet(np.ones(3), size=256)
    loads = {m: np.zeros(n) for m in ("constant", "linear_x", "linear_y", "radial_quadratic")}
    for face in faces:
        pts = mesh.coords[face]
        if len(face) != 3:
            raise NotImplementedError("quad top face")
        smp = bary @ pts
        inside = square_radius(smp) <= A_TES
        if not inside.any():
            continue
        area = 0.5 * np.linalg.norm(np.cross(pts[1] - pts[0], pts[2] - pts[0]))
        md = modes(smp[inside])
        d = mesh.node_to_dof[face]
        for m, w in md.items():
            np.add.at(loads[m], d, (bary[inside] * w[:, None]).sum(axis=0) * area / len(bary))
    total = loads["constant"].sum()
    mask = np.zeros(n, dtype=bool); mask[net_dofs] = True; mask[dir_rows] = False
    free = np.where(mask)[0]
    lu = spla.splu(K_net.tocsc()[free][:, free], permc_spec="COLAMD")
    out: dict[str, Any] = {"loaded_area_m2": float(total), "reference_area_m2": (2 * A_TES) ** 2}
    for m, q in loads.items():
        q = q / total
        u = lu.solve(q[free])
        out[f"mode_{m}"] = {"energy_K_per_W": float(q[free] @ u), "net_load": float(q.sum())}
    out["G_uniform_flux_W_per_K"] = 1.0 / out["mode_constant"]["energy_K_per_W"]
    return out


# --------------------------------------------------------------------------
# per-case analysis


def analyze_case(key: str) -> dict[str, Any]:
    spec = CASES[key]
    mesh = load_mesh(spec)
    caps = {f: load_capture(spec.capture(f)) for f in FRACTIONS}
    cap1 = caps[1.0]
    n = cap1.n_primal
    names = mesh.body_names
    xyz = node_xyz_by_dof(mesh, n)
    s_node = square_radius(xyz)
    kes = element_matrices(mesh)
    dir_rows = dirichlet_rows(cap1.K)
    dmask = np.zeros(n, dtype=bool); dmask[dir_rows] = True
    ver = verify(mesh, kes, cap1.K, dir_rows)
    pairs = mortar_pairs(mesh, cap1)
    top_mask = np.zeros(n, dtype=bool)
    for name in TOP:
        top_mask[body_dofs(mesh, name)] = True
    network_mask = np.zeros(n, dtype=bool)
    for b, nm in names.items():
        if nm not in TOP:
            network_mask[body_dofs(mesh, nm)] = True
    tes_mem_rows = pairs.get(("Membrane_SiNx", "TES"), np.zeros(0, dtype=np.int64))
    # element geometry
    cent = np.array([mesh.coords[nodes].mean(axis=0) for nodes in mesh.elem_nodes])
    s_elem = square_radius(cent)
    b = mesh.name_to_body
    faces = {
        "TES_bottom": plane_faces(mesh, b("TES"), 192 * UM),
        "Membrane_top_under_TES": plane_faces(mesh, b("Membrane_SiNx"), 192 * UM, A_TES),
        "MembraneSiNx_sidewall_to_SiNx": [f for f in body_faces(mesh, b("SiNx")) if np.all(np.abs(square_radius(mesh.coords[f]) - A_WIN) < STOL) and np.ptp(mesh.coords[f, 2]) > ZTOL],
        "MembraneSi1_sidewall_to_Si_1": [f for f in body_faces(mesh, b("Si_1")) if np.all(np.abs(square_radius(mesh.coords[f]) - A_WIN) < STOL) and np.ptp(mesh.coords[f, 2]) > ZTOL],
        "frame_top_SiNx": plane_faces(mesh, b("SiNx"), 192 * UM),
        "bath": plane_faces(mesh, b("SiO2_2"), -177 * UM),
    }

    def sub_assembly(elems: list[int]) -> sp.csr_matrix:
        rows, cols, vals = [], [], []
        for i in elems:
            d = mesh.node_to_dof[mesh.elem_nodes[i]]
            rows.append(np.repeat(d, len(d))); cols.append(np.tile(d, len(d))); vals.append(kes[i].ravel())
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))

    # parallel branches of the actual graph: the three suspended window sheets
    sheet_elems = {sheet: [i for i in range(len(kes)) if names[int(mesh.elem_body[i])] == sheet and s_elem[i] < A_WIN] for sheet in SHEETS}
    sheet_K = {sheet: sub_assembly(el) for sheet, el in sheet_elems.items()}
    sheet_dofs = {sheet: np.unique(np.concatenate([mesh.node_to_dof[mesh.elem_nodes[i]] for i in el])) for sheet, el in sheet_elems.items()}
    tes_dofs_mask = np.zeros(n, dtype=bool); tes_dofs_mask[body_dofs(mesh, "TES")] = True

    points: dict[str, Any] = {}
    for f, cap in caps.items():
        T = cap.T
        src = cap.b[:n].copy(); src[dmask] = 0.0
        P = float(src.sum())
        lam = cap.lam
        C = cap.C.tocsr()
        m_pair = {k: np.asarray(C[rows].T @ lam[rows]).ravel() for k, rows in pairs.items()}
        m_tm = m_pair.get(("Membrane_SiNx", "TES"), np.zeros(n))
        # TES -> Membrane: closed cut around the TES/Stycast/abs stack
        q_tm_total, _ = cut_flow(cap, dmask, top_mask, m_tm)
        q_tm_mortar = float(m_tm[top_mask].sum())
        # TES -> Stycast (dead-end branch): cut around Stycast+abs only
        sa = np.zeros(n, dtype=bool); sa[body_dofs(mesh, "Stycast")] = True; sa[body_dofs(mesh, "abs")] = True
        sa &= ~np.isin(np.arange(n), body_dofs(mesh, "TES"))
        m_ts = m_pair.get(("Stycast", "TES"), np.zeros(n))
        q_into_stycast, _ = cut_flow(cap, dmask, sa, m_ts)
        # bath: flow into Dirichlet rows from their neighbours
        q_bath, _ = cut_flow(cap, dmask, ~dmask)
        # lateral closed cuts through the suspended stack
        cuts = {}
        for label, radius in (("TES_footprint_edge", A_TES), ("window_perimeter", A_WIN)):
            for own in ("perimeter_inside", "perimeter_outside"):
                lim = radius + STOL if own == "perimeter_inside" else radius - STOL
                inside = top_mask | (network_mask & (s_node <= lim) & (xyz[:, 2] >= 175 * UM - ZTOL))
                inside &= ~dmask
                cuts[f"{label}|{own}"] = layer_attributed_cut(cap, dmask, inside, xyz)
        # element dissipation, exact partition of P (T_src - T_bath)
        D = dissipation(mesh, kes, T, TBATH)
        t_src = float(np.dot(src, T) / P)
        dT_src = t_src - TBATH
        t_tes = body_mean_temperature_like_harness(mesh, T, b("TES"))
        body_D = {names[int(bb)]: float(D[mesh.elem_body == bb].sum()) for bb in np.unique(mesh.elem_body)}
        zone_D = {}
        for name in names.values():
            sel = mesh.elem_body == b(name)
            for lo, hi in ZONES:
                m = sel & (s_elem >= lo * UM) & (s_elem < hi * UM)
                zone_D[f"{name}|{lo}-{hi if hi < 1e8 else 'inf'}"] = float(D[m].sum())
        # native-verified element reactions per branch sheet
        t_shift = T - TBATH
        master_mask = np.abs(m_tm) > 0
        branches = {}
        for sheet in SHEETS:
            r = np.asarray(sheet_K[sheet] @ t_shift).ravel()
            dofs = sheet_dofs[sheet]
            cat = np.full(n, "", dtype=object)
            cat[dofs] = "residual"
            for other in SHEETS:
                if other != sheet:
                    shared = np.intersect1d(dofs, sheet_dofs[other])
                    cat[shared] = f"exchange_from_{other}"
            tes_if = dofs[tes_dofs_mask[dofs] | master_mask[dofs]]
            cat[tes_if] = "TES_interface"
            cat[dofs[s_node[dofs] >= A_WIN - STOL]] = "perimeter"
            row = {"Q_in_from_TES_W": float(r[cat == "TES_interface"].sum()),
                   "Q_exit_perimeter_W": float(-r[cat == "perimeter"].sum())}
            for other in SHEETS:
                if other != sheet:
                    row[f"Q_in_from_{other}_W"] = float(r[cat == f"exchange_from_{other}"].sum())
            row["balance_residual_W"] = float(r[cat == "residual"].sum())
            hot = dofs[s_node[dofs] <= A_TES + STOL]
            cold = dofs[s_node[dofs] >= A_WIN - STOL]
            row["T_hot_footprint_mean_K"] = float(T[hot].mean())
            row["T_cold_perimeter_mean_K"] = float(T[cold].mean())
            row["dT_branch_K"] = row["T_hot_footprint_mean_K"] - row["T_cold_perimeter_mean_K"]
            row["dissipation_WK"] = float(sum(D[i] for i in sheet_elems[sheet]))
            branches[sheet] = row
        traces = {k: aw(mesh, T, v) for k, v in faces.items()}
        traces["TES_average_harness"] = t_tes
        traces["TES_source_weighted"] = t_src
        traces["bath"] = TBATH
        resid = cap.A @ cap.x - cap.b
        points[tag(f)] = {
            "P_source_W": P, "P_nominal_W": P0 * f,
            "Q_TES_to_Membrane_W": q_tm_total, "Q_TES_to_Membrane_mortar_part_W": q_tm_mortar,
            "Q_TES_to_Membrane_conformal_part_W": q_tm_total - q_tm_mortar,
            "Q_into_Stycast_abs_W": q_into_stycast, "Q_bath_W": q_bath,
            "lateral_cuts_edge_attribution": cuts,
            "branch_sheets": branches,
            "dissipation_total_WK": float(D.sum()), "P_dT_src_WK": P * dT_src,
            "dissipation_closure_rel": float(D.sum() / (P * dT_src) - 1.0),
            "body_dissipation_WK": body_D, "zone_dissipation_WK": zone_D,
            "traces_K": traces,
            "G_eff_secant_harness_W_per_K": P / (t_tes - TBATH),
            "G_src_W_per_K": P / dT_src,
            "native_residual_primal_inf": float(np.abs(resid[:n]).max()),
            "native_residual_constraint_inf": float(np.abs(resid[n:]).max()) if len(resid) > n else 0.0,
            "native_balance_rel": float((q_bath - P) / P),
        }
    tlo, thi = (body_mean_temperature_like_harness(mesh, caps[x].T, b("TES")) for x in (0.95, 1.05))
    g_eff = 0.10 * P0 / (thi - tlo)

    # ---------------- operators (P independent: K identical at all points)
    net_elems = [i for i in range(len(kes)) if names[int(mesh.elem_body[i])] not in TOP]

    K_net = sub_assembly(net_elems)
    C_tm = cap1.C.tocsr()[tes_mem_rows] if len(tes_mem_rows) else None
    ops: dict[str, Any] = {}
    mode_rows: list[dict[str, Any]] = []
    net_dofs = np.where(network_mask)[0]
    t_actual = cap1.T - TBATH
    P1 = float(cap1.b[:n][~dmask].sum())
    # Membrane-side trace of the TES/Membrane coupling: membrane top nodes under
    # the TES footprint plus every mortar master DOF (historical mortar masters
    # extend to s~264um and carry ~35% of the heat).
    mem_nodes, mem_area = lumped(mesh, faces["Membrane_top_under_TES"])
    f_mem = mesh.dofs(mem_nodes)
    area_by_dof = dict(zip(f_mem.tolist(), mem_area.tolist()))
    if C_tm is not None:
        masters = np.unique(C_tm.indices)
        masters = masters[network_mask[masters]]
        f_mem = np.unique(np.concatenate([f_mem, masters]))
    area_mem = np.array([area_by_dof.get(int(d), 0.0) for d in f_mem])
    # (1) network DtN seen from the membrane-side trace (coupling excluded)
    op = LocalOperator(K_net, f_mem, dir_rows, net_dofs, None)
    ops["Lambda_membrane_trace_to_bath"], rows = operator_report("Lambda_membrane_trace_to_bath", op, xyz[f_mem], area_mem)
    ops["Lambda_membrane_trace_to_bath"]["actual_trace_reproduction_Q_over_P"] = float(op.apply(t_actual[op.f]).sum() / P1)
    mode_rows += rows
    # (1b) same network, but loaded by a *uniform heat flux* over the TES
    #      footprint on the Membrane_SiNx top face (no isothermal TES edge, so
    #      no edge singularity; face-sampled consistent loads, snapping-free).
    ops["Stack_NtD_uniform_footprint_flux"] = stack_ntd(mesh, K_net, dir_rows, net_dofs, b("Membrane_SiNx"))
    # (2) TES/Membrane local operator on the same trace: Membrane_SiNx layer only,
    #     grounded where it meets Membrane_Si1 / SiNx / Si_1
    mem_elems = [i for i in range(len(kes)) if names[int(mesh.elem_body[i])] == "Membrane_SiNx"]
    K_mem = sub_assembly(mem_elems)
    mem_dofs = body_dofs(mesh, "Membrane_SiNx")
    others = np.zeros(n, dtype=bool)
    for nm in ("Membrane_Si1", "SiNx", "Si_1"):
        others[body_dofs(mesh, nm)] = True
    ground = mem_dofs[others[mem_dofs]]
    op = LocalOperator(K_mem, f_mem, ground, mem_dofs, None)
    ops["TES_Membrane_local_operator"], rows = operator_report("TES_Membrane_local_operator", op, xyz[f_mem], area_mem)
    mode_rows += rows
    # (3) TES-side response of the complete native saddle operator (NtD):
    #     native source pattern modulated by the low-order modes; well posed
    #     for mortar and conformal couplings alike.
    ops["TES_side_native_NtD"] = native_ntd(cap1, mesh, dmask, xyz)
    # slave-side Dirichlet well-posedness of the mortar coupling
    if C_tm is not None:
        Ci = C_tm[:, np.where(~top_mask)[0]]
        Cd = Ci.toarray()
        Cd = Cd[:, np.abs(Cd).sum(axis=0) > 0]
        sv = np.linalg.svd(Cd, compute_uv=False)
        ops["TES_Membrane_mortar_slave_Dirichlet_wellposedness"] = {
            "constraint_rows": int(C_tm.shape[0]), "master_dofs": int(Cd.shape[1]),
            "sigma_min_over_max": float(sv.min() / sv.max()),
            "near_null_count_rel_1e-6": int((sv < 1e-6 * sv.max()).sum()),
            "conclusion": "slave-side (TES) Dirichlet trace-to-flux map is ill-posed" if sv.min() / sv.max() < 1e-8 else "well posed",
        }
    # (4) per-sheet lateral annulus operators (parallel branch operators)
    sheet_ops = {}
    for sheet in SHEETS:
        elems = [i for i in range(len(kes)) if names[int(mesh.elem_body[i])] == sheet and s_elem[i] < A_WIN]
        Ks = sub_assembly(elems)
        dofs = np.unique(np.concatenate([mesh.node_to_dof[mesh.elem_nodes[i]] for i in elems]))
        ds = s_node[dofs]
        hot = dofs[ds <= A_TES + STOL]
        cold = dofs[ds >= A_WIN - STOL]
        op = LocalOperator(Ks, hot, cold, dofs, None)
        G = float(op.apply(np.ones(len(hot))).sum())
        kt = conductivity(sheet) * {"Membrane_SiNx": 1e-6, "Membrane_Si1": 15e-6, "SiO2_1": 1e-6}[sheet]
        ring = ds[(ds > A_TES + STOL)]
        # NtD version (insensitive to Dirichlet node snapping): uniform heat
        # density over the TES footprint, perimeter isothermal; the lumped
        # load uses the sampled element volume fraction inside s<=250um.
        load = np.zeros(n)
        rng = np.random.default_rng(0)
        bary = rng.dirichlet(np.ones(4), size=64)
        for i in elems:
            nodes = mesh.elem_nodes[i]
            pts = mesh.coords[nodes]
            if int(mesh.elem_type[i]) == 504:
                sample = bary @ pts
                vol = abs(np.linalg.det(np.c_[np.ones(4), pts])) / 6.0
            else:
                sample = np.vstack([bary[:, :3] / bary[:, :3].sum(1, keepdims=True) @ pts[:3], bary[:, :3] / bary[:, :3].sum(1, keepdims=True) @ pts[3:]])
                vol = 0.5 * np.linalg.norm(np.cross(pts[1] - pts[0], pts[2] - pts[0])) * abs(pts[3:, 2].mean() - pts[:3, 2].mean())
            frac = float((square_radius(sample) <= A_TES).mean())
            if frac:
                np.add.at(load, mesh.node_to_dof[nodes], vol * frac / len(nodes))
        load /= load.sum()
        freeN = dofs[ds < A_WIN - STOL]
        lu = spla.splu(Ks.tocsc()[freeN][:, freeN], permc_spec="COLAMD")
        u = lu.solve(load[freeN])
        r_ntd = float(load[freeN] @ u)          # K/W for 1 W
        sheet_ops[sheet] = {
            "G_sheet_W_per_K": G, "k_times_t_W_per_K": kt, "shape_factor_mesh": G / kt,
            "NtD_uniform_footprint_G_W_per_K": 1.0 / r_ntd, "NtD_shape_factor_mesh": 1.0 / (r_ntd * kt),
            "hot_dofs": int(len(hot)), "cold_dofs": int(len(cold)), "free_dofs": int(len(op.free)),
            "hot_nodes_exactly_on_s250": int((np.abs(ds - A_TES) < STOL).sum()),
            "cold_nodes_exactly_on_s350": int((np.abs(ds - A_WIN) < STOL).sum()),
            "first_free_node_s_um": float(ring.min() / UM) if len(ring) else float("nan"),
            "annulus_elements": int(sum(1 for i in elems if A_TES < s_elem[i] < A_WIN)),
            "annulus_mean_edge_um": float(np.mean([np.mean([np.linalg.norm(mesh.coords[a] - mesh.coords[c]) for a in mesh.elem_nodes[i] for c in mesh.elem_nodes[i] if a < c]) for i in elems if A_TES < s_elem[i] < A_WIN]) / UM),
        }
    ops["sheet_annulus_operators"] = sheet_ops
    # (5) frame -> bath network (downstream of the window perimeter)
    frame_dofs = np.where(network_mask & (s_node >= A_WIN - STOL))[0]
    perim = frame_dofs[np.abs(s_node[frame_dofs] - A_WIN) < STOL]
    frame_elems = [i for i in net_elems if s_elem[i] > A_WIN]
    op = LocalOperator(sub_assembly(frame_elems), perim, dir_rows, frame_dofs, None)
    ops["Lambda_window_perimeter_to_bath"] = {"trace_dof_count": int(len(perim)),
                                              "constant_mode_total_reaction_W_per_K": float(op.apply(np.ones(len(perim))).sum())}
    # mesh statistics near the TES edge
    stats = {}
    for sheet in ("TES",) + SHEETS:
        sel = [i for i in range(len(kes)) if names[int(mesh.elem_body[i])] == sheet]
        stats[sheet] = {"elements": len(sel), "nodes": int(len(mesh.body_node_sets[b(sheet)]))}
    # in-plane resolution of the TES-footprint edge band (s = 230..270 um)
    band = {}
    for sheet in ("TES",) + SHEETS:
        lengths: dict[int, list[float]] = {}
        for i in np.where(mesh.elem_body == b(sheet))[0]:
            if not (230 * UM <= s_elem[i] <= 270 * UM):
                continue
            nodes = mesh.elem_nodes[i]
            for x in range(len(nodes)):
                for y in range(x + 1, len(nodes)):
                    pa, pb = mesh.coords[nodes[x]], mesh.coords[nodes[y]]
                    if abs(pa[2] - pb[2]) < ZTOL:
                        lengths.setdefault(int(round(pa[2] / UM)), []).append(float(np.linalg.norm(pa - pb) / UM))
        on_edge = {}
        for zc in sorted(lengths):
            sel = (np.abs(xyz[:, 2] - zc * UM) < ZTOL) & (np.abs(s_node - A_TES) < STOL)
            on_edge[str(zc)] = int(sel.sum())
        band[sheet] = {"inplane_edge_mean_um_by_plane": {str(z): float(np.mean(v)) for z, v in sorted(lengths.items())},
                       "nodes_on_TES_edge_line_by_plane": on_edge,
                       "elements_in_band": int(sum(1 for i in np.where(mesh.elem_body == b(sheet))[0] if 230 * UM <= s_elem[i] <= 270 * UM))}
    stats["TES_edge_band_230_270um"] = band
    stats["TES_bottom_trace_faces"] = len(faces["TES_bottom"])
    stats["Membrane_top_under_TES_faces"] = len(faces["Membrane_top_under_TES"])
    return {
        "case": key, "label": spec.label, "G_eff_harness_symmetric_W_per_K": g_eff,
        "element_reassembly_verification": ver, "mortar_rows": {"|".join(k): int(len(v)) for k, v in pairs.items()},
        "points": points, "operators": ops, "mode_rows": mode_rows, "mesh_stats": stats,
        "trace_face_counts": {k: len(v) for k, v in faces.items()},
        "trace_areas_m2": {k: float(sum(face_area(mesh, x) for x in v)) for k, v in faces.items()},
        "graph": actual_graph(mesh, pairs, dir_rows, faces),
        "provenance": provenance(key, spec, mesh, caps),
    }


def shared_faces(mesh: Mesh, a: int, bb: int) -> int:
    na = set(mesh.body_node_sets[a].tolist())
    nb = set(mesh.body_node_sets[bb].tolist())
    fa = {tuple(sorted(f.tolist())) for f in body_faces(mesh, a) if all(int(x) in nb for x in f)}
    fb = {tuple(sorted(f.tolist())) for f in body_faces(mesh, bb) if all(int(x) in na for x in f)}
    return len(fa & fb)


def actual_graph(mesh: Mesh, pairs, dir_rows, faces) -> dict[str, Any]:
    names = mesh.body_names
    bodies = sorted(mesh.body_node_sets)
    edges = []
    for i, a in enumerate(bodies):
        for bb in bodies[i + 1:]:
            shared = np.intersect1d(mesh.body_node_sets[a], mesh.body_node_sets[bb])
            if len(shared):
                nf = shared_faces(mesh, a, bb)
                edges.append({"bodies": [names[a], names[bb]], "coupling_type": "conformal shared-face" if nf else "conformal shared-node line only (rim)",
                              "shared_nodes": int(len(shared)), "shared_faces": nf})
    for pair, rows in pairs.items():
        edges.append({"bodies": list(pair), "coupling_type": "mortar Lagrange multiplier (Galerkin plane projector)", "constraint_rows": int(len(rows))})
    return {"bodies": {str(k): v for k, v in names.items()}, "edges": edges,
            "bath": {"dirichlet_rows": int(len(dir_rows)), "T_K": TBATH, "location": "SiO2_2 bottom z=-177um"}}


def provenance(key: str, spec, mesh: Mesh, caps: dict[float, Capture]) -> dict[str, Any]:
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    rows = {}
    for f, cap in caps.items():
        sif, log = spec.sif(f), spec.log(f)
        text = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
        rows[tag(f)] = {
            "P_W": P0 * f, "capture": str(cap.path.relative_to(ROOT)).replace("\\", "/"),
            "matrix_sha256": sha256(cap.path / "full_A_before.dat"), "rhs_sha256": sha256(cap.path / "full_b_before.dat"),
            "solution_sha256": sha256(cap.path / "full_x_after.dat"), "metadata_matrix_sha_before": cap.meta.get("matrix_sha_before"),
            "sif": str(sif.relative_to(ROOT)).replace("\\", "/") if sif.is_file() else None, "sif_sha256": sha256(sif) if sif.is_file() else None,
            "solver_log": str(log.relative_to(ROOT)).replace("\\", "/") if log.is_file() else None,
            "solver_ALL_DONE": "ALL DONE" in text,
            "backend": f"CPU {cap.meta.get('provenance', {}).get('solver_name')} {cap.meta.get('provenance', {}).get('direct_method')}",
            "capture_kind": cap.meta.get("capture_kind"),
            "capture_flags": {"Phase24 Full Restriction Capture": True, "TES Inner Circuit Freeze Power": True, "Nonlinear System Max Iterations": 1},
        }
    return {"case": key, "mesh": f"work/meshes/{spec.mesh}", "mesh_hash_sha256": mesh.hash, "audit_commit_parent_sha": head,
            "heatsolve_binary": "tools/elmer-hypre/install-steady-full-capture/bin/ElmerSolver.exe (CPU MUMPS; capture hook opt-in, default OFF)",
            "points": rows}


def main(argv: list[str]) -> int:
    (OUT / "raw").mkdir(parents=True, exist_ok=True)
    keys = argv or list(CASES)
    for key in keys:
        data = analyze_case(key)
        (OUT / "raw" / f"{key}.json").write_text(json.dumps(jsonable(data), indent=1) + "\n", encoding="utf-8")
        print(key, "G_eff", data["G_eff_harness_symmetric_W_per_K"], "verify", data["element_reassembly_verification"]["max_entry_diff_relative_to_row_diagonal"], flush=True)
    if not argv or argv == ["reference"]:
        ref = {"dirichlet": {f"h{h:g}um": annulus_shape_factor(h) for h in (4.0, 2.0, 1.0)},
               "ntd_uniform_footprint": {f"h{h:g}um": annulus_ntd_shape_factor(h) for h in (4.0, 2.0, 1.0)}}
        (OUT / "raw" / "annulus_reference.json").write_text(json.dumps(ref, indent=1) + "\n", encoding="utf-8")
        print(ref)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
