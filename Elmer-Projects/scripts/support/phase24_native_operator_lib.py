"""Shared loaders for the Phase24 native branch/operator audit.

Everything here reads *existing* solver-native artifacts only:

* ``full_A_before.dat`` / ``full_b_before.dat`` / ``full_x_after.dat`` written by
  the opt-in ``Phase24 Full Restriction Capture`` hook inside
  ``SolveWithLinearRestriction`` (native HeatSolve CPU/MUMPS), i.e. the exact
  saddle-point system ``[[K, C^T], [C, D]] [T; lambda] = [f; g]`` that MUMPS
  factorized, and the solution MUMPS returned;
* the Elmer mesh files and the ``.result`` permutation (node -> DOF).

No solver is re-run and no physics is evaluated in Python.  Heat flows are
derived from the assembled operator only:

* conductive edge flow  ``q_ij = K_ij (T_j - T_i)``  (out of node i into j),
  so that for every non-Dirichlet primal row ``f_i = sum_j q_ij + m_i``;
* mortar reaction       ``m = C^T lambda`` restricted to the constraint rows of
  one mortar pair (flow leaving the node through that interface).
"""
from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
P0 = 3.203004762115138e-10
TBATH = 0.15
FRACTIONS = (0.95, 1.0, 1.05)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def tag(fraction: float) -> str:
    return f"{fraction:.2f}P".replace(".", "p")


@dataclass(frozen=True)
class CaseSpec:
    key: str
    label: str
    mesh: str
    capture_root: Path
    capture_prefix: str
    result_name: str
    sif_dir: Path

    def capture(self, fraction: float) -> Path:
        return self.capture_root / f"{self.capture_prefix}_{tag(fraction)}" / "ts0001_nl0001"

    def sif(self, fraction: float) -> Path:
        return self.sif_dir / f"{self.capture_prefix}_{tag(fraction)}.sif"

    def log(self, fraction: float) -> Path:
        return self.sif_dir / f"{self.capture_prefix}_{tag(fraction)}.solver.log"


CASES = {
    "historical": CaseSpec(
        "historical", "Historical reference", "mesh_singlepixel_prod_v2",
        ROOT / "artifacts/phase24_thermal_network_localization/capture", "historical",
        "historical_1p00p.result", ROOT / "artifacts/phase24_thermal_network_localization",
    ),
    "best_phase24": CaseSpec(
        "best_phase24", "Best Phase24 (Membrane/substrate h=10um control)",
        "mesh_phase24_trace_membrane_substrate_historical",
        ROOT / "artifacts/p24trace/ms/capture", "ms_hist", "ms_hist_1p00p.result",
        ROOT / "artifacts/p24trace/ms",
    ),
    "refined_parent": CaseSpec(
        "refined_parent", "Refined mortar parent (Stycast 10um)", "mesh_phase24_stycast_density_10um",
        ROOT / "artifacts/p24d10b/capture", "phase24_historical_density",
        "phase24_historical_density_1p00p.result", ROOT / "artifacts/p24d10b",
    ),
}


def _read_table(path: Path) -> np.ndarray:
    text = path.read_text(encoding="utf-8", errors="replace").replace("D", "E")
    return np.array(text.split(), dtype=float)


@dataclass
class Capture:
    path: Path
    A: sp.csr_matrix
    b: np.ndarray
    x: np.ndarray
    n_primal: int
    meta: dict

    @property
    def K(self) -> sp.csr_matrix:
        return self.A[: self.n_primal, : self.n_primal]

    @property
    def C(self) -> sp.csr_matrix:
        return self.A[self.n_primal:, : self.n_primal]

    @property
    def T(self) -> np.ndarray:
        return self.x[: self.n_primal]

    @property
    def lam(self) -> np.ndarray:
        return self.x[self.n_primal:]


def load_capture(path: Path) -> Capture:
    meta = json.loads((path / "metadata.json").read_text(encoding="utf-8"))
    n = int(meta["runtime"]["total_rows"])
    npr = int(meta["runtime"]["primal_rows"])
    a = _read_table(path / "full_A_before.dat").reshape(-1, 3)
    A = sp.csr_matrix((a[:, 2], (a[:, 0].astype(np.int64) - 1, a[:, 1].astype(np.int64) - 1)), shape=(n, n))
    A.sum_duplicates()

    def vec(name: str) -> np.ndarray:
        v = _read_table(path / name).reshape(-1, 2)
        out = np.zeros(n)
        out[v[:, 0].astype(np.int64) - 1] = v[:, 1]
        return out

    return Capture(path, A, vec("full_b_before.dat"), vec("full_x_after.dat"), npr, meta)


@dataclass
class Mesh:
    path: Path
    coords: np.ndarray                      # (n_nodes, 3), node id = index + 1
    elem_body: np.ndarray                   # (n_elem,)
    elem_type: np.ndarray                   # Elmer element code (504 tet, 706 wedge, ...)
    elem_nodes: list[np.ndarray]
    boundary_faces: dict[int, list[np.ndarray]]
    body_names: dict[int, str]
    boundary_names: dict[int, str]
    node_to_dof: np.ndarray                 # 0-based dof per 0-based node, -1 if absent
    body_node_sets: dict[int, np.ndarray] = field(default_factory=dict)

    def name_to_body(self, name: str) -> int:
        return next(k for k, v in self.body_names.items() if v == name)

    def dofs(self, nodes0: np.ndarray) -> np.ndarray:
        d = self.node_to_dof[nodes0]
        if (d < 0).any():
            raise RuntimeError("node without temperature DOF")
        return d

    def boundary_nodes(self, bid: int) -> np.ndarray:
        faces = self.boundary_faces.get(bid, [])
        return np.unique(np.concatenate(faces)) if faces else np.zeros(0, dtype=np.int64)

    def face_area_weights(self, bid: int) -> tuple[np.ndarray, np.ndarray]:
        """Lumped nodal area weights of a named boundary (node0, weight)."""
        weights: dict[int, float] = {}
        for face in self.boundary_faces.get(bid, []):
            p = self.coords[face]
            area = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]))
            if len(face) == 4:
                area += 0.5 * np.linalg.norm(np.cross(p[3] - p[0], p[2] - p[0]))
            for node in face:
                weights[int(node)] = weights.get(int(node), 0.0) + area / len(face)
        keys = np.array(sorted(weights), dtype=np.int64)
        return keys, np.array([weights[k] for k in keys])

    @cached_property
    def hash(self) -> str:
        return ":".join(sha256(self.path / f) for f in ("mesh.nodes", "mesh.elements", "mesh.boundary"))


def _names(mesh: Path) -> tuple[dict[int, str], dict[int, str]]:
    bodies: dict[int, str] = {}
    boundaries: dict[int, str] = {}
    target: dict[int, str] | None = None
    for line in (mesh / "mesh.names").read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if "names for bodies" in s:
            target = bodies
        elif "names for boundaries" in s:
            target = boundaries
        elif target is not None:
            m = re.match(r"\$\s*(\S+)\s*=\s*(-?\d+)", s)
            if m:
                target[int(m.group(2))] = m.group(1)
    return bodies, boundaries


def _perm(result: Path, n_nodes: int) -> np.ndarray:
    lines = result.read_text(encoding="utf-8", errors="replace").splitlines()
    marker = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
    count = int(lines[marker + 1].split()[1])
    out = -np.ones(n_nodes, dtype=np.int64)
    for offset in range(count):
        node, dof = lines[marker + 2 + offset].split()[:2]
        out[int(node) - 1] = int(dof) - 1
    return out


def load_mesh(spec: CaseSpec) -> Mesh:
    mesh = ROOT / "work/meshes" / spec.mesh
    nodes = _read_table(mesh / "mesh.nodes").reshape(-1, 5)
    order = np.argsort(nodes[:, 0])
    coords = nodes[order, 2:5]
    elem_body: list[int] = []
    elem_type: list[int] = []
    elem_nodes: list[np.ndarray] = []
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8").splitlines():
        f = line.split()
        if len(f) >= 4:
            elem_body.append(int(f[1]))
            elem_type.append(int(f[2]))
            elem_nodes.append(np.array(f[3:], dtype=np.int64) - 1)
    faces: dict[int, list[np.ndarray]] = {}
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8").splitlines():
        f = line.split()
        if len(f) >= 8:
            code = int(f[4])
            n = {303: 3, 404: 4}.get(code, len(f) - 5)
            faces.setdefault(int(f[1]), []).append(np.array(f[5:5 + n], dtype=np.int64) - 1)
    bodies, boundaries = _names(mesh)
    perm = _perm(mesh / spec.result_name, len(coords))
    m = Mesh(mesh, coords, np.array(elem_body), np.array(elem_type), elem_nodes, faces, bodies, boundaries, perm)
    for body in np.unique(m.elem_body):
        idx = np.where(m.elem_body == body)[0]
        m.body_node_sets[int(body)] = np.unique(np.concatenate([elem_nodes[i] for i in idx]))
    return m


def body_mean_temperature_like_harness(mesh: Mesh, T: np.ndarray, body: int) -> float:
    """TES average exactly as the frozen-power harness defines it (element-node mean)."""
    idx = np.where(mesh.elem_body == body)[0]
    nodes = np.concatenate([mesh.elem_nodes[i] for i in idx])
    return float(np.mean(T[mesh.dofs(nodes)]))


def dirichlet_rows(K: sp.csr_matrix) -> np.ndarray:
    """Rows whose only stored nonzero is the diagonal (Elmer Dirichlet rows)."""
    K = K.tocsr()
    diag = K.diagonal()
    offabs = np.asarray(abs(K).sum(axis=1)).ravel() - np.abs(diag)
    return np.where(offabs == 0.0)[0]
