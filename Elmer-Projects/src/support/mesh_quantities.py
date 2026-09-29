"""Mesh-derived quantities needed at SIF build time.

Computes, for a given Elmer mesh directory, the absorber-body centroid (used
as the default pulse center) and the FE integral of a nodal-sampled Gaussian
over the absorber (the `Pulse Discrete Norm` that makes the deposited pulse
energy exact regardless of how coarsely the mesh resolves the profile).
"""
from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path
_src_path = str(_Path(__file__).resolve().parents[1])
if _src_path not in _sys.path:
    _sys.path.insert(0, _src_path)

from pathlib import Path

from support.mesh_names import parse_mesh_names

TET_TYPE = "504"
PRISM_TYPE = "706"
# Bodies that receive the pulse (their union is "the absorber" for the pulse
# centre and the discrete norm).  build_cases sets it from the mesh role table;
# the default is the conventional body name.
PULSE_BODY_NAMES: tuple[str, ...] = ("abs",)


def set_pulse_bodies(names) -> None:
    global PULSE_BODY_NAMES
    PULSE_BODY_NAMES = tuple(names) if names else ("abs",)


def _load_abs_tets(mesh_dir: Path):
    """Nodes, linear tetrahedra and their volumes of the pulse bodies.

    Wedges (706) are split into three tetrahedra, so a prism absorber works
    too (the nodal-interpolated integrals stay exact for linear fields)."""
    import numpy as np

    # Body IDs are assigned by the mesh converter: resolve them by name.
    bodies = parse_mesh_names(mesh_dir / "mesh.names").bodies
    missing = [n for n in PULSE_BODY_NAMES if n not in bodies]
    if missing:
        raise ValueError(f"pulse bodies {missing} not in {mesh_dir / 'mesh.names'}")
    body_ids = {str(bodies[n]) for n in PULSE_BODY_NAMES}
    nodes = np.loadtxt(mesh_dir / "mesh.nodes", usecols=(2, 3, 4))
    tets = []
    with (mesh_dir / "mesh.elements").open() as f:
        for line in f:
            parts = line.split()
            if len(parts) < 7 or parts[1] not in body_ids:
                continue
            if parts[2] == TET_TYPE:
                tets.append([int(x) - 1 for x in parts[3:7]])
            elif parts[2] == PRISM_TYPE and len(parts) >= 9:
                a, b, c, d, e, g = (int(x) - 1 for x in parts[3:9])
                tets += [[a, b, c, d], [b, c, d, e], [c, d, e, g]]
    if not tets:
        raise ValueError(f"No tetrahedra/wedges in pulse bodies {PULSE_BODY_NAMES} of {mesh_dir}")
    tets = np.array(tets)
    p1, p2, p3, p4 = (nodes[tets[:, i]] for i in range(4))
    volumes = np.abs(np.einsum("ij,ij->i", p2 - p1, np.cross(p3 - p1, p4 - p1))) / 6.0
    return nodes, tets, volumes


def absorber_centroid(mesh_dir: Path) -> tuple[float, float, float]:
    """Volume-weighted centroid of the absorber body."""
    import numpy as np

    nodes, tets, volumes = _load_abs_tets(mesh_dir)
    tet_centroids = nodes[tets].mean(axis=1)
    centroid = (volumes[:, None] * tet_centroids).sum(axis=0) / volumes.sum()
    return tuple(float(c) for c in centroid)


def gaussian_discrete_norm(
    mesh_dir: Path, center: tuple[float, float, float], sigma: float
) -> float:
    """FE integral over the absorber of the nodal-interpolated Gaussian
    exp(-r^2 / 2 sigma^2) centred at *center* (exact for linear tets)."""
    import numpy as np

    nodes, tets, volumes = _load_abs_tets(mesh_dir)
    r2 = ((nodes - np.asarray(center)) ** 2).sum(axis=1)
    weights = np.exp(-r2 / (2.0 * sigma * sigma))
    return float((volumes * weights[tets].mean(axis=1)).sum())


def sphere_discrete_norm(
    mesh_dir: Path, center: tuple[float, float, float], radius: float
) -> float:
    """FE integral of a nodal indicator for a uniform spherical source."""
    import numpy as np

    nodes, tets, volumes = _load_abs_tets(mesh_dir)
    r2 = ((nodes - np.asarray(center)) ** 2).sum(axis=1)
    weights = (r2 <= radius * radius).astype(float)
    norm = float((volumes * weights[tets].mean(axis=1)).sum())
    if norm <= 0.0:
        raise ValueError("uniform-sphere pulse has zero discrete absorber volume")
    return norm
