"""Mesh-derived quantities needed at SIF build time.

Computes, for a given Elmer mesh directory, the absorber-body centroid (used
as the default pulse center) and the FE integral of a nodal-sampled Gaussian
over the absorber (the `Pulse Discrete Norm` that makes the deposited pulse
energy exact regardless of how coarsely the mesh resolves the profile).
"""
from __future__ import annotations

from pathlib import Path

from scripts.support.mesh_names import parse_mesh_names

TET_TYPE = "504"


def _load_abs_tets(mesh_dir: Path):
    import numpy as np

    # Body IDs are assigned by the mesh converter.  The original tetra mesh
    # happened to use 100 for `abs`; the hybrid prism mesh correctly uses a
    # compact physical-group ID instead, so resolve the name from mesh.names.
    absorber_body_id = parse_mesh_names(mesh_dir / "mesh.names").bodies.get("abs")
    if absorber_body_id is None:
        raise ValueError(f"No absorber body named 'abs' in {mesh_dir / 'mesh.names'}")
    nodes = np.loadtxt(mesh_dir / "mesh.nodes", usecols=(2, 3, 4))
    tets = []
    with (mesh_dir / "mesh.elements").open() as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 7 and parts[1] == str(absorber_body_id) and parts[2] == TET_TYPE:
                tets.append([int(x) - 1 for x in parts[3:7]])
    if not tets:
        raise ValueError(f"No absorber (body {ABS_BODY_ID}) tetrahedra found in {mesh_dir}")
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
