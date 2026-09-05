"""Validate shared-node contact topology in an Elmer ASCII mesh.

This is deliberately a post-ElmerGrid gate.  CAD coincidence or a successful
Gmsh ``setPeriodic`` call is not enough: the converted Elmer mesh must contain
the same node IDs and the same surface-element partition on both sides of
each physical contact patch.

Example::

    python scripts/analysis/check_conformal_interfaces.py \
      work/meshes/mesh_singlepixel_conformal_gpu \
      --output artifacts/phase20_conformal/interface_connectivity.json
"""
from __future__ import annotations

import argparse
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


DEFAULT_PAIRS = (
    ("Membrane_SiNx__zmax", "TES__zmin", "Membrane_TES", "Membrane_SiNx", "TES"),
    ("TES__zmax", "Stycast__zmin", "TES_Stycast", "TES", "Stycast"),
    ("Stycast__zmax", "abs__zmin", "Stycast_absorber", "Stycast", "abs"),
)


def read_names(path: Path) -> tuple[dict[str, int], dict[str, int]]:
    bodies: dict[str, int] = {}
    boundaries: dict[str, int] = {}
    section = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "names for bodies" in line.lower():
            section = bodies
        elif "names for boundaries" in line.lower():
            section = boundaries
        match = re.match(r"^\$\s+(.+?)\s*=\s*(-?\d+)\s*$", line)
        if match and section is not None:
            section[match.group(1)] = int(match.group(2))
    return bodies, boundaries


def read_nodes(path: Path) -> dict[int, tuple[float, float, float]]:
    nodes: dict[int, tuple[float, float, float]] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5:
            nodes[int(fields[0])] = tuple(float(value) for value in fields[2:5])
    return nodes


def read_boundary(path: Path) -> dict[int, list[frozenset[int]]]:
    result: dict[int, list[frozenset[int]]] = defaultdict(list)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 6:
            continue
        result[int(fields[1])].append(frozenset(int(value) for value in fields[5:]))
    return dict(result)


def read_elements(path: Path) -> dict[int, list[frozenset[int]]]:
    result: dict[int, list[frozenset[int]]] = defaultdict(list)
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 5:
            continue
        result[int(fields[1])].append(frozenset(int(value) for value in fields[3:]))
    return dict(result)


def nearest_gap(
    left: Iterable[tuple[float, float, float]],
    right: Iterable[tuple[float, float, float]],
) -> float:
    right_points = list(right)
    if not right_points:
        return float("inf")
    maximum = 0.0
    for point in left:
        maximum = max(
            maximum,
            min(math.dist(point, candidate) for candidate in right_points),
        )
    return maximum


def mesh_quality(
    mesh: Path,
    nodes: dict[int, tuple[float, float, float]],
    elements: dict[int, list[frozenset[int]]],
) -> dict[str, object]:
    zero = 0
    nonfinite = 0
    tetrahedra = 0
    volumes: list[float] = []
    for body_elements in elements.values():
        for element in body_elements:
            if len(element) != 4:
                continue
            tetrahedra += 1
            points = [nodes[node] for node in element]
            a, b, c, d = points
            cross = (
                (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
                (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
                (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]),
            )
            signed = (
                cross[0] * (d[0] - a[0])
                + cross[1] * (d[1] - a[1])
                + cross[2] * (d[2] - a[2])
            ) / 6.0
            volume = abs(signed)
            if not math.isfinite(volume):
                nonfinite += 1
            elif volume == 0.0:
                zero += 1
            else:
                volumes.append(volume)
    return {
        "tetrahedra": tetrahedra,
        "min_abs_tet_volume_m3": min(volumes) if volumes else None,
        "max_abs_tet_volume_m3": max(volumes) if volumes else None,
        "zero_volume_elements": zero,
        "nonfinite_volume_elements": nonfinite,
        "duplicate_volume_connectivity": duplicate_connectivity(elements),
        "status": "PASS" if zero == 0 and nonfinite == 0 else "FAIL",
    }


def duplicate_connectivity(elements: dict[int, list[frozenset[int]]]) -> int:
    count = 0
    for body_elements in elements.values():
        frequencies = Counter(body_elements)
        count += sum(value - 1 for value in frequencies.values() if value > 1)
    return count


def check_pair(
    pair: tuple[str, str, str, str, str],
    boundaries: dict[str, int],
    boundary_elements: dict[int, list[frozenset[int]]],
    nodes: dict[int, tuple[float, float, float]],
    body_ids: dict[str, int],
) -> dict[str, object]:
    left_name, right_name, label, left_body, right_body = pair
    left_elements = boundary_elements.get(boundaries.get(left_name, -1), [])
    right_elements = boundary_elements.get(boundaries.get(right_name, -1), [])
    left_nodes = set().union(*left_elements) if left_elements else set()
    right_nodes = set().union(*right_elements) if right_elements else set()
    left_coords = [nodes[node] for node in left_nodes]
    right_coords = [nodes[node] for node in right_nodes]
    common_surface_elements = len(set(left_elements).intersection(right_elements))
    coord_gap = nearest_gap(left_coords, right_coords)
    body_presence = left_body in body_ids and right_body in body_ids
    status = (
        bool(left_elements)
        and bool(right_elements)
        and left_nodes == right_nodes
        and set(left_elements) == set(right_elements)
        and coord_gap <= 1.0e-12
        and body_presence
    )
    return {
        "interface": label,
        "left_surface": left_name,
        "right_surface": right_name,
        "left_boundary_id": boundaries.get(left_name),
        "right_boundary_id": boundaries.get(right_name),
        "shared_nodes": len(left_nodes.intersection(right_nodes)),
        "left_only_nodes": len(left_nodes - right_nodes),
        "right_only_nodes": len(right_nodes - left_nodes),
        "left_surface_elements": len(left_elements),
        "right_surface_elements": len(right_elements),
        "matched_surface_elements": common_surface_elements,
        "max_coordinate_gap_m": coord_gap,
        "connected_element_adjacency": common_surface_elements > 0,
        "status": "PASS" if status else "FAIL",
        "failure_reasons": [
            reason
            for condition, reason in (
                (not left_elements, "left_surface_missing"),
                (not right_elements, "right_surface_missing"),
                (left_nodes != right_nodes, "node_id_sets_differ"),
                (set(left_elements) != set(right_elements), "surface_partition_differs"),
                (coord_gap > 1.0e-12, "coordinate_gap_exceeds_tolerance"),
                (not body_presence, "body_missing"),
            )
            if condition
        ],
    }


def inspect(mesh: Path) -> dict[str, object]:
    bodies, boundaries = read_names(mesh / "mesh.names")
    nodes = read_nodes(mesh / "mesh.nodes")
    boundary_elements = read_boundary(mesh / "mesh.boundary")
    elements = read_elements(mesh / "mesh.elements")
    pairs = [check_pair(pair, boundaries, boundary_elements, nodes, bodies) for pair in DEFAULT_PAIRS]
    return {
        "mesh": str(mesh.resolve()),
        "node_count": len(nodes),
        "volume_element_count": sum(len(value) for value in elements.values()),
        "body_ids": bodies,
        "boundary_ids": boundaries,
        "interfaces": pairs,
        "mesh_quality": mesh_quality(mesh, nodes, elements),
        "temperature_jump": {entry["interface"]: None for entry in pairs},
        "heat_flux_consistency": {entry["interface"]: None for entry in pairs},
        "status": "PASS"
        if all(entry["status"] == "PASS" for entry in pairs)
        else "FAIL",
        "limitations": [
            "temperature_jump and heat_flux_consistency require a result-field reader and are reported as null",
            "coordinate and topology checks are post-ElmerGrid; they do not infer connectivity from CAD coincidence",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = inspect(args.mesh)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "mesh": report["mesh"]}, indent=2))
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
