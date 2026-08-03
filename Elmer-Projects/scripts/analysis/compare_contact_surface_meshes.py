"""Summarize named contact-surface discretizations in Elmer ASCII meshes."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


DEFAULT_SURFACES = (
    "TES__zmin", "TES__zmax", "Stycast__zmin", "Stycast__zmax",
    "abs__zmin", "Membrane_SiNx__zmax",
)


def names(path: Path) -> dict[str, int]:
    return {
        name: int(value)
        for name, value in re.findall(r"^\$\s+(.+?)\s+=\s+(\d+)\s*$", path.read_text(), re.M)
    }


def nodes(path: Path) -> dict[int, tuple[float, float, float]]:
    result: dict[int, tuple[float, float, float]] = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        result[int(fields[0])] = tuple(map(float, fields[2:5]))
    return result


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def triangle_area(a, b, c) -> float:
    u = [b[i] - a[i] for i in range(3)]
    v = [c[i] - a[i] for i in range(3)]
    cross = (
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    )
    return 0.5 * math.sqrt(sum(x * x for x in cross))


def summarize(mesh: Path, requested: tuple[str, ...]) -> dict[str, dict]:
    ids = names(mesh / "mesh.names")
    xyz = nodes(mesh / "mesh.nodes")
    requested_ids = {ids[name]: name for name in requested if name in ids}
    groups: dict[str, dict] = {
        name: {"boundary_id": boundary_id, "elements": 0, "area_m2": 0.0, "edges": []}
        for boundary_id, name in requested_ids.items()
    }
    for line in (mesh / "mesh.boundary").read_text().splitlines():
        fields = line.split()
        boundary_id = int(fields[1])
        name = requested_ids.get(boundary_id)
        if name is None:
            continue
        kind = int(fields[4])
        node_ids = [int(value) for value in fields[5:]]
        points = [xyz[node_id] for node_id in node_ids]
        if kind == 303:
            area = triangle_area(*points)
        elif kind == 404:
            area = triangle_area(points[0], points[1], points[2]) + triangle_area(points[0], points[2], points[3])
        else:
            raise ValueError(f"{mesh}: unsupported boundary type {kind} on {name}")
        entry = groups[name]
        entry["elements"] += 1
        entry["area_m2"] += area
        entry["edges"].extend(distance(points[i], points[(i + 1) % len(points)]) for i in range(len(points)))
    for entry in groups.values():
        edges = entry.pop("edges")
        entry["mean_edge_um"] = math.fsum(edges) / len(edges) * 1.0e6
        entry["max_edge_um"] = max(edges) * 1.0e6
    return groups


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reference = summarize(args.reference, DEFAULT_SURFACES)
    candidate = summarize(args.candidate, DEFAULT_SURFACES)
    comparison = {}
    for name, ref in reference.items():
        if name in candidate:
            cand = candidate[name]
            comparison[name] = {
                "area_relative_difference": (cand["area_m2"] - ref["area_m2"]) / ref["area_m2"],
                "element_ratio": cand["elements"] / ref["elements"],
                "mean_edge_ratio": cand["mean_edge_um"] / ref["mean_edge_um"],
            }
    args.output.write_text(json.dumps({"reference": reference, "candidate": candidate, "comparison": comparison}, indent=2) + "\n")


if __name__ == "__main__":
    main()
