"""Compare direct and downstream trace operators for historical/refined meshes."""

from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
import sys
sys.path.insert(0, str(ROOT / "scripts" / "support"))
from analyze_phase24_mortar_element_stiffness import element_stiffness  # noqa: E402

OUT = ROOT / "artifacts" / "phase24_trace_resistance_controlled"
ELEMENT_NODES = {504: 4, 706: 6, 808: 8}
FACE_PATTERNS = {
    4: ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
    6: ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (0, 2, 5, 3)),
    8: ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
}

MESHES = {
    "historical": (ROOT / "work/meshes/mesh_singlepixel_prod_v2", 24, 23),
    "refined_mortar": (ROOT / "work/meshes/mesh_phase24_stycast_density_10um", 1104, 1305),
    "tes_membrane_control": (ROOT / "work/meshes/mesh_phase24_trace_tes_membrane_historical", 1104, 1305),
    "membrane_substrate_control": (ROOT / "work/meshes/mesh_phase24_trace_membrane_substrate_historical", 1104, 1305),
}


def parse_names(mesh: Path) -> dict[int, str]:
    section = False
    result: dict[int, str] = {}
    for line in (mesh / "mesh.names").read_text(encoding="utf-8", errors="replace").splitlines():
        if "names for bodies" in line:
            section = True
        elif "names for boundaries" in line:
            section = False
        if section:
            match = re.match(r"\s*\$\s*(.*?)\s*=\s*(-?\d+)\s*$", line)
            if match:
                result[int(match.group(2))] = match.group(1).strip()
    return result


def face_area(points: list[np.ndarray]) -> float:
    area = 0.5 * float(np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])))
    if len(points) == 4:
        area += 0.5 * float(np.linalg.norm(np.cross(points[3] - points[0], points[2] - points[0])))
    return area


def local_condensed_trace(points: list[np.ndarray], element_type: int, face_indices: list[int]) -> float:
    # Unit conductivity is intentional: this is a geometry/trace operator
    # comparison, not a material-property change.
    try:
        k = element_stiffness(points, element_type, 1.0)
    except (np.linalg.LinAlgError, ValueError):
        return float("nan")
    face = np.asarray(face_indices, dtype=int)
    interior = np.asarray([i for i in range(4) if i not in face], dtype=int)
    kff = k[np.ix_(face, face)]
    if interior.size:
        kii = k[np.ix_(interior, interior)]
        kfi = k[np.ix_(face, interior)]
        condensed = kff - kfi @ np.linalg.solve(kii, kfi.T)
    else:
        condensed = kff
    return float(np.trace(condensed))


def load_mesh(mesh: Path) -> tuple[dict[int, np.ndarray], dict[int, tuple[int, int, tuple[int, ...]]], dict[tuple[int, ...], list[tuple[int, int, tuple[int, ...]]]]]:
    nodes: dict[int, np.ndarray] = {}
    for fields in (line.split() for line in (mesh / "mesh.nodes").read_text(encoding="utf-8", errors="replace").splitlines()):
        if len(fields) >= 5 and fields[0].lstrip("-").isdigit():
            nodes[int(fields[0])] = np.asarray(list(map(float, fields[2:5])), dtype=float)
    elements: dict[int, tuple[int, int, tuple[int, ...]]] = {}
    faces: dict[tuple[int, ...], list[tuple[int, int, tuple[int, ...]]]] = defaultdict(list)
    for fields in (line.split() for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines()):
        if len(fields) < 4 or not fields[0].lstrip("-").isdigit():
            continue
        eid, body, kind = map(int, fields[:3])
        count = ELEMENT_NODES.get(kind)
        if count is None or len(fields) < 3 + count:
            continue
        ids = tuple(map(int, fields[3:3 + count]))
        elements[eid] = (body, kind, ids)
        for pattern in FACE_PATTERNS[count]:
            face_nodes = tuple(ids[index] for index in pattern)
            faces[tuple(sorted(face_nodes))].append((body, eid, face_nodes))
    return nodes, elements, faces


def boundary_rows(mesh: Path, boundary_id: int) -> list[tuple[int, tuple[int, ...]]]:
    rows = []
    for fields in (line.split() for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines()):
        if len(fields) < 6 or not fields[0].lstrip("-").isdigit() or int(fields[1]) != boundary_id:
            continue
        count = 3 if int(fields[4]) == 303 else 4 if int(fields[4]) == 404 else 0
        if count:
            rows.append((int(fields[2]), tuple(map(int, fields[5:5 + count]))))
    return rows


def direct_trace(label: str, mesh: Path, slave_id: int, master_id: int) -> list[dict[str, Any]]:
    nodes, elements, _ = load_mesh(mesh)
    output = []
    for side, boundary_id in (("slave", slave_id), ("master", master_id)):
        rows = boundary_rows(mesh, boundary_id)
        areas = []
        edges = []
        condensed = []
        for parent, face_nodes in rows:
            body, kind, element_nodes = elements[parent]
            points = [nodes[node] for node in face_nodes]
            areas.append(face_area(points))
            edges.extend(float(np.linalg.norm(points[i] - points[j])) for i in range(len(points)) for j in range(i))
            indices = [element_nodes.index(node) for node in face_nodes]
            condensed.append(local_condensed_trace([nodes[node] for node in element_nodes], kind, indices))
        output.append({
            "case": label,
            "trace": f"TES_to_Membrane_{side}",
            "boundary_id": boundary_id,
            "face_count": len(rows),
            "area_m2": sum(areas),
            "mean_edge_m": float(np.mean(edges)) if edges else float("nan"),
            "condensed_trace_K_unit_sum": float(np.nansum(condensed)),
            "condensed_trace_K_unit_median": float(np.nanmedian(condensed)) if condensed else float("nan"),
        })
    return output


def shared_trace_summary(label: str, mesh: Path) -> list[dict[str, Any]]:
    nodes, elements, faces = load_mesh(mesh)
    names = parse_names(mesh)
    rows: dict[tuple[int, int], dict[str, Any]] = {}
    for owners in faces.values():
        by_body: dict[int, tuple[int, int, tuple[int, ...]]] = {}
        for owner in owners:
            by_body.setdefault(owner[0], owner)
        if len(by_body) < 2:
            continue
        body_ids = sorted(by_body)
        for left_index, left in enumerate(body_ids):
            for right in body_ids[left_index + 1:]:
                if left in (100, 101, 102) or right in (100, 101, 102):
                    continue
                key = (left, right)
                row = rows.setdefault(key, {"case": label, "trace": f"{names.get(left, left)}->{names.get(right, right)}", "body_a": left, "body_b": right, "face_count": 0, "area_m2": 0.0, "condensed_a_sum": 0.0, "condensed_b_sum": 0.0})
                row["face_count"] += 1
                for body, target in ((left, "condensed_a_sum"), (right, "condensed_b_sum")):
                    _, element_id, face_nodes = by_body[body]
                    _, kind, element_nodes = elements[element_id]
                    points = [nodes[node] for node in face_nodes]
                    row["area_m2"] += face_area(points) if body == left else 0.0
                    row[target] += local_condensed_trace([nodes[node] for node in element_nodes], kind, [element_nodes.index(node) for node in face_nodes])
    return list(rows.values())


def main() -> int:
    direct_rows: list[dict[str, Any]] = []
    shared_rows: list[dict[str, Any]] = []
    for label, (mesh, slave, master) in MESHES.items():
        direct_rows.extend(direct_trace(label, mesh, slave, master))
        shared_rows.extend(shared_trace_summary(label, mesh))
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "direct_trace_operator.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(direct_rows[0])); writer.writeheader(); writer.writerows(direct_rows)
    with (OUT / "shared_trace_operator.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(shared_rows[0])); writer.writeheader(); writer.writerows(shared_rows)
    print("direct TES/Membrane trace")
    for row in direct_rows:
        print(row)
    print("downstream shared traces")
    for row in shared_rows:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
