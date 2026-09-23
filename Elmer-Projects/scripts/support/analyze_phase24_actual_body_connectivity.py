#!/usr/bin/env python3
"""Reconstruct the thermal body graph from Elmer mesh connectivity.

This diagnostic deliberately does not infer the physical network from boundary
names.  It reads volume elements, node ownership, boundary facets, and the
body/material declarations used by the SIF.  The output is intended for the
Phase24 abs/mortar investigation and is safe to run against existing meshes.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


ELEMENT_NODE_COUNTS = {504: 4, 706: 6, 808: 8, 409: 10, 510: 20}
FACE_PATTERNS = {
    4: ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
    6: ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (0, 2, 5, 3)),
    8: ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
}


def parse_names(path: Path) -> tuple[dict[int, str], dict[int, str]]:
    section = None
    bodies: dict[int, str] = {}
    boundaries: dict[int, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "names for bodies" in line:
            section = "body"
        elif "names for boundaries" in line:
            section = "boundary"
        match = re.match(r"\s*\$\s+(.+?)\s*=\s*(-?\d+)\s*$", line)
        if not match or section is None:
            continue
        name, number = match.group(1), int(match.group(2))
        (bodies if section == "body" else boundaries)[number] = name
    return bodies, boundaries


def parse_sif(path: Path) -> tuple[dict[int, dict], dict[int, dict], list[dict], list[int]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    materials: dict[int, dict] = {}
    for match in re.finditer(r"(?ms)^Material\s+(\d+)\s*\n(.*?)^End\s*$", text):
        block = match.group(2)
        name = re.search(r'^\s*Name\s*=\s*"([^"]+)"', block, re.M)
        density = re.search(r"^\s*Density\s*=\s*(.+)$", block, re.M)
        cp = re.search(r"^\s*Heat Capacity\s*=\s*(.+)$", block, re.M)
        conductivity = re.search(r"^\s*Heat Conductivity\s*=\s*(.+)$", block, re.M)
        materials[int(match.group(1))] = {
            "name": name.group(1) if name else "",
            "density": density.group(1).strip() if density else "",
            "cp": cp.group(1).strip() if cp else "",
            "conductivity": conductivity.group(1).strip() if conductivity else "",
            "block_normalized": " ".join(block.split()),
        }
    bodies: dict[int, dict] = {}
    for match in re.finditer(r"(?ms)^Body\s+(\d+)\s*\n(.*?)^End\s*$", text):
        block = match.group(2)
        target = re.search(r"Target Bodies\(1\)\s*=\s*(\d+)", block)
        name = re.search(r'^\s*Name\s*=\s*"([^"]+)"', block, re.M)
        material = re.search(r"Material\s*=\s*(\d+)", block)
        if target:
            material_id = int(material.group(1)) if material else None
            bodies[int(target.group(1))] = {
                "body_id": int(target.group(1)),
                "body_name": name.group(1) if name else "",
                "material_id": material_id,
                "material": materials.get(material_id, {}).get("name", "") if material_id else "",
            }
    mortar: list[dict] = []
    for match in re.finditer(r"(?ms)^Boundary Condition\s+(\d+)\s*\n(.*?)^End\s*$", text):
        block = match.group(2)
        target = re.search(r"Target Boundaries\(1\)\s*=\s*(\d+)", block)
        mortar_bc = re.search(r"Mortar BC\s*=\s*(\d+)", block)
        if target and mortar_bc:
            mortar.append({
                "bc_number": int(match.group(1)),
                "slave_boundary_id": int(target.group(1)),
                "master_bc_number": int(mortar_bc.group(1)),
            })
    bath_boundaries: list[int] = []
    for match in re.finditer(r"(?ms)^Boundary Condition\s+(\d+)\s*\n(.*?)^End\s*$", text):
        block = match.group(2)
        target = re.search(r"Target Boundaries\(1\)\s*=\s*(\d+)", block)
        if target and re.search(r"Temperature\s*=\s*0\.15", block):
            bath_boundaries.append(int(target.group(1)))
    return bodies, materials, mortar, bath_boundaries


def parse_mesh(mesh: Path) -> dict:
    body_names, boundary_names = parse_names(mesh / "mesh.names")
    nodes: dict[int, tuple[float, float, float]] = {}
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[0].lstrip("-").isdigit() and fields[1] == "-1":
            nodes[int(fields[0])] = tuple(float(v) for v in fields[2:5])

    elements: dict[int, tuple[int, int, tuple[int, ...]]] = {}
    body_nodes: defaultdict[int, set[int]] = defaultdict(set)
    element_counts: defaultdict[int, Counter] = defaultdict(Counter)
    face_owners: defaultdict[tuple[int, ...], list[tuple[int, int]]] = defaultdict(list)
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 4:
            continue
        try:
            element_id, body_id, element_type = map(int, fields[:3])
        except ValueError:
            continue
        count = ELEMENT_NODE_COUNTS.get(element_type)
        if count is None or len(fields) < 3 + count:
            continue
        node_ids = tuple(int(value) for value in fields[3:3 + count])
        elements[element_id] = (body_id, element_type, node_ids)
        body_nodes[body_id].update(node_ids)
        element_counts[body_id][element_type] += 1
        for face in FACE_PATTERNS.get(count, ()):
            key = tuple(sorted(node_ids[index] for index in face))
            face_owners[key].append((body_id, element_id))

    boundaries: list[dict] = []
    boundary_by_id: defaultdict[int, list[dict]] = defaultdict(list)
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 6:
            continue
        try:
            boundary_id, condition_id, parent_element, boundary_type = map(int, fields[:4])
        except ValueError:
            continue
        # Elmer boundary records are: id, condition, parent element, local
        # side, element type, node ids...  The element type is not a node.
        node_ids = tuple(int(value) for value in fields[5:])
        body_id = elements.get(parent_element, (None, None, ()))[0]
        row = {
            "boundary_id": boundary_id,
            "condition_id": condition_id,
            "parent_element": parent_element,
            "boundary_type": boundary_type,
            "node_ids": node_ids,
            "body_id": body_id,
        }
        boundaries.append(row)
        boundary_by_id[condition_id].append(row)

    def triangle_area(node_ids: tuple[int, ...]) -> float:
        if len(node_ids) < 3:
            return 0.0
        a, b, c = (nodes[node_ids[i]] for i in range(3))
        ab = tuple(b[i] - a[i] for i in range(3))
        ac = tuple(c[i] - a[i] for i in range(3))
        cross = (
            ab[1] * ac[2] - ab[2] * ac[1],
            ab[2] * ac[0] - ab[0] * ac[2],
            ab[0] * ac[1] - ab[1] * ac[0],
        )
        area = 0.5 * math.sqrt(sum(value * value for value in cross))
        if len(node_ids) == 4:
            a2, b2, c2 = (nodes[node_ids[i]] for i in (0, 2, 3))
            ab2 = tuple(b2[i] - a2[i] for i in range(3))
            ac2 = tuple(c2[i] - a2[i] for i in range(3))
            cross2 = (
                ab2[1] * ac2[2] - ab2[2] * ac2[1],
                ab2[2] * ac2[0] - ab2[0] * ac2[2],
                ab2[0] * ac2[1] - ab2[1] * ac2[0],
            )
            area += 0.5 * math.sqrt(sum(value * value for value in cross2))
        return area

    def element_volume(element_type: int, node_ids: tuple[int, ...]) -> float:
        points = [nodes[node_id] for node_id in node_ids]
        if element_type == 504:
            a, b, c, d = points
            mat = [[b[i] - a[i], c[i] - a[i], d[i] - a[i]] for i in range(3)]
            det = (
                mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
                - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
                + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0])
            )
            return abs(det) / 6.0
        if element_type == 706:
            # Elmer prism ordering is two corresponding triangles.
            a, b, c, d, e, f = points
            v1 = tuple(b[i] - a[i] for i in range(3))
            v2 = tuple(c[i] - a[i] for i in range(3))
            normal = (
                v1[1] * v2[2] - v1[2] * v2[1],
                v1[2] * v2[0] - v1[0] * v2[2],
                v1[0] * v2[1] - v1[1] * v2[0],
            )
            height = tuple(((d[i] + e[i] + f[i]) - (a[i] + b[i] + c[i])) / 3.0 for i in range(3))
            return abs(sum(normal[i] * height[i] for i in range(3))) / 2.0
        return 0.0

    body_stats: dict[int, dict] = {}
    for body_id, node_ids in body_nodes.items():
        coords = [nodes[node_id] for node_id in node_ids]
        body_elements = [row for row in elements.values() if row[0] == body_id]
        body_stats[body_id] = {
            "body_id": body_id,
            "body_name": body_names.get(body_id, f"body_{body_id}"),
            "element_count": len(body_elements),
            "element_types": dict(element_counts[body_id]),
            "node_count": len(node_ids),
            "volume_m3": sum(element_volume(kind, ids) for _, kind, ids in body_elements),
            "bbox_m": {
                axis: [min(point[index] for point in coords), max(point[index] for point in coords)]
                for index, axis in enumerate(("x", "y", "z"))
            },
        }
        body_stats[body_id]["thickness_m"] = body_stats[body_id]["bbox_m"]["z"][1] - body_stats[body_id]["bbox_m"]["z"][0]
        body_stats[body_id]["footprint_bbox_area_m2"] = (
            body_stats[body_id]["bbox_m"]["x"][1] - body_stats[body_id]["bbox_m"]["x"][0]
        ) * (
            body_stats[body_id]["bbox_m"]["y"][1] - body_stats[body_id]["bbox_m"]["y"][0]
        )

    shared_node: dict[tuple[int, int], int] = {}
    shared_face: dict[tuple[int, int], int] = {}
    body_ids = sorted(body_nodes)
    for index, body_a in enumerate(body_ids):
        for body_b in body_ids[index + 1:]:
            shared_node[(body_a, body_b)] = len(body_nodes[body_a] & body_nodes[body_b])
    for owners in face_owners.values():
        unique_bodies = sorted({body_id for body_id, _ in owners})
        if len(unique_bodies) < 2:
            continue
        for index, body_a in enumerate(unique_bodies):
            for body_b in unique_bodies[index + 1:]:
                shared_face[(body_a, body_b)] = shared_face.get((body_a, body_b), 0) + 1

    boundary_stats = {}
    for boundary_id, rows in boundary_by_id.items():
        body_counter = Counter(row["body_id"] for row in rows)
        boundary_stats[boundary_id] = {
            "boundary_id": boundary_id,
            "boundary_name": boundary_names.get(boundary_id, f"boundary_{boundary_id}"),
            "facet_count": len(rows),
            "physical_area_m2": sum(triangle_area(row["node_ids"]) for row in rows),
            "body_ids": dict(body_counter),
        }

    return {
        "mesh": str(mesh),
        "body_names": body_names,
        "boundary_names": boundary_names,
        "body_stats": body_stats,
        "shared_node_counts": {f"{a}:{b}": count for (a, b), count in shared_node.items()},
        "conforming_face_counts": {f"{a}:{b}": count for (a, b), count in shared_face.items()},
        "boundary_stats": boundary_stats,
        "element_count": len(elements),
        "node_count": len(nodes),
        "boundary_facet_count": len(boundaries),
    }


def key_pair(a: int, b: int) -> str:
    return f"{min(a, b)}:{max(a, b)}"


def materialize_graph(label: str, mesh_data: dict, sif_bodies: dict[int, dict], mortar: list[dict]) -> tuple[list[dict], dict]:
    body_stats = mesh_data["body_stats"]
    by_name = {row["body_name"]: row for row in body_stats.values()}
    body_id_by_name = {row["body_name"]: int(body_id) for body_id, row in body_stats.items()}
    sif_by_target = {int(body_id): row for body_id, row in sif_bodies.items()}
    boundary_stats = mesh_data["boundary_stats"]
    boundary_id_by_name = {row["boundary_name"]: int(boundary_id) for boundary_id, row in boundary_stats.items()}

    contacts = [
        ("Membrane_SiNx", "TES", ("Membrane_SiNx__zmax", "TES__zmin"), "TES bottom mortar"),
        ("TES", "Stycast", ("TES__zmax", "Stycast__zmin"), "TES top mortar"),
        ("Stycast", "abs", ("Stycast__zmax", "abs__zmin"), "Pb bottom mortar"),
    ]
    rows: list[dict] = []
    for name_a, name_b, boundary_names, role in contacts:
        if name_a not in body_id_by_name or name_b not in body_id_by_name:
            continue
        body_a, body_b = body_id_by_name[name_a], body_id_by_name[name_b]
        pair_key = key_pair(body_a, body_b)
        resolved_boundary_names = list(boundary_names)
        if resolved_boundary_names[0] not in boundary_id_by_name and name_a == "Membrane_SiNx":
            alternatives = sorted(
                name for name in boundary_id_by_name
                if name.startswith("Membrane_SiNx") and name.endswith("__zmax")
                and "free" not in name
            )
            if alternatives:
                resolved_boundary_names[0] = alternatives[0]
        ids = [boundary_id_by_name[name] for name in resolved_boundary_names if name in boundary_id_by_name]
        areas = [boundary_stats[boundary_id]["physical_area_m2"] for boundary_id in ids]
        mortar_rows = [entry for entry in mortar if entry["slave_boundary_id"] in ids]
        coupling = "mortar" if mortar_rows else ("shared-node" if mesh_data["shared_node_counts"].get(pair_key, 0) else "nonconforming-unclassified")
        material_a = sif_by_target.get(body_a, {}).get("material", "")
        material_b = sif_by_target.get(body_b, {}).get("material", "")
        master_area = boundary_stats[ids[0]]["physical_area_m2"] if ids else 0.0
        slave_area = boundary_stats[ids[1]]["physical_area_m2"] if len(ids) > 1 else 0.0
        rows.append({
            "case": label,
            "route": f"{name_a} -> {name_b}",
            "body_a_id": body_a,
            "body_a_name": name_a,
            "body_a_material": material_a,
            "body_b_id": body_b,
            "body_b_name": name_b,
            "body_b_material": material_b,
            "interface_boundary_ids": ";".join(str(value) for value in ids),
            "physical_area_m2": min(areas) if areas else 0.0,
            "master_boundary_area_m2": master_area,
            "slave_boundary_area_m2": slave_area,
            "master_to_slave_area_ratio": master_area / slave_area if slave_area else 0.0,
            "shared_node_count": mesh_data["shared_node_counts"].get(pair_key, 0),
            "conforming_face_count": mesh_data["conforming_face_counts"].get(pair_key, 0),
            "mortar_constraint_count": "from_capture" if mortar_rows else 0,
            "coupling_type": coupling,
            "master": resolved_boundary_names[0] if ids and len(ids) > 1 else "unknown",
            "slave": resolved_boundary_names[1] if ids and len(ids) > 1 else "unknown",
            "physical_role": role,
        })
    # Add all actual shared-node body edges, including the substrate stack.
    for pair_key, count in mesh_data["shared_node_counts"].items():
        if count <= 0:
            continue
        body_a, body_b = (int(value) for value in pair_key.split(":"))
        if any({row["body_a_id"], row["body_b_id"]} == {body_a, body_b} for row in rows):
            continue
        name_a = body_stats[body_a]["body_name"]
        name_b = body_stats[body_b]["body_name"]
        rows.append({
            "case": label,
            "route": f"{name_a} -> {name_b}",
            "body_a_id": body_a,
            "body_a_name": name_a,
            "body_a_material": sif_by_target.get(body_a, {}).get("material", ""),
            "body_b_id": body_b,
            "body_b_name": name_b,
            "body_b_material": sif_by_target.get(body_b, {}).get("material", ""),
            "interface_boundary_ids": "",
            "physical_area_m2": 0.0,
            "master_boundary_area_m2": 0.0,
            "slave_boundary_area_m2": 0.0,
            "master_to_slave_area_ratio": 0.0,
            "shared_node_count": count,
            "conforming_face_count": mesh_data["conforming_face_counts"].get(pair_key, 0),
            "mortar_constraint_count": 0,
            "coupling_type": "shared-node",
            "master": "",
            "slave": "",
            "physical_role": "internal conforming body adjacency",
        })
    return rows, {name: body_stats[body_id_by_name[name]] for name in body_id_by_name}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--case", action="append", nargs=3, metavar=("LABEL", "MESH", "SIF"), required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict] = []
    cases: dict[str, dict] = {}
    for label, mesh_text, sif_text in args.case:
        mesh_path, sif_path = Path(mesh_text), Path(sif_text)
        sif_bodies, materials, mortar, bath_boundaries = parse_sif(sif_path)
        data = parse_mesh(mesh_path)
        rows, named_stats = materialize_graph(label, data, sif_bodies, mortar)
        for row in rows:
            if row["mortar_constraint_count"] == "from_capture":
                row["mortar_constraint_count"] = None
        cases[label] = {
            "mesh": data,
            "sif": {"bodies": sif_bodies, "materials": materials, "mortar_blocks": mortar, "bath_boundary_ids": bath_boundaries},
            "named_body_stats": named_stats,
            "graph_rows": rows,
            "bath_audit": [data["boundary_stats"].get(str(boundary_id), data["boundary_stats"].get(boundary_id, {})) for boundary_id in bath_boundaries],
        }
        all_rows.extend(rows)
    (args.out / "actual_body_connectivity.json").write_text(json.dumps(cases, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    fields = list(all_rows[0]) if all_rows else []
    with (args.out / "actual_body_connectivity.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(all_rows)
    material_rows: list[dict] = []
    for label, value in cases.items():
        stats_by_name = {row["body_name"]: row for row in value["mesh"]["body_stats"].values()}
        for body_id, body in value["sif"]["bodies"].items():
            material = value["sif"]["materials"].get(body["material_id"], {})
            stats = stats_by_name.get(body["body_name"], {})
            material_rows.append({
                "case": label,
                "body_id": body_id,
                "body_name": body["body_name"],
                "material_id": body["material_id"],
                "material_name": material.get("name", ""),
                "conductivity": material.get("conductivity", ""),
                "density": material.get("density", ""),
                "cp": material.get("cp", ""),
                "element_count": stats.get("element_count", ""),
                "volume_m3": stats.get("volume_m3", ""),
            })
    with (args.out / "material_network_comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(material_rows[0]))
        writer.writeheader()
        writer.writerows(material_rows)
    print(json.dumps({label: {"nodes": value["mesh"]["node_count"], "elements": value["mesh"]["element_count"], "bodies": len(value["mesh"]["body_stats"])} for label, value in cases.items()}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
