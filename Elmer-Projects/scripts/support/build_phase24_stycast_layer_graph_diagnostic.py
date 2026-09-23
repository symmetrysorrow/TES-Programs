#!/usr/bin/env python3
"""Materialize and audit the nominal 32-layer Stycast connectivity.

This is a topology-only diagnostic. It copies the already converted Phase24
mesh, changes the body label of the elements that belong to Stycast body 102,
and writes no production mesh or physics input. The purpose is to expose
shortcuts, non-neighbor shared-node links, and unintended external faces.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from collections import defaultdict
from pathlib import Path


SOURCE_DEFAULT = Path("work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar")
TARGET_DEFAULT = Path("work/meshes/mesh_phase24_stycast32_explicit_layers_diag")
ARTIFACT_DEFAULT = Path("artifacts/phase24_stycast_layer_graph_diagnostic")
STYCAST_BODY = 102
LAYER_FIRST_BODY = 200
LAYER_COUNT = 32


def read_nodes(path: Path) -> dict[int, tuple[float, float, float]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    nodes: dict[int, tuple[float, float, float]] = {}
    for line in lines[1:]:
        fields = line.split()
        if not fields:
            continue
        # Elmer mesh.nodes: node_id, partition/marker, x, y, z.
        nodes[int(fields[0])] = tuple(map(float, fields[2:5]))
    return nodes


def read_elements(path: Path) -> tuple[str, list[dict]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    records = []
    for line in lines[1:]:
        fields = line.split()
        if not fields:
            continue
        records.append(
            {
                "id": int(fields[0]),
                "body": int(fields[1]),
                "type": int(fields[2]),
                "nodes": [int(value) for value in fields[3:]],
                "line": line,
            }
        )
    return lines[0], records


def read_names(path: Path) -> tuple[dict[int, str], dict[int, str]]:
    bodies: dict[int, str] = {}
    boundaries: dict[int, str] = {}
    section = ""
    pattern = re.compile(r"^\$\s*(.*?)\s*=\s*(-?\d+)\s*$")
    for line in path.read_text(encoding="utf-8").splitlines():
        lower = line.lower()
        if "names for bodies" in lower:
            section = "bodies"
            continue
        if "names for boundaries" in lower:
            section = "boundaries"
            continue
        match = pattern.match(line)
        if not match:
            continue
        name, number = match.group(1), int(match.group(2))
        if section == "bodies":
            bodies[number] = name
        elif section == "boundaries":
            boundaries[number] = name
    return bodies, boundaries


def layer_for_z(z: float, zmin: float, dz: float) -> int:
    layer = int((z - zmin) / dz)
    return min(LAYER_COUNT - 1, max(0, layer))


def csv_write(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE_DEFAULT)
    parser.add_argument("--target", type=Path, default=TARGET_DEFAULT)
    parser.add_argument("--artifacts", type=Path, default=ARTIFACT_DEFAULT)
    args = parser.parse_args()

    source = args.source
    target = args.target
    artifacts = args.artifacts
    target.mkdir(parents=True, exist_ok=True)
    artifacts.mkdir(parents=True, exist_ok=True)

    nodes = read_nodes(source / "mesh.nodes")
    element_header, elements = read_elements(source / "mesh.elements")
    bodies, boundaries = read_names(source / "mesh.names")
    stycast_elements = [record for record in elements if record["body"] == STYCAST_BODY]
    if not stycast_elements:
        raise RuntimeError("no Stycast elements found in source mesh")

    stycast_nodes = {
        node_id
        for record in stycast_elements
        for node_id in record["nodes"]
    }
    z_values = [nodes[node_id][2] for node_id in stycast_nodes]
    zmin, zmax = min(z_values), max(z_values)
    dz = (zmax - zmin) / LAYER_COUNT
    tolerance = max(1.0e-12, dz * 1.0e-8)

    layer_records: dict[int, list[dict]] = defaultdict(list)
    crossing_elements = []
    for record in stycast_elements:
        element_z = [nodes[node_id][2] for node_id in record["nodes"]]
        centroid = sum(element_z) / len(element_z)
        layer = layer_for_z(centroid, zmin, dz)
        lower = zmin + layer * dz
        upper = zmin + (layer + 1) * dz
        if min(element_z) < lower - tolerance or max(element_z) > upper + tolerance:
            crossing_elements.append(
                {
                    "element_id": record["id"],
                    "layer_by_centroid": layer + 1,
                    "z_min": min(element_z),
                    "z_max": max(element_z),
                }
            )
        record["layer"] = layer
        layer_records[layer].append(record)

    # Copy the mesh files before rewriting only mesh.elements and mesh.names.
    for name in ("mesh.nodes", "mesh.boundary", "mesh.header"):
        shutil.copy2(source / name, target / name)
    original_names = source.joinpath("mesh.names").read_text(encoding="utf-8")
    shutil.copy2(source / "mesh.names", target / "mesh.names")

    element_by_id = {record["id"]: record for record in elements}
    with (target / "mesh.elements").open("w", encoding="utf-8", newline="") as handle:
        handle.write(element_header + "\n")
        for record in elements:
            if record["body"] == STYCAST_BODY:
                body = LAYER_FIRST_BODY + record["layer"]
                fields = [
                    str(record["id"]),
                    str(body),
                    str(record["type"]),
                    *map(str, record["nodes"]),
                ]
                handle.write(" ".join(fields) + "\n")
            else:
                handle.write(record["line"] + "\n")

    body_lines = ["! ----- names for bodies -----", "$ abs = 100"]
    for layer in range(LAYER_COUNT):
        body_lines.append(
            f"$ Stycast_layer_{layer + 1:02d} = {LAYER_FIRST_BODY + layer}"
        )
    for body_id, name in sorted(bodies.items()):
        if body_id not in (100, STYCAST_BODY):
            body_lines.append(f"$ {name} = {body_id}")
    boundary_lines = ["! ----- names for boundaries -----"]
    for boundary_id, name in sorted(boundaries.items()):
        boundary_lines.append(f"$ {name} = {boundary_id}")
    (target / "mesh.names").write_text(
        "\n".join(body_lines + boundary_lines) + "\n",
        encoding="utf-8",
    )

    layer_node_sets = {
        layer: {
            node_id
            for record in records
            for node_id in record["nodes"]
        }
        for layer, records in layer_records.items()
    }
    shared_edges = []
    element_pair_counts: dict[tuple[int, int], int] = defaultdict(int)
    node_to_layer_elements: dict[int, dict[int, list[int]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for layer, records in layer_records.items():
        for record in records:
            for node_id in set(record["nodes"]):
                node_to_layer_elements[node_id][layer].append(record["id"])
    for layer_map in node_to_layer_elements.values():
        layers_at_node = sorted(layer_map)
        for left_index, left in enumerate(layers_at_node):
            for right in layers_at_node[left_index + 1:]:
                element_pair_counts[(left, right)] += (
                    len(layer_map[left]) * len(layer_map[right])
                )

    element_edges = []
    for left in range(LAYER_COUNT):
        for right in range(left + 1, LAYER_COUNT):
            shared_nodes = layer_node_sets[left] & layer_node_sets[right]
            if shared_nodes:
                shared_edges.append(
                    {
                        "layer_a": left + 1,
                        "layer_b": right + 1,
                        "layer_gap": right - left,
                        "shared_node_count": len(shared_nodes),
                        "neighbor_expected": right == left + 1,
                        "edge_class": "neighbor" if right == left + 1 else "skip",
                    }
                )
            connected_elements = element_pair_counts.get((left, right), 0)
            if connected_elements:
                element_edges.append(
                    {
                        "layer_a": left + 1,
                        "layer_b": right + 1,
                        "layer_gap": right - left,
                        "element_pair_count": connected_elements,
                        "neighbor_expected": right == left + 1,
                        "edge_class": "neighbor" if right == left + 1 else "skip",
                    }
                )

    body_node_sets: dict[int, set[int]] = defaultdict(set)
    for record in elements:
        body_node_sets[record["body"]].update(record["nodes"])
    cross_body_rows = []
    for layer in range(LAYER_COUNT):
        for body_id, body_nodes in sorted(body_node_sets.items()):
            if body_id == STYCAST_BODY:
                continue
            shared_count = len(layer_node_sets[layer] & body_nodes)
            if shared_count:
                cross_body_rows.append(
                    {
                        "layer": layer + 1,
                        "layer_body_id": LAYER_FIRST_BODY + layer,
                        "other_body_id": body_id,
                        "other_body_name": bodies.get(body_id, f"body_{body_id}"),
                        "shared_node_count": shared_count,
                    }
                )

    boundary_rows = []
    boundary_path = source / "mesh.boundary"
    for line in boundary_path.read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if len(fields) < 6:
            continue
        boundary_id = int(fields[1])
        parent_element = int(fields[2])
        parent = element_by_id.get(parent_element)
        if not parent or parent["body"] != STYCAST_BODY:
            continue
        layer = parent["layer"]
        name = boundaries.get(boundary_id, f"boundary_{boundary_id}")
        allowed = name in {"Stycast__zmin", "Stycast__zmax"}
        boundary_rows.append(
            {
                "boundary_id": boundary_id,
                "boundary_name": name,
                "parent_element": parent_element,
                "layer": layer + 1,
                "allowed_outer_face": allowed,
                "classification": "allowed_outer_face" if allowed else "unexpected_lateral_or_other",
            }
        )

    layer_rows = []
    for layer in range(LAYER_COUNT):
        layer_zmin = zmin + layer * dz
        layer_zmax = zmin + (layer + 1) * dz
        layer_rows.append(
            {
                "layer": layer + 1,
                "body_id": LAYER_FIRST_BODY + layer,
                "z_min_nominal": layer_zmin,
                "z_max_nominal": layer_zmax,
                "element_count": len(layer_records[layer]),
                "node_count": len(layer_node_sets[layer]),
                "element_z_min": min(
                    nodes[node_id][2]
                    for record in layer_records[layer]
                    for node_id in record["nodes"]
                ),
                "element_z_max": max(
                    nodes[node_id][2]
                    for record in layer_records[layer]
                    for node_id in record["nodes"]
                ),
            }
        )

    skip_shared = [row for row in shared_edges if row["edge_class"] == "skip"]
    skip_elements = [row for row in element_edges if row["edge_class"] == "skip"]
    unexpected_boundary_rows = [
        row for row in boundary_rows if not row["allowed_outer_face"]
    ]
    all_layers_present = all(layer_records[layer] for layer in range(LAYER_COUNT))
    series_graph = (
        all_layers_present
        and not crossing_elements
        and not skip_shared
        and not skip_elements
        and not unexpected_boundary_rows
    )

    audit = {
        "source_mesh": str(source),
        "diagnostic_mesh": str(target),
        "stycast_source_body": STYCAST_BODY,
        "layer_count": LAYER_COUNT,
        "diagnostic_body_id_range": [
            LAYER_FIRST_BODY,
            LAYER_FIRST_BODY + LAYER_COUNT - 1,
        ],
        "stycast_element_count": len(stycast_elements),
        "stycast_node_count": len(stycast_nodes),
        "z_min": zmin,
        "z_max": zmax,
        "nominal_layer_thickness": dz,
        "crossing_element_count": len(crossing_elements),
        "crossing_elements": crossing_elements[:100],
        "all_layers_present": all_layers_present,
        "shared_node_edges": shared_edges,
        "element_connectivity_edges": element_edges,
        "cross_body_shared_nodes": cross_body_rows,
        "boundary_face_count_on_stycast_parent_elements": len(boundary_rows),
        "unexpected_lateral_or_other_boundary_count": len(unexpected_boundary_rows),
        "unexpected_lateral_or_other_boundaries": unexpected_boundary_rows[:100],
        "series_graph": series_graph,
    }
    (artifacts / "stycast_layer_graph_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    csv_write(
        artifacts / "stycast_layer_bounds.csv",
        list(layer_rows[0]),
        layer_rows,
    )
    csv_write(
        artifacts / "stycast_layer_shared_edges.csv",
        list(shared_edges[0]) if shared_edges else [
            "layer_a",
            "layer_b",
            "layer_gap",
            "shared_node_count",
            "neighbor_expected",
            "edge_class",
        ],
        shared_edges,
    )
    csv_write(
        artifacts / "stycast_layer_element_edges.csv",
        list(element_edges[0]) if element_edges else [
            "layer_a",
            "layer_b",
            "layer_gap",
            "element_pair_count",
            "neighbor_expected",
            "edge_class",
        ],
        element_edges,
    )
    csv_write(
        artifacts / "stycast_layer_lateral_faces.csv",
        list(boundary_rows[0]) if boundary_rows else [
            "boundary_id",
            "boundary_name",
            "parent_element",
            "layer",
            "allowed_outer_face",
            "classification",
        ],
        boundary_rows,
    )
    csv_write(
        artifacts / "stycast_layer_cross_body_shared_nodes.csv",
        list(cross_body_rows[0]) if cross_body_rows else [
            "layer",
            "layer_body_id",
            "other_body_id",
            "other_body_name",
            "shared_node_count",
        ],
        cross_body_rows,
    )

    summary = [
        "# Phase24 Stycast 32-layer graph diagnostic",
        "",
        f"- Source mesh: {source}",
        f"- Diagnostic mesh: {target}",
        f"- Stycast elements relabelled: {len(stycast_elements):,}",
        f"- Stycast nodes: {len(stycast_nodes):,}",
        f"- z range: {zmin:.16g} .. {zmax:.16g} m",
        f"- Nominal layer thickness: {dz:.16g} m",
        f"- Element-plane crossings: {len(crossing_elements)} "
        "(nonzero means existing elements are not layer-resolved)",
        f"- All 32 layers populated: {all_layers_present}",
        f"- Shared-node edges: {len(shared_edges)} total, {len(skip_shared)} non-neighbor skip edges",
        f"- Element-connectivity edges: {len(element_edges)} total, {len(skip_elements)} non-neighbor skip edges",
        f"- Cross-body shared-node rows: {len(cross_body_rows)}",
        f"- Stycast-parent external faces: {len(boundary_rows)}",
        f"- Unexpected lateral/other faces: {len(unexpected_boundary_rows)}",
        "",
        "## Interpretation",
        "",
        (
            "The layer graph is a clean series graph under this diagnostic relabelling."
            if series_graph
            else "The existing mesh is not a clean explicit-layer series graph; "
            "inspect crossing elements, edge CSVs, and boundary CSVs."
        ),
        "",
        "The copied diagnostic mesh is topology-only; no material or solver input was generated.",
        "",
        "Source mesh.names was preserved in memory only for name discovery: "
        + str(bool(original_names)),
    ]
    (artifacts / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")

    # Keep the generated diagnostic mesh self-describing without changing any
    # solver setup. This file is intentionally not used as a production SIF.
    (target / "entities.sif").write_text(
        "! Diagnostic only: body IDs 200..231 are Stycast layers.\n"
        "! No material, equation, solver, or boundary conditions are defined.\n",
        encoding="utf-8",
    )

    print(json.dumps(
        {
            "series_graph": series_graph,
            "stycast_elements": len(stycast_elements),
            "stycast_nodes": len(stycast_nodes),
            "shared_edges": len(shared_edges),
            "skip_shared_edges": len(skip_shared),
            "element_edges": len(element_edges),
            "skip_element_edges": len(skip_elements),
            "cross_body_rows": len(cross_body_rows),
            "stycast_boundary_faces": len(boundary_rows),
            "unexpected_boundary_faces": len(unexpected_boundary_rows),
            "target": str(target),
            "artifacts": str(artifacts),
        }, indent=2
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
