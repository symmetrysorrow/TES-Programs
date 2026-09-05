"""Derive a non-conforming Mortar control mesh from one conformal Elmer mesh.

Only node IDs on the right-hand body of each contact are duplicated. Element
connectivity, coordinates, body volumes, surface geometry, and element counts
remain unchanged; the resulting mesh differs only by the interface coupling
topology.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.check_conformal_interfaces import (
    DEFAULT_PAIRS,
    read_boundary,
    read_names,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_mesh(mesh: Path):
    nodes = mesh.joinpath("mesh.nodes").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    elements = mesh.joinpath("mesh.elements").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    boundary = mesh.joinpath("mesh.boundary").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    return nodes, elements, boundary


def derive(source: Path, output: Path) -> dict[str, object]:
    bodies, boundaries = read_names(source / "mesh.names")
    boundary_sets = read_boundary(source / "mesh.boundary")
    nodes, elements, boundary = parse_mesh(source)

    body_by_element: dict[int, int] = {}
    element_fields: list[list[str]] = []
    for line in elements:
        fields = line.split()
        if len(fields) >= 5:
            element_fields.append(fields)
            body_by_element[int(fields[0])] = int(fields[1])

    duplicate_by_body: dict[int, set[int]] = {}
    interface_records: list[dict[str, object]] = []
    for left_name, right_name, label, left_body, right_body in DEFAULT_PAIRS:
        left_id = boundaries[left_name]
        right_id = boundaries[right_name]
        left_nodes = set().union(*boundary_sets[left_id])
        right_nodes = set().union(*boundary_sets[right_id])
        duplicate_by_body.setdefault(bodies[right_body], set()).update(right_nodes)
        interface_records.append(
            {
                "interface": label,
                "left_surface": left_name,
                "right_surface": right_name,
                "duplicated_body": right_body,
                "duplicated_nodes": len(right_nodes),
                "source_shared_nodes": len(left_nodes & right_nodes),
            }
        )

    max_node = max(int(line.split()[0]) for line in nodes if line.split())
    duplicate_id: dict[tuple[int, int], int] = {}
    output_nodes = list(nodes)
    node_lines = {int(line.split()[0]): line for line in nodes if line.split()}
    for body_id, node_ids in sorted(duplicate_by_body.items()):
        for node_id in sorted(node_ids):
            max_node += 1
            duplicate_id[(body_id, node_id)] = max_node
            fields = node_lines[node_id].split()
            fields[0] = str(max_node)
            output_nodes.append(" ".join(fields))

    def remap(body_id: int, node_id: int) -> int:
        return duplicate_id.get((body_id, node_id), node_id)

    output_elements: list[str] = []
    for fields in element_fields:
        body_id = int(fields[1])
        fields[3:] = [str(remap(body_id, int(node))) for node in fields[3:]]
        output_elements.append(" ".join(fields))

    output_boundary: list[str] = []
    for line in boundary:
        fields = line.split()
        if len(fields) < 6:
            output_boundary.append(line)
            continue
        parent = int(fields[2])
        body_id = body_by_element.get(parent)
        if body_id is not None:
            fields[5:] = [str(remap(body_id, int(node))) for node in fields[5:]]
        output_boundary.append(" ".join(fields))

    output.mkdir(parents=True, exist_ok=True)
    for name in ("mesh.names", "entities.sif"):
        shutil.copy2(source / name, output / name)
    output.joinpath("mesh.nodes").write_text(
        chr(10).join(output_nodes) + chr(10), encoding="utf-8"
    )
    output.joinpath("mesh.elements").write_text(
        chr(10).join(output_elements) + chr(10), encoding="utf-8"
    )
    output.joinpath("mesh.boundary").write_text(
        chr(10).join(output_boundary) + chr(10), encoding="utf-8"
    )
    header_lines = source.joinpath("mesh.header").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()
    header_fields = header_lines[0].split()
    header_fields[0] = str(len(output_nodes))
    header_lines[0] = " ".join(header_fields)
    output.joinpath("mesh.header").write_text(
        chr(10).join(header_lines) + chr(10), encoding="utf-8"
    )
    provenance = {
        "route": "mortar_control",
        "source_mesh": str(source.resolve()),
        "source_mesh_sha256": {
            name: sha256(source / name)
            for name in ("mesh.nodes", "mesh.elements", "mesh.boundary", "mesh.header")
        },
        "method": "duplicate right-side interface node IDs; preserve all coordinates and volume connectivity",
        "interfaces": interface_records,
        "node_count": len(output_nodes),
        "volume_element_count": len(output_elements),
    }
    output.joinpath("PROVENANCE.json").write_text(
        json.dumps(provenance, indent=2) + chr(10), encoding="utf-8"
    )
    return provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = derive(args.source, args.output)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + chr(10), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
