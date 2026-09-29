"""Make a diagnostic Elmer mesh nonconforming at the TES--membrane contact.

ElmerGrid's coordinate merge can collapse coincident contact nodes even when
the CAD surfaces were intended to be independent.  This post-process is
restricted to the isolated diagnostic mesh: it duplicates only nodes shared
by body 101 (TES) and body 103 (Membrane_SiNx), rewrites the membrane volume
and membrane boundary connectivity, and leaves coordinates/materials/
geometry unchanged.
"""
from __future__ import annotations

import argparse
from pathlib import Path


TES_BODY = 101
MEMBRANE_BODY = 103
MEMBRANE_BOUNDARIES = set(range(1300, 1307))


def fields(path: Path) -> list[list[str]]:
    return [line.split() for line in path.read_text(encoding="utf-8").splitlines()]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    args = parser.parse_args()
    mesh = args.mesh
    node_path = mesh / "mesh.nodes"
    element_path = mesh / "mesh.elements"
    boundary_path = mesh / "mesh.boundary"
    node_rows = fields(node_path)
    element_rows = fields(element_path)
    boundary_rows = fields(boundary_path)

    tes_nodes: set[int] = set()
    membrane_nodes: set[int] = set()
    for row in element_rows:
        if len(row) < 4:
            continue
        body = int(row[1])
        nodes = [int(value) for value in row[3:]]
        if body == TES_BODY:
            tes_nodes.update(nodes)
        elif body == MEMBRANE_BODY:
            membrane_nodes.update(nodes)
    shared = sorted(tes_nodes & membrane_nodes)
    if not shared:
        raise SystemExit("TES--membrane interface has no shared nodes")

    max_node = max(int(row[0]) for row in node_rows)
    duplicate = {old: max_node + index for index, old in enumerate(shared, start=1)}
    node_by_id = {int(row[0]): row for row in node_rows}
    node_rows.extend([[str(new), *node_by_id[old][1:]] for old, new in duplicate.items()])

    changed_elements = 0
    for row in element_rows:
        if len(row) >= 4 and int(row[1]) == MEMBRANE_BODY:
            for index in range(3, len(row)):
                row[index] = str(duplicate.get(int(row[index]), int(row[index])))
            changed_elements += 1

    changed_boundaries = 0
    for row in boundary_rows:
        if len(row) < 6 or int(row[1]) not in MEMBRANE_BOUNDARIES:
            continue
        code = int(row[4])
        count = {303: 3, 404: 4}.get(code, max(3, len(row) - 5))
        for index in range(5, min(5 + count, len(row))):
            row[index] = str(duplicate.get(int(row[index]), int(row[index])))
        changed_boundaries += 1

    node_path.write_text("\n".join(" ".join(row) for row in node_rows) + "\n", encoding="utf-8")
    element_path.write_text("\n".join(" ".join(row) for row in element_rows) + "\n", encoding="utf-8")
    boundary_path.write_text("\n".join(" ".join(row) for row in boundary_rows) + "\n", encoding="utf-8")
    header = mesh / "mesh.header"
    header_lines = header.read_text(encoding="utf-8").splitlines()
    counts = header_lines[0].split()
    counts[0] = str(len(node_rows))
    header_lines[0] = "  ".join(counts)
    header.write_text("\n".join(header_lines) + "\n", encoding="utf-8")
    print(f"duplicated {len(shared)} TES--membrane nodes; membrane elements={changed_elements}; boundary facets={changed_boundaries}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
