"""Remove Elmer boundary facets that are internal to a conforming mesh.

Elmer's no-mortar path must not expose a conforming body interface as an
ordinary boundary on both sides: doing so can leave the heat equation with a
zero-flux boundary instead of using the shared volume-face connectivity.
This post-process keeps all nodes and volume elements unchanged and removes
only boundary facets whose node set is owned by volume elements from two
different bodies.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path


FACE_NODES = {
    504: ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
    706: ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
}


def _volume_face_owners(path: Path) -> set[tuple[int, ...]]:
    owners: defaultdict[tuple[int, ...], set[int]] = defaultdict(set)
    for line in (path / "mesh.elements").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 4:
            continue
        body, element_type = int(fields[1]), int(fields[2])
        local_faces = FACE_NODES.get(element_type)
        if local_faces is None:
            continue
        nodes = tuple(map(int, fields[3:]))
        for face in local_faces:
            owners[tuple(sorted(nodes[index] for index in face))].add(body)
    return {face for face, body_ids in owners.items() if len(body_ids) >= 2}


def strip(mesh_dir: Path) -> dict[str, int]:
    internal_faces = _volume_face_owners(mesh_dir)
    source = mesh_dir / "mesh.boundary"
    lines = source.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"invalid Elmer boundary file: {source}")
    kept = []
    removed = 0
    for line in lines:
        fields = line.split()
        if len(fields) >= 6:
            node_ids = tuple(sorted(map(int, fields[5:])))
            if node_ids in internal_faces:
                removed += 1
                continue
        kept.append(line)
    boundary_count = len(kept)
    boundary_type_counts = Counter(
        int(fields[4]) for line in kept if (fields := line.split()) and len(fields) >= 5
    )
    # ElmerGrid writes a headerless mesh.boundary; the element count is stored
    # in mesh.header.  Keep the file as a pure boundary-element list.
    source.write_text("\n".join(kept) + "\n", encoding="utf-8")

    header_path = mesh_dir / "mesh.header"
    header_lines = header_path.read_text(encoding="utf-8").splitlines()
    counts = header_lines[0].split()
    if len(counts) >= 3:
        counts[2] = str(boundary_count)
        header_lines[0] = "  ".join(counts)
        for index, line in enumerate(header_lines[2:], start=2):
            fields = line.split()
            if len(fields) >= 2 and fields[0].isdigit():
                fields[1] = str(boundary_type_counts.get(int(fields[0]), 0))
                header_lines[index] = "  ".join(fields)
        header_path.write_text("\n".join(header_lines) + "\n", encoding="utf-8")
    return {"internal_faces": len(internal_faces), "removed_boundary_facets": removed,
            "remaining_boundary_facets": boundary_count}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh_dir", type=Path)
    args = parser.parse_args()
    result = strip(args.mesh_dir)
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
