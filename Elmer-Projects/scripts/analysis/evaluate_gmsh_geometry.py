"""Compare ideal, pre-ElmerGrid Gmsh, and Elmer FEM volumes."""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from evaluate_physical_parity import expected_geometry, mesh_data, tetra_volume


def occ_entities(path: Path) -> list[dict[str, object]]:
    try:
        import gmsh
    except ImportError:
        return []
    gmsh.initialize(["-"])
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(str(path))
        entries = []
        for _, tag in gmsh.model.getEntities(3):
            entries.append(
                {
                    "tag": tag,
                    "volume_m3": gmsh.model.occ.getMass(3, tag),
                    "bounding_box_m": list(gmsh.model.getBoundingBox(3, tag)),
                }
            )
        return entries
    finally:
        gmsh.finalize()


def occ_body_volume(name: str, entries: list[dict[str, object]], ideal: float) -> float | None:
    # The current single-pixel OCC tree has one volume each for absorber and
    # Stycast.  TES is split into the Au and Ti films; identify that pair by
    # its thin z-range and sum it.  Other stack bodies may be split into many
    # boolean fragments and remain explicitly unresolved rather than guessed.
    if name in {"abs", "Stycast"}:
        candidates = [float(entry["volume_m3"]) for entry in entries]
        return min(candidates, key=lambda value: abs(value - ideal), default=None)
    if name == "TES":
        candidates = [
            float(entry["volume_m3"])
            for entry in entries
            if float(entry["bounding_box_m"][2]) > 1.918e-4
            and float(entry["bounding_box_m"][5]) < 1.923e-4
        ]
        if candidates and abs(sum(candidates) - ideal) / ideal < 1.0e-8:
            return sum(candidates)
    return None


def gmsh_tetra_volumes(path: Path) -> dict[int, dict[str, float | int]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    nodes: dict[int, tuple[float, float, float]] = {}
    elements: list[tuple[int, tuple[int, ...]]] = []
    index = lines.index("$Nodes") + 1
    count = int(lines[index])
    for line in lines[index + 1:index + 1 + count]:
        fields = line.split()
        nodes[int(fields[0])] = tuple(float(value) for value in fields[1:4])
    index = lines.index("$Elements") + 1
    count = int(lines[index])
    for line in lines[index + 1:index + 1 + count]:
        fields = line.split()
        if int(fields[1]) != 4:
            continue
        tag_count = int(fields[2])
        tags = [int(value) for value in fields[3:3 + tag_count]]
        conn = tuple(int(value) for value in fields[3 + tag_count:])
        elements.append((tags[0], conn))
    result: dict[int, dict[str, float | int]] = defaultdict(lambda: {"volume_m3": 0.0, "element_count": 0})
    for physical_id, conn in elements:
        volume = tetra_volume([nodes[node] for node in conn])
        result[physical_id]["volume_m3"] += volume
        result[physical_id]["element_count"] += 1
    return dict(result)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--gmsh", type=Path, required=True)
    parser.add_argument("--brep", type=Path)
    parser.add_argument("--elmer-mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    expected, _ = expected_geometry(args.project)
    gmsh = gmsh_tetra_volumes(args.gmsh)
    occ = occ_entities(args.brep) if args.brep else []
    bodies, _, nodes, _, elements_by_body, _, _ = mesh_data(args.elmer_mesh)
    reverse = {value: key for key, value in bodies.items()}
    elmer = {}
    for body_id, elements in elements_by_body.items():
        volume = 0.0
        for conn in elements:
            volume += tetra_volume([nodes[node] for node in conn])
        elmer[reverse[body_id]] = {"volume_m3": volume, "element_count": len(elements)}
    physical_ids = {"abs": 100, "TES": 101, "Stycast": 102}
    bodies_report = {}
    for name, ideal in expected.items():
        gmsh_entry = gmsh.get(physical_ids.get(name, -1), {})
        elmer_entry = elmer.get(name, {})
        bodies_report[name] = {
            "ideal_volume_m3": ideal,
            "gmsh_pre_elmer_volume_m3": gmsh_entry.get("volume_m3"),
            "elmer_fem_volume_m3": elmer_entry.get("volume_m3"),
            "gmsh_element_count": gmsh_entry.get("element_count"),
            "elmer_element_count": elmer_entry.get("element_count"),
            "occ_kernel_volume_m3": occ_body_volume(name, occ, ideal),
            "gmsh_vs_ideal_relative_error": ((float(gmsh_entry["volume_m3"]) - ideal) / ideal) if gmsh_entry else None,
            "elmer_vs_ideal_relative_error": ((float(elmer_entry["volume_m3"]) - ideal) / ideal) if elmer_entry else None,
            "elmer_vs_gmsh_relative_error": ((float(elmer_entry["volume_m3"]) - float(gmsh_entry["volume_m3"])) / float(gmsh_entry["volume_m3"])) if gmsh_entry and elmer_entry else None,
        }
    report = {
        "mesh": str(args.elmer_mesh.resolve()),
        "gmsh": str(args.gmsh.resolve()),
        "definition": "Gmsh pre-ElmerGrid tetra integration is the CAD-discretized proxy; OCC kernel volumes are read from the generated BREP when identifiable.",
        "occ_entities": occ,
        "bodies": bodies_report,
        "status": "REPORT_ONLY",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
