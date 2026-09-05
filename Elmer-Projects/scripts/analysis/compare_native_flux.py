"""Compare Elmer FluxSolver output with the raw tetra-gradient diagnostic."""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.evaluate_physical_parity import (
    PAIR_DEFS,
    heat_flux_consistency,
    mesh_data,
)


def scalar_field(path: Path, label: str) -> dict[int, float]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    matches = [i for i, line in enumerate(lines) if line.strip().lower() == label.lower()]
    if not matches:
        raise ValueError(f"native result has no scalar field {label!r}")
    label_index = matches[-1]
    perm_index = next(i for i in range(label_index + 1, len(lines)) if lines[i].startswith("Perm:"))
    perm_fields = lines[perm_index].split()
    if len(perm_fields) >= 3 and perm_fields[1].lower() == "use" and perm_fields[2].lower() == "previous":
        temperature_label = next(i for i, line in enumerate(lines) if line.strip().lower() == "temperature")
        temperature_perm = next(i for i in range(temperature_label + 1, len(lines)) if lines[i].startswith("Perm:"))
        count = int(lines[temperature_perm].split()[1])
        permutation = [int(lines[temperature_perm + 1 + i].split()[0]) for i in range(count)]
        start = perm_index + 1
    else:
        count = int(perm_fields[1])
        permutation = [int(lines[perm_index + 1 + i].split()[0]) for i in range(count)]
        start = perm_index + 1 + count
    values = [float(lines[start + i].replace("D", "E")) for i in range(count)]
    return {permutation[position]: values[position] for position in range(count)}


def normal_and_area(points: list[tuple[float, float, float]], face: tuple[int, ...], conn: tuple[int, ...]) -> tuple[tuple[float, float, float], float]:
    indexes = tuple(conn.index(node) for node in face)
    a, b, c = (points[index] for index in indexes)
    ab = tuple(b[index] - a[index] for index in range(3))
    ac = tuple(c[index] - a[index] for index in range(3))
    raw = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    length = math.sqrt(sum(value * value for value in raw))
    face_centroid = tuple(sum(points[index][axis] for index in indexes) / 3.0 for axis in range(3))
    tetra_centroid = tuple(sum(point[axis] for point in points) / 4.0 for axis in range(3))
    sign = 1.0 if sum(raw[axis] * (face_centroid[axis] - tetra_centroid[axis]) for axis in range(3)) >= 0.0 else -1.0
    return tuple(sign * value / length for value in raw), 0.5 * length


def native_interface_flux(mesh: Path, native_result: Path) -> dict[str, dict[str, float]]:
    bodies, _, nodes, elements, _, faces, parents = mesh_data(mesh)
    del bodies
    if native_result.suffix.lower() == ".vtu":
        import meshio

        native_mesh = meshio.read(native_result)
        vector = native_mesh.point_data["temperature flux"]
        by_coordinate = {
            tuple(round(float(value), 12) for value in point): index + 1
            for index, point in enumerate(native_mesh.points)
        }
        flux_components = [dict() for _ in range(3)]
        for index, point in enumerate(native_mesh.points):
            node = by_coordinate[tuple(round(float(value), 12) for value in point)]
            for axis in range(3):
                flux_components[axis][node] = float(vector[index, axis])
        native_source = "VTK PointData 'temperature flux'"
    else:
        flux_components = [scalar_field(native_result, f"temperature flux {axis}") for axis in (1, 2, 3)]
        native_source = "Elmer ASCII result scalar fields"
    report = {}
    for left, right, label, *_ in PAIR_DEFS:
        sides = []
        for surface in (left, right):
            total = 0.0
            area = 0.0
            for face, parent in zip(faces[surface], parents[surface]):
                body_id, conn = elements[parent]
                del body_id
                points = [nodes[node] for node in conn]
                normal, facet_area = normal_and_area(points, face, conn)
                face_flux = tuple(
                    sum(flux_components[axis][node] for node in face) / len(face)
                    for axis in range(3)
                )
                # FluxSolver stores +k*grad; physical heat flux is its negative.
                total += -sum(face_flux[axis] * normal[axis] for axis in range(3)) * facet_area
                area += facet_area
            sides.append({"integrated_outward_flux_W": total, "surface_area_m2": area})
        q_left = sides[0]["integrated_outward_flux_W"]
        q_right = sides[1]["integrated_outward_flux_W"]
        report[label] = {
            "left": sides[0],
            "right": sides[1],
            "absolute_imbalance_W": abs(q_left + q_right),
            "normalized_imbalance": abs(q_left + q_right) / max(abs(q_left), abs(q_right), 1.0e-30),
        }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--raw-result", type=Path, required=True)
    parser.add_argument("--native-result", type=Path, required=True)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = {
        "native_solver": {
            "available": True,
            "implementation": "Elmer 26.1 FluxSolver; Calculate Flux + Calculate Grad",
            "stored_fields": ["temperature flux 1", "temperature flux 2", "temperature flux 3", "temperature grad 1", "temperature grad 2", "temperature grad 3"],
            "sign_conversion": "FluxSolver stores +k*grad; comparison converts to physical -k*grad.",
        },
        "raw_gradient": heat_flux_consistency(args.mesh, args.raw_result, args.project),
        "native_flux": native_interface_flux(args.mesh, args.native_result),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
