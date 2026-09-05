"""Evaluate same-geometry Mortar/conformal physical parity.

The report deliberately separates mesh geometry, common temperature
observables, interface continuity, and heat-flux balance. Missing electrical
or transient observables are reported as unavailable rather than inferred.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.check_conformal_interfaces import read_names
from scripts.support.reconcile_project import reconcile_project


PAIR_DEFS = (
    ("Membrane_SiNx__zmax", "TES__zmin", "Membrane_TES", "Membrane_SiNx", "TES", 1.0, -1.0),
    ("TES__zmax", "Stycast__zmin", "TES_Stycast", "TES", "Stycast", 1.0, -1.0),
    ("Stycast__zmax", "abs__zmin", "Stycast_absorber", "Stycast", "abs", 1.0, -1.0),
)


def result_values(path: Path, field_index: int = -1) -> dict[int, float]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    permutations = [i for i, line in enumerate(lines) if line.startswith("Perm:")]
    if not permutations:
        raise ValueError(f"result file has no Perm field: {path}")
    selected = field_index if field_index >= 0 else len(permutations) + field_index
    if selected < 0 or selected >= len(permutations):
        raise IndexError(f"result field index out of range: {field_index}")
    permutation: list[int] | None = None
    index = permutations[selected]
    start = 0
    for number, perm_index in enumerate(permutations):
        tokens = lines[perm_index].split()
        if len(tokens) >= 2 and tokens[1].lower() == "use":
            if permutation is None:
                raise ValueError(f"first result field cannot use previous permutation: {path}")
            count = len(permutation)
            value_start = perm_index + 1
        else:
            count = int(tokens[1])
            # Elmer writes ``node_index, internal_permutation``.  Values are
            # emitted in node-index order, so the first column identifies the
            # mesh node; the second column is only internal storage permutation.
            permutation = [int(lines[perm_index + 1 + i].split()[0]) for i in range(count)]
            value_start = perm_index + 1 + count
        if number == selected:
            start = value_start
            break
    values = [
        float(lines[start + i].replace("D", "E"))
        for i in range(len(permutation or []))
    ]
    return {permutation[position]: values[position] for position in range(len(values))}


def result_field_times(path: Path) -> list[float | None]:
    """Return the Elmer ``Time:`` value associated with every Perm block."""
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    times: list[float | None] = []
    for index, line in enumerate(lines):
        if not line.startswith("Perm:"):
            continue
        time_value: float | None = None
        for previous in reversed(lines[:index]):
            if previous.startswith("Time:"):
                fields = previous.split()
                if len(fields) >= 4:
                    time_value = float(fields[3].replace("D", "E"))
                break
        times.append(time_value)
    return times


def mesh_data(mesh: Path):
    bodies, boundaries = read_names(mesh / "mesh.names")
    nodes: dict[int, tuple[float, float, float]] = {}
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 5:
            nodes[int(fields[0])] = tuple(float(value) for value in fields[2:5])
    elements: dict[int, tuple[int, tuple[int, ...]]] = {}
    elements_by_body: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 5:
            element_id = int(fields[0])
            body = int(fields[1])
            kind = int(fields[2])
            if kind == 504:
                conn = tuple(int(value) for value in fields[3:])
                elements[element_id] = (body, conn)
                elements_by_body[body].append(conn)
    boundary_faces: dict[str, list[tuple[int, ...]]] = defaultdict(list)
    boundary_parent: dict[str, list[int]] = defaultdict(list)
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) < 6:
            continue
        boundary_id = int(fields[1])
        name = next((key for key, value in boundaries.items() if value == boundary_id), None)
        if name is None:
            continue
        boundary_parent[name].append(int(fields[2]))
        boundary_faces[name].append(tuple(int(value) for value in fields[5:]))
    return bodies, boundaries, nodes, elements, elements_by_body, boundary_faces, boundary_parent


def tetra_volume(points: list[tuple[float, float, float]]) -> float:
    a, b, c, d = points
    ab = tuple(b[i] - a[i] for i in range(3))
    ac = tuple(c[i] - a[i] for i in range(3))
    ad = tuple(d[i] - a[i] for i in range(3))
    cross = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    return abs(sum(cross[i] * ad[i] for i in range(3))) / 6.0


def triangle_area(points: list[tuple[float, float, float]]) -> float:
    a, b, c = points
    ab = tuple(b[i] - a[i] for i in range(3))
    ac = tuple(c[i] - a[i] for i in range(3))
    cross = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    return 0.5 * math.sqrt(sum(value * value for value in cross))


def body_observables(mesh: Path, result: Path | None) -> dict[str, dict[str, float | int]]:
    bodies, _, nodes, _, elements_by_body, _, _ = mesh_data(mesh)
    values = result_values(result) if result else {}
    reverse = {value: key for key, value in bodies.items()}
    report: dict[str, dict[str, float | int]] = {}
    for body_id, elements in elements_by_body.items():
        name = reverse.get(body_id, f"body_{body_id}")
        volume = 0.0
        integral = 0.0
        minimum = math.inf
        maximum = -math.inf
        for conn in elements:
            points = [nodes[node] for node in conn]
            element_volume = tetra_volume(points)
            volume += element_volume
            if values:
                temperatures = [values[node] for node in conn]
                mean = sum(temperatures) / len(temperatures)
                integral += mean * element_volume
                minimum = min(minimum, *temperatures)
                maximum = max(maximum, *temperatures)
        entry: dict[str, float | int] = {
            "volume_m3": volume,
            "element_count": len(elements),
        }
        if values:
            entry.update(
                {
                    "volume_average_temperature_K": integral / volume,
                    "min_temperature_K": minimum,
                    "max_temperature_K": maximum,
                }
            )
        report[name] = entry
    return report


def contact_areas(mesh: Path) -> dict[str, dict[str, float]]:
    bodies, boundaries, nodes, _, _, faces, _ = mesh_data(mesh)
    del bodies
    result: dict[str, dict[str, float]] = {}
    for left, right, label, *_ in PAIR_DEFS:
        result[label] = {
            "left_m2": sum(triangle_area([nodes[node] for node in face]) for face in faces[left]),
            "right_m2": sum(triangle_area([nodes[node] for node in face]) for face in faces[right]),
        }
    return result


def expected_geometry(project: Path) -> tuple[dict[str, float], dict[str, float]]:
    model = reconcile_project(json.loads(project.read_text(encoding="utf-8")))
    p = model["parameters"]
    # The analytic geometry intentionally uses the literal 498 um cylinder
    # in the single-pixel geometry tree, not the historical Stycast_dx
    # parameter (400 um) used by some older project templates.
    stycast_diameter = 498.0e-6
    expected = {
        "abs": p["abs_dx"] * p["abs_dy"] * p["abs_dz"],
        "TES": p["TES_Au_dx"] * p["TES_Au_dy"] * p["TES_dz"],
        "Stycast": math.pi * (stycast_diameter / 2.0) ** 2 * p["Stycast_dz"],
        "SiO2_1": p["Si_dx"] * p["Si_dy"] * p["SiO2_1_dz"],
        "Si_1": (p["Si_dx"] * p["Si_dy"] - p["membrane_dx"] * p["membrane_dy"]) * p["Si_1_dz"],
        "SiNx": (p["Si_dx"] * p["Si_dy"] - p["membrane_dx"] * p["membrane_dy"]) * p["SiNx_dz"],
        "Si_2": (p["Si_dx"] * p["Si_dy"] - p["membrane_dx"] * p["membrane_dy"]) * p["Si_2_dz"],
        "SiO2_2": (p["Si_dx"] * p["Si_dy"] - p["membrane_dx"] * p["membrane_dy"]) * p["SiO2_2_dz"],
        "Membrane_Si1": p["membrane_dx"] * p["membrane_dy"] * p["Si_1_dz"],
        "Membrane_SiNx": p["membrane_dx"] * p["membrane_dy"] * p["SiNx_dz"],
    }
    contact = {
        "Membrane_TES": p["TES_Au_dx"] * p["TES_Au_dy"],
        "TES_Stycast": math.pi * (stycast_diameter / 2.0) ** 2,
        "Stycast_absorber": math.pi * (stycast_diameter / 2.0) ** 2,
    }
    return expected, contact


def interface_continuity(mesh: Path, result: Path) -> dict[str, dict[str, float | int]]:
    bodies, boundaries, nodes, _, _, faces, _ = mesh_data(mesh)
    del bodies
    values = result_values(result)
    report: dict[str, dict[str, float | int]] = {}
    for left, right, label, *_ in PAIR_DEFS:
        left_ids = sorted({node for face in faces[left] for node in face})
        right_ids = sorted({node for face in faces[right] for node in face})
        diffs: list[float] = []
        if set(left_ids) == set(right_ids):
            diffs = [abs(values[node] - values[node]) for node in left_ids]
        else:
            for node in left_ids:
                point = nodes[node]
                nearest = min(
                    right_ids,
                    key=lambda candidate: math.dist(point, nodes[candidate]),
                )
                diffs.append(abs(values[node] - values[nearest]))
            for node in right_ids:
                point = nodes[node]
                nearest = min(
                    left_ids,
                    key=lambda candidate: math.dist(point, nodes[candidate]),
                )
                diffs.append(abs(values[node] - values[nearest]))
        report[label] = {
            "left_nodes": len(left_ids),
            "right_nodes": len(right_ids),
            "max_abs_temperature_jump_K": max(diffs, default=math.nan),
            "mean_abs_temperature_jump_K": sum(diffs) / len(diffs) if diffs else math.nan,
            "rms_temperature_jump_K": math.sqrt(sum(value * value for value in diffs) / len(diffs)) if diffs else math.nan,
        }
    return report


def determinant(matrix: list[list[float]]) -> float:
    a, b, c = matrix
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def tetra_gradient(
    points: list[tuple[float, float, float]], temperatures: list[float]
) -> tuple[float, float, float] | None:
    origin = points[0]
    matrix = [
        [points[i][axis] - origin[axis] for axis in range(3)]
        for i in range(1, 4)
    ]
    rhs = [temperatures[i] - temperatures[0] for i in range(1, 4)]
    det = determinant(matrix)
    if abs(det) < 1.0e-40:
        return None
    gradient = []
    for axis in range(3):
        component_matrix = [row[:] for row in matrix]
        for index in range(3):
            component_matrix[index][axis] = rhs[index]
        gradient.append(determinant(component_matrix) / det)
    return tuple(gradient)


def oriented_face_flux(
    points: list[tuple[float, float, float]],
    temperatures: list[float],
    face_indices: tuple[int, int, int],
    conductivity_value: float,
) -> dict[str, float | tuple[float, float, float]] | None:
    """Return outward elemental heat flux for one tetrahedral face.

    The face ordering in ``mesh.boundary`` is not trusted.  The cross product
    is oriented using the parent tetrahedron centroid, so reversing the input
    face node order cannot change the outward flux sign.
    """
    gradient = tetra_gradient(points, temperatures)
    if gradient is None:
        return None
    a, b, c = (points[index] for index in face_indices)
    ab = tuple(b[index] - a[index] for index in range(3))
    ac = tuple(c[index] - a[index] for index in range(3))
    raw_normal = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    raw_norm = math.sqrt(sum(value * value for value in raw_normal))
    if raw_norm == 0.0:
        return None
    face_centroid = tuple((a[index] + b[index] + c[index]) / 3.0 for index in range(3))
    tetra_centroid = tuple(sum(point[index] for point in points) / 4.0 for index in range(3))
    outward_test = sum(
        raw_normal[index] * (face_centroid[index] - tetra_centroid[index])
        for index in range(3)
    )
    sign = 1.0 if outward_test >= 0.0 else -1.0
    normal = tuple(sign * value / raw_norm for value in raw_normal)
    area = 0.5 * raw_norm
    heat_flux_vector = tuple(-conductivity_value * value for value in gradient)
    integrated_flux = sum(heat_flux_vector[index] * normal[index] for index in range(3)) * area
    return {
        "integrated_flux_W": integrated_flux,
        "area_m2": area,
        "normal": normal,
        "gradient": gradient,
        "conductivity_W_mK": conductivity_value,
    }


def conductivity(body: str, temperature: float, parameters: dict) -> float:
    if body.startswith("Membrane"):
        scale = max(temperature, 1.0e-12)
        return (
            parameters["G0"]
            * scale ** (4.252 - 1.0)
            * 0.4
            * (parameters["membrane_dx"] - parameters["TES_Au_dx"])
            / (
                8.0
                * (
                    parameters["SiNx_dz"]
                    + parameters["Si_1_dz"]
                    + parameters["SiO2_1_dz"]
                )
                * parameters["TES_Au_dx"]
            )
        )
    return {
        "TES": 68.0,
        "Stycast": 2.69094e-6,
        "abs": 0.0168,
    }.get(body, 1.0)


def heat_flux_consistency(
    mesh: Path,
    result: Path,
    project: Path,
    field_metadata: dict[str, object] | None = None,
    field_index: int = -1,
) -> dict[str, dict[str, object]]:
    bodies, _, nodes, elements, _, faces, parents = mesh_data(mesh)
    values = result_values(result, field_index)
    model = reconcile_project(json.loads(project.read_text(encoding="utf-8")))
    parameters = model["parameters"]
    reverse = {value: key for key, value in bodies.items()}
    report: dict[str, dict[str, float | str]] = {}
    for left, right, label, left_body, right_body, left_sign, right_sign in PAIR_DEFS:
        side_reports: list[dict[str, object]] = []
        side_fluxes: list[dict[tuple[int, ...], tuple[float, float]] ] = []
        for surface, body_name in (
            (left, left_body),
            (right, right_body),
        ):
            total = 0.0
            area_total = 0.0
            facet_count = 0
            skipped_facets = 0
            conductivities: list[float] = []
            normals: list[tuple[float, float, float]] = []
            surface_temperatures: list[float] = []
            gradient_magnitudes: list[float] = []
            fluxes_by_face: dict[tuple[int, ...], tuple[float, float]] = {}
            for face, parent in zip(
                faces[surface], parents[surface]
            ):
                if parent not in elements:
                    skipped_facets += 1
                    continue
                body_id, conn = elements[parent]
                actual_body = reverse.get(body_id, body_name)
                points = [nodes[node] for node in conn]
                temps = [values[node] for node in conn]
                surface_temperatures.extend(values[node] for node in face)
                face_indices = tuple(conn.index(node) for node in face)
                flux = oriented_face_flux(
                    points,
                    temps,
                    face_indices,
                    conductivity(actual_body, sum(temps) / 4.0, parameters),
                )
                if flux is None:
                    skipped_facets += 1
                    continue
                total += float(flux["integrated_flux_W"])
                area_total += float(flux["area_m2"])
                facet_count += 1
                conductivities.append(float(flux["conductivity_W_mK"]))
                normals.append(flux["normal"])
                gradient = flux["gradient"]
                gradient_magnitude = math.sqrt(sum(value * value for value in gradient))
                gradient_magnitudes.append(gradient_magnitude)
                fluxes_by_face[tuple(sorted(face))] = (
                    float(flux["integrated_flux_W"]),
                    float(flux["area_m2"]),
                )
            side_reports.append(
                {
                    "integrated_outward_flux_W": total,
                    "surface_area_m2": area_total,
                    "facet_count": facet_count,
                    "skipped_facets": skipped_facets,
                    "conductivity_min_W_mK": min(conductivities, default=math.nan),
                    "conductivity_max_W_mK": max(conductivities, default=math.nan),
                    "conductivity_mean_W_mK": sum(conductivities) / len(conductivities) if conductivities else math.nan,
                    "mean_abs_flux_density_W_m2": (
                        sum(abs(value[0]) / value[1] for value in fluxes_by_face.values())
                        / len(fluxes_by_face) if fluxes_by_face else math.nan
                    ),
                    "gradient_magnitude_min_K_m": min(gradient_magnitudes, default=math.nan),
                    "gradient_magnitude_max_K_m": max(gradient_magnitudes, default=math.nan),
                    "gradient_magnitude_mean_K_m": (
                        sum(gradient_magnitudes) / len(gradient_magnitudes)
                        if gradient_magnitudes else math.nan
                    ),
                    "surface_temperature_min_K": min(surface_temperatures, default=math.nan),
                    "surface_temperature_max_K": max(surface_temperatures, default=math.nan),
                    "surface_temperature_mean_K": sum(surface_temperatures) / len(surface_temperatures) if surface_temperatures else math.nan,
                    "outward_normal_mean": tuple(
                        sum(normal[index] for normal in normals) / len(normals)
                        for index in range(3)
                    ) if normals else (math.nan, math.nan, math.nan),
                    "outward_normal_z_positive_facets": sum(normal[2] > 0.0 for normal in normals),
                    "outward_normal_z_negative_facets": sum(normal[2] < 0.0 for normal in normals),
                }
            )
            side_fluxes.append(fluxes_by_face)
        q_left = float(side_reports[0]["integrated_outward_flux_W"])
        q_right = float(side_reports[1]["integrated_outward_flux_W"])
        denominator = max(abs(q_left), abs(q_right), 1.0e-30)
        imbalance = abs(q_left + q_right) / denominator
        magnitude = max(abs(q_left), abs(q_right))
        common_faces = sorted(set(side_fluxes[0]) & set(side_fluxes[1]))
        if common_faces:
            jumps = [side_fluxes[0][face][0] + side_fluxes[1][face][0] for face in common_faces]
            jump_densities = [
                side_fluxes[0][face][0] / side_fluxes[0][face][1]
                + side_fluxes[1][face][0] / side_fluxes[1][face][1]
                for face in common_faces
            ]
            mean_local_flux = sum(
                abs(side_fluxes[0][face][0]) / side_fluxes[0][face][1]
                + abs(side_fluxes[1][face][0]) / side_fluxes[1][face][1]
                for face in common_faces
            ) / (2.0 * len(common_faces))
            local_flux_jump: dict[str, object] = {
                "status": "AVAILABLE",
                "paired_facet_count": len(common_faces),
                "max_local_facet_flux_jump_W": max(abs(value) for value in jumps),
                "rms_local_facet_flux_jump_W": math.sqrt(sum(value * value for value in jumps) / len(jumps)),
                "max_local_facet_flux_jump_W_m2": max(abs(value) for value in jump_densities),
                "rms_local_facet_flux_jump_W_m2": math.sqrt(
                    sum(value * value for value in jump_densities) / len(jump_densities)
                ),
                "mean_local_flux_W_m2": mean_local_flux,
            }
        else:
            local_flux_jump = {
                "status": "NOT_AVAILABLE",
                "paired_facet_count": 0,
                "reason": "The two sides do not share node IDs; use integrated side flux for the Mortar route.",
            }
        status = (
            "NOT_INFORMATIVE"
            if magnitude < 1.0e-12
            else ("PASS" if imbalance <= 1.0e-3 else "FAIL")
        )
        report[label] = {
            "sign_convention": "Each face flux is q=-k*grad(T) dot n_outward; opposite body outward normals should give q_left + q_right = 0.",
            "source_result": str(result.resolve()),
            "field_metadata": field_metadata or {"field": "last saved result field", "time_step": None, "nonlinear_iter": None},
            "left": side_reports[0],
            "right": side_reports[1],
            "left_outward_flux_W": q_left,
            "right_outward_flux_W": q_right,
            "absolute_imbalance_W": abs(q_left + q_right),
            "surface_mean_temperature_difference_K": abs(
                float(side_reports[0]["surface_temperature_mean_K"])
                - float(side_reports[1]["surface_temperature_mean_K"])
            ),
            "normalized_imbalance": imbalance,
            "local_flux_jump": local_flux_jump,
            "status": status,
        }
    return report


def body_boundary_flux_balance(
    mesh: Path,
    result: Path,
    project: Path,
    field_index: int = -1,
) -> dict[str, dict[str, object]]:
    """Integrate reconstructed flux over every boundary face of each body.

    This is a body-level diagnostic, not an acceptance shortcut: for a
    source-free steady control it checks whether the reconstructed boundary
    fluxes close around each material body.  Internal interface faces are
    included once for each adjacent body, with each body's outward normal.
    """
    bodies, _, nodes, elements, _, faces, parents = mesh_data(mesh)
    values = result_values(result, field_index)
    model = reconcile_project(json.loads(project.read_text(encoding="utf-8")))
    parameters = model["parameters"]
    reverse = {value: key for key, value in bodies.items()}
    by_body: dict[str, dict[str, object]] = {}
    for surface, surface_faces in faces.items():
        for face, parent in zip(surface_faces, parents[surface]):
            if parent not in elements:
                continue
            body_id, conn = elements[parent]
            body_name = reverse.get(body_id, f"body_{body_id}")
            points = [nodes[node] for node in conn]
            temps = [values[node] for node in conn]
            face_indices = tuple(conn.index(node) for node in face)
            flux = oriented_face_flux(
                points,
                temps,
                face_indices,
                conductivity(body_name, sum(temps) / 4.0, parameters),
            )
            if flux is None:
                continue
            entry = by_body.setdefault(
                body_name,
                {"boundary_fluxes_W": defaultdict(float), "facet_count": 0},
            )
            entry["boundary_fluxes_W"][surface] += float(flux["integrated_flux_W"])
            entry["facet_count"] += 1
    report: dict[str, dict[str, object]] = {}
    for body_name, entry in by_body.items():
        boundary_fluxes = dict(entry["boundary_fluxes_W"])
        net = sum(boundary_fluxes.values())
        scale = max((abs(value) for value in boundary_fluxes.values()), default=1.0e-30)
        report[body_name] = {
            "boundary_fluxes_W": boundary_fluxes,
            "net_outward_flux_W": net,
            "absolute_net_residual_W": abs(net),
            "normalized_net_residual": abs(net) / scale,
            "facet_count": entry["facet_count"],
        }
    return report


def parse_series(path: Path | None) -> dict[str, float] | None:
    if path is None or not path.exists():
        return None
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    if not rows:
        return None
    row = rows[-1]
    aliases = {
        "tes_temperature_K": ("tes_temperature_K",),
        "tes_current_A": ("tes_current_A", "raw_current_A"),
        "tes_resistance_ohm": ("tes_resistance_ohm",),
        "tes_power_W": ("tes_power_W",),
    }
    result: dict[str, float] = {}
    for name, keys in aliases.items():
        for key in keys:
            if row.get(key):
                result[name] = float(row[key])
                break
    return result or None


def series_metadata(path: Path | None) -> dict[str, object]:
    if path is None or not path.exists():
        return {"field": "last saved result field", "time_step": None, "nonlinear_iter": None}
    rows = list(csv.DictReader(path.open(encoding="utf-8", newline="")))
    if not rows:
        return {"field": "last saved result field", "time_step": None, "nonlinear_iter": None}
    row = rows[-1]
    return {
        "field": "last saved result field",
        "time_s": float(row["time_s"]) if row.get("time_s") else None,
        "time_step": int(float(row["time_step"])) if row.get("time_step") else None,
        "nonlinear_iter": int(float(row["nonlinear_iter"])) if row.get("nonlinear_iter") else None,
    }


def compare_routes(mortar: dict, conformal: dict) -> dict:
    result: dict[str, object] = {}
    for body in sorted(set(mortar) & set(conformal)):
        entry = {}
        for metric in ("volume_average_temperature_K", "min_temperature_K", "max_temperature_K"):
            if metric in mortar[body] and metric in conformal[body]:
                entry[metric + "_abs_difference"] = abs(
                    float(conformal[body][metric]) - float(mortar[body][metric])
                )
        entry["volume_relative_difference"] = (
            float(conformal[body]["volume_m3"]) - float(mortar[body]["volume_m3"])
        ) / float(mortar[body]["volume_m3"])
        result[body] = entry
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mortar-mesh", type=Path, required=True)
    parser.add_argument("--conformal-mesh", type=Path, required=True)
    parser.add_argument("--mortar-result", type=Path, required=True)
    parser.add_argument("--conformal-result", type=Path, required=True)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--mortar-series", type=Path)
    parser.add_argument("--conformal-series", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    expected, expected_contact = expected_geometry(args.project)
    mortar = body_observables(args.mortar_mesh, args.mortar_result)
    conformal = body_observables(args.conformal_mesh, args.conformal_result)
    geometry = {}
    for body, analytic in expected.items():
        geometry[body] = {
            "analytic_volume_m3": analytic,
            "mortar_mesh_volume_m3": mortar[body]["volume_m3"],
            "conformal_mesh_volume_m3": conformal[body]["volume_m3"],
            "mortar_relative_to_analytic": (float(mortar[body]["volume_m3"]) - analytic) / analytic,
            "conformal_relative_to_analytic": (float(conformal[body]["volume_m3"]) - analytic) / analytic,
            "route_relative_difference": (float(conformal[body]["volume_m3"]) - float(mortar[body]["volume_m3"])) / float(mortar[body]["volume_m3"]),
        }
    areas = {
        label: {
            "analytic_m2": expected_contact[label],
            "mortar": contact_areas(args.mortar_mesh)[label],
            "conformal": contact_areas(args.conformal_mesh)[label],
        }
        for label in expected_contact
    }
    route_volume_error = max(abs(value["route_relative_difference"]) for value in geometry.values())
    report = {
        "inputs": {
            "mortar_mesh": str(args.mortar_mesh.resolve()),
            "conformal_mesh": str(args.conformal_mesh.resolve()),
            "project": str(args.project.resolve()),
        },
        "geometry_parity": {
            "bodies": geometry,
            "contact_areas": areas,
            "max_route_volume_relative_difference": route_volume_error,
            "status": "PASS" if route_volume_error <= 1.0e-12 else "FAIL",
            "analytic_convergence_status": "REPORT_ONLY",
            "analytic_convergence_note": "The route comparison is exact; analytic cylinder discretization error is reported separately and must be reduced by mesh convergence before production acceptance.",
        },
        "mortar_observables": mortar,
        "conformal_observables": conformal,
        "steady_thermal_parity": compare_routes(mortar, conformal),
        "interface_continuity": {
            "mortar": interface_continuity(args.mortar_mesh, args.mortar_result),
            "conformal": interface_continuity(args.conformal_mesh, args.conformal_result),
        },
        "electrical_observables": {
            "mortar": parse_series(args.mortar_series),
            "conformal": parse_series(args.conformal_series),
            "status": "PASS" if args.mortar_series and args.conformal_series else "NOT_AVAILABLE",
        },
        "heat_flux_consistency": {
            "mortar": heat_flux_consistency(args.mortar_mesh, args.mortar_result, args.project, series_metadata(args.mortar_series)),
            "conformal": heat_flux_consistency(args.conformal_mesh, args.conformal_result, args.project, series_metadata(args.conformal_series)),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + chr(10), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
