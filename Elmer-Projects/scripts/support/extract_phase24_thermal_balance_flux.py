"""Extract FE body balances and integrated interface fluxes.

The extraction uses the saved temperature result and the existing mesh only;
it does not alter a solver input or rerun production physics.  Values on
nonconforming mortar interfaces are reported as the two one-sided FE flux
integrals, with the interface association taken from ``mesh.names``.
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.support.run_phase24_historical_mesh_strict_validation import result_field


HIST_MESH = ROOT / "work" / "meshes" / "mesh_singlepixel_prod_v2"
HIST_RESULT = HIST_MESH / "case_phase24_historical_one_shot.result"
CURRENT_MESH = ROOT / "work" / "meshes" / "mesh_singlepixel_gpu_fine_stycast32_mortar"
CURRENT_RESULT = CURRENT_MESH / "case_phase24_restart_refine_fine_stycast32_mortar_mumps_mortar.result"
OUT = ROOT / "artifacts" / "phase24_historical_mesh_strict_steady_validation"


def names(mesh: Path) -> tuple[dict[int, str], dict[int, str]]:
    bodies: dict[int, str] = {}
    boundaries: dict[int, str] = {}
    in_boundaries = False
    for line in (mesh / "mesh.names").read_text(encoding="utf-8").splitlines():
        if "names for boundaries" in line.lower():
            in_boundaries = True
        elif "names for bodies" in line.lower():
            in_boundaries = False
        hit = re.match(r"\$\s+(.+?)\s+=\s+(\d+)$", line)
        if not hit:
            continue
        label, number = hit.group(1), int(hit.group(2))
        if in_boundaries:
            boundaries[number] = label
        else:
            bodies[number] = label
    return bodies, boundaries


def load_nodes(mesh: Path) -> np.ndarray:
    return np.loadtxt(mesh / "mesh.nodes", usecols=(2, 3, 4), dtype=np.float64)


def load_elements(mesh: Path) -> tuple[dict[int, tuple[int, int, list[int]]], dict[int, list[int]], dict[tuple[int, int], list[int]]]:
    elements: dict[int, tuple[int, int, list[int]]] = {}
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 7:
            elements[int(fields[0])] = (int(fields[1]), int(fields[2]), [int(value) for value in fields[3:]])
    boundary_faces: dict[int, list[int]] = defaultdict(list)
    boundary_records: dict[tuple[int, int], list[int]] = {}
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 8:
            boundary_faces[int(fields[1])].append(int(fields[2]))
            boundary_records[(int(fields[1]), int(fields[2]))] = [int(value) for value in fields[5:]]
    return elements, boundary_faces, boundary_records


def material_k(body: str, temperature: float) -> float:
    if body == "TES":
        return 68.0
    if body == "Stycast":
        return 2.69094e-6
    if "Membrane" in body:
        return (7.854e-6) * max(temperature, 1.0e-12) ** (4.252 - 1.0) * 0.4 * (0.0007 - 0.0005) / (8.0 * (1.0e-6 + 1.5e-5 + 1.0e-6) * 0.0005)
    if "SiNx" in body:
        return 1.02e-4
    if "SiO2" in body:
        return 4.94e-4
    if body.startswith("Si"):
        return 0.37
    return 0.0168


def tes_average_local(mesh: Path, field: np.ndarray, permutation: np.ndarray, tes_body_id: int) -> float:
    weights = np.zeros(permutation.size, dtype=np.float64)
    count = 0
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 4 or int(fields[1]) != tes_body_id:
            continue
        node_ids = [int(value) for value in fields[3:]]
        for node in node_ids:
            dof = int(permutation[node - 1])
            if dof > 0:
                weights[dof - 1] += 1.0 / len(node_ids)
        count += 1
    return float(np.dot(weights / count, field))


def tet_data(points: np.ndarray, temperatures: np.ndarray, conductivity: float) -> tuple[float, np.ndarray]:
    matrix = np.column_stack((np.ones(4), points))
    coefficients = np.linalg.solve(matrix, temperatures)
    gradient = coefficients[1:]
    volume = abs(np.linalg.det(points[1:] - points[0])) / 6.0
    return volume, -conductivity * gradient


def split_element(kind: int, node_ids: list[int]) -> list[list[int]]:
    if kind == 504:
        return [node_ids[:4]]
    if kind in {303, 706} and len(node_ids) >= 6:
        return [[node_ids[i] for i in tet] for tet in ((0, 1, 2, 3), (1, 2, 4, 3), (2, 4, 5, 3))]
    return []


def outward_face_flux(points: np.ndarray, field_by_node: np.ndarray, face_nodes: list[int], element_nodes: list[int], conductivity: float) -> float:
    selected = [node for node in element_nodes if node in face_nodes]
    if len(selected) < 3:
        return 0.0
    tet_nodes = next((tet for tet in split_element(706, element_nodes) if set(face_nodes).issubset(tet)), None)
    if tet_nodes is None:
        tet_nodes = next((tet for tet in split_element(504, element_nodes) if set(face_nodes).issubset(tet)), None)
    if tet_nodes is None:
        tet_nodes = element_nodes[:4]
    tet_points = points[np.asarray(tet_nodes) - 1]
    tet_values = field_by_node[np.asarray(tet_nodes) - 1]
    _, heat_flux = tet_data(tet_points, tet_values, conductivity)
    face_points = points[np.asarray(face_nodes) - 1]
    normal = np.cross(face_points[1] - face_points[0], face_points[2] - face_points[0])
    area_vector = 0.5 * normal
    face_centroid = face_points.mean(axis=0)
    if np.dot(area_vector, face_centroid - tet_points.mean(axis=0)) < 0.0:
        area_vector = -area_vector
    return float(np.dot(heat_flux, area_vector))


def analyze(mesh: Path, result: Path, route: str, joule_power: float) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    body_names, boundary_names = names(mesh)
    points = load_nodes(mesh)
    elements, boundary_faces, boundary_records = load_elements(mesh)
    field, permutation = result_field(result)
    field_by_node = np.zeros(len(points), dtype=np.float64)
    field_by_node[:] = field[permutation - 1]
    body_flux: dict[int, float] = defaultdict(float)
    boundary_flux: dict[int, float] = defaultdict(float)
    body_volume: dict[int, float] = defaultdict(float)
    for body_id, kind, node_ids in elements.values():
        body = body_names.get(body_id, str(body_id))
        for tet in split_element(kind, node_ids):
            tet_points = points[np.asarray(tet) - 1]
            tet_values = field_by_node[np.asarray(tet) - 1]
            volume, _ = tet_data(tet_points, tet_values, material_k(body, float(tet_values.mean())))
            body_volume[body_id] += volume
    for boundary_id, parents in boundary_faces.items():
        for parent_id in parents:
            if parent_id not in elements:
                continue
            body_id, kind, node_ids = elements[parent_id]
            body = body_names.get(body_id, str(body_id))
            face_nodes = boundary_records.get((boundary_id, parent_id), [])
            if face_nodes:
                values = field_by_node[np.asarray(face_nodes) - 1]
                flux = outward_face_flux(points, field_by_node, face_nodes, node_ids, material_k(body, float(values.mean())))
                boundary_flux[boundary_id] += flux
                body_flux[body_id] += flux
    balance_rows: list[dict[str, Any]] = []
    for body_id, body in sorted(body_names.items()):
        source = joule_power if body == "TES" else 0.0
        outgoing = body_flux[body_id]
        balance_rows.append({
            "route": route, "body_id": body_id, "body": body,
            "volume_m3": body_volume[body_id], "joule_input_W": source,
            "total_boundary_outgoing_W": outgoing,
            "energy_balance_error_W": outgoing - source,
            "evidence": "FE boundary-face flux from saved result; mortar faces are one-sided",
        })
    interface_rows: list[dict[str, Any]] = []
    by_name = {label: number for number, label in boundary_names.items()}
    pairs = [
        ("TES_to_membrane", ("TES__zmin", "Membrane_SiNx__zmax")),
        ("TES_to_Stycast", ("TES__zmax", "Stycast__zmin")),
        ("Stycast_to_substrate_or_bath", ("Stycast__zmax", "abs__zmin")),
    ]
    for interface, (left, right) in pairs:
        left_id, right_id = by_name.get(left), by_name.get(right)
        interface_rows.append({
            "route": route, "interface": interface,
            "left_boundary": left, "right_boundary": right,
            "left_flux_outgoing_W": boundary_flux[left_id] if left_id is not None else None,
            "right_flux_outgoing_W": boundary_flux[right_id] if right_id is not None else None,
            "signed_mismatch_W": (boundary_flux[left_id] + boundary_flux[right_id]) if left_id is not None and right_id is not None else None,
            "evidence": "integrated FE normal flux; sign is outward from each named body; mortar association from boundary names",
        })
    bath_id = by_name.get("bath")
    bath_flux = boundary_flux[bath_id] if bath_id is not None else None
    summary = {
        "route": route, "mesh": str(mesh.relative_to(ROOT)), "result": str(result.relative_to(ROOT)),
        "tes_temperature_K": tes_average_local(mesh, field, permutation, next(body_id for body_id, name in body_names.items() if name == "TES")), "joule_power_W": joule_power,
        "bath_boundary_id": bath_id, "bath_outgoing_W": bath_flux,
        "body_count": len(body_names), "boundary_count": len(boundary_names),
        "method_note": "Wedge elements are split into three linear tetrahedra for the flux diagnostic; this is a diagnostic FE integral, not a solver balance report.",
    }
    return balance_rows, interface_rows, summary


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    historical = analyze(HIST_MESH, HIST_RESULT, "historical_refined_or_one_shot", 3.2029879304735996e-10)
    current = analyze(CURRENT_MESH, CURRENT_RESULT, "phase24_refined", 4.2324271156482167e-10)
    balance = historical[0] + current[0]
    interface = historical[1] + current[1]
    for path, rows in ((OUT / "thermal_balance.csv", balance), (OUT / "interface_flux.csv", interface)):
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    summaries = {"historical": historical[2], "phase24": current[2]}
    for value in summaries.values():
        temp_delta = value["tes_temperature_K"] - 0.15
        value["secant_G_eff_W_per_K"] = value["bath_outgoing_W"] / temp_delta if value["bath_outgoing_W"] is not None and temp_delta else None
        value["G_eff_definition"] = "bath integrated outgoing flux / (TES average temperature - 150 mK); secant diagnostic, not dQ/dT perturbation"
    (OUT / "effective_conductance.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")
    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
