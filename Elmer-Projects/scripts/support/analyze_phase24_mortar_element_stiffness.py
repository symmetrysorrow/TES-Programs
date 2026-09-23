#!/usr/bin/env python3
"""Audit element-level thermal stiffness proxies at the two mortar interfaces.

No solver or production input is changed. This reconstructs linear-element
temperature gradients from saved steady result files and combines them with
the existing mesh boundary faces. The k*A/h quantity is a local normal
stiffness scale, not a replacement for the assembled mortar operator.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HISTORICAL_MESH = ROOT / "work/meshes/mesh_singlepixel_prod_v2"
PHASE24_MESH = ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar"
HISTORICAL_RESULT = HISTORICAL_MESH / "historical_1p00p.result"
PHASE24_RESULT = PHASE24_MESH / "phase24_stycast_substrate_coupling_only_1p00P.result"
HISTORICAL_CAPTURE = ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p00P/ts0001_nl0001"
PHASE24_CAPTURE = ROOT / "artifacts/phase24_thermal_network_localization/capture/original_phase24_1p00P/ts0001_nl0001"
SIF = ROOT / "generated/cases/case_phase24_historical_one_shot.sif"
OUTPUT = ROOT / "artifacts/phase24_mortar_element_stiffness"


def read_nodes(path: Path) -> dict[int, np.ndarray]:
    nodes: dict[int, np.ndarray] = {}
    for line in path.joinpath("mesh.nodes").read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if fields:
            nodes[int(fields[0])] = np.asarray(
                list(map(float, fields[2:5])), dtype=np.float64
            )
    return nodes


def read_names(path: Path) -> tuple[dict[int, str], dict[int, str]]:
    bodies: dict[int, str] = {}
    boundaries: dict[int, str] = {}
    section = ""
    pattern = re.compile(r"^\$\s*(.*?)\s*=\s*(-?\d+)\s*$")
    for line in path.joinpath("mesh.names").read_text(encoding="utf-8").splitlines():
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
        name, value = match.group(1), int(match.group(2))
        if section == "bodies":
            bodies[value] = name
        elif section == "boundaries":
            boundaries[value] = name
    return bodies, boundaries


def read_elements(path: Path) -> tuple[dict[int, dict], dict[int, set[int]]]:
    elements: dict[int, dict] = {}
    body_nodes: dict[int, set[int]] = defaultdict(set)
    for line in path.joinpath("mesh.elements").read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if not fields:
            continue
        record = {
            "id": int(fields[0]),
            "body": int(fields[1]),
            "type": int(fields[2]),
            "nodes": [int(value) for value in fields[3:]],
        }
        elements[record["id"]] = record
        body_nodes[record["body"]].update(record["nodes"])
    return elements, body_nodes


def read_material(sif: Path) -> dict[str, float]:
    text = sif.read_text(encoding="utf-8")
    for _, block in re.findall(r"Material\s+(\d+)\s*\n(.*?)(?=\nEnd\s*)", text, re.S):
        if re.search(r'Name\s*=\s*"Stycast"', block):
            match = re.search(r"Heat Conductivity\s*=\s*([0-9.eE+-]+)", block)
            if not match:
                raise RuntimeError("Stycast conductivity is not constant in SIF")
            return {"k_W_mK": float(match.group(1))}
    raise RuntimeError(f"Stycast material not found in {sif}")


def load_temperature(result: Path, capture: Path | None = None) -> dict[int, float]:
    lines = result.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip().lower() != "temperature":
            continue
        cursor = index + 1
        while cursor < len(lines) and not lines[cursor].strip().lower().startswith("perm:"):
            cursor += 1
        if cursor >= len(lines):
            continue
        header = lines[cursor].split()
        count = int(header[1])
        permutation: dict[int, int] = {}
        for row in range(cursor + 1, cursor + 1 + count):
            fields = lines[row].split()
            permutation[int(fields[0])] = int(fields[1])
        if capture is None:
            values_start = cursor + 1 + count
            values = [
                float(lines[row].strip())
                for row in range(values_start, values_start + count)
            ]
        else:
            values_array = np.zeros(count, dtype=np.float64)
            with (capture / "full_x_after.dat").open(
                encoding="utf-8", errors="replace"
            ) as stream:
                for value_line in stream:
                    fields = value_line.split()
                    if len(fields) >= 2:
                        row_index = int(fields[0])
                        if 1 <= row_index <= count:
                            values_array[row_index - 1] = float(fields[1])
            values = values_array.tolist()
        if len(values) != count or not np.all(np.isfinite(values)):
            raise RuntimeError(f"temperature vector length mismatch in {result}")
        return {
            node_id: values[dof - 1]
            for node_id, dof in permutation.items()
            if dof > 0
        }
    raise RuntimeError(f"Temperature Perm table not found in {result}")


def triangle_area(points: list[np.ndarray]) -> float:
    return 0.5 * float(np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])))


def face_area(points: list[np.ndarray]) -> float:
    if len(points) == 3:
        return triangle_area(points)
    if len(points) == 4:
        return triangle_area([points[0], points[1], points[2]]) + triangle_area(
            [points[0], points[2], points[3]]
        )
    raise RuntimeError(f"unsupported boundary face with {len(points)} nodes")


def element_volume(points: list[np.ndarray], element_type: int) -> float:
    def tetra(ps: list[np.ndarray]) -> float:
        return abs(float(np.linalg.det(np.column_stack(
            (ps[1] - ps[0], ps[2] - ps[0], ps[3] - ps[0])
        )))) / 6.0

    if element_type == 504:
        return tetra(points)
    if element_type == 706:
        return (
            tetra([points[0], points[1], points[2], points[3]])
            + tetra([points[1], points[2], points[3], points[4]])
            + tetra([points[2], points[3], points[4], points[5]])
        )
    raise RuntimeError(f"unsupported volume element type {element_type}")


def element_stiffness(
    points: list[np.ndarray], element_type: int, k: float
) -> np.ndarray:
    """Return the exact/degree-2 integrated local conduction matrix."""
    if element_type == 504:
        matrix = np.asarray([[1.0, *point] for point in points], dtype=np.float64)
        inverse = np.linalg.inv(matrix)
        gradients = inverse[1:, :].T
        return k * element_volume(points, element_type) * (gradients @ gradients.T)
    if element_type != 706:
        raise RuntimeError(f"unsupported stiffness element type {element_type}")

    # Wedge shape functions on a triangular reference base and t in [0,1].
    triangle_points = [
        (2.0 / 3.0, 1.0 / 6.0),
        (1.0 / 6.0, 2.0 / 3.0),
        (1.0 / 6.0, 1.0 / 6.0),
    ]
    triangle_weights = [1.0 / 3.0] * 3
    line_points = [(0.2113248654051871, 0.5), (0.7886751345948129, 0.5)]
    local = np.zeros((6, 6), dtype=np.float64)
    for (xi, eta), triangle_weight in zip(triangle_points, triangle_weights):
        barycentric = np.asarray([1.0 - xi - eta, xi, eta])
        for t, line_weight in line_points:
            bottom = np.asarray(points[:3])
            top = np.asarray(points[3:6])
            dxi = (1.0 - t) * (bottom[1] - bottom[0]) + t * (top[1] - top[0])
            deta = (1.0 - t) * (bottom[2] - bottom[0]) + t * (top[2] - top[0])
            dt = (
                barycentric[0] * (top[0] - bottom[0])
                + barycentric[1] * (top[1] - bottom[1])
                + barycentric[2] * (top[2] - bottom[2])
            )
            jacobian = np.column_stack((dxi, deta, dt))
            determinant = abs(float(np.linalg.det(jacobian)))
            if determinant <= 0.0:
                raise RuntimeError("degenerate wedge Jacobian")
            dref = np.asarray(
                [
                    [-1.0 + t, 1.0 - t, 0.0, -t, t, 0.0],
                    [-1.0 + t, 0.0, 1.0 - t, -t, 0.0, t],
                    [-barycentric[0], -barycentric[1], -barycentric[2],
                     barycentric[0], barycentric[1], barycentric[2]],
                ],
                dtype=np.float64,
            )
            gradients = np.linalg.inv(jacobian).T @ dref
            weight = determinant * triangle_weight * line_weight
            local += k * weight * (gradients.T @ gradients)
    return local


def affine_gradient(points: list[np.ndarray], temperatures: list[float]) -> np.ndarray:
    design = np.asarray([[1.0, *point] for point in points], dtype=np.float64)
    values = np.asarray(temperatures, dtype=np.float64)
    coefficients, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    return coefficients[1:]


def percentile_stats(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "count": int(array.size),
        "sum": float(array.sum()),
        "mean": float(array.mean()),
        "p50": float(np.percentile(array, 50)),
        "p95": float(np.percentile(array, 95)),
        "max": float(array.max()),
        "min": float(array.min()),
    }


def parse_boundary_faces(
    mesh: Path,
    elements: dict[int, dict],
    nodes: dict[int, np.ndarray],
    temperatures: dict[int, float],
    boundary_ids: set[int],
    interface: str,
    side: str,
    k: float,
) -> list[dict]:
    rows = []
    boundary_text = mesh.joinpath("mesh.boundary").read_text(
        encoding="utf-8", errors="replace"
    )
    for line in boundary_text.splitlines()[1:]:
        fields = line.split()
        if len(fields) < 8 or int(fields[1]) not in boundary_ids:
            continue
        boundary_id = int(fields[1])
        parent_id = int(fields[2])
        parent = elements.get(parent_id)
        if parent is None:
            continue
        face_type = int(fields[4])
        face_node_count = 3 if face_type == 303 else 4 if face_type == 404 else 0
        if not face_node_count:
            continue
        face_node_ids = [int(value) for value in fields[5:5 + face_node_count]]
        face_points = [nodes[node_id] for node_id in face_node_ids]
        element_points = [nodes[node_id] for node_id in parent["nodes"]]
        face_center = np.mean(face_points, axis=0)
        element_center = np.mean(element_points, axis=0)
        normal_raw = np.cross(face_points[1] - face_points[0], face_points[2] - face_points[0])
        normal_norm = float(np.linalg.norm(normal_raw))
        normal = normal_raw / max(normal_norm, 1.0e-300)
        area = face_area(face_points)
        h = abs(float(np.dot(element_center - face_center, normal)))
        face_temperature = float(np.mean([temperatures[node_id] for node_id in face_node_ids]))
        element_temperature = float(np.mean([
            temperatures[node_id] for node_id in parent["nodes"]
        ]))
        gradient = affine_gradient(
            element_points,
            [temperatures[node_id] for node_id in parent["nodes"]],
        )
        normal_gradient = abs(float(np.dot(gradient, normal)))
        face_indices = [
            parent["nodes"].index(node_id) for node_id in face_node_ids
        ]
        local_matrix = element_stiffness(element_points, parent["type"], k)
        face_block_diagonal = float(np.sum(np.diag(local_matrix)[face_indices]))
        matrix_trace = float(np.trace(local_matrix))
        rows.append(
            {
                "interface": interface,
                "side": side,
                "boundary_id": boundary_id,
                "parent_element": parent_id,
                "element_type": parent["type"],
                "face_type": face_type,
                "face_node_count": face_node_count,
                "face_area_m2": area,
                "element_volume_m3": element_volume(element_points, parent["type"]),
                "centroid_to_face_h_m": h,
                "local_kA_over_h_W_K": k * area / max(h, 1.0e-300),
                "local_face_block_diagonal_W_K": face_block_diagonal,
                "local_matrix_trace_W_K": matrix_trace,
                "local_face_block_trace_fraction": face_block_diagonal
                / max(matrix_trace, 1.0e-300),
                "face_temperature_K": face_temperature,
                "element_temperature_K": element_temperature,
                "abs_element_face_deltaT_K": abs(element_temperature - face_temperature),
                "abs_normal_gradient_K_m": normal_gradient,
                "abs_normal_flux_W_m2": k * normal_gradient,
                "abs_face_heatflow_proxy_W": k * area * normal_gradient,
            }
        )
    return rows


def analyze_case(
    label: str,
    mesh: Path,
    result: Path,
    capture: Path,
    material: dict[str, float],
) -> tuple[dict, list[dict]]:
    nodes = read_nodes(mesh)
    elements, _ = read_elements(mesh)
    bodies, _ = read_names(mesh)
    temperatures = load_temperature(result, capture)
    k = material["k_W_mK"]
    if label == "historical":
        interfaces = [
            ("TES_to_Stycast", "TES", "Stycast", {25}, {26}),
            ("Stycast_to_substrate", "Stycast", "abs", {27}, {28}),
        ]
    else:
        interfaces = [
            ("TES_to_Stycast", "TES", "Stycast", {1105}, {1204}),
            ("Stycast_to_substrate", "Stycast", "abs", {1205}, {1004}),
        ]
    rows = []
    interface_summary = {}
    for interface, side_a, side_b, boundaries_a, boundaries_b in interfaces:
        for side, side_boundaries, body_name in (
            ("side_a", boundaries_a, side_a),
            ("side_b", boundaries_b, side_b),
        ):
            side_rows = parse_boundary_faces(
                mesh, elements, nodes, temperatures, side_boundaries, interface, side, k
            )
            for row in side_rows:
                row["case"] = label
            rows.extend(side_rows)
            areas = [row["face_area_m2"] for row in side_rows]
            stiffness = [row["local_kA_over_h_W_K"] for row in side_rows]
            heights = [row["centroid_to_face_h_m"] for row in side_rows]
            gradients = [row["abs_normal_gradient_K_m"] for row in side_rows]
            heatflows = [row["abs_face_heatflow_proxy_W"] for row in side_rows]
            face_stiffness = [
                row["local_face_block_diagonal_W_K"] for row in side_rows
            ]
            if not side_rows:
                raise RuntimeError(f"no faces found for {label} {interface} {side}")
            interface_summary[f"{interface}:{side}"] = {
                "body": body_name,
                "body_id": next(
                    (body_id for body_id, name in bodies.items() if name == body_name), None
                ),
                "boundary_ids": sorted(side_boundaries),
                "face_count": len(side_rows),
                "face_area": percentile_stats(areas),
                "centroid_to_face_h": percentile_stats(heights),
                "local_kA_over_h": percentile_stats(stiffness),
                "local_face_block_diagonal": percentile_stats(face_stiffness),
                "normal_gradient": percentile_stats(gradients),
                "heatflow_proxy": percentile_stats(heatflows),
                "temperature_face_min_K": min(
                    row["face_temperature_K"] for row in side_rows
                ),
                "temperature_face_max_K": max(
                    row["face_temperature_K"] for row in side_rows
                ),
            }
    return {
        "mesh": str(mesh),
        "result": str(result),
        "capture": str(capture),
        "element_count": len(elements),
        "node_count": len(nodes),
        "material_k_W_mK": k,
        "body_element_counts": dict(Counter(
            record["body"] for record in elements.values()
        )),
        "interface_summary": interface_summary,
    }, rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical-mesh", type=Path, default=HISTORICAL_MESH)
    parser.add_argument("--phase24-mesh", type=Path, default=PHASE24_MESH)
    parser.add_argument("--historical-result", type=Path, default=HISTORICAL_RESULT)
    parser.add_argument("--phase24-result", type=Path, default=PHASE24_RESULT)
    parser.add_argument("--historical-capture", type=Path, default=HISTORICAL_CAPTURE)
    parser.add_argument("--phase24-capture", type=Path, default=PHASE24_CAPTURE)
    parser.add_argument("--sif", type=Path, default=SIF)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    material = read_material(args.sif)
    historical, historical_rows = analyze_case(
        "historical",
        args.historical_mesh,
        args.historical_result,
        args.historical_capture,
        material,
    )
    phase24, phase24_rows = analyze_case(
        "phase24",
        args.phase24_mesh,
        args.phase24_result,
        args.phase24_capture,
        material,
    )
    all_rows = historical_rows + phase24_rows
    write_csv(args.output / "interface_element_local_metrics.csv", all_rows)

    comparison = {}
    for key in historical["interface_summary"]:
        h = historical["interface_summary"][key]
        p = phase24["interface_summary"][key]
        comparison[key] = {
            "face_count_phase24_over_historical": p["face_count"] / h["face_count"],
            "area_sum_phase24_over_historical": p["face_area"]["sum"] / h["face_area"]["sum"],
            "centroid_to_face_h_mean_phase24_over_historical": p[
                "centroid_to_face_h"
            ]["mean"] / h["centroid_to_face_h"]["mean"],
            "kA_over_h_sum_phase24_over_historical": p["local_kA_over_h"]["sum"]
            / h["local_kA_over_h"]["sum"],
            "face_block_diagonal_sum_phase24_over_historical": p[
                "local_face_block_diagonal"
            ]["sum"] / h["local_face_block_diagonal"]["sum"],
            "gradient_mean_phase24_over_historical": p["normal_gradient"]["mean"]
            / max(h["normal_gradient"]["mean"], 1.0e-300),
            "heatflow_proxy_sum_phase24_over_historical": p["heatflow_proxy"]["sum"]
            / max(h["heatflow_proxy"]["sum"], 1.0e-300),
        }

    report = {
        "method": {
            "temperature_source": "saved 1.00P native full_x_after capture, mapped through result Temperature Perm",
            "local_stiffness_proxy": "k * face_area / centroid_to_face_distance",
            "local_gradient": "affine least-squares gradient over parent element nodes",
            "mortar_caveat": "face tessellations are nonconforming; side sums are not paired flux conservation",
        },
        "material": material,
        "historical": historical,
        "phase24": phase24,
        "comparison": comparison,
    }
    (args.output / "mortar_element_stiffness.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# Phase24 mortar-neighborhood element stiffness audit",
        "",
        f"- Stycast conductivity used: {material['k_W_mK']:.12g} W/(m K)",
        "- k*A/h is a local parent-element normal stiffness scale, not the assembled mortar conductance.",
        "- Temperature gradients come from saved 1.00P native full_x_after captures mapped through the result Temperature Perm table.",
        "",
        "| interface side | faces historical | faces Phase24 | area ratio | mean h ratio | sum kA/h ratio | sum exact face-block K ratio | mean gradient ratio | heatflow proxy ratio |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for key, value in comparison.items():
        h = historical["interface_summary"][key]
        p = phase24["interface_summary"][key]
        lines.append(
            f"| {key} | {h['face_count']} | {p['face_count']} | "
            f"{value['area_sum_phase24_over_historical']:.9f} | "
            f"{value['centroid_to_face_h_mean_phase24_over_historical']:.9f} | "
            f"{value['kA_over_h_sum_phase24_over_historical']:.9f} | "
            f"{value['face_block_diagonal_sum_phase24_over_historical']:.9f} | "
            f"{value['gradient_mean_phase24_over_historical']:.9f} | "
            f"{value['heatflow_proxy_sum_phase24_over_historical']:.9f} |"
        )
    lines.extend(
        [
            "",
            "The face-area and local-stiffness comparisons are split by side because "
            "the mortar interfaces are nonconforming. A large side-to-side gradient "
            "change with near-equal area would point to local element geometry or "
            "temperature interpolation; a large area or kA/h change would point to "
            "face/height integration. The exact face-block K metric is computed "
            "from the element conduction matrix and is the stronger comparison.",
        ]
    )
    (args.output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
