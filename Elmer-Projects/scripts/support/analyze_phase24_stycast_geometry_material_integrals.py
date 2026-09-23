#!/usr/bin/env python3
"""Compare historical and Phase24 Stycast geometry/material integrals.

The analysis is deliberately mesh-only. It uses the existing converted meshes
and the existing SIF material definition; it does not run Elmer or modify any
production input. For the nominal 32 layers it reports V, V/dz, k*V, and the
1-D series resistance implied by each layer's effective area.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
HISTORICAL_DEFAULT = ROOT / "work/meshes/mesh_singlepixel_prod_v2"
PHASE24_DEFAULT = ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar"
SIF_DEFAULT = ROOT / "generated/cases/case_phase24_historical_one_shot.sif"
ARTIFACT_DEFAULT = ROOT / "artifacts/phase24_stycast_geometry_material_integrals"
LAYER_COUNT = 32


def read_nodes(path: Path) -> dict[int, tuple[float, float, float]]:
    nodes: dict[int, tuple[float, float, float]] = {}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if fields:
            # Elmer mesh.nodes: node_id, partition/marker, x, y, z.
            nodes[int(fields[0])] = tuple(map(float, fields[2:5]))
    return nodes


def read_body_id(path: Path, body_name: str) -> int:
    in_bodies = False
    pattern = re.compile(r"^\$\s*(.*?)\s*=\s*(-?\d+)\s*$")
    for line in path.joinpath("mesh.names").read_text(encoding="utf-8").splitlines():
        lower = line.lower()
        if "names for bodies" in lower:
            in_bodies = True
            continue
        if "names for boundaries" in lower:
            in_bodies = False
        if not in_bodies:
            continue
        match = pattern.match(line)
        if match and match.group(1) == body_name:
            return int(match.group(2))
    raise RuntimeError(f"body {body_name!r} not found in {path / 'mesh.names'}")


def read_elements(path: Path, body_id: int) -> list[dict]:
    rows = []
    for line in path.joinpath("mesh.elements").read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if fields and int(fields[1]) == body_id:
            rows.append(
                {
                    "id": int(fields[0]),
                    "type": int(fields[2]),
                    "nodes": [int(value) for value in fields[3:]],
                }
            )
    return rows


def sub(a: tuple[float, float, float], b: tuple[float, float, float]):
    return tuple(a[index] - b[index] for index in range(3))


def determinant(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> float:
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def tetra_volume(points: list[tuple[float, float, float]]) -> float:
    return abs(
        determinant(
            sub(points[1], points[0]),
            sub(points[2], points[0]),
            sub(points[3], points[0]),
        )
    ) / 6.0


def element_volume(nodes: dict[int, tuple[float, float, float]], record: dict) -> float:
    points = [nodes[node_id] for node_id in record["nodes"]]
    if record["type"] == 504:
        return tetra_volume(points)
    if record["type"] == 706:
        # Elmer wedge ordering: triangular base [0,1,2], translated base
        # [3,4,5]. This three-tetra decomposition is exact for a prism.
        return (
            tetra_volume([points[0], points[1], points[2], points[3]])
            + tetra_volume([points[1], points[2], points[3], points[4]])
            + tetra_volume([points[2], points[3], points[4], points[5]])
        )
    raise RuntimeError(f"unsupported Stycast element type {record['type']}")


def read_stycast_material(sif: Path) -> dict[str, float | str]:
    text = sif.read_text(encoding="utf-8")
    blocks = re.findall(r"Material\s+(\d+)\s*\n(.*?)(?=\nEnd\s*)", text, re.S)
    for material_id, block in blocks:
        if re.search(r'Name\s*=\s*"Stycast"', block):
            conductivity = re.search(r"Heat Conductivity\s*=\s*([0-9.eE+-]+)", block)
            density = re.search(r"Density\s*=\s*([0-9.eE+-]+)", block)
            capacity = re.search(r"Heat Capacity\s*=\s*([0-9.eE+-]+)", block)
            if not (conductivity and density and capacity):
                raise RuntimeError("Stycast material is not a constant-property block")
            return {
                "material_id": int(material_id),
                "name": "Stycast",
                "k_W_mK": float(conductivity.group(1)),
                "rho_kg_m3": float(density.group(1)),
                "cp_J_kgK": float(capacity.group(1)),
            }
    raise RuntimeError(f"Stycast material not found in {sif}")


def analyze_mesh(
    label: str,
    mesh: Path,
    body_name: str,
    material: dict[str, float | str],
) -> tuple[dict, list[dict]]:
    nodes = read_nodes(mesh / "mesh.nodes")
    body_id = read_body_id(mesh, body_name)
    elements = read_elements(mesh, body_id)
    if not elements:
        raise RuntimeError(f"no {body_name} elements in {mesh}")

    body_node_ids = {node_id for record in elements for node_id in record["nodes"]}
    body_z = [nodes[node_id][2] for node_id in body_node_ids]
    z_min, z_max = min(body_z), max(body_z)
    dz = (z_max - z_min) / LAYER_COUNT
    tolerance = max(1.0e-14, dz * 1.0e-8)
    layer_data: dict[int, dict] = defaultdict(
        lambda: {"volume_m3": 0.0, "element_count": 0, "node_ids": set()}
    )
    crossing_elements = []
    volume_by_type = Counter()

    for record in elements:
        z_values = [nodes[node_id][2] for node_id in record["nodes"]]
        centroid_z = sum(z_values) / len(z_values)
        layer = min(LAYER_COUNT - 1, max(0, int((centroid_z - z_min) / dz)))
        lower = z_min + layer * dz
        upper = z_min + (layer + 1) * dz
        if min(z_values) < lower - tolerance or max(z_values) > upper + tolerance:
            crossing_elements.append(
                {
                    "element_id": record["id"],
                    "layer": layer + 1,
                    "z_min": min(z_values),
                    "z_max": max(z_values),
                }
            )
        volume = element_volume(nodes, record)
        volume_by_type[record["type"]] += 1
        layer_data[layer]["volume_m3"] += volume
        layer_data[layer]["element_count"] += 1
        layer_data[layer]["node_ids"].update(record["nodes"])

    k = float(material["k_W_mK"])
    rho_cp = float(material["rho_kg_m3"]) * float(material["cp_J_kgK"])
    rows = []
    for layer in range(LAYER_COUNT):
        data = layer_data[layer]
        area = data["volume_m3"] / dz
        rows.append(
            {
                "mesh": label,
                "layer": layer + 1,
                "z_min_m": z_min + layer * dz,
                "z_max_m": z_min + (layer + 1) * dz,
                "thickness_m": dz,
                "element_count": data["element_count"],
                "node_count": len(data["node_ids"]),
                "volume_m3": data["volume_m3"],
                "effective_area_m2": area,
                "k_W_mK": k,
                "integral_k_dV_W_m": k * data["volume_m3"],
                "integral_rho_cp_dV_J_K": rho_cp * data["volume_m3"],
                "one_dimensional_layer_R_K_W": dz / (k * area),
            }
        )

    total_volume = sum(row["volume_m3"] for row in rows)
    total_area = total_volume / (z_max - z_min)
    total_resistance = sum(row["one_dimensional_layer_R_K_W"] for row in rows)
    summary = {
        "mesh": label,
        "path": str(mesh),
        "body_name": body_name,
        "body_id": body_id,
        "element_count": len(elements),
        "element_types": dict(volume_by_type),
        "node_count": len(body_node_ids),
        "x_min_m": min(nodes[node_id][0] for node_id in body_node_ids),
        "x_max_m": max(nodes[node_id][0] for node_id in body_node_ids),
        "y_min_m": min(nodes[node_id][1] for node_id in body_node_ids),
        "y_max_m": max(nodes[node_id][1] for node_id in body_node_ids),
        "z_min_m": z_min,
        "z_max_m": z_max,
        "thickness_m": z_max - z_min,
        "nominal_layer_thickness_m": dz,
        "crossing_element_count": len(crossing_elements),
        "total_volume_m3": total_volume,
        "volume_over_thickness_effective_area_m2": total_area,
        "integral_k_dV_W_m": k * total_volume,
        "integral_rho_cp_dV_J_K": rho_cp * total_volume,
        "one_dimensional_series_R_K_W": total_resistance,
        "layer_area_min_m2": min(row["effective_area_m2"] for row in rows),
        "layer_area_max_m2": max(row["effective_area_m2"] for row in rows),
        "layer_area_minmax_ratio": min(row["effective_area_m2"] for row in rows)
        / max(row["effective_area_m2"] for row in rows),
        "crossing_elements": crossing_elements[:20],
    }
    return summary, rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical", type=Path, default=HISTORICAL_DEFAULT)
    parser.add_argument("--phase24", type=Path, default=PHASE24_DEFAULT)
    parser.add_argument("--sif", type=Path, default=SIF_DEFAULT)
    parser.add_argument("--artifacts", type=Path, default=ARTIFACT_DEFAULT)
    args = parser.parse_args()
    args.artifacts.mkdir(parents=True, exist_ok=True)

    material = read_stycast_material(args.sif)
    historical, historical_rows = analyze_mesh(
        "historical", args.historical, "Stycast", material
    )
    phase24, phase24_rows = analyze_mesh(
        "phase24", args.phase24, "Stycast", material
    )
    comparison = {
        "phase24_over_historical_total_volume": phase24["total_volume_m3"]
        / historical["total_volume_m3"],
        "phase24_over_historical_effective_area": phase24[
            "volume_over_thickness_effective_area_m2"
        ]
        / historical["volume_over_thickness_effective_area_m2"],
        "phase24_over_historical_integral_k_dV": phase24["integral_k_dV_W_m"]
        / historical["integral_k_dV_W_m"],
        "phase24_over_historical_integral_rho_cp_dV": phase24["integral_rho_cp_dV_J_K"]
        / historical["integral_rho_cp_dV_J_K"],
        "phase24_over_historical_1d_series_R": phase24["one_dimensional_series_R_K_W"]
        / historical["one_dimensional_series_R_K_W"],
    }
    audit = {
        "material": material,
        "historical": historical,
        "phase24": phase24,
        "comparison": comparison,
    }
    (args.artifacts / "geometry_material_integrals.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    write_csv(args.artifacts / "stycast_layer_geometry_integrals.csv", historical_rows + phase24_rows)

    summary = [
        "# Phase24 Stycast geometry and material integrals",
        "",
        f"- Stycast conductivity: {material['k_W_mK']:.12g} W/(m K)",
        f"- Stycast density: {material['rho_kg_m3']:.12g} kg/m3",
        f"- Stycast heat capacity: {material['cp_J_kgK']:.12g} J/(kg K)",
        "",
        "| quantity | historical | Phase24 | Phase24 / historical |",
        "|---|---:|---:|---:|",
        f"| total volume (m3) | {historical['total_volume_m3']:.12e} | {phase24['total_volume_m3']:.12e} | {comparison['phase24_over_historical_total_volume']:.9f} |",
        f"| V / thickness effective area (m2) | {historical['volume_over_thickness_effective_area_m2']:.12e} | {phase24['volume_over_thickness_effective_area_m2']:.12e} | {comparison['phase24_over_historical_effective_area']:.9f} |",
        f"| integral k dV (W m) | {historical['integral_k_dV_W_m']:.12e} | {phase24['integral_k_dV_W_m']:.12e} | {comparison['phase24_over_historical_integral_k_dV']:.9f} |",
        f"| integral rho cp dV (J/K) | {historical['integral_rho_cp_dV_J_K']:.12e} | {phase24['integral_rho_cp_dV_J_K']:.12e} | {comparison['phase24_over_historical_integral_rho_cp_dV']:.9f} |",
        f"| 1-D series R (K/W) | {historical['one_dimensional_series_R_K_W']:.12e} | {phase24['one_dimensional_series_R_K_W']:.12e} | {comparison['phase24_over_historical_1d_series_R']:.9f} |",
        "",
        "## Layer uniformity",
        "",
        f"- Historical layer effective-area min/max ratio: {historical['layer_area_minmax_ratio']:.12f}",
        f"- Phase24 layer effective-area min/max ratio: {phase24['layer_area_minmax_ratio']:.12f}",
        f"- Historical element-plane crossings: {historical['crossing_element_count']}",
        f"- Phase24 element-plane crossings: {phase24['crossing_element_count']}",
        "",
        "The Stycast material is constant-k in the SIF, so the material conductivity "
        "integral is exactly k times mesh volume. The volume and effective-area ratios "
        "are approximately 0.995, far from the approximately 1.48 conductance excess.",
    ]
    (args.artifacts / "summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
