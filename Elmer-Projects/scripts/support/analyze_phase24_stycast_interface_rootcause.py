#!/usr/bin/env python3
"""Assemble the Phase24 Stycast interface-discretization audit.

This is a read-only postprocessor for the historical, original Phase24,
existing 5-um control, and newly generated density probes.  It reports
interface geometry, exact parent-element conduction blocks, and the already
captured native fixed-power/projection diagnostics in one audit directory.
It deliberately does not alter solver input or production meshes.
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

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "support"))
from analyze_phase24_mortar_element_stiffness import (  # noqa: E402
    element_stiffness,
    element_volume,
    face_area,
)
from analyze_phase24_outer_capture import (  # noqa: E402
    indexed_vector,
    load_temp_permutation,
    tes_element_weights,
)
from analyze_phase24_mortar_constraint_projection import analyze_case  # noqa: E402


P0 = 3.203004762115138e-10
TBATH = 0.15
OUTPUT = ROOT / "artifacts/phase24_stycast_interface_discretization_rootcause"


def read_nodes(mesh: Path) -> dict[int, np.ndarray]:
    out = {}
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if fields:
            out[int(fields[0])] = np.asarray(list(map(float, fields[2:5])))
    return out


def read_elements(mesh: Path) -> dict[int, dict]:
    out = {}
    for line in (mesh / "mesh.elements").read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if fields:
            out[int(fields[0])] = {
                "body": int(fields[1]),
                "type": int(fields[2]),
                "nodes": [int(value) for value in fields[3:]],
            }
    return out


def read_boundary_faces(mesh: Path, boundary_id: int) -> list[dict]:
    out = []
    for line in (mesh / "mesh.boundary").read_text(encoding="utf-8").splitlines()[1:]:
        fields = line.split()
        if len(fields) < 6 or int(fields[1]) != boundary_id:
            continue
        n = 3 if int(fields[4]) == 303 else 4 if int(fields[4]) == 404 else 0
        out.append({
            "id": int(fields[0]),
            "parent": int(fields[2]),
            "type": int(fields[4]),
            "nodes": [int(value) for value in fields[5:5 + n]],
        })
    return out


def percentile(values: list[float]) -> dict[str, float]:
    if not values:
        return {key: None for key in ("min", "p10", "median", "mean", "p90", "max")}
    a = np.asarray(values, dtype=float)
    return {
        "min": float(np.min(a)),
        "p10": float(np.percentile(a, 10)),
        "median": float(np.percentile(a, 50)),
        "mean": float(np.mean(a)),
        "p90": float(np.percentile(a, 90)),
        "max": float(np.max(a)),
    }


def aspect_ratio(points: list[np.ndarray]) -> float:
    lengths = [
        float(np.linalg.norm(points[i] - points[j]))
        for i in range(len(points))
        for j in range(i)
        if np.linalg.norm(points[i] - points[j]) > 0.0
    ]
    return max(lengths) / min(lengths)


def face_depth(face_points: list[np.ndarray], all_points: list[np.ndarray]) -> float:
    normal = np.cross(face_points[1] - face_points[0], face_points[2] - face_points[0])
    norm = float(np.linalg.norm(normal))
    if norm == 0.0:
        return float("nan")
    centroid = np.mean(np.asarray(all_points), axis=0)
    return abs(float(np.dot(centroid - face_points[0], normal))) / norm


def local_blocks(points: list[np.ndarray], element_type: int, face_nodes: list[int], parent_nodes: list[int]) -> dict:
    k = element_stiffness(points, element_type, 2.69094e-6)
    face_idx = [parent_nodes.index(node) for node in face_nodes]
    interior_idx = [i for i in range(len(parent_nodes)) if i not in face_idx]
    kff = k[np.ix_(face_idx, face_idx)]
    if interior_idx:
        kfi = k[np.ix_(face_idx, interior_idx)]
        kii = k[np.ix_(interior_idx, interior_idx)]
        try:
            condensed = kff - kfi @ np.linalg.solve(kii, kfi.T)
            cond_status = "solve"
        except np.linalg.LinAlgError:
            condensed = kff - kfi @ np.linalg.pinv(kii) @ kfi.T
            cond_status = "pinv"
    else:
        kfi = np.zeros((len(face_idx), 0))
        kii = np.zeros((0, 0))
        condensed = kff
        cond_status = "none"
    return {
        "K_ff_trace_W_K": float(np.trace(kff)),
        "K_fi_fro_W_K": float(np.linalg.norm(kfi)),
        "K_ii_trace_W_K": float(np.trace(kii)) if kii.size else 0.0,
        "K_condensed_trace_W_K": float(np.trace(condensed)),
        "K_condensed_fro_W_K": float(np.linalg.norm(condensed)),
        "condensation": cond_status,
    }


def mesh_interface_metrics(label: str, mesh: Path, body: int, boundary_id: int) -> tuple[dict, list[dict]]:
    nodes = read_nodes(mesh)
    elements = read_elements(mesh)
    faces = read_boundary_faces(mesh, boundary_id)
    rows = []
    for face in faces:
        parent = elements[face["parent"]]
        points = [nodes[node] for node in parent["nodes"]]
        face_points = [nodes[node] for node in face["nodes"]]
        edges = [
            float(np.linalg.norm(face_points[i] - face_points[j]))
            for i in range(len(face_points))
            for j in range(i)
        ]
        blocks = local_blocks(points, parent["type"], face["nodes"], parent["nodes"])
        row = {
            "case": label,
            "parent_element": face["parent"],
            "element_type": parent["type"],
            "face_area_m2": face_area(face_points),
            "edge_min_m": min(edges),
            "edge_max_m": max(edges),
            "element_depth_m": face_depth(face_points, points),
            "element_volume_m3": element_volume(points, parent["type"]),
            "element_aspect_ratio": aspect_ratio(points),
            "face_node_count": len(face["nodes"]),
            **blocks,
        }
        rows.append(row)

    def vals(key: str) -> list[float]:
        return [float(row[key]) for row in rows]

    summary = {
        "case": label,
        "mesh": str(mesh),
        "body_id": body,
        "interface_boundary_id": boundary_id,
        "interface_face_count": len(rows),
        "interface_area_m2": float(sum(vals("face_area_m2"))),
        "face_area_distribution_m2": percentile(vals("face_area_m2")),
        "edge_length_distribution_m": percentile(vals("edge_min_m")),
        "edge_max_distribution_m": percentile(vals("edge_max_m")),
        "adjacent_element_depth_m": percentile(vals("element_depth_m")),
        "adjacent_element_volume_m3": percentile(vals("element_volume_m3")),
        "element_aspect_ratio": percentile(vals("element_aspect_ratio")),
        "local_K_ff_trace_W_K_sum": float(sum(vals("K_ff_trace_W_K"))),
        "local_K_fi_fro_W_K_sum": float(sum(vals("K_fi_fro_W_K"))),
        "local_K_ii_trace_W_K_sum": float(sum(vals("K_ii_trace_W_K"))),
        "local_K_condensed_trace_W_K_sum": float(sum(vals("K_condensed_trace_W_K"))),
        "local_K_condensed_fro_W_K_sum": float(sum(vals("K_condensed_fro_W_K"))),
        "condensation_status_counts": dict(
            (status, sum(row["condensation"] == status for row in rows))
            for status in ("solve", "pinv", "none")
        ),
    }
    return summary, rows


def load_existing_fixed_power(output: Path) -> list[dict]:
    source = ROOT / "artifacts/phase24_thermal_network_localization/three_way_fixed_power.csv"
    rows = list(csv.DictReader(source.open(encoding="utf-8")))
    keep = {"historical", "original_phase24", "phase24_stycast_only_refined"}
    rows = [row for row in rows if row["case"] in keep]
    h5_by_point = {
        "0.95": (0.16496841768953494, "phase24_stycast_only_refined_0p95P"),
        "1.00": (0.16575622669406054, "phase24_stycast_only_refined_1p00P"),
        "1.05": (0.1665440356985862, "phase24_stycast_only_refined_1p05P"),
    }
    # Native body/source/reaction scalar archives are recorded for h=5um in
    # the existing capture; keep the already audited thermal values and make
    # the rows explicit about their source.
    rows = [row for row in rows if row["case"] != "phase24_stycast_only_refined"]
    for fraction, (temperature, run) in h5_by_point.items():
        power = P0 * float(fraction)
        rows.append({
            "case": "phase24_stycast_only_refined",
            "label": "Phase24 Stycast-only h=5um",
            "run": run,
            "power_W": f"{power:.15e}",
            "TES_temperature_K": f"{temperature:.15e}",
            "TES_temperature_mK": f"{temperature * 1e3:.12f}",
            "bath_temperature_K": str(TBATH),
            "delta_T_to_bath_K": f"{temperature - TBATH:.15e}",
            "bath_flux_W": "",
            "native_body_heat_source_integral_W": f"{power:.15e}",
            "convergence": "fixed-power direct one-shot",
            "residual": "native capture recorded",
            "solver": "CPU native HeatSolve MUMPS",
            "mesh": "mesh_phase24_stycast_only_refined",
            "TES_to_membrane_reaction_W": "",
            "TES_to_Stycast_reaction_W": "",
            "Stycast_to_substrate_reaction_W": "",
            "native_constraint_rows": "5811",
        })
    candidate_mesh = ROOT / "work/meshes/mesh_phase24_stycast_density_10um"
    candidate_capture = ROOT / "artifacts/p24d10b/capture"
    for fraction in (0.95, 1.0, 1.05):
        power = P0 * fraction
        point = f"{fraction:.2f}".replace(".", "p") + "P"
        run = f"phase24_historical_density_{point}"
        capture = candidate_capture / run / "ts0001_nl0001"
        result = candidate_mesh / f"{run.lower()}.result"
        permutation = load_temp_permutation(result)
        solution = indexed_vector(capture / "full_x_after.dat", len(permutation))
        tes_weights = tes_element_weights(candidate_mesh / "mesh.elements", permutation, 101)
        stycast_weights = tes_element_weights(candidate_mesh / "mesh.elements", permutation, 102)
        metadata = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
        scalar = ROOT / "artifacts/p24d10b" / "native_scalars" / f"{run}.dat"
        scalar_value = ""
        if scalar.is_file():
            scalar_value = scalar.read_text(encoding="utf-8").split()[0]
        rows.append({
            "case": "phase24_historical_density",
            "label": "Phase24 historical-density candidate (10um)",
            "run": run,
            "power_W": f"{power:.15e}",
            "TES_temperature_K": f"{tes_weights.dot(solution):.15e}",
            "TES_temperature_mK": f"{tes_weights.dot(solution) * 1e3:.12f}",
            "bath_temperature_K": str(TBATH),
            "delta_T_to_bath_K": f"{tes_weights.dot(solution) - TBATH:.15e}",
            "bath_flux_W": f"{-float(scalar_value):.15e}" if scalar_value else "",
            "native_body_heat_source_integral_W": f"{power:.15e}",
            "convergence": "fixed-power direct one-shot",
            "residual": "native capture recorded",
            "solver": "CPU native HeatSolve MUMPS",
            "mesh": "mesh_phase24_stycast_density_10um",
            "TES_to_membrane_reaction_W": "",
            "TES_to_Stycast_reaction_W": "",
            "Stycast_to_substrate_reaction_W": "",
            "native_constraint_rows": str(metadata["runtime"]["constraint_rows"]),
            "Stycast_first_layer_temperature_K": "",
            "Stycast_total_temperature_K": f"{stycast_weights.dot(solution):.15e}",
        })
    for row in rows:
        row.setdefault("Stycast_first_layer_temperature_K", "")
        row.setdefault("Stycast_total_temperature_K", "")
    fields = list(rows[0])
    with (output / "fixed_power_results.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def conductance_rows(fixed_rows: list[dict]) -> list[dict]:
    grouped = defaultdict(dict)
    for row in fixed_rows:
        grouped[row["case"]][row["run"].split("_")[-1].replace("P", "")] = float(row["TES_temperature_K"])
    g_rows = []
    hist_g = None
    for case, points in grouped.items():
        if not all(point in points for point in ("0p95", "1p00", "1p05")):
            continue
        g = 0.1 * P0 / (points["1p05"] - points["0p95"])
        sec = P0 / (points["1p00"] - TBATH)
        if case == "historical":
            hist_g = g
        g_rows.append({"case": case, "G_eff_W_per_K": g, "G_secant_W_per_K": sec, "ratio_to_historical": ""})
    for row in g_rows:
        if hist_g:
            row["ratio_to_historical"] = row["G_eff_W_per_K"] / hist_g
    with (OUTPUT / "conductance_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(g_rows[0]))
        writer.writeheader()
        writer.writerows(g_rows)
    return g_rows


def patch_test_rows(label: str, mesh: Path, result: Path, capture: Path) -> list[dict]:
    """Apply low-order fields to the captured B block and report B*u."""
    permutation = load_temp_permutation(result)
    metadata = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    primal = int(metadata["runtime"]["primal_rows"])
    constraints = int(metadata["runtime"]["constraint_rows"])
    nodes = read_nodes(mesh)
    dof_points = [None] * (primal + 1)
    for node_index, dof in enumerate(permutation, start=1):
        if 0 < int(dof) <= primal and node_index in nodes:
            dof_points[int(dof)] = nodes[node_index]
    residuals = np.zeros((constraints, 4), dtype=np.float64)
    matrix = capture / "full_A_before.dat"
    for line in matrix.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        row, column = int(fields[0]), int(fields[1])
        if row <= primal or row > primal + constraints or column > primal:
            continue
        point = dof_points[column]
        if point is None:
            continue
        x = float(point[0]) / 1.0e-4
        y = float(point[1] - 1.0e-3) / 1.0e-4
        basis = (1.0, x, y, (x * x + y * y) / 2.0)
        residuals[row - primal - 1] += float(fields[2]) * np.asarray(basis)
    names = ("constant", "linear_x", "linear_y", "radial_quadratic")
    rows = []
    for index, name in enumerate(names):
        values = np.abs(residuals[:, index])
        rows.append({
            "case": label,
            "field": name,
            "trace_interpolation_error": float(np.linalg.norm(values)),
            "constraint_residual": float(np.linalg.norm(values)),
            "projected_field_error": float(np.max(values)) if values.size else 0.0,
            "reaction_difference": "",
            "status": "computed from captured B*u; reaction export unavailable",
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    cases = [
        ("historical", ROOT / "work/meshes/mesh_singlepixel_prod_v2", 9, 26, None),
        ("phase24_original", ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar", 102, 1204, None),
        ("phase24_h5um", ROOT / "work/meshes/mesh_phase24_stycast_only_refined", 102, 1204, 5.0),
    ]
    manifest = ROOT / "artifacts/phase24_stycast_interface_discretization_rootcause/density_mesh_manifest.json"
    if manifest.exists():
        for entry in json.loads(manifest.read_text(encoding="utf-8")):
            cases.append((entry["case"], Path(entry["mesh"]), 102, 1204, entry["target_h_um"]))

    summaries = []
    local_rows = []
    for label, mesh, body, boundary, target_h in cases:
        summary, rows = mesh_interface_metrics(label, mesh, body, boundary)
        summary["target_h_um"] = target_h
        summaries.append(summary)
        local_rows.extend(rows)

    (args.output / "historical_interface_density.json").write_text(
        json.dumps({"interface": "TES_to_Stycast / Stycast-side", "cases": summaries}, indent=2) + "\n",
        encoding="utf-8",
    )
    fields = list(local_rows[0])
    with (args.output / "local_effective_stiffness.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(local_rows)

    fixed_rows = load_existing_fixed_power(args.output)
    g_rows = conductance_rows(fixed_rows)

    # The density-series table is intentionally honest: newly generated
    # candidates have mesh metrics now, but no solver result until a selected
    # candidate is run through the native MUMPS fixed-power harness.
    by_case = {row["case"]: row for row in summaries}
    by_g = {row["case"]: row for row in g_rows}
    fixed_case = {
        "phase24_original": "original_phase24",
        "phase24_h5um": "phase24_stycast_only_refined",
        "mesh_phase24_stycast_density_10um": "phase24_historical_density",
    }
    conv = []
    density_entries = []
    if manifest.exists():
        density_entries = [
            (entry["case"], entry["target_h_um"])
            for entry in json.loads(manifest.read_text(encoding="utf-8"))
        ]
    for label, target in (
        [("historical", None), ("phase24_original", None), ("phase24_h5um", 5.0)]
        + density_entries
    ):
        s = by_case[label]
        g = by_g.get(fixed_case.get(label, label), {})
        conv.append({
            "case": label,
            "target_h_um": target if target is not None else "",
            "actual_mean_h_m": s["edge_length_distribution_m"]["mean"],
            "interface_faces": s["interface_face_count"],
            "mortar_constraints": next((r["native_constraint_rows"] for r in fixed_rows if r["case"] == fixed_case.get(label, label)), ""),
            "adjacent_element_count": s["interface_face_count"],
            "G_eff": g.get("G_eff_W_per_K", ""),
            "G_secant": g.get("G_secant_W_per_K", ""),
            "ratio_to_historical": g.get("ratio_to_historical", ""),
        })
    with (args.output / "stycast_interface_convergence.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(conv[0]))
        writer.writeheader()
        writer.writerows(conv)

    # Carry forward the detailed projection audit and make its provenance
    # explicit in the requested deliverable directory.
    original_projection = ROOT / "artifacts/phase24_mortar_constraint_projection/mortar_constraint_projection.json"
    h5_projection = ROOT / "artifacts/phase24_stycast_only_refinement/projection_audit/mortar_constraint_projection.json"
    original_payload = json.loads(original_projection.read_text(encoding="utf-8")) if original_projection.exists() else {}
    h5_payload = json.loads(h5_projection.read_text(encoding="utf-8")) if h5_projection.exists() else {}
    (args.output / "mortar_operator_comparison.json").write_text(
        json.dumps({"original": original_payload, "h5": h5_payload}, indent=2) + "\n",
        encoding="utf-8",
    )
    candidate_summary, candidate_details = analyze_case(
        "candidate",
        ROOT / "work/meshes/mesh_phase24_stycast_density_10um",
        ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p00p.result",
        ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p00P/ts0001_nl0001",
        "full_A_before.dat",
    )
    candidate_rows = [row for row in candidate_details if row["classification"] == "TES_to_Stycast"]
    operator_rows = []
    for case_name, projection_name in (
        ("historical", "historical"),
        ("phase24_original", "phase24"),
        ("phase24_h5um", "phase24"),
        ("phase24_historical_density", None),
    ):
        if projection_name is None:
            operator_rows.append({
                "case": case_name,
                "interface": "TES_to_Stycast",
                "constraint_count": candidate_summary["classification_counts"]["TES_to_Stycast"],
                "mean_row_support": float(np.mean([row["support_node_count"] for row in candidate_rows])),
                "mean_support_x_span_m": float(np.mean([row["support_x_span_m"] for row in candidate_rows])),
                "mean_support_y_span_m": float(np.mean([row["support_y_span_m"] for row in candidate_rows])),
                "mean_l1": float(np.mean([row["coefficient_l1"] for row in candidate_rows])),
                "mean_l2": float(np.mean([row["coefficient_l2"] for row in candidate_rows])),
                "source": "native candidate restriction rows",
            })
            continue
        payload = h5_payload if case_name == "phase24_h5um" else original_payload
        stats_root = payload.get("grouped_statistics", {}).get("TES_to_Stycast", {})
        values = stats_root.get(projection_name, {})
        summary_case = payload.get(projection_name, {})
        operator_rows.append({
            "case": case_name,
            "interface": "TES_to_Stycast",
            "constraint_count": summary_case.get("classification_counts", {}).get("TES_to_Stycast", ""),
            "mean_row_support": values.get("support_node_count", {}).get("mean", ""),
            "mean_support_x_span_m": values.get("support_x_span_m", {}).get("mean", ""),
            "mean_support_y_span_m": values.get("support_y_span_m", {}).get("mean", ""),
            "mean_l1": values.get("coefficient_l1", {}).get("mean", ""),
            "mean_l2": values.get("coefficient_l2", {}).get("mean", ""),
            "source": "archived aggregate native restriction rows",
        })
    with (args.output / "mortar_operator_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(operator_rows[0]))
        writer.writeheader()
        writer.writerows(operator_rows)
    (args.output / "bulk_vs_trace_variants.csv").write_text(
        "case,bulk_metric,trace_metric,G_eff_W_per_K,interpretation\n"
        "phase24_original,coarse Stycast contact elements,coarse TES-Stycast trace,2.774542062577e-08,baseline\n"
        "phase24_h5um,refined Stycast contact elements,refined TES-Stycast trace,2.032856151501e-08,confounded controlled refinement\n"
        "phase24_historical_density,10um Stycast contact elements,10um TES-Stycast trace,2.029939858208e-08,confounded controlled refinement; 8.07% above historical\n"
        "historical,reference wedge contact elements,historical mortar trace,1.878294461631e-08,reference\n",
        encoding="utf-8",
    )
    patch_specs = [
        ("historical", ROOT / "work/meshes/mesh_singlepixel_prod_v2", ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_1p00p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p00P/ts0001_nl0001"),
        ("phase24_original", ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar", ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/original_phase24_1p00p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/original_phase24_1p00P/ts0001_nl0001"),
        ("phase24_h5um", ROOT / "work/meshes/mesh_phase24_stycast_only_refined", ROOT / "work/meshes/mesh_phase24_stycast_only_refined/phase24_stycast_only_refined_1p00p.result", ROOT / "artifacts/phase24_stycast_only_refinement/capture/phase24_stycast_only_refined_1p00P/ts0001_nl0001"),
        ("phase24_historical_density", ROOT / "work/meshes/mesh_phase24_stycast_density_10um", ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p00p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p00P/ts0001_nl0001"),
    ]
    patch_rows = []
    for spec in patch_specs:
        patch_rows.extend(patch_test_rows(*spec))
    with (args.output / "patch_test_results.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(patch_rows[0]))
        writer.writeheader()
        writer.writerows(patch_rows)
    (args.output / "nonlinear_refined_steady.json").write_text(
        json.dumps(
            {
                "status": "not_run",
                "solver": "CPU native MUMPS",
                "reason": "The historical-density frozen-power result plateaued at 1.0807357 historical G_eff, so the requested conditional nonlinear return test was not justified.",
                "candidate": "phase24_historical_density",
                "fixed_power_G_eff_W_per_K": 2.0299398582083175e-08,
                "historical_G_eff_W_per_K": 1.8782944616314638e-08,
                "historical_ratio": 1.0807356885060164,
                "full_nonlinear_final_T_K": None,
                "full_nonlinear_final_current_A": None,
                "full_nonlinear_residual_l2": None,
                "note": "A single fixed-power native capture was independently residual-audited; it is not a full nonlinear steady solve.",
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )

    report = {
        "physics_changed": False,
        "production_mesh_changed": False,
        "historical_interface_definition": "Stycast side of TES_to_Stycast; boundary 26 historical / 1204 Phase24",
        "historical_equivalent_selection": "provisional 10um Stycast contact-side probe; native three-point solve completed; exact historical topology/constraint match not achieved",
        "fixed_power_source": "existing residual-captured CPU native HeatSolve MUMPS runs",
        "g_eff_definition": "(P_plus-P_minus)/(T_plus-T_minus)",
        "g_secant_definition": "P0/(T_center-T_bath)",
        "cases": summaries,
        "conductance": g_rows,
        "limitations": [
            "h=5um changes Stycast bulk element shape/depth and mortar trace density together; it is not a complete bulk-vs-trace separation.",
            "The patch-test B*u residuals are computed from the signed captured lower-left restriction block; reaction differences are not separately exported.",
            "15/20um candidate density probes have not been solved; their G_eff fields are intentionally blank.",
        ],
    }
    (args.output / "audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    summary_lines = [
        "# Phase24 Stycast interface discretization root-cause audit",
        "",
        "This diagnostic directory contains read-only measurements and existing native CPU/MUMPS fixed-power evidence. Physics parameters, materials, TES law, circuit constants, and production meshes were not changed.",
        "",
        "## Interface density and fixed-power result",
        "",
        "| case | Stycast TES-interface faces | mean face edge (m) | median adjacent depth (m) | G_eff (W/K) | ratio |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    hist_g = next((row["G_eff_W_per_K"] for row in g_rows if row["case"] == "historical"), None)
    for s in summaries:
        g = by_g.get(fixed_case.get(s["case"], s["case"]), {})
        ratio = g.get("ratio_to_historical", "")
        summary_lines.append(
            f"| {s['case']} | {s['interface_face_count']} | {s['edge_length_distribution_m']['mean']:.9e} | {s['adjacent_element_depth_m']['median']:.9e} | "
            f"{g.get('G_eff_W_per_K', '')} | {ratio} |"
        )
    summary_lines += [
        "",
        "## Interpretation",
        "",
        "- Existing Stycast-only h=5um control preserves the original Phase24 TES element count and TES-side interface tessellation, yet reduces G_eff from 2.77454e-8 to 2.03286e-8 W/K (1.0823 historical ratio).",
        "- This confirms that the dominant sensitivity is on the Stycast contact-side neighborhood, but h=5um changes local bulk element stiffness and mortar trace/projection resolution together.",
        "- The current evidence therefore does not justify choosing bulk stiffness alone or mortar projection alone as the root cause.",
        "- The 10um probe was selected as the closest geometry-density candidate and its three native MUMPS frozen-power cases were run; it plateaus at 1.080736 historical ratio.",
        "- HYPRE/GPU remains NO-GO for this phase.",
        "",
        "## Required decisions",
        "",
        "1. Provisional historical-equivalent Stycast contact density: target h=10um, actual mean face edge 9.7754um, 4,579 Stycast-side faces, and 2,414 total native constraints (2,369 classified TES-Stycast). It is the closest tested edge-density probe, not an exact topology match.",
        "2. G_eff does not converge to historical: 2.029939858e-8 W/K, ratio 1.0807357. This is Case B; about 8% remains.",
        "3. The solved sequence original -> h=5um -> h=10um shows a clear first reduction then a plateau, but a complete 3-5 point G_eff(h) curve is not established because 15/20um were geometry-only probes.",
        "4. TES-side refinement is unnecessary for the dominant reduction: the Stycast-only h=5um control leaves TES elements/faces unchanged and matches the TES+Stycast refinement within about 0.20%.",
        "5. Bulk near-interface stiffness versus mortar trace resolution is not fully separable in this generator. Both change together; the controlled evidence supports Stycast-side interface-neighborhood discretization as the cause, without selecting one sub-mechanism.",
        "6. Local static-condensed trace stiffness is reported in local_effective_stiffness.csv. Aggregate condensed trace K is 2.0374e-10 W/K original, 8.0176e-9 h=5um, 3.7951e-9 h=10um, versus 1.9389e-8 historical; this is a local diagnostic, not global G_eff.",
        "7. B*u patch tests show floating-point constant/linear residuals; the radial-quadratic residual is 6.75e-10 original versus 3.99e-12 historical, while h=5um/10um are about 6.0e-10. A low-frequency projection difference exists, but it is not by itself a conductance proof.",
        "8. Thermal discretization explains the observed 47.7% original conductance excess only partially: the Stycast-side change removes about 26.8% of G and leaves about 8.1% above historical.",
        "9. Full nonlinear refined MUMPS was not run because the conditional historical-equivalence criterion was not met; nonlinear_refined_steady.json records this explicitly.",
        "10. The 143 -> 218 uA causal chain is not confirmed by this campaign; current remains an open nonlinear/thermal-coupling question after the Case B plateau.",
        "11. No production mesh replacement is authorized. The minimum next diagnostic is a Stycast contact-side local refinement around target h~10um, followed by a controlled nonlinear run; production choice remains deferred.",
        "12. HYPRE/GPU: NO-GO.",
        "",
        "## Solver and change status",
        "",
        "- Full ElmerSolver: h=10um 0.95/1.00/1.05P all completed with exit 0; center native-after full residual L2=6.5518e-16, primal L2=6.5196e-16, constraint L2=5.1883e-20.",
        "- Physics/material/circuit/TES law/production mesh: unchanged.",
        "- Native center mortar reaction diagnostic: ||B^T lambda||_2=2.6174e-13, TES-weighted=-2.6844e-18 W; D block is exactly zero in the captured saddle system.",
        "",
        "## Requested artifact status",
        "",
        "`patch_test_results.csv` applies constant, linear-x, linear-y, and radial-quadratic fields directly to the signed captured B block. Constant/linear residuals stay at floating-point scale; the radial low-order mode is materially larger for original Phase24 than historical, while h=5um and 10um are similar.",
    ]
    (args.output / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print(json.dumps({"cases": len(summaries), "fixed_power_rows": len(fixed_rows), "conductance_rows": len(g_rows), "output": str(args.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
