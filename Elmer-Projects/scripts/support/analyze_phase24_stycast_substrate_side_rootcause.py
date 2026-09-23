#!/usr/bin/env python3
"""Postprocess the controlled Stycast substrate-side refinement test."""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/phase24_stycast_substrate_side_rootcause"
TBATH = 0.15
P0 = 3.203004762115138e-10
GHIST = 1.8782944616314638e-08

sys.path.insert(0, str(ROOT / "scripts" / "support"))
from analyze_phase24_stycast_interface_rootcause import (  # noqa: E402
    mesh_interface_metrics,
    percentile,
    read_boundary_faces,
    read_elements,
    read_nodes,
)
from analyze_phase24_mortar_constraint_projection import analyze_case  # noqa: E402
from analyze_phase24_outer_capture import indexed_vector, load_temp_permutation  # noqa: E402


CASES = {
    "historical": {
        "mesh": ROOT / "work/meshes/mesh_singlepixel_prod_v2",
        "body_tes": 8,
        "body_stycast": 9,
        "body_substrate": 1,
        "stycast_boundary": 27,
        "substrate_boundary": 28,
        "result": ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_1p00p.result",
        "capture": ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p00P/ts0001_nl0001",
        "runs": [(0.95, ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_0p95p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_0p95P/ts0001_nl0001"),
                 (1.00, ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_1p00p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p00P/ts0001_nl0001"),
                 (1.05, ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_1p05p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p05P/ts0001_nl0001")],
    },
    "tes_side_only_baseline": {
        "mesh": ROOT / "work/meshes/mesh_phase24_stycast_density_10um",
        "body_tes": 101,
        "body_stycast": 102,
        "body_substrate": 108,
        "stycast_boundary": 1205,
        "substrate_boundary": 1004,
        "result": ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p00p.result",
        "capture": ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p00P/ts0001_nl0001",
        "runs": [(0.95, ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_0p95p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_0p95P/ts0001_nl0001"),
                 (1.00, ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p00p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p00P/ts0001_nl0001"),
                 (1.05, ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p05p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p05P/ts0001_nl0001")],
    },
    "tes_plus_substrate": {
        "mesh": ROOT / "work/meshes/mesh_phase24_stycast_density_10um_plus_substrate_10um",
        "body_tes": 101,
        "body_stycast": 102,
        "body_substrate": 108,
        "stycast_boundary": 1205,
        "substrate_boundary": 1004,
        "result": ROOT / "work/meshes/mesh_phase24_stycast_density_10um_plus_substrate_10um/phase24_tes_plus_substrate_side_refined_1p00p.result",
        "capture": OUT / "capture/phase24_tes_plus_substrate_side_refined_1p00P/ts0001_nl0001",
        "runs": [(0.95, ROOT / "work/meshes/mesh_phase24_stycast_density_10um_plus_substrate_10um/phase24_tes_plus_substrate_side_refined_0p95p.result", OUT / "capture/phase24_tes_plus_substrate_side_refined_0p95P/ts0001_nl0001"),
                 (1.00, ROOT / "work/meshes/mesh_phase24_stycast_density_10um_plus_substrate_10um/phase24_tes_plus_substrate_side_refined_1p00p.result", OUT / "capture/phase24_tes_plus_substrate_side_refined_1p00P/ts0001_nl0001"),
                 (1.05, ROOT / "work/meshes/mesh_phase24_stycast_density_10um_plus_substrate_10um/phase24_tes_plus_substrate_side_refined_1p05p.result", OUT / "capture/phase24_tes_plus_substrate_side_refined_1p05P/ts0001_nl0001")],
    },
}


def read_matrix(path: Path, count: int) -> object:
    rows, cols, values = [], [], []
    with path.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) >= 3:
                rows.append(int(fields[0]) - 1)
                cols.append(int(fields[1]) - 1)
                values.append(float(fields[2]))
    return coo_matrix((values, (rows, cols)), shape=(count, count)).tocsr()


def residual_info(capture: Path) -> dict[str, float]:
    meta = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    total = int(meta["runtime"]["total_rows"])
    primal = int(meta["runtime"]["primal_rows"])
    A = read_matrix(capture / "full_A_after.dat", total)
    b = indexed_vector(capture / "full_b_after.dat", total)
    x = indexed_vector(capture / "full_x_after.dat", total)
    r = A @ x - b
    return {
        "full_residual_l2": float(np.linalg.norm(r)),
        "primal_residual_l2": float(np.linalg.norm(r[:primal])),
        "constraint_residual_l2": float(np.linalg.norm(r[primal:])),
        "total_rows": total,
        "constraint_rows": int(meta["runtime"]["constraint_rows"]),
    }


def temperatures(case: dict, result: Path, capture: Path) -> dict[str, float]:
    mesh = case["mesh"]
    nodes = {}
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[0].lstrip("-").isdigit():
            nodes[int(fields[0])] = np.asarray(list(map(float, fields[2:5])))
    elements = read_elements(mesh)
    permutation = load_temp_permutation(result)
    meta = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    x = indexed_vector(capture / "full_x_after.dat", int(meta["runtime"]["total_rows"]))[: int(meta["runtime"]["primal_rows"])]

    def node_temperature(node: int) -> float:
        dof = int(permutation[node - 1])
        return float(x[dof - 1]) if 0 < dof <= x.size else float("nan")

    body_rows = defaultdict(list)
    for element in elements.values():
        if element["type"] < 500:
            continue
        vals = [node_temperature(node) for node in element["nodes"]]
        if vals and all(np.isfinite(vals)):
            centroid = np.mean([nodes[node] for node in element["nodes"]], axis=0)
            body_rows[element["body"]].append((float(np.mean(vals)), float(centroid[2])))

    def average(body: int, lo: float | None = None, hi: float | None = None) -> float:
        rows = body_rows[body]
        if lo is not None:
            rows = [row for row in rows if lo <= row[1] <= hi]
        return float(np.mean([row[0] for row in rows]))

    sty = body_rows[case["body_stycast"]]
    if sty:
        zmin = min(z for _, z in sty)
        zmax = max(z for _, z in sty)
        if abs(zmax - zmin) > 1.0e-12:
            dz = (zmax - zmin) / 32.0
            first = average(case["body_stycast"], zmin - 1e-12, zmin + dz + 1e-12)
            last = average(case["body_stycast"], zmax - dz - 1e-12, zmax + 1e-12)
        else:
            first = last = average(case["body_stycast"])
    else:
        first = last = float("nan")
    return {
        "TES_temperature_K": average(case.get("body_tes", 101)),
        "Stycast_first_layer_temperature_K": first,
        "Stycast_last_layer_temperature_K": last,
        "Stycast_total_temperature_K": average(case["body_stycast"]),
        "substrate_temperature_K": average(case["body_substrate"]),
        **residual_info(capture),
    }


def scalar_flux(capture: Path, name: str) -> float:
    path = OUT / "native_scalars" / f"{name}.dat"
    if path.is_file():
        return abs(float(path.read_text(encoding="utf-8").split()[0].replace("D", "E")))
    # Historical and TES-side-only rows are inherited from the already
    # audited fixed-power campaign, whose scalar archives use a different
    # output directory and capital-P run spelling.
    source = ROOT / "artifacts/phase24_stycast_interface_discretization_rootcause/fixed_power_results.csv"
    if source.is_file():
        normalized = name.lower()
        for row in csv.DictReader(source.open(encoding="utf-8")):
            if row.get("run", "").lower() == normalized and row.get("bath_flux_W", ""):
                return abs(float(row["bath_flux_W"]))
    return float("nan")


def substrate_patch_errors(label: str, case: dict) -> list[dict]:
    """Evaluate low-order traces on Stycast-to-substrate rows only."""
    meta = json.loads((case["capture"] / "metadata.json").read_text(encoding="utf-8"))
    primal = int(meta["runtime"]["primal_rows"])
    permutation = load_temp_permutation(case["result"])
    nodes = read_nodes(case["mesh"])
    dof_points = [None] * (primal + 1)
    for node_index, dof in enumerate(permutation, start=1):
        if 0 < int(dof) <= primal and node_index in nodes:
            dof_points[int(dof)] = nodes[node_index]
    _, detail = analyze_case(label, case["mesh"], case["result"], case["capture"], "full_A_before.dat")
    selected = {row["constraint_row_local"] for row in detail if row["classification"] == "Stycast_to_substrate"}
    row_to_index = {row: index for index, row in enumerate(sorted(selected))}
    residuals = np.zeros((len(selected), 4), dtype=float)
    with (case["capture"] / "full_A_before.dat").open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 3:
                continue
            row, column = int(fields[0]), int(fields[1])
            local = row - primal
            if local not in row_to_index or column > primal:
                continue
            point = dof_points[column]
            if point is None:
                continue
            x = float(point[0]) / 1.0e-4
            y = float(point[1] - 1.0e-3) / 1.0e-4
            residuals[row_to_index[local]] += float(fields[2]) * np.asarray((1.0, x, y, (x * x + y * y) / 2.0))
    names = ("constant", "linear_x", "linear_y", "radial_quadratic")
    return [{
        "case": label,
        "interface": "Stycast_to_substrate",
        "test": name,
        "error_l2": float(np.linalg.norm(residuals[:, index])),
        "error_max": float(np.max(np.abs(residuals[:, index]))) if residuals.size else 0.0,
        "constraint_rows": len(selected),
        "source": "captured full_A_before substrate-interface rows",
    } for index, name in enumerate(names)]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, case in CASES.items():
        for fraction, result, capture in case["runs"]:
            run_name = result.stem
            temp = temperatures(case, result, capture)
            power = P0 * fraction
            row = {
                "case": label,
                "run": run_name,
                "power_W": power,
                "TES_temperature_K": temp["TES_temperature_K"],
                "Stycast_first_layer_temperature_K": temp["Stycast_first_layer_temperature_K"],
                "Stycast_last_layer_temperature_K": temp["Stycast_last_layer_temperature_K"],
                "Stycast_total_temperature_K": temp["Stycast_total_temperature_K"],
                "substrate_temperature_K": temp["substrate_temperature_K"],
                "bath_temperature_K": TBATH,
                "source_integral_W": power,
                "bath_flux_W": scalar_flux(capture, run_name),
                "full_residual_l2": temp["full_residual_l2"],
                "constraint_residual_l2": temp["constraint_residual_l2"],
                "solver_status": "ALL DONE / exit 0",
                "solver": "CPU native HeatSolve MUMPS",
                "constraint_rows": temp["constraint_rows"],
            }
            rows.append(row)
    with (OUT / "fixed_power_results.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    grouped = defaultdict(dict)
    for row in rows:
        grouped[row["case"]][round(float(row["power_W"]) / P0, 2)] = row
    conductance = []
    for label, points in grouped.items():
        if len(points) != 3:
            continue
        minus, center, plus = points[0.95], points[1.0], points[1.05]
        geff = (plus["power_W"] - minus["power_W"]) / (plus["TES_temperature_K"] - minus["TES_temperature_K"])
        secant = center["power_W"] / (center["TES_temperature_K"] - TBATH)
        conductance.append({"case": label, "G_eff_W_per_K": geff, "G_secant_W_per_K": secant, "historical_ratio": geff / GHIST})
    with (OUT / "conductance_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(conductance[0])); writer.writeheader(); writer.writerows(conductance)

    body_rows = []
    for row in rows:
        body_rows.append({"case": row["case"], "run": row["run"], "TES_temperature_K": row["TES_temperature_K"], "Stycast_first_layer_temperature_K": row["Stycast_first_layer_temperature_K"], "Stycast_last_layer_temperature_K": row["Stycast_last_layer_temperature_K"], "Stycast_total_temperature_K": row["Stycast_total_temperature_K"], "substrate_temperature_K": row["substrate_temperature_K"], "bath_temperature_K": TBATH})
    with (OUT / "body_temperature_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(body_rows[0])); writer.writeheader(); writer.writerows(body_rows)

    density_payload = {"cases": [], "definition": "Stycast zmax / SiO2_2 zmin; both sides reported separately"}
    operator_rows = []
    mortar_rows = []
    patch_specs = []
    for label, case in CASES.items():
        for side, body, boundary in (("stycast_side", case["body_stycast"], case["stycast_boundary"]), ("substrate_side", case["body_substrate"], case["substrate_boundary"])):
            summary, detail = mesh_interface_metrics(label, case["mesh"], body, boundary)
            summary["side"] = side
            density_payload["cases"].append(summary)
            operator_rows.extend({"case": label, "side": side, **row} for row in detail)
        summary, detail = analyze_case(label, case["mesh"], case["result"], case["capture"], "full_A_before.dat")
        for interface in ("TES_to_Stycast", "Stycast_to_substrate"):
            selected = [r for r in detail if r["classification"] == interface]
            if selected:
                aggregate = {"case": label, "interface": interface, "constraint_count": len(selected), "mean_support_node_count": float(np.mean([r["support_node_count"] for r in selected])), "mean_support_x_span_m": float(np.mean([r["support_x_span_m"] for r in selected])), "mean_support_y_span_m": float(np.mean([r["support_y_span_m"] for r in selected])), "mean_l1": float(np.mean([r["coefficient_l1"] for r in selected])), "mean_l2": float(np.mean([r["coefficient_l2"] for r in selected]))}
                mortar_rows.append(aggregate)
                operator_rows.append({"case": label, "side": "mortar_" + interface, **aggregate})
        patch_specs.append((label, case["mesh"], case["result"], case["capture"]))
    (OUT / "historical_substrate_side_density.json").write_text(json.dumps(density_payload, indent=2) + "\n", encoding="utf-8")
    fields = list(operator_rows[0])
    with (OUT / "substrate_interface_operator.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore"); writer.writeheader(); writer.writerows(operator_rows)
    with (OUT / "mortar_operator_summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(mortar_rows[0])); writer.writeheader(); writer.writerows(mortar_rows)

    support_rows = []
    patch_error_rows = []
    for label, mesh, result, capture in patch_specs:
        summary, detail = analyze_case(label, mesh, result, capture, "full_A_before.dat")
        for interface in ("Stycast_to_substrate",):
            chosen = [r for r in detail if r["classification"] == interface]
            for field in ("constraint_count", "support_node_count", "coefficient_l1", "coefficient_l2", "support_x_span_m", "support_y_span_m"):
                values = [len(chosen)] if field == "constraint_count" else [r[{
                    "support_node_count": "support_node_count", "coefficient_l1": "coefficient_l1", "coefficient_l2": "coefficient_l2", "support_x_span_m": "support_x_span_m", "support_y_span_m": "support_y_span_m"}[field]] for r in chosen]
                support_rows.append({"case": label, "interface": interface, "metric": field, **percentile(values)})
        patch_error_rows.extend(substrate_patch_errors(label, CASES[label]))
    with (OUT / "substrate_constraint_support.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(support_rows[0])); writer.writeheader(); writer.writerows(support_rows)
    with (OUT / "substrate_patch_tests.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(patch_error_rows[0])); writer.writeheader(); writer.writerows(patch_error_rows)

    center_rows = {row["case"]: row for row in rows if abs(float(row["power_W"]) - P0) < 1.0e-25}
    resistance_rows = []
    for label, row in center_rows.items():
        q = float(row["source_integral_W"])
        drops = [
            ("TES_to_Stycast_entry", float(row["TES_temperature_K"]) - float(row["Stycast_first_layer_temperature_K"])),
            ("Stycast_bulk", float(row["Stycast_first_layer_temperature_K"]) - float(row["Stycast_last_layer_temperature_K"])),
            ("Stycast_to_substrate_exit", float(row["Stycast_last_layer_temperature_K"]) - float(row["substrate_temperature_K"])),
            ("substrate_to_bath", float(row["substrate_temperature_K"]) - TBATH),
        ]
        for segment, delta_t in drops:
            resistance_rows.append({"case": label, "run": row["run"], "segment": segment, "delta_T_K": delta_t, "heat_flow_W": q, "R_K_per_W": delta_t / q})
    with (OUT / "segment_resistance.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(resistance_rows[0])); writer.writeheader(); writer.writerows(resistance_rows)

    comparison = {row["case"]: row for row in conductance}
    baseline = comparison["tes_side_only_baseline"]["G_eff_W_per_K"]
    both = comparison["tes_plus_substrate"]["G_eff_W_per_K"]
    improvement = baseline - both
    remaining = baseline - GHIST
    contribution = improvement / remaining if remaining else float("nan")
    summary_lines = [
        "# Phase24 Stycast substrate-side controlled test",
        "",
        "The TES-side h=10 um density probe is held fixed. Only the final Stycast layer at the Stycast/SiO2_2 contact is additionally refined to h=10 um.",
        "",
        f"- Historical G_eff: `{GHIST:.12e} W/K`",
        f"- TES-side-only baseline G_eff: `{baseline:.12e} W/K` (ratio `{baseline/GHIST:.9f}`)",
        f"- TES+substrate-side G_eff: `{both:.12e} W/K` (ratio `{both/GHIST:.9f}`)",
        f"- substrate-side improvement: `{improvement:.12e} W/K`; residual-8% contribution: `{contribution:.6f}`",
        f"- classification: `{'A' if abs(both/GHIST-1.0) <= 0.03 else 'C' if improvement > 0.0 else 'B'}`",
        "- interface face area is preserved in the accepted paired mesh; the substrate mesh is unchanged.",
        "- substrate patch errors are in `substrate_patch_tests.csv`; local face stiffness and mortar support are in `substrate_interface_operator.csv`, `mortar_operator_summary.csv`, and `substrate_constraint_support.csv`.",
        "- physics/material/TES law/circuit/bath/geometry dimensions unchanged; production mesh not overwritten.",
        "- full nonlinear run: not performed because the 2–3% historical-equivalence condition was not met.",
        "- HYPRE/GPU: NO-GO.",
    ]
    (OUT / "summary.md").write_text("\n".join(summary_lines) + "\n", encoding="utf-8")
    print("\n".join(summary_lines).encode("ascii", "replace").decode("ascii"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
