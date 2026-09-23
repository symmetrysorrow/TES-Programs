#!/usr/bin/env python3
"""Audit native mortar constraint projection rows from saved full matrices.

The script does not rerun Elmer and does not modify solver input. It maps each
restriction row back to mesh boundary node sets using the saved Temperature
permutation, then compares support widths and coefficient scales between the
historical and Phase24 interfaces.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
HIST_MESH = ROOT / "work/meshes/mesh_singlepixel_prod_v2"
P24_MESH = ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar"
HIST_RESULT = HIST_MESH / "historical_1p00p.result"
P24_RESULT = P24_MESH / "original_phase24_1p00p.result"
HIST_CAPTURE = ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p00P/ts0001_nl0001"
P24_CAPTURE = ROOT / "artifacts/phase24_thermal_network_localization/capture/original_phase24_1p00P/ts0001_nl0001"
OUTPUT = ROOT / "artifacts/phase24_mortar_constraint_projection"


def load_permutation(result: Path) -> tuple[dict[int, int], int]:
    lines = result.read_text(encoding="utf-8", errors="replace").splitlines()
    for index, line in enumerate(lines):
        if line.strip().lower() != "temperature":
            continue
        cursor = index + 1
        while cursor < len(lines) and not lines[cursor].strip().lower().startswith("perm:"):
            cursor += 1
        if cursor == len(lines):
            continue
        count = int(lines[cursor].split()[1])
        permutation = {}
        for row in range(cursor + 1, cursor + 1 + count):
            fields = lines[row].split()
            permutation[int(fields[0])] = int(fields[1])
        return permutation, count
    raise RuntimeError(f"Temperature Perm table not found in {result}")


def boundary_nodes(mesh: Path, boundary_id: int) -> set[int]:
    nodes: set[int] = set()
    for line in mesh.joinpath("mesh.boundary").read_text(
        encoding="utf-8", errors="replace"
    ).splitlines()[1:]:
        fields = line.split()
        if len(fields) < 6 or int(fields[1]) != boundary_id:
            continue
        face_type = int(fields[4])
        count = 3 if face_type == 303 else 4 if face_type == 404 else 0
        nodes.update(map(int, fields[5:5 + count]))
    return nodes


def read_runtime(capture: Path) -> tuple[int, int]:
    metadata = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    runtime = metadata["runtime"]
    return int(runtime["primal_rows"]), int(runtime["constraint_rows"])


def read_constraint_rows(
    matrix: Path,
    primal_rows: int,
    constraint_rows: int,
) -> tuple[dict[int, dict[int, float]], int, int]:
    rows: dict[int, dict[int, float]] = defaultdict(dict)
    records = 0
    constraint_constraint_records = 0
    with matrix.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 3:
                continue
            row = int(fields[0])
            column = int(fields[1])
            if row <= primal_rows or row > primal_rows + constraint_rows:
                continue
            value = float(fields[2])
            records += 1
            if column > primal_rows:
                constraint_constraint_records += 1
                continue
            row_map = rows[row - primal_rows]
            row_map[column] = row_map.get(column, 0.0) + value
    return rows, records, constraint_constraint_records


def classify_row(
    row: dict[int, float],
    inverse_permutation: dict[int, int],
    interface_nodes: dict[str, tuple[set[int], set[int]]],
) -> tuple[str, dict[str, int]]:
    node_ids = {
        inverse_permutation[column]
        for column in row
        if column in inverse_permutation
    }
    scores = {}
    side_counts = {}
    for name, (side_a, side_b) in interface_nodes.items():
        count_a = len(node_ids & side_a)
        count_b = len(node_ids & side_b)
        side_counts[name] = (count_a, count_b)
        scores[name] = count_a + count_b
    best = max(scores, key=scores.get)
    if scores[best] == 0:
        return "unknown", {"side_a_nodes": 0, "side_b_nodes": 0}
    return best, {
        "side_a_nodes": side_counts[best][0],
        "side_b_nodes": side_counts[best][1],
    }


def analyze_case(
    label: str,
    mesh: Path,
    result: Path,
    capture: Path,
    matrix_name: str,
) -> tuple[dict, list[dict]]:
    permutation, primal_rows = load_permutation(result)
    inverse = {dof: node for node, dof in permutation.items() if dof > 0}
    runtime_primal, constraint_rows = read_runtime(capture)
    if runtime_primal != primal_rows:
        raise RuntimeError("primal row metadata mismatch")

    if label == "historical":
        interface_specs = {
            "TES_to_Stycast": ({25}, {26}),
            "Stycast_to_substrate": ({27}, {28}),
        }
    else:
        interface_specs = {
            "TES_to_Stycast": ({1105}, {1204}),
            "Stycast_to_substrate": ({1205}, {1004}),
        }
    interface_nodes = {
        name: (
            set().union(*(boundary_nodes(mesh, boundary) for boundary in side_a)),
            set().union(*(boundary_nodes(mesh, boundary) for boundary in side_b)),
        )
        for name, (side_a, side_b) in interface_specs.items()
    }

    matrix_path = capture / matrix_name
    rows, matrix_records, cc_records = read_constraint_rows(
        matrix_path, primal_rows, constraint_rows
    )
    detail_rows = []
    for row_index in range(1, constraint_rows + 1):
        row = rows.get(row_index, {})
        classification, side_counts = classify_row(row, inverse, interface_nodes)
        node_ids = {
            inverse[column]
            for column in row
            if column in inverse
        }
        values = np.asarray(list(row.values()), dtype=np.float64)
        points = []
        if node_ids:
            node_file = mesh / "mesh.nodes"
            if not hasattr(analyze_case, "_node_cache"):
                analyze_case._node_cache = {}
            if str(mesh) not in analyze_case._node_cache:
                cache = {}
                for line in node_file.read_text(
                    encoding="utf-8", errors="replace"
                ).splitlines()[1:]:
                    fields = line.split()
                    if fields:
                        cache[int(fields[0])] = tuple(map(float, fields[2:5]))
                analyze_case._node_cache[str(mesh)] = cache
            coordinates = analyze_case._node_cache[str(mesh)]
            points = [coordinates[node_id] for node_id in node_ids if node_id in coordinates]
        x_span = max(point[0] for point in points) - min(point[0] for point in points) if points else 0.0
        y_span = max(point[1] for point in points) - min(point[1] for point in points) if points else 0.0
        detail_rows.append(
            {
                "case": label,
                "constraint_row_local": row_index,
                "classification": classification,
                "coefficient_count": len(row),
                "support_node_count": len(node_ids),
                "side_a_node_count": side_counts["side_a_nodes"],
                "side_b_node_count": side_counts["side_b_nodes"],
                "coefficient_l1": float(np.sum(np.abs(values))) if values.size else 0.0,
                "coefficient_l2": float(np.linalg.norm(values)) if values.size else 0.0,
                "coefficient_max_abs": float(np.max(np.abs(values))) if values.size else 0.0,
                "support_x_span_m": x_span,
                "support_y_span_m": y_span,
            }
        )

    summary = {
        "case": label,
        "mesh": str(mesh),
        "result": str(result),
        "capture": str(capture),
        "matrix": str(matrix_path),
        "primal_rows": primal_rows,
        "constraint_rows": constraint_rows,
        "matrix_constraint_records": matrix_records,
        "constraint_constraint_records": cc_records,
        "interface_boundary_node_counts": {
            name: {
                "side_a": len(side_a),
                "side_b": len(side_b),
            }
            for name, (side_a, side_b) in interface_nodes.items()
        },
        "classification_counts": {
            name: sum(row["classification"] == name for row in detail_rows)
            for name in interface_nodes
        }
        | {"unknown": sum(row["classification"] == "unknown" for row in detail_rows)},
    }
    return summary, detail_rows


def stats(rows: list[dict], key: str) -> dict[str, float | int]:
    values = np.asarray([row[key] for row in rows], dtype=np.float64)
    return {
        "count": int(values.size),
        "sum": float(values.sum()),
        "mean": float(values.mean()),
        "p50": float(np.percentile(values, 50)),
        "p95": float(np.percentile(values, 95)),
        "max": float(values.max()),
        "min": float(values.min()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--historical-mesh", type=Path, default=HIST_MESH)
    parser.add_argument("--phase24-mesh", type=Path, default=P24_MESH)
    parser.add_argument("--historical-result", type=Path, default=HIST_RESULT)
    parser.add_argument("--phase24-result", type=Path, default=P24_RESULT)
    parser.add_argument("--historical-capture", type=Path, default=HIST_CAPTURE)
    parser.add_argument("--phase24-capture", type=Path, default=P24_CAPTURE)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    historical, historical_rows = analyze_case(
        "historical",
        args.historical_mesh,
        args.historical_result,
        args.historical_capture,
        "full_A_before.dat",
    )
    phase24, phase24_rows = analyze_case(
        "phase24",
        args.phase24_mesh,
        args.phase24_result,
        args.phase24_capture,
        "full_A_before.dat",
    )
    all_rows = historical_rows + phase24_rows
    with (args.output / "constraint_row_projection.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(all_rows[0]))
        writer.writeheader()
        writer.writerows(all_rows)

    grouped = {}
    for interface in ("TES_to_Stycast", "Stycast_to_substrate", "unknown"):
        hrows = [row for row in historical_rows if row["classification"] == interface]
        prows = [row for row in phase24_rows if row["classification"] == interface]
        if not hrows or not prows:
            continue
        grouped[interface] = {
            "historical": {
                key: stats(hrows, key)
                for key in (
                    "coefficient_count",
                    "support_node_count",
                    "side_a_node_count",
                    "side_b_node_count",
                    "coefficient_l1",
                    "coefficient_l2",
                    "coefficient_max_abs",
                    "support_x_span_m",
                    "support_y_span_m",
                )
            },
            "phase24": {
                key: stats(prows, key)
                for key in (
                    "coefficient_count",
                    "support_node_count",
                    "side_a_node_count",
                    "side_b_node_count",
                    "coefficient_l1",
                    "coefficient_l2",
                    "coefficient_max_abs",
                    "support_x_span_m",
                    "support_y_span_m",
                )
            },
            "sum_ratios_phase24_over_historical": {
                key: (
                    stats(prows, key)["sum"]
                    / max(stats(hrows, key)["sum"], 1.0e-300)
                )
                for key in ("coefficient_l1", "coefficient_l2", "coefficient_max_abs")
            },
        }

    report = {
        "method": {
            "matrix": "full_A_before.dat lower-left restriction block",
            "classification": "constraint support node intersection with each interface boundary node set",
            "coefficient_statistics": "duplicate row/column entries aggregated before statistics",
        },
        "historical": historical,
        "phase24": phase24,
        "grouped_statistics": grouped,
    }
    (args.output / "mortar_constraint_projection.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# Phase24 mortar constraint projection audit",
        "",
        "| case | total constraint rows | TES–Stycast | Stycast–substrate | unknown |",
        "|---|---:|---:|---:|---:|",
        f"| historical | {historical['constraint_rows']} | {historical['classification_counts']['TES_to_Stycast']} | {historical['classification_counts']['Stycast_to_substrate']} | {historical['classification_counts']['unknown']} |",
        f"| Phase24 | {phase24['constraint_rows']} | {phase24['classification_counts']['TES_to_Stycast']} | {phase24['classification_counts']['Stycast_to_substrate']} | {phase24['classification_counts']['unknown']} |",
        "",
        "| interface | metric | historical mean | Phase24 mean | Phase24 / historical |",
        "|---|---|---:|---:|---:|",
    ]
    for interface, values in grouped.items():
        if interface == "unknown":
            continue
        for key in ("coefficient_count", "support_node_count", "coefficient_l2", "support_x_span_m", "support_y_span_m"):
            h = values["historical"][key]["mean"]
            p = values["phase24"][key]["mean"]
            lines.append(f"| {interface} | {key} | {h:.9e} | {p:.9e} | {p / max(h, 1.0e-300):.9e} |")
    lines.extend(
        [
            "",
            "Constraint rows are classified from their primal support nodes. "
            "The raw constraint-row count is not itself a conductance; it diagnoses "
            "mortar tessellation and projection stencil changes. Phase24 has fewer "
            "but much larger support patches; the total coefficient L1/L2 mass is "
            "reported in the JSON to distinguish area aggregation from an actual "
            "projection-strength change.",
        ]
    )
    (args.output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
