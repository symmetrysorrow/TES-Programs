#!/usr/bin/env python3
"""Record the valid saddle-system interpretation and the interface-refinement control.

The captured Elmer restriction systems have a zero lower-right D block, so the
literal B^T D^-1 B expression is undefined. This report verifies that fact and
compares the existing historical, original Phase24, and TES/Stycast-interface
refined Phase24 fixed-power controls.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "support"))
from analyze_phase24_outer_capture import (  # noqa: E402
    indexed_vector,
    load_temp_permutation,
    tes_element_weights,
)


P0 = 3.203004762115138e-10
OUTPUT = ROOT / "artifacts/phase24_mortar_schur_control"


CASES = {
    "historical": {
        "mesh": ROOT / "work/meshes/mesh_singlepixel_prod_v2",
        "result_pattern": "historical_{point}.result",
        "capture_pattern": "historical_{power}P",
        "body": 8,
        "points": {
            "minus": "0p95p",
            "center": "1p00p",
            "plus": "1p05p",
        },
    },
    "original_phase24": {
        "mesh": ROOT / "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar",
        "result_pattern": "original_phase24_{point}.result",
        "capture_pattern": "original_phase24_{power}P",
        "body": 101,
        "points": {
            "minus": "0p95p",
            "center": "1p00p",
            "plus": "1p05p",
        },
    },
    "phase24_interface_refined": {
        "mesh": ROOT / "work/meshes/mesh_phase24_membrane_stycast_diag_v2",
        "capture_root": ROOT / "artifacts/phase24_thermal_network_localization/capture",
        "result_pattern": "phase24_membrane_stycast_variant_v2_{point}.result",
        "capture_pattern": "phase24_membrane_stycast_variant_v2_{power}P",
        "body": 101,
        "points": {
            "minus": "0p95p",
            "center": "1p00p",
            "plus": "1p05p",
        },
    },
    "phase24_stycast_only_refined": {
        "mesh": ROOT / "work/meshes/mesh_phase24_stycast_only_refined",
        "capture_root": ROOT / "artifacts/phase24_stycast_only_refinement/capture",
        "result_pattern": "phase24_stycast_only_refined_{point}.result",
        "capture_pattern": "phase24_stycast_only_refined_{power}P",
        "body": 101,
        "points": {
            "minus": "0p95p",
            "center": "1p00p",
            "plus": "1p05p",
        },
    },
}


def scan_d_block(matrix: Path, primal: int, constraints: int) -> dict:
    records = 0
    nonzero = 0
    max_abs = 0.0
    with matrix.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 3:
                continue
            row, column = int(fields[0]), int(fields[1])
            if row <= primal or row > primal + constraints or column <= primal:
                continue
            records += 1
            value = abs(float(fields[2]))
            max_abs = max(max_abs, value)
            if value != 0.0:
                nonzero += 1
    return {
        "records": records,
        "nonzero_records": nonzero,
        "max_abs": max_abs,
        "is_zero": nonzero == 0,
    }


def read_dimensions(capture: Path) -> tuple[int, int]:
    metadata = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    runtime = metadata["runtime"]
    return int(runtime["primal_rows"]), int(runtime["constraint_rows"])


def tes_temperature(mesh: Path, result: Path, capture: Path, body: int) -> float:
    permutation = load_temp_permutation(result)
    solution = indexed_vector(capture / "full_x_after.dat", len(permutation))
    weights = tes_element_weights(mesh / "mesh.elements", permutation, body)
    return float(weights.dot(solution))


def case_result(case_name: str, case: dict) -> dict:
    mesh = case["mesh"]
    temperatures = {}
    d_blocks = {}
    for label, point in case["points"].items():
        result_point = point
        # Result files use a trailing lower-case ``p`` while capture folders
        # use the numeric power label followed by one upper-case ``P``.
        power = point[:-1]
        result = mesh / case["result_pattern"].format(point=result_point)
        capture = (
            case.get("capture_root", ROOT / "artifacts/phase24_thermal_network_localization/capture")
            / case["capture_pattern"].format(power=power)
            / "ts0001_nl0001"
        )
        temperatures[label] = tes_temperature(mesh, result, capture, case["body"])
        primal, constraints = read_dimensions(capture)
        d_blocks[label] = scan_d_block(
            capture / "full_A_before.dat", primal, constraints
        )
    delta_power = 0.10 * P0
    derivative = delta_power / (temperatures["plus"] - temperatures["minus"])
    center_secant = P0 / (temperatures["center"] - 0.15)
    return {
        "mesh": str(mesh),
        "temperatures_K": temperatures,
        "temperatures_mK": {key: value * 1000.0 for key, value in temperatures.items()},
        "G_derivative_W_K": derivative,
        "G_secant_to_bath_W_K": center_secant,
        "D_blocks": d_blocks,
        "constraint_rows_center": read_dimensions(
            case.get("capture_root", ROOT / "artifacts/phase24_thermal_network_localization/capture")
            / case["capture_pattern"].format(power="1p00")
            / "ts0001_nl0001"
        )[1],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    results = {name: case_result(name, case) for name, case in CASES.items()}
    historical_g = results["historical"]["G_derivative_W_K"]
    original_g = results["original_phase24"]["G_derivative_W_K"]
    refined_g = results["phase24_interface_refined"]["G_derivative_W_K"]
    stycast_only_g = results["phase24_stycast_only_refined"]["G_derivative_W_K"]
    comparison = {
        "original_over_historical_G": original_g / historical_g,
        "refined_over_historical_G": refined_g / historical_g,
        "refined_over_original_G": refined_g / original_g,
        "original_to_refined_reduction": 1.0 - refined_g / original_g,
        "stycast_only_over_historical_G": stycast_only_g / historical_g,
        "stycast_only_over_original_G": stycast_only_g / original_g,
        "stycast_only_over_interface_refined_G": stycast_only_g / refined_g,
        "stycast_only_vs_interface_refined_relative_difference": stycast_only_g / refined_g - 1.0,
    }
    report = {
        "literal_expression": "B^T D^-1 B",
        "system_form": "[K B^T; B D] [x; lambda] = [f; g]",
        "D_block_result": "D is exactly zero in all captured cases; D^-1 does not exist.",
        "valid_saddle_system_statement": "For D=0, eliminate x to obtain the constraint Schur complement -B K^-1 B^T. A primal additive B^T D^-1 B term is not defined.",
        "cases": results,
        "comparison": comparison,
        "control": {
            "mesh": str(CASES["phase24_interface_refined"]["mesh"]),
            "notes": "Existing diagnostic mesh: TES top and first 1 um of Stycast locally refined; outer Stycast-substrate interface remains original Phase24 coarse topology.",
            "constraint_rows_center_original_phase24": results["original_phase24"]["constraint_rows_center"],
            "constraint_rows_center_interface_refined": results["phase24_interface_refined"]["constraint_rows_center"],
            "constraint_rows_center_stycast_only_refined": results["phase24_stycast_only_refined"]["constraint_rows_center"],
        },
    }
    (args.output / "mortar_schur_control.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# Phase24 mortar Schur/control audit",
        "",
        "## D-block check",
        "",
        "- The captured saddle systems have an exactly zero lower-right D block.",
        "- Therefore the literal B^T D^-1 B expression is undefined.",
        "- The valid constraint Schur complement after eliminating primal x is -B K^-1 B^T.",
        "",
        "| case | center constraints | T(0.95P) mK | T(1.00P) mK | T(1.05P) mK | G derivative W/K | G / historical |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, result in results.items():
        lines.append(
            f"| {name} | {result['constraint_rows_center']} | "
            f"{result['temperatures_mK']['minus']:.9f} | "
            f"{result['temperatures_mK']['center']:.9f} | "
            f"{result['temperatures_mK']['plus']:.9f} | "
            f"{result['G_derivative_W_K']:.12e} | "
            f"{result['G_derivative_W_K'] / historical_g:.9f} |"
        )
    lines.extend(
        [
            "",
            "## Controlled interface refinement",
            "",
            f"- Original Phase24 G: {original_g:.12e} W/K",
            f"- TES/Stycast interface-refined G: {refined_g:.12e} W/K",
            f"- Reduction from original Phase24: {comparison['original_to_refined_reduction']:.9f}",
            f"- Refined/historical G ratio: {comparison['refined_over_historical_G']:.9f}",
            f"- Stycast-only refined G: {stycast_only_g:.12e} W/K",
            f"- Stycast-only/historical G ratio: {comparison['stycast_only_over_historical_G']:.9f}",
            f"- Stycast-only vs TES+Stycast refined: {comparison['stycast_only_vs_interface_refined_relative_difference']:.9f}",
            "",
            "The refinement changes only the local TES top / first 1 um Stycast mesh "
            "in this existing diagnostic case. The outer Stycast-substrate interface "
            "remains coarse. The large G reduction therefore strongly implicates the "
            "TES-Stycast interface discretization or its local FE/mortar coupling, "
            "although it is not a pure algebra-only mortar test because the local "
            "element mesh is refined as well.",
            "The separate Stycast-side-only control leaves the TES element mesh at the original Phase24 density and still gives essentially the same G as the TES+Stycast refinement. This isolates the dominant sensitivity to the Stycast-side interface neighborhood, while the 5 um control is somewhat denser than the historical Stycast face density.",
        ]
    )
    (args.output / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
