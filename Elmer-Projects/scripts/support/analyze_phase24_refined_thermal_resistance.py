"""Decompose the historical-vs-refined Phase24 thermal resistance budget.

This is a read-only postprocessor for already captured CPU/MUMPS fixed-power
solutions.  It does not rerun ElmerSolver and does not modify meshes or SIFs.

The endpoint budget is deliberately explicit about the unresolved middle
network:

    TES -> Membrane_SiNx -> bath-side SiO2_2 -> bath

The Membrane_SiNx-to-SiO2_2 term contains the intervening substrate stack and
possible parallel paths.  Per-body temperatures are emitted separately so
that this aggregate is not mistaken for a unique single-interface resistance.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_refined_thermal_resistance"
P0 = 3.203004762115138e-10
TBATH = 0.150

sys.path.insert(0, str(ROOT / "scripts" / "support"))
from analyze_phase24_mortar_element_stiffness import element_volume, read_elements  # noqa: E402
from analyze_phase24_outer_capture import indexed_vector, load_temp_permutation  # noqa: E402


CASES = {
    "historical": {
        "mesh": ROOT / "work/meshes/mesh_singlepixel_prod_v2",
        "body_ids": {"TES": 8, "Membrane_SiNx": 7, "Stycast": 9, "SiO2_1": 3, "Si_1": 4, "SiNx": 6, "Si_2": 2, "SiO2_2": 1, "Membrane_Si1": 5},
        "runs": [
            (0.95, ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_0p95p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_0p95P/ts0001_nl0001"),
            (1.00, ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_1p00p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p00P/ts0001_nl0001"),
            (1.05, ROOT / "work/meshes/mesh_singlepixel_prod_v2/historical_1p05p.result", ROOT / "artifacts/phase24_thermal_network_localization/capture/historical_1p05P/ts0001_nl0001"),
        ],
    },
    "refined_mortar": {
        "mesh": ROOT / "work/meshes/mesh_phase24_stycast_density_10um",
        "body_ids": {"TES": 101, "Membrane_SiNx": 103, "Stycast": 102, "SiO2_1": 104, "Si_1": 105, "SiNx": 106, "Si_2": 107, "SiO2_2": 108, "Membrane_Si1": 109},
        "runs": [
            (0.95, ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_0p95p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_0p95P/ts0001_nl0001"),
            (1.00, ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p00p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p00P/ts0001_nl0001"),
            (1.05, ROOT / "work/meshes/mesh_phase24_stycast_density_10um/phase24_historical_density_1p05p.result", ROOT / "artifacts/p24d10b/capture/phase24_historical_density_1p05P/ts0001_nl0001"),
        ],
    },
}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def body_temperatures(mesh: Path, result: Path, capture: Path, body_ids: dict[str, int]) -> dict[str, float]:
    meta = json.loads((capture / "metadata.json").read_text(encoding="utf-8"))
    total_rows = int(meta["runtime"]["total_rows"])
    primal_rows = int(meta["runtime"]["primal_rows"])
    x = indexed_vector(capture / "full_x_after.dat", total_rows)[:primal_rows]
    permutation = load_temp_permutation(result)
    nodes: dict[int, np.ndarray] = {}
    # ElmerGrid mesh.nodes has no header in these generated meshes.
    for line in (mesh / "mesh.nodes").read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 5 and fields[0].lstrip("-").isdigit():
            nodes[int(fields[0])] = np.asarray(list(map(float, fields[2:5])), dtype=float)
    elements, _ = read_elements(mesh)
    target_ids = set(body_ids.values())
    sums: dict[int, float] = defaultdict(float)
    volumes: dict[int, float] = defaultdict(float)
    for element in elements.values():
        body = int(element["body"])
        if body not in target_ids or int(element["type"]) < 500:
            continue
        dofs = [int(permutation[node - 1]) for node in element["nodes"]]
        if not dofs or any(dof <= 0 or dof > x.size for dof in dofs):
            continue
        points = [nodes[node] for node in element["nodes"]]
        volume = float(element_volume(points, int(element["type"])))
        temperature = float(np.mean([x[dof - 1] for dof in dofs]))
        sums[body] += temperature * volume
        volumes[body] += volume
    reverse = {body: name for name, body in body_ids.items()}
    return {name: sums[body] / volumes[body] for body, name in reverse.items() if volumes[body] > 0.0}


def slope(minus: float, plus: float) -> float:
    return (plus - minus) / (0.10 * P0)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    profile_rows: list[dict[str, Any]] = []
    grouped: dict[str, dict[float, dict[str, float]]] = defaultdict(dict)
    for case_name, case in CASES.items():
        for fraction, result, capture in case["runs"]:
            temps = body_temperatures(case["mesh"], result, capture, case["body_ids"])
            grouped[case_name][fraction] = temps
            for body, temperature in temps.items():
                profile_rows.append({
                    "case": case_name,
                    "power_fraction": fraction,
                    "power_W": fraction * P0,
                    "body": body,
                    "temperature_K": temperature,
                    "temperature_mK": temperature * 1.0e3,
                    "bath_temperature_K": TBATH,
                })
    write_csv(OUT / "body_temperature_profile.csv", profile_rows)

    body_slope_rows: list[dict[str, Any]] = []
    for case_name, points in grouped.items():
        for body in sorted(points[1.00]):
            dt_dp = slope(points[0.95][body], points[1.05][body])
            body_slope_rows.append({
                "case": case_name,
                "body": body,
                "dT_dP_K_per_W": dt_dp,
                "body_endpoint_dT_dP_K_per_W": dt_dp,
            })
    write_csv(OUT / "body_temperature_slopes.csv", body_slope_rows)

    budget_rows: list[dict[str, Any]] = []
    for case_name, points in grouped.items():
        minus, center, plus = points[0.95], points[1.00], points[1.05]
        segments = [
            ("TES_to_Membrane_SiNx", "TES", "Membrane_SiNx"),
            ("Membrane_SiNx_to_bath_side_SiO2_2_including_stack", "Membrane_SiNx", "SiO2_2"),
            ("bath_side_SiO2_2_to_bath", "SiO2_2", None),
            ("total_TES_to_bath", "TES", None),
        ]
        for segment, first, second in segments:
            def drop(point: dict[str, float]) -> float:
                left = point[first]
                right = TBATH if second is None else point[second]
                return left - right
            delta_minus, delta_center, delta_plus = drop(minus), drop(center), drop(plus)
            resistance = slope(delta_minus, delta_plus)
            budget_rows.append({
                "case": case_name,
                "segment": segment,
                "delta_T_0p95P_K": delta_minus,
                "delta_T_1p00P_K": delta_center,
                "delta_T_1p05P_K": delta_plus,
                "R_differential_K_per_W": resistance,
                "R_secant_K_per_W": delta_center / P0,
            })
    write_csv(OUT / "resistance_budget.csv", budget_rows)

    budget_by_case = defaultdict(dict)
    for row in budget_rows:
        budget_by_case[row["case"]][row["segment"]] = row
    comparison_rows: list[dict[str, Any]] = []
    historical = budget_by_case["historical"]
    refined = budget_by_case["refined_mortar"]
    total_delta = float(refined["total_TES_to_bath"]["R_differential_K_per_W"]) - float(historical["total_TES_to_bath"]["R_differential_K_per_W"])
    for segment in ("TES_to_Membrane_SiNx", "Membrane_SiNx_to_bath_side_SiO2_2_including_stack", "bath_side_SiO2_2_to_bath", "total_TES_to_bath"):
        delta = float(refined[segment]["R_differential_K_per_W"]) - float(historical[segment]["R_differential_K_per_W"])
        comparison_rows.append({
            "segment": segment,
            "historical_R_K_per_W": historical[segment]["R_differential_K_per_W"],
            "refined_mortar_R_K_per_W": refined[segment]["R_differential_K_per_W"],
            "refined_minus_historical_R_K_per_W": delta,
            "share_of_total_delta": delta / total_delta if total_delta else float("nan"),
        })
    write_csv(OUT / "resistance_comparison.csv", comparison_rows)

    summary = [
        "# Refined Phase24 thermal-resistance decomposition",
        "",
        "This read-only analysis reuses the existing historical and refined h=10 um mortar native MUMPS captures. No solver, production mesh, material, or TES law was changed.",
        "",
        "The endpoint budget is TES -> Membrane_SiNx -> bath-side SiO2_2 -> bath. The middle term includes the intervening SiO2/Si/SiNx/Membrane stack and any parallel network paths; it is not a unique single-interface resistance.",
        "",
    ]
    for row in comparison_rows:
        summary.append(f"- {row['segment']}: historical `{float(row['historical_R_K_per_W']):.9e} K/W`, refined mortar `{float(row['refined_mortar_R_K_per_W']):.9e} K/W`, delta `{float(row['refined_minus_historical_R_K_per_W']):.9e} K/W`, share `{float(row['share_of_total_delta']):.3f}`")
    summary += [
        "",
        "Per-body temperature slopes are in `body_temperature_slopes.csv`; the three-point temperature profiles are in `body_temperature_profile.csv`.",
        "The TES->Membrane and bath-side substrate->bath endpoint terms are separated directly. Any dominant residual in the middle term points to the membrane/substrate network realization, not TES-Stycast mortar formulation.",
        "",
    ]
    (OUT / "summary.md").write_text("\n".join(summary), encoding="utf-8")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
