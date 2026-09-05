"""Measure steady body flux closure and transient thermal storage."""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from scripts.analysis.evaluate_physical_parity import (
    body_boundary_flux_balance,
    heat_flux_consistency,
    mesh_data,
    result_field_times,
    result_values,
    tetra_volume,
)
from scripts.support.reconcile_project import reconcile_project


BODY_MATERIAL = {
    "abs": "Pb",
    "TES": "TES",
    "Stycast": "Stycast",
    "Membrane_SiNx": "Membrane",
    "Membrane_Si1": "Membrane",
    "SiO2_1": "SiO2",
    "Si_1": "Si",
    "SiNx": "SiNx",
    "Si_2": "Si",
    "SiO2_2": "SiO2",
}


def body_thermal_energy(mesh: Path, result: Path, project: Path, field_index: int) -> dict[str, float]:
    bodies, _, nodes, _, elements_by_body, _, _ = mesh_data(mesh)
    values = result_values(result, field_index)
    model = reconcile_project(json.loads(project.read_text(encoding="utf-8")))
    reverse = {value: key for key, value in bodies.items()}
    energies: dict[str, float] = {}
    for body_id, elements in elements_by_body.items():
        body_name = reverse.get(body_id, f"body_{body_id}")
        material_name = BODY_MATERIAL.get(body_name)
        if material_name is None:
            continue
        material = model["materials"][material_name]
        rho = float(material["rho"]["nominal"])
        cp = float(material["cp"]["nominal"])
        integral_temperature = 0.0
        volume = 0.0
        for conn in elements:
            element_volume = tetra_volume([nodes[node] for node in conn])
            volume += element_volume
            integral_temperature += (
                sum(values[node] for node in conn) / len(conn) * element_volume
            )
        energies[body_name] = rho * cp * integral_temperature
    return energies


def parse_level(value: str) -> tuple[str, Path, Path]:
    label, mesh, result = value.split("=", 2)
    return label, Path(mesh), Path(result)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--steady-level", action="append", required=True, help="label=mesh_dir=result_file")
    parser.add_argument("--transient-mesh", type=Path)
    parser.add_argument("--transient-result", type=Path)
    parser.add_argument("--transient-project", type=Path)
    parser.add_argument("--transient-series", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    steady = []
    for raw in args.steady_level:
        label, mesh, result = parse_level(raw)
        steady.append(
            {
                "label": label,
                "mesh": str(mesh.resolve()),
                "result": str(result.resolve()),
                "body_boundary_flux_balance": body_boundary_flux_balance(mesh, result, args.project),
            }
        )

    report: dict[str, object] = {
        "steady_control": {
            "source_term": "none",
            "levels": steady,
            "interpretation": "For the source-free steady control, body net outward reconstructed flux should approach zero. This is a raw-gradient diagnostic, not a relaxed interface gate.",
        }
    }
    if args.transient_mesh and args.transient_result and args.transient_project:
        times = result_field_times(args.transient_result)
        series_rows = []
        if args.transient_series and args.transient_series.exists():
            with args.transient_series.open(encoding="utf-8", newline="") as handle:
                series_rows = list(csv.DictReader(handle))
        fields = []
        for index, time_s in enumerate(times):
            body_flux = body_boundary_flux_balance(
                args.transient_mesh,
                args.transient_result,
                args.transient_project,
                index,
            )
            interface_flux = heat_flux_consistency(
                args.transient_mesh,
                args.transient_result,
                args.transient_project,
                field_index=index,
            )
            source_power = None
            if index < len(series_rows) and series_rows[index].get("tes_power_W"):
                source_power = float(series_rows[index]["tes_power_W"])
            fields.append(
                {
                    "field_index": index,
                    "time_s": time_s,
                    "internal_energy_J": body_thermal_energy(
                        args.transient_mesh,
                        args.transient_result,
                        args.transient_project,
                        index,
                    ),
                    "body_boundary_flux_balance": body_flux,
                    "interface_flux": {
                        name: {
                            "left_outward_flux_W": value["left_outward_flux_W"],
                            "right_outward_flux_W": value["right_outward_flux_W"],
                            "absolute_imbalance_W": value["absolute_imbalance_W"],
                            "normalized_imbalance": value["normalized_imbalance"],
                        }
                        for name, value in interface_flux.items()
                    },
                    "tes_source_power_W": source_power,
                }
            )
        storage = []
        for previous, current in zip(fields, fields[1:]):
            dt = (current["time_s"] or 0.0) - (previous["time_s"] or 0.0)
            delta = {
                body: current["internal_energy_J"].get(body, math.nan)
                - previous["internal_energy_J"].get(body, math.nan)
                for body in current["internal_energy_J"]
            }
            storage_power = {body: value / dt for body, value in delta.items()}
            tes_source = current.get("tes_source_power_W")
            tes_net_flux = current["body_boundary_flux_balance"].get("TES", {}).get("net_outward_flux_W")
            stycast_net_flux = current["body_boundary_flux_balance"].get("Stycast", {}).get("net_outward_flux_W")
            residuals = {}
            if tes_net_flux is not None and tes_source is not None:
                residuals["TES"] = storage_power.get("TES", math.nan) + tes_net_flux - tes_source
            if stycast_net_flux is not None:
                residuals["Stycast"] = storage_power.get("Stycast", math.nan) + stycast_net_flux
            storage.append({"from_field": previous["field_index"], "to_field": current["field_index"], "dt_s": dt, "delta_internal_energy_J": delta, "storage_power_W": storage_power, "tes_source_power_W": tes_source, "discrete_balance_residual_W": residuals})
        report["transient_storage"] = {
            "mesh": str(args.transient_mesh.resolve()),
            "result": str(args.transient_result.resolve()),
            "fields": fields,
            "intervals": storage,
            "note": "Internal energy is rho*cp*integral(T dV) from the saved nodal fields. Source/interface terms are reported separately and are not inferred from storage alone.",
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
