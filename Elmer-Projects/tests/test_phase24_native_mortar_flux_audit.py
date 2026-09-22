from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "artifacts" / "phase24_native_mortar_flux_audit"


def rows(name: str) -> list[dict[str, str]]:
    with (AUDIT / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_native_reaction_is_partitioned_by_interface() -> None:
    values = rows("native_mortar_reaction_integrals.csv")
    historical = {row["interface"]: row for row in values if row["case"] == "historical"}
    phase24 = {row["interface"]: row for row in values if row["case"] == "Phase24"}
    assert int(historical["TES_to_membrane"]["constraint_count"]) == 2218
    assert abs(float(historical["TES_to_membrane"]["weighted_reaction_W"]) - 3.203008493e-10) < 1.0e-18
    assert int(phase24["TES_to_membrane"]["constraint_count"]) == 0


def test_native_partition_closes_against_bath_with_explicit_residual() -> None:
    values = {row["case"]: row for row in rows("native_heatflow_partition.csv")}
    assert abs(float(values["historical"]["global_balance_relative"])) < 2.0e-6
    assert abs(float(values["Phase24"]["global_balance_relative"])) < 5.0e-4
    assert float(values["Phase24"]["other_or_conformal_path_W"]) > 4.0e-10


def test_matched_point_outputs_are_marked_as_secant_estimates() -> None:
    values = rows("matched_temperature_comparison.csv")
    assert len(values) == 2
    assert all(row["diagnostic_status"] == "estimated_from_native_operating_point" for row in values)
    assert all("no ±delta-T solve" in row["method"] for row in values)
