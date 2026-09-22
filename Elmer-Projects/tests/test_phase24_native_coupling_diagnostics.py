from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "artifacts" / "phase24_native_coupling_diagnostics"


def test_fixed_power_native_cases_have_identical_tes_membrane_constraint_topology() -> None:
    with (OUT / "native_mortar_reaction_integrals.csv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    values = [row for row in rows if row["interface"] == "TES_to_membrane"]
    assert {row["constraint_count"] for row in values} == {"179"}
    assert all(abs(float(row["reaction_balance_W"])) < 1.0e-22 for row in values)


def test_symmetric_delta_t_is_an_actual_fixed_power_native_run() -> None:
    payload = json.loads((OUT / "symmetric_delta_T_native.json").read_text(encoding="utf-8"))
    assert payload["delta_T_mK"] > 1.0
    assert payload["G_native_bath_W_per_K"] > 0.0
    assert "actual three-point fixed-power native MUMPS solves" in payload["method"]
