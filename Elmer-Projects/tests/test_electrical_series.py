from __future__ import annotations

import csv
from pathlib import Path

from run import ensure_electrical_series


def test_iteration_csv_is_normalized_to_one_canonical_row_per_timestep(tmp_path: Path) -> None:
    iteration = tmp_path / "iterations.csv"
    iteration.write_text(
        "time_s,time_step,nonlinear_iter,tes_temperature_K,previous_current_A,"
        "raw_current_A,tes_resistance_ohm,raw_power_W,residual_W,omega,omega_cap,relaxed_power_W\n"
        "0.1,1,1,0.16,1e-4,2e-4,0.01,4e-10,1e-11,0.5,0.5,4e-10\n"
        "0.1,1,2,0.17,1e-4,1.9e-4,0.011,3.9e-10,1e-12,0.5,0.5,3.99e-10\n"
        "0.2,2,1,0.18,1.9e-4,1.8e-4,0.012,3.8e-10,1e-12,0.5,0.5,3.8e-10\n",
        encoding="utf-8",
    )
    out = tmp_path / "results"
    out.mkdir()
    source = ensure_electrical_series(
        tmp_path,
        out,
        {"state_file": "missing.state"},
        {"parameters": {"I_bias": 7.15e-4, "R_sh": 3.9e-3}},
        "series.csv",
        iteration,
    )
    assert source == "normalized_from_iteration_or_state"
    with (out / "series.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[0]["nonlinear_iter"] == "2"
    assert rows[0]["tes_current_A"] == "0.00019"
    assert rows[0]["bias_current_A"] == "0.000715"
    assert rows[1]["time_step"] == "2"


def test_state_file_is_a_steady_series_fallback(tmp_path: Path) -> None:
    state = tmp_path / "steady.state"
    state.write_text("0.17 0.00019 0.011 3.96e-10 0.00019\n", encoding="utf-8")
    out = tmp_path / "results"
    out.mkdir()
    source = ensure_electrical_series(
        tmp_path,
        out,
        {"state_file": "steady.state"},
        {"parameters": {"I_bias": 7.15e-4, "R_sh": 3.9e-3}},
        "series.csv",
        None,
    )
    assert source == "normalized_from_iteration_or_state"
    text = (out / "series.csv").read_text(encoding="utf-8")
    assert "tes_temperature_K" in text
    assert "0.17" in text
