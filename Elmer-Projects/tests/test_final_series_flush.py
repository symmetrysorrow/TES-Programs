from __future__ import annotations

import csv
from pathlib import Path

from scripts.support.final_series_flush import flush_final_series_row


HEADER = "time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W\n"


def test_flushes_only_a_new_final_row(tmp_path: Path) -> None:
    series = tmp_path / "series.csv"
    iterations = tmp_path / "iterations.csv"
    series.write_text(HEADER + "1.0,0.1,0.2,0.3,0.4\n", encoding="utf-8")
    iterations.write_text(
        "time_s,time_step,nonlinear_iter,tes_temperature_K,previous_current_A,raw_current_A,"
        "tes_resistance_ohm,raw_power_W,residual_W,omega,omega_cap,relaxed_power_W\n"
        "1.5,2,2,0.11,0.2,0.21,0.31,0.41,0,0.5,0.5,0.42\n", encoding="utf-8")
    result = flush_final_series_row(series, iterations, True)
    assert result["appended"] is True
    assert len(series.read_text(encoding="utf-8").splitlines()) == 3
    assert flush_final_series_row(series, iterations, True)["appended"] is False
    assert flush_final_series_row(series, iterations, False)["reason"] == "solver_not_completed_normally"
