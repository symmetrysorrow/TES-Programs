from __future__ import annotations

import csv
from pathlib import Path

from scripts.analysis.summarize_iteration_series import summarize


def test_summarize_counts_iterations_per_timestep(tmp_path: Path) -> None:
    path = tmp_path / "iterations.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["time_s", "time_step", "nonlinear_iter"]
        )
        writer.writeheader()
        writer.writerows(
            [
                {"time_s": 1, "time_step": 1, "nonlinear_iter": 1},
                {"time_s": 1, "time_step": 1, "nonlinear_iter": 2},
                {"time_s": 2, "time_step": 2, "nonlinear_iter": 1},
                {"time_s": 2, "time_step": 2, "nonlinear_iter": 2},
                {"time_s": 2, "time_step": 2, "nonlinear_iter": 3},
            ]
        )
    result = summarize(path)
    assert result["rows"] == 5
    assert result["timesteps"] == 2
    assert result["iterations_per_step_histogram"] == {"2": 1, "3": 1}
    assert result["iterations_by_timestep"] == {"1": 2, "2": 3}
