"""Flush the accepted final circuit row after a normal solver completion."""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Any


def _last(path: Path) -> dict[str, str] | None:
    if not path.is_file():
        return None
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return rows[-1] if rows else None


def flush_final_series_row(series: Path, iteration_series: Path, normal_completion: bool) -> dict[str, Any]:
    """Append one missing final row, returning an auditable decision."""
    result: dict[str, Any] = {"series": str(series), "iteration_series": str(iteration_series),
                              "normal_completion": normal_completion, "appended": False,
                              "reason": ""}
    if not normal_completion:
        result["reason"] = "solver_not_completed_normally"
        return result
    final = _last(iteration_series)
    existing = _last(series)
    if final is None or existing is None:
        result["reason"] = "missing_series_input"
        return result
    candidate_time = float(final["time_s"])
    existing_time = float(existing["time_s"])
    result.update({"candidate_time_s": candidate_time, "existing_time_s": existing_time})
    if candidate_time <= existing_time + 1.0e-18:
        result["reason"] = "already_flushed_or_not_newer"
        return result
    with series.open("a", encoding="utf-8", newline="") as handle:
        handle.write("{:.16E},{:.16E},{:.16E},{:.16E},{:.16E}\n".format(
            candidate_time, float(final["tes_temperature_K"]), float(final["raw_current_A"]),
            float(final["tes_resistance_ohm"]), float(final["relaxed_power_W"])))
    result["appended"] = True
    result["reason"] = "final_accepted_iteration_row"
    return result
