from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts.analysis.compare_mpi_series import compare


def _series(path: Path, currents: list[float]) -> None:
    times = [0.020018, 0.020019, 0.020020, 0.020021, 0.020030, 0.020040]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["time_s", "tes_current_A"])
        writer.writeheader()
        writer.writerows(
            {"time_s": time, "tes_current_A": current}
            for time, current in zip(times, currents, strict=True)
        )


def test_identical_series_pass_all_metrics(tmp_path: Path) -> None:
    reference = tmp_path / "reference.csv"
    candidate = tmp_path / "candidate.csv"
    currents = [10.0, 10.0, 10.0, 9.0, 8.0, 9.0]
    _series(reference, currents)
    _series(candidate, currents)
    result = compare(reference, candidate)
    assert result["overall_passed"]
    assert result["matched_samples"] == 6
    assert all(metric["value"] == 0.0 for metric in result["metrics"].values())
    assert result["inputs"]["reference"]["sha256"]


def test_shape_difference_fails_pointwise_metrics(tmp_path: Path) -> None:
    reference = tmp_path / "reference.csv"
    candidate = tmp_path / "candidate.csv"
    _series(reference, [10.0, 10.0, 10.0, 9.0, 8.0, 9.0])
    _series(candidate, [10.0, 10.0, 10.0, 9.2, 8.0, 9.0])
    result = compare(reference, candidate)
    assert not result["overall_passed"]
    assert not result["metrics"]["raw_max_current"]["passed"]
    assert not result["metrics"]["baseline_corrected_max_current"]["passed"]


def test_missing_internal_timestamp_is_rejected(tmp_path: Path) -> None:
    reference = tmp_path / "reference.csv"
    candidate = tmp_path / "candidate.csv"
    _series(reference, [10.0, 10.0, 10.0, 9.0, 8.0, 9.0])
    _series(candidate, [10.0, 10.0, 10.0, 9.0, 8.0, 9.0])
    rows = candidate.read_text(encoding="utf-8").splitlines()
    candidate.write_text("\n".join([*rows[:4], *rows[5:]]) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing reference timestamps"):
        compare(reference, candidate)
