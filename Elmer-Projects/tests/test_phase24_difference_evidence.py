from pathlib import Path

from scripts.analysis.phase24_difference_evidence import load_series, mean_baseline


def test_load_series_converts_units_and_collapses_duplicate_times(tmp_path: Path) -> None:
    path = tmp_path / "series.csv"
    path.write_text(
        "time_s,tes_current_A\n"
        "1.0e-3,2.0e-4\n"
        "1.0e-3,2.5e-4\n"
        "1.1e-3,1.5e-4\n",
        encoding="utf-8",
    )

    times, currents = load_series(path)

    assert times == [1000.0, 1100.0]
    assert currents == [250.0, 150.0]


def test_mean_baseline_uses_pre_event_window() -> None:
    times = [-2.0, -1.0, 0.0, 1.0]
    currents = [10.0, 12.0, 99.0, 8.0]

    assert mean_baseline(times, currents, -2.0, 0.0) == 11.0
