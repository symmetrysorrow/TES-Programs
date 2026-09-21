from pathlib import Path

from scripts.analysis import phase24_absolute_offset_diagnostic as diagnostic


def test_pair_report_separates_baseline_offset(tmp_path: Path) -> None:
    def make(label: str, offset: float) -> diagnostic.Series:
        times = [-1.0, 0.0, 0.1, 0.5, 0.9]
        values = {field: [offset, offset, offset - 1.0, offset - 2.0, offset - 3.0] for field in diagnostic.FIELDS}
        baseline = {field: offset for field in diagnostic.FIELDS}
        return diagnostic.Series(label, times, values, baseline)

    pair = diagnostic.pair_report(make("left", 10.0), make("right", 8.0))
    current = pair["metrics"]["tes_current_A"]
    assert current["baseline_offset"] == 2.0
    assert current["max_corrected_abs"] == 0.0


def test_load_reads_all_state_columns(tmp_path: Path) -> None:
    path = tmp_path / "series.csv"
    path.write_text(
        "time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W\n"
        "0.020019,1,0.0001,2,3\n"
        "0.020020,1.1,0.00009,2.1,2.9\n",
        encoding="utf-8",
    )
    series = diagnostic.load(path, "test")
    assert len(series.time_us) == 2
    assert set(series.baseline) == set(diagnostic.FIELDS)
