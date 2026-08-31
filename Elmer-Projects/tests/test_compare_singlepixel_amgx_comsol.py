from __future__ import annotations

import unittest
import tempfile
from pathlib import Path

import numpy as np

from scripts.analysis.compare_singlepixel_amgx_comsol import (
    Series,
    compare,
    crossing,
    formatted_time,
    normalized_series,
    read_elmer,
)


class SinglePixelAmgxComsolTests(unittest.TestCase):
    def test_normalization_uses_pre_pulse_mean(self) -> None:
        time_s = np.asarray([19.5e-3, 20.0e-3, 20.02e-3, 20.03e-3])
        series = normalized_series(time_s, np.asarray([10.0, 12.0, 9.0, 8.0]))
        self.assertEqual(series.baseline_uA, 11.0)
        np.testing.assert_allclose(series.drop_uA, [1.0, -1.0, 2.0, 3.0])

    def test_crossing_interpolates_linearly(self) -> None:
        value = crossing(
            np.asarray([-1.0, 0.0, 10.0]),
            np.asarray([0.0, 0.0, 2.0]),
            1.0,
        )
        self.assertEqual(value, 5.0)

    def test_missing_crossing_formats_as_na(self) -> None:
        self.assertEqual(formatted_time(None), "n/a")

    def test_elmer_reader_sorts_and_deduplicates_equal_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "series.csv"
            path.write_text(
                "time_s,tes_current_A\n"
                "0.020020,1.0e-4\n"
                "0.019500,2.0e-4\n"
                "0.020020,1.0e-4\n",
                encoding="utf-8",
            )
            series = read_elmer(path)
            np.testing.assert_allclose(series.time_us, [-520.0, 0.0])
            np.testing.assert_allclose(series.current_uA, [200.0, 100.0])

    def test_identical_series_have_zero_waveform_error(self) -> None:
        series = Series(
            np.asarray([-1.0, 0.0, 50.0, 100.0]),
            np.asarray([10.0, 10.0, 8.0, 6.0]),
            10.0,
        )
        metrics, aligned = compare(series, series, 100.0)
        self.assertEqual(metrics["max_abs_difference_uA"], 0.0)
        self.assertEqual(metrics["rmse_uA"], 0.0)
        np.testing.assert_allclose(aligned[:, 3], 0.0)

        with self.assertRaisesRegex(ValueError, "ends at 100.000 us"):
            compare(series, series, 120.0)


if __name__ == "__main__":
    unittest.main()
