import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))

from Analyze_Results import (
    ReadPulse,
    _energy_column,
    bessel_filter,
    reference_peak_indices,
)


def _write_input(path):
    (path / "input.json").write_text(
        json.dumps({"rate": 1_000_000, "cutoff": 10_000, "SettlingTime": 20}),
        encoding="utf-8",
    )


def _pulse():
    pulse = np.zeros(120)
    pulse[10:21] = np.linspace(0.0, 1.0, 11)
    pulse[21:35] = np.linspace(1.0, 0.2, 14)
    return pulse


def test_readpulse_preserves_adaptive_height_and_fixed_energy_window(tmp_path):
    _write_input(tmp_path)
    values = ReadPulse(
        tmp_path,
        _pulse(),
        target="Pulse_ms",
        fixed_peak_index=15,
    )

    assert len(values) == 5
    assert values[1] == 20
    # The adaptive peak feature remains the peak-centered average.
    assert values[0] > values[4]
    # The fixed feature is the mean of samples 10..20, including the rising edge.
    np.testing.assert_allclose(values[4], np.mean(_pulse()[10:21]))


def test_readpulse_failure_and_invalid_fixed_index_have_five_nan_fields(tmp_path):
    _write_input(tmp_path)
    too_early = np.zeros(40)
    too_early[3] = 1
    assert len(ReadPulse(tmp_path, too_early)) == 5
    assert np.all(np.isnan(ReadPulse(tmp_path, too_early)))

    invalid = ReadPulse(tmp_path, _pulse(), fixed_peak_index=999)
    assert len(invalid) == 5
    assert np.all(np.isnan(invalid))


def test_reference_peak_indices_uses_requested_filter(tmp_path):
    _write_input(tmp_path)
    reference = np.exp(-0.5 * ((np.arange(80) - 30) / 5) ** 2)
    with h5py.File(tmp_path / "pulses.h5", "w") as file:
        file.attrs["format"] = "tes-pulses"
        file.create_dataset("event_id", data=np.asarray(["1"], dtype="S1"))
        file.create_dataset("ch0", data=reference[None, :])
        file.create_dataset("ch1", data=(reference * 0.5)[None, :])

    assert reference_peak_indices(tmp_path, 1, "Pulse_ms", "none") == (30, 30)
    filtered = reference_peak_indices(tmp_path, 1, "Pulse_noise", "bessel")
    expected = int(np.argmax(bessel_filter(reference, 1_000_000, 10_000)))
    assert filtered == (expected, expected)


def test_energy_column_prefers_new_csv_field_and_supports_legacy():
    legacy = pd.DataFrame({"height": [1.0, 2.0]})
    current = pd.DataFrame({"height": [1.0, 2.0], "energy_height": [3.0, 4.0]})
    pd.testing.assert_series_equal(_energy_column(legacy), legacy["height"])
    pd.testing.assert_series_equal(_energy_column(current), current["energy_height"])
