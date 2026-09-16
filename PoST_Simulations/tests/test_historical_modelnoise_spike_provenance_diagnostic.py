from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import historical_modelnoise_spike_provenance_diagnostic as diag


def test_discover_historical_settings_prefers_setting_json(tmp_path):
    payload = {
        "Config": {
            "threshold": 5e-5,
            "presamples": 5000,
            "eta_uA_per_V": 2.5,
        },
        "main": {"cutoff": 10000},
    }
    (tmp_path / "setting.json").write_text(
        json.dumps(payload),
        encoding="utf-8",
    )
    result = diag.discover_historical_settings(
        tmp_path,
        acquisition_cutoff=12000,
    )
    assert result["threshold"] == 5e-5
    assert result["presamples"] == 5000
    assert result["cutoff_Hz"] == 10000
    assert result["eta_uA_per_V"] == 2.5
    assert result["historical_peak_mask_available"]
    assert result["historical_bessel_available"]


def test_discover_historical_settings_uses_cli_override(tmp_path):
    result = diag.discover_historical_settings(
        tmp_path,
        threshold_override=1.2e-4,
        presamples_override=6000,
        cutoff_override=9000,
    )
    assert result["threshold"] == 1.2e-4
    assert result["presamples"] == 6000
    assert result["cutoff_Hz"] == 9000
    assert result["sources"]["threshold"] == "CLI --peak-threshold"


def test_historical_peak_accept_matches_noise_main_semantics():
    values = np.zeros(10000)
    values[7000] = 4e-5
    accepted, baseline, peak = diag.historical_peak_accept(
        values,
        presamples=5000,
        threshold=5e-5,
    )
    assert accepted
    assert baseline == 0.0
    assert peak == 4e-5

    values[7000] = 6e-5
    accepted, _baseline, peak = diag.historical_peak_accept(
        values,
        presamples=5000,
        threshold=5e-5,
    )
    assert not accepted
    assert peak == 6e-5


def test_mask_overlap_counts_and_jaccard():
    result = diag.mask_overlap({"a", "b", "c"}, {"b", "c", "d"})
    assert result["intersection_count"] == 2
    assert result["a_only_count"] == 1
    assert result["b_only_count"] == 1
    assert result["union_count"] == 4
    assert result["jaccard"] == 0.5


def test_local_line_excess_is_scale_invariant():
    frequency = np.arange(0.0, 20000.0, 5.0)
    asd = np.ones_like(frequency)
    index = int(np.argmin(np.abs(frequency - 10000.0)))
    asd[index] = 10.0

    first = diag.local_line_metric(
        frequency,
        asd,
        10000.0,
    )
    second = diag.local_line_metric(
        frequency,
        asd * 123.0,
        10000.0,
    )
    assert np.isclose(
        first["excess_over_local_baseline_dB"],
        20.0,
        atol=1e-12,
    )
    assert np.isclose(
        first["excess_over_local_baseline_dB"],
        second["excess_over_local_baseline_dB"],
        atol=1e-12,
    )


def test_theoretical_bessel_filtfilt_suppresses_100k_more_than_10k():
    values = diag.theoretical_filtfilt_transfer_db(
        500000.0,
        10000.0,
        [10000.0, 100000.0],
    )
    assert values[1] < values[0] - 50.0
    assert values[1] < -70.0


def test_discover_stored_modelnoise_uses_setting_output(tmp_path):
    output = tmp_path / "CH0_noise" / "output" / "run01"
    output.mkdir(parents=True)
    modelnoise = output / "modelnoise.txt"
    modelnoise.write_text("1\n2\n", encoding="utf-8")
    path, meta = diag.discover_stored_modelnoise_path(
        tmp_path,
        {"output_name": "run01"},
        None,
    )
    assert path == modelnoise
    assert meta["source"] == "setting.json Config.output"
