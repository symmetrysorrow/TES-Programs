from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import cross_day_spike_provenance_diagnostic as diag


def _line(frequency_hz, excess_db):
    return {
        "frequency_Hz": float(frequency_hz),
        "asd": 1.0,
        "excess_over_local_baseline_dB": float(excess_db),
        "offset_from_anchor_Hz": 0.0,
    }


def test_cluster_frequencies_groups_nearby_day_shift():
    clusters = diag.cluster_frequencies(
        [59640.0, 59720.0, 80000.0, 80100.0, 100000.0],
        tolerance_hz=150.0,
    )
    assert len(clusters) == 3
    assert clusters[0]["members_Hz"] == [59640.0, 59720.0]
    assert clusters[1]["members_Hz"] == [80000.0, 80100.0]
    assert clusters[2]["members_Hz"] == [100000.0]


def test_localize_line_finds_shifted_peak_inside_window():
    frequency = np.arange(90000.0, 91005.0, 5.0)
    asd = np.ones_like(frequency)
    shifted = int(np.argmin(np.abs(frequency - 90425.0)))
    asd[shifted] = 10.0
    result = diag.localize_line(
        frequency,
        asd,
        anchor_hz=90500.0,
        search_half_width_hz=500.0,
        baseline_width_hz=200.0,
    )
    assert result is not None
    assert result["frequency_Hz"] == 90425.0
    assert result["excess_over_local_baseline_dB"] > 19.0


def test_classify_day_variation_when_stored_matches_reference():
    repeat = _line(119290.0, 15.0)
    reference = _line(119300.0, 2.0)
    stored = _line(119300.0, 2.4)
    result = diag.classify_line(
        repeat,
        reference,
        stored,
        material_db=3.0,
        stored_match_db=1.5,
        line_present_db=3.0,
    )
    assert result == "day_or_acquisition_variation_supported"


def test_classify_stored_generation_difference():
    repeat = _line(100000.0, 17.0)
    reference = _line(100000.0, 17.2)
    stored = _line(100000.0, 3.0)
    result = diag.classify_line(
        repeat,
        reference,
        stored,
        material_db=3.0,
        stored_match_db=1.5,
        line_present_db=3.0,
    )
    assert result == "stored_generation_difference_supported"


def test_classify_visual_concealment_when_line_persists():
    repeat = _line(100000.0, 17.0)
    reference = _line(100000.0, 17.2)
    stored = _line(100000.0, 16.8)
    result = diag.classify_line(
        repeat,
        reference,
        stored,
        material_db=3.0,
        stored_match_db=1.5,
        line_present_db=3.0,
    )
    assert result == "line_persists_in_stored_target_visual_concealment_candidate"


def test_classify_mixed_when_day_and_stored_both_change():
    repeat = _line(96000.0, 12.0)
    reference = _line(96000.0, 4.0)
    stored = _line(96000.0, 10.0)
    result = diag.classify_line(
        repeat,
        reference,
        stored,
        material_db=3.0,
        stored_match_db=1.5,
        line_present_db=3.0,
    )
    assert result == "mixed_day_and_stored_difference"


def test_cleaned_result_converts_path_and_numpy_types(tmp_path):
    payload = {
        "comparison": {
            "repeat": {
                "experiment_path": tmp_path / "run",
                "count": np.int64(3),
                "values": np.asarray([1.0, 2.0]),
            }
        },
        "_plot": {"large": np.asarray([1.0])},
    }
    cleaned = diag.cleaned_result(payload)
    assert cleaned["comparison"]["repeat"]["experiment_path"] == str(
        tmp_path / "run"
    )
    assert cleaned["comparison"]["repeat"]["count"] == 3
    assert cleaned["comparison"]["repeat"]["values"] == [1.0, 2.0]
    assert "_plot" not in cleaned
