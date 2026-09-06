"""Unit regressions for the v5 domain, edge, and surrogate contracts."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from PoST_Simulations.subScript import pulse_contamination_v5 as v5


def test_phase_randomized_surrogate_preserves_non_dc_periodogram():
    rng = np.random.default_rng(12)
    x = rng.normal(size=256)
    y = v5._phase_randomized(x, np.random.default_rng(13))
    np.testing.assert_allclose(
        np.abs(np.fft.rfft(x - x.mean()))[1:],
        np.abs(np.fft.rfft(y))[1:],
        rtol=1e-12,
        atol=1e-12,
    )


def test_edge_score_reports_negative_lag_and_available_overlap():
    kernel = np.zeros(16)
    kernel[3:8] = (1, 2, 3, 2, 1)
    record = np.zeros(32)
    record[:12] = kernel[4:]
    result = v5._edge_score(record, kernel, 1.0)
    assert result["lag_samples"] == -4
    assert result["overlap_samples"] == 12
    assert result["overlap_fraction"] == 12 / 16


def test_tail_age_bank_offsets_have_measured_source_support():
    source = np.arange(v5.TAIL_SOURCE_LENGTH + 100, dtype=float)
    rows = [{"raw_tail": source, "processed_tail": source}]
    for age in v5.TAIL_AGES_MS:
        offset = int(age * 1e-3 * v5.RATE_HZ)
        assert offset + v5.TAIL_LENGTH <= rows[0]["raw_tail"].size


def test_classifier_uses_one_preprocessing_pass_and_exclusive_class():
    template = np.zeros(8)
    template[2:5] = (1, 2, 1)
    library = {
        "CH0": {
            "templates": {
                "full_pulse": {"processed_values": template.tolist()},
                "tail_age_0ms": {"processed_values": (template / 2).tolist()},
            }
        }
    }
    null = {"family_rho": [0.0, 0.1, 0.2]}
    detector = {
        "library": library,
        "models": {"CH0": {"sigma": 1.0}},
        "short_null": {"CH0": null},
        "long_null": {"CH0": null},
    }
    hit = v5.classify_record_v5(np.zeros(32), "CH0", detector)
    assert hit["preprocessing_passes"] == 1
    assert hit["primary_class"] == "pulse_free_candidate"
    assert sum(hit["primary_class"] == name for name in ("full_pulse", "long_tail", "ambiguous", "pulse_free_candidate")) == 1
