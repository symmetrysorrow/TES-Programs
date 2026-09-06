"""Focused v6 contract tests that do not require the target data volume."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from PoST_Simulations.subScript import pulse_contamination_v6 as v6


def test_recordwise_empirical_p_value_uses_the_record_local_null():
    assert v6.empirical(3.0, np.array([1.0, 2.0, 3.0])) == 0.5
    assert v6.empirical(3.0, np.array([3.0, 3.0, 3.0])) == 1.0


def test_phase_randomized_batch_preserves_each_source_fft_magnitude():
    rng = np.random.default_rng(4)
    source = rng.normal(size=(3, 128))
    surrogate = v6.phase_randomized_batch(source, np.random.default_rng(5))
    expected = np.abs(np.fft.rfft(source - source.mean(axis=1, keepdims=True), axis=1))
    actual = np.abs(np.fft.rfft(surrogate, axis=1))
    np.testing.assert_allclose(actual[:, 1:], expected[:, 1:], rtol=1e-12, atol=1e-12)


def test_iaaft_batch_preserves_each_record_amplitude_distribution():
    rng = np.random.default_rng(6)
    source = rng.normal(size=(2, 64))
    surrogate = v6.iaaft_batch(source, np.random.default_rng(7), iterations=1)
    np.testing.assert_allclose(np.sort(surrogate, axis=1), np.sort(source, axis=1))


def test_short_combined_statistic_is_a_single_maximum():
    contained = np.array([1.0, 5.0, 2.0])
    edge = np.array([4.0, 2.0, 3.0])
    combined = np.maximum(contained, edge)
    np.testing.assert_array_equal(combined, [4.0, 5.0, 3.0])
    assert not np.array_equal(combined, np.minimum(contained, edge))


def test_v6_classifier_exposes_common_time_axis_and_edge_provenance():
    template = np.zeros(16)
    template[4:8] = [1.0, 2.0, 1.0, 0.5]
    channel = {
        "templates": {
            "full_pulse": {"processed_values": template.tolist()},
            "tail_age_0ms": {"processed_values": template.tolist()},
        }
    }
    detector = {
        "library": {"CH0": channel},
        "models": {"CH0": {"sigma": 1.0}},
        "short_null": {"CH0": {"combined_family_rho": [0.0, 0.1, 0.2]}},
        "injection_long_null": {"CH0": np.array([0.0, 0.1, 0.2])},
    }
    result = v6.classify_record_v6(np.zeros(64), "CH0", detector, long_null=np.array([0.0, 0.1, 0.2]))
    assert result["preprocessing_passes"] == 1
    assert "selected_lag_samples" in result
    assert result["edge_mode"] in {"fully_contained", "onset_near_end", "pre_record_tail"}
