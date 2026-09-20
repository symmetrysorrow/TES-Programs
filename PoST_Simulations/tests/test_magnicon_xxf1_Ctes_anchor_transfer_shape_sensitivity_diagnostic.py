from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_Ctes_anchor_transfer_shape_sensitivity_diagnostic as diag


def test_geometric_mean_transfer():
    a = np.array([1.0, 4.0, 9.0])
    b = np.array([1.0, 1.0, 4.0])
    np.testing.assert_allclose(
        diag.geometric_mean_transfer(a, b),
        np.array([1.0, 2.0, 6.0]),
    )


def test_geometric_mean_transfer_rejects_nonpositive():
    with pytest.raises(ValueError):
        diag.geometric_mean_transfer(
            np.array([1.0, 0.0]),
            np.array([1.0, 1.0]),
        )


def test_c2_is_censored():
    row = {
        "c2_free": {
            "readout_boundary_hits": {
                "c2": {
                    "at_lower": False,
                    "at_upper": True,
                }
            }
        }
    }
    assert diag.c2_is_censored(row)


def test_c2_is_not_censored():
    row = {
        "c2_free": {
            "readout_boundary_hits": {
                "c2": {
                    "at_lower": False,
                    "at_upper": False,
                }
            }
        }
    }
    assert not diag.c2_is_censored(row)


def test_spread_summary_identical_transfers():
    rows = diag.spread_summary(
        {
            "a": np.array([1.0, 2.0, 3.0]),
            "b": np.array([1.0, 2.0, 3.0]),
        }
    )
    assert rows["rms_peak_to_peak_spread_dB"] == pytest.approx(0.0)
    assert rows["max_peak_to_peak_spread_dB"] == pytest.approx(0.0)


def test_difference_summary_identity():
    frequency = np.array([1_000.0, 2_000.0, 10_000.0])
    values = np.array([1.0, 2.0, 4.0])
    row = diag.difference_summary(
        frequency,
        values,
        values,
        [{"name": "all", "min": 1_000.0, "max": 10_000.0}],
        [1_000.0, 10_000.0],
    )
    assert row["full_1_200k"]["rms_difference_dB"] == pytest.approx(0.0)
    assert row["full_1_200k"]["max_abs_difference_dB"] == pytest.approx(0.0)


def test_classify_pairwise_stable():
    pairwise = {
        "a_vs_b": {
            "full_1_200k": {
                "rms_difference_dB": 0.3,
                "max_abs_difference_dB": 0.8,
            }
        }
    }
    screen = {
        "max_pairwise_full_band_rms_difference_dB": 0.5,
        "max_pairwise_abs_difference_dB": 1.5,
    }
    result = diag.classify_pairwise(pairwise, screen)
    assert result["passes_anchor_sensitivity_screen"]
    assert result["classification"] == (
        "total_readout_transfer_stable_across_Ctes_anchors"
    )


def test_classify_pairwise_sensitive():
    pairwise = {
        "a_vs_b": {
            "full_1_200k": {
                "rms_difference_dB": 1.2,
                "max_abs_difference_dB": 2.0,
            }
        }
    }
    screen = {
        "max_pairwise_full_band_rms_difference_dB": 0.5,
        "max_pairwise_abs_difference_dB": 1.5,
    }
    result = diag.classify_pairwise(pairwise, screen)
    assert not result["passes_anchor_sensitivity_screen"]
    assert result["classification"] == (
        "required_total_readout_transfer_is_Ctes_anchor_sensitive"
    )


def test_classify_pairwise_empty():
    result = diag.classify_pairwise(
        {},
        {
            "max_pairwise_full_band_rms_difference_dB": 0.5,
            "max_pairwise_abs_difference_dB": 1.5,
        },
    )
    assert result["classification"] == "insufficient_uncensored_anchor_pairs"
    assert result["passes_anchor_sensitivity_screen"] is None
