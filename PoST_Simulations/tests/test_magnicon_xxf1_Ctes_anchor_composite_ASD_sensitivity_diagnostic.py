from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_Ctes_anchor_composite_ASD_sensitivity_diagnostic as diag


def test_normalized_shape_reference_is_unity():
    frequency = np.array([1_000.0, 2_000.0, 4_000.0])
    values = np.array([2.0, 4.0, 8.0])
    result = diag.normalized_shape(
        frequency,
        values,
        reference_hz=1_000.0,
    )
    np.testing.assert_allclose(result, np.array([1.0, 2.0, 4.0]))


def test_raw_c2_is_censored():
    row = {
        "readout_boundary_hits": {
            "c2": {
                "at_lower": False,
                "at_upper": True,
            }
        }
    }
    assert diag.raw_c2_is_censored(row)


def test_raw_c2_not_censored_without_hit():
    row = {
        "readout_boundary_hits": {
            "c2": {
                "at_lower": False,
                "at_upper": False,
            }
        }
    }
    assert not diag.raw_c2_is_censored(row)


def test_sensitivity_ratio():
    assert diag.sensitivity_ratio(2.0, 8.0) == pytest.approx(0.25)
    assert diag.sensitivity_ratio(2.0, 0.0) is None


def test_classify_decomposition_nonidentifiable():
    observable = {"passes_anchor_sensitivity_screen": True}
    readout = {"passes_anchor_sensitivity_screen": False}
    assert diag.classify_decomposition(
        observable_classification=observable,
        readout_classification=readout,
    ) == (
        "Ctes_readout_decomposition_nonidentifiable_"
        "but_observable_ASD_stable"
    )


def test_classify_decomposition_observable_sensitive():
    observable = {"passes_anchor_sensitivity_screen": False}
    readout = {"passes_anchor_sensitivity_screen": False}
    assert diag.classify_decomposition(
        observable_classification=observable,
        readout_classification=readout,
    ) == "observable_ASD_is_Ctes_anchor_sensitive"


def test_classify_decomposition_both_stable():
    observable = {"passes_anchor_sensitivity_screen": True}
    readout = {"passes_anchor_sensitivity_screen": True}
    assert diag.classify_decomposition(
        observable_classification=observable,
        readout_classification=readout,
    ) == "observable_and_readout_shapes_stable_across_Ctes_anchors"


def test_pairwise_summaries_identity_and_difference():
    frequency = np.array([1_000.0, 2_000.0, 10_000.0])
    shapes = {
        "a": np.array([1.0, 1.0, 1.0]),
        "b": np.array([1.0, 2.0, 2.0]),
        "c": np.array([1.0, 1.0, 1.0]),
    }
    rows = diag.pairwise_summaries(
        frequency=frequency,
        shapes=shapes,
        included_keys=["a", "b", "c"],
        bands=[
            {
                "name": "all",
                "min": 1_000.0,
                "max": 10_000.0,
            }
        ],
        anchors=[1_000.0, 10_000.0],
    )
    assert set(rows) == {"a_vs_b", "a_vs_c", "b_vs_c"}
    assert rows["a_vs_c"]["full_1_200k"][
        "rms_difference_dB"
    ] == pytest.approx(0.0)
    assert rows["a_vs_b"]["full_1_200k"][
        "rms_difference_dB"
    ] > 0.0


def test_geometric_mean_shapes():
    reference = {
        "full_fitted_ASD": np.array([1.0, 4.0]),
        "readout_transfer": np.array([1.0, 9.0]),
    }
    repeat = {
        "full_fitted_ASD": np.array([1.0, 1.0]),
        "readout_transfer": np.array([1.0, 4.0]),
    }
    result = diag._geometric_mean_shapes(reference, repeat)
    np.testing.assert_allclose(
        result["full_fitted_ASD"],
        np.array([1.0, 2.0]),
    )
    np.testing.assert_allclose(
        result["readout_transfer"],
        np.array([1.0, 6.0]),
    )
