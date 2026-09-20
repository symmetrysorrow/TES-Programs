from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_c2_on_off_fullband_plot as diag


def test_residual_db_identity():
    values = np.array([1.0, 2.0, 4.0])
    np.testing.assert_allclose(
        diag.residual_db(values, values),
        np.zeros(3),
    )


def test_residual_db_known_ratio():
    target = np.ones(2)
    model = np.array([1.0, 10.0])
    result = diag.residual_db(model, target)
    np.testing.assert_allclose(result, np.array([0.0, 20.0]))


def test_residual_db_rejects_nonpositive():
    with pytest.raises(ValueError):
        diag.residual_db(
            np.array([1.0, 0.0]),
            np.array([1.0, 1.0]),
        )


def test_multiplier_slug():
    assert diag.multiplier_slug(0.25) == "0p25x"
    assert diag.multiplier_slug(1.0) == "1x"


def test_figure_path(tmp_path):
    result = diag.figure_path(tmp_path, 0.25)
    assert result == (
        tmp_path
        / "magnicon_xxf1_c2_on_off_fullband_Ctes_0p25x.png"
    )


def test_display_mask():
    frequency = np.array([500.0, 1000.0, 5000.0, 200000.0, 250000.0])
    mask = diag.display_mask(
        frequency,
        {"min": 1000.0, "max": 200000.0},
    )
    np.testing.assert_array_equal(
        mask,
        np.array([False, True, True, True, False]),
    )


def test_display_mask_rejects_invalid_range():
    with pytest.raises(ValueError):
        diag.display_mask(
            np.array([1000.0, 2000.0]),
            {"min": 2000.0, "max": 1000.0},
        )


def test_day_plot_payload():
    problem = {
        "target": np.array([1.0, 2.0]),
        "raw_target": np.array([1.1, 2.2]),
    }
    row = {
        "c2_zero": {
            "_model_full": np.array([1.0, 4.0]),
        },
        "c2_free": {
            "_model_full": np.array([1.0, 2.0]),
        },
    }
    result = diag._day_plot_payload(problem, row)
    np.testing.assert_allclose(
        result["c2_zero_residual_dB"],
        np.array([0.0, 20.0 * np.log10(2.0)]),
    )
    np.testing.assert_allclose(
        result["c2_free_residual_dB"],
        np.zeros(2),
    )
