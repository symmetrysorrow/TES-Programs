from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import preanalysis_readout_biquad_diagnostic as diag  # noqa: E402


def test_biquad_is_unity_when_pole_and_zero_match():
    frequency = np.geomspace(100.0, 500_000.0, 301)
    response = diag.biquad_magnitude(
        frequency,
        pole_hz=12_000.0,
        pole_q=2.5,
        zero_hz=12_000.0,
        zero_q=2.5,
    )
    np.testing.assert_allclose(response, np.ones_like(response), rtol=1e-13)


def test_biquad_requires_complex_conjugate_q():
    with pytest.raises(ValueError):
        diag.biquad_magnitude(
            np.array([1_000.0]),
            pole_hz=10_000.0,
            pole_q=0.5,
            zero_hz=20_000.0,
            zero_q=1.0,
        )


def test_pre_analysis_model_preserves_post_filter_white(monkeypatch):
    monkeypatch.setattr(
        diag.opt,
        "hardware_filter_magnitude",
        lambda frequency, **kwargs: np.ones_like(
            np.asarray(frequency, dtype=float)
        ),
    )
    context = {
        "frequency_Hz": np.array([1_000.0, 2_000.0]),
        "alias_frequency_Hz": np.array([99_000.0, 98_000.0]),
        "main_intrinsic_asd": np.array([3.0, 4.0]),
        "alias_intrinsic_asd": np.array([4.0, 3.0]),
        "same_bin": np.array([False, False]),
        "post_filter_white_asd_A_rtHz": 5.0,
    }
    result = diag.pre_analysis_model(context, [])
    absolute = np.sqrt(
        context["main_intrinsic_asd"] ** 2
        + context["alias_intrinsic_asd"] ** 2
        + 5.0**2
    )
    expected = absolute / absolute[0]
    np.testing.assert_allclose(result, expected)


def test_biquad_is_applied_before_alias_fold(monkeypatch):
    monkeypatch.setattr(
        diag.opt,
        "hardware_filter_magnitude",
        lambda frequency, **kwargs: np.ones_like(
            np.asarray(frequency, dtype=float)
        ),
    )
    context = {
        "frequency_Hz": np.array([1_000.0, 20_000.0]),
        "alias_frequency_Hz": np.array([99_000.0, 80_000.0]),
        "main_intrinsic_asd": np.array([2.0, 2.0]),
        "alias_intrinsic_asd": np.array([1.0, 1.0]),
        "same_bin": np.array([False, False]),
        "post_filter_white_asd_A_rtHz": 0.0,
    }
    section = {
        "pole_Hz": 12_000.0,
        "pole_Q": 1.5,
        "zero_Hz": 80_000.0,
        "zero_Q": 1.2,
    }
    result = diag.pre_analysis_model(context, [section])

    main_h = diag.biquad_magnitude(
        context["frequency_Hz"], 12_000.0, 1.5, 80_000.0, 1.2
    )
    alias_h = diag.biquad_magnitude(
        context["alias_frequency_Hz"], 12_000.0, 1.5, 80_000.0, 1.2
    )
    absolute = np.sqrt(
        (context["main_intrinsic_asd"] * main_h) ** 2
        + (context["alias_intrinsic_asd"] * alias_h) ** 2
    )
    expected = absolute / absolute[0]
    np.testing.assert_allclose(result, expected)


def test_decode_returns_physical_positive_sections_sorted_by_pole():
    vector = np.log10(
        np.array(
            [
                80_000.0,
                2.0,
                20_000.0,
                1.0,
                10_000.0,
                3.0,
                90_000.0,
                4.0,
            ]
        )
    )
    sections = diag.decode(vector, 2)
    assert len(sections) == 2
    assert sections[0]["pole_Hz"] < sections[1]["pole_Hz"]
    for section in sections:
        assert section["pole_Hz"] > 0.0
        assert section["zero_Hz"] > 0.0
        assert section["pole_Q"] > 0.5
        assert section["zero_Q"] > 0.5
