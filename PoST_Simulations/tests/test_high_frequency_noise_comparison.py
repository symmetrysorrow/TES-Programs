"""Lightweight contract tests for the frozen high-frequency comparison."""
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))
from PoST_Simulations.subScript.high_frequency_noise_comparison import log_interp, intervals


def test_log_interp_is_log_log():
    x = np.array([1.0, 10.0, 100.0]); y = x ** 2
    np.testing.assert_allclose(log_interp(x, y, np.array([3.16227766])), [10.0], rtol=1e-7)


def test_intervals_are_deterministic_and_contiguous():
    freq = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert intervals(freq, np.array([False, True, True, False, True])) == [[2.0, 3.0], [5.0, 5.0]]


def test_primary_policy_is_frozen_shape_only():
    source = Path(__file__).parents[1] / "subScript" / "high_frequency_noise_comparison.py"
    text = source.read_text(encoding="utf-8")
    assert "pulse_consistent_scenarios" in text
    assert "noise_residual_fit" in text
    assert "parameter_optimization" in text
    assert "pulse_free_candidate mask not used" in text


def test_primary_band_excludes_low_frequency_from_score():
    source = Path(__file__).parents[1] / "subScript" / "high_frequency_noise_comparison.py"
    text = source.read_text(encoding="utf-8")
    assert "1000.0, 10000.0" in text
    assert "eval_freq = np.logspace" in text


def test_display_band_covers_full_positive_frequency_range():
    source = Path(__file__).parents[1] / "subScript" / "high_frequency_noise_comparison.py"
    text = source.read_text(encoding="utf-8")
    assert "DISPLAY_MIN_HZ = RATE_HZ / SAMPLES" in text
    assert "DISPLAY_MAX_HZ = RATE_HZ / 2.0" in text
    assert "plot_freq = np.logspace(np.log10(DISPLAY_MIN_HZ), np.log10(DISPLAY_MAX_HZ), 1200)" in text
    assert '"display_band_Hz": [float(DISPLAY_MIN_HZ), float(DISPLAY_MAX_HZ)]' in text
    assert '"display_excludes_DC": True' in text
