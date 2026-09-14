from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))

from subScript import thermal_path_partition_diagnostic as diagnostic  # noqa: E402


def _args() -> SimpleNamespace:
    return SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )


def test_grid_with_baseline_includes_exact_baseline() -> None:
    grid = diagnostic.grid_with_baseline(1.0, 100.0, 3.0, 5)
    assert grid[0] == pytest.approx(1.0)
    assert grid[-1] == pytest.approx(100.0)
    assert 3.0 in grid


def test_grid_with_baseline_rejects_out_of_bounds_baseline() -> None:
    with pytest.raises(ValueError, match="outside physical bounds"):
        diagnostic.grid_with_baseline(1.0, 10.0, 20.0, 5)


def test_thermal_path_scan_finds_joint_correct_direction(monkeypatch) -> None:
    fit_freq = np.geomspace(1_000.0, 200_000.0, 121)
    target = np.ones_like(fit_freq)
    target[(fit_freq >= 5_000.0) & (fit_freq <= 15_000.0)] = 1.25
    target[(fit_freq >= 40_000.0) & (fit_freq <= 100_000.0)] = 1.0 / 1.2

    candidate = {
        "thermal_link_model": "stycast_node",
        "G_tes-stycast": 1.0,
        "G_stycast-abs": 1.0,
        "G_tes-bath": 10.0,
    }

    def fake_operating_point(_candidate):
        return {
            "valid": True,
            "stable": True,
            "reason": None,
            "current_A": 1.0,
            "joule_power_W": 2.0,
        }

    def fake_spectrum(trial, frequency):
        model = np.ones_like(frequency)
        if trial["G_tes-stycast"] > 1.0:
            model[(frequency >= 5_000.0) & (frequency <= 15_000.0)] = 1.1
            model[(frequency >= 40_000.0) & (frequency <= 100_000.0)] = 0.9
        return model, 1.0

    monkeypatch.setattr(
        diagnostic.optimizer,
        "tes_operating_point",
        fake_operating_point,
    )
    monkeypatch.setattr(
        diagnostic.optimizer,
        "deterministic_simulated_spectrum",
        fake_spectrum,
    )

    result = diagnostic.run_diagnostic(
        candidate,
        target,
        fit_freq,
        _args(),
        g_tes_stycast_grid=(1.0, 2.0),
        g_stycast_abs_grid=(1.0,),
    )

    assert result["diagnostic_only"] is True
    assert result["production_optimizer_unchanged"] is True
    assert result["added_noise_source"] is False
    assert result[
        "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction"
    ] is True
    assert result["can_improve_5_15k_and_40_100k_absolute_error"] is True
    best = result["best_correct_direction_row"]
    assert best is not None
    assert best["G_tes_stycast_W_per_K"] == pytest.approx(2.0)
    assert best["current_ratio_to_baseline"] == pytest.approx(1.0)
    assert best["joule_power_ratio_to_baseline"] == pytest.approx(1.0)
    assert candidate["G_tes-stycast"] == pytest.approx(1.0)


def test_fit_args_from_summary_preserves_objective_settings() -> None:
    summary = {
        "fit": {
            "min_hz": 1_000.0,
            "max_hz": 200_000.0,
            "weight_start_hz": 40_000.0,
            "high_frequency_weight": 4.0,
            "robust_delta_dex": 0.3,
            "points": 601,
            "absolute_asd_weight": 0.0,
        }
    }
    args = diagnostic.fit_args_from_summary(summary)
    assert args.fit_points == 601
    assert args.high_frequency_weight == pytest.approx(4.0)
    assert args.absolute_asd_weight == pytest.approx(0.0)
