from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import etf_ctes_inductance_diagnostic as diag  # noqa: E402


def args():
    return SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )


def test_alpha_inward_grid_never_exceeds_baseline():
    values = diag.alpha_inward_grid(200.0, 7, 0.05)
    assert values[0] == pytest.approx(10.0)
    assert values[-1] == pytest.approx(200.0)
    assert max(values) <= 200.0


def test_grid_includes_exact_baseline():
    values = diag.grid(1.0e-13, 1.0e-11, 2.1e-13, 5, log=True)
    assert 2.1e-13 in values


def test_run_detects_dynamical_joint_improvement(monkeypatch):
    freq = np.geomspace(1_000.0, 200_000.0, 121)
    target = np.ones_like(freq)
    target[(freq >= 5_000.0) & (freq <= 15_000.0)] = 1.25
    target[(freq >= 40_000.0) & (freq <= 100_000.0)] = 0.8

    candidate = {
        "alpha": 200.0,
        "C_tes": 2.0e-13,
        "L": 10.0e-9,
    }

    def fake_point(_trial):
        return {
            "valid": True,
            "stable": True,
            "reason": None,
            "current_A": 1.0,
            "joule_power_W": 2.0,
        }

    def fake_spectrum(trial, frequency):
        model = np.ones_like(frequency)
        if (
            trial["alpha"] < 200.0
            and trial["C_tes"] > 2.0e-13
            and trial["L"] < 10.0e-9
        ):
            model[(frequency >= 5_000.0) & (frequency <= 15_000.0)] = 1.15
            model[(frequency >= 40_000.0) & (frequency <= 100_000.0)] = 0.90
        return model, 1.0

    monkeypatch.setattr(diag.opt, "tes_operating_point", fake_point)
    monkeypatch.setattr(
        diag.opt, "deterministic_simulated_spectrum", fake_spectrum
    )
    monkeypatch.setattr(diag.opt, "C_TES_MATERIAL_J_PER_K", 1.0e-12)

    result = diag.run(
        candidate,
        target,
        freq,
        args(),
        (100.0, 200.0),
        (2.0e-13, 4.0e-13),
        (5.0e-9, 10.0e-9),
    )

    assert result["diagnostic_only"] is True
    assert result[
        "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction"
    ] is True
    assert result["can_improve_5_15k_and_40_100k_absolute_error"] is True
    assert result[
        "can_improve_5_15k_and_40_100k_without_worsening_100_200k"
    ] is True
    best = result["best_correct_direction_row"]
    assert best is not None
    assert best["current_ratio_to_baseline"] == pytest.approx(1.0)
    assert best["joule_power_ratio_to_baseline"] == pytest.approx(1.0)
    assert candidate["alpha"] == pytest.approx(200.0)
    assert candidate["C_tes"] == pytest.approx(2.0e-13)
    assert candidate["L"] == pytest.approx(10.0e-9)
