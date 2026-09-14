from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import iv_power_thermal_path_diagnostic as diag  # noqa: E402


def args():
    return SimpleNamespace(
        fit_min_hz=1_000.0,
        fit_max_hz=200_000.0,
        fit_weight_start_hz=40_000.0,
        high_frequency_weight=4.0,
        robust_delta_dex=0.3,
    )


def test_target_iv_power_selects_best_rsh_case():
    summary = {
        "best_case_R_SH_ohm": 0.00385,
        "cases": [
            {
                "R_SH_ohm": 0.00380,
                "iv_operating_point": {"P_J_W": 1.0e-9},
            },
            {
                "R_SH_ohm": 0.00385,
                "iv_operating_point": {"P_J_W": 1.2e-9},
            },
        ],
    }
    assert diag.target_iv_power(summary) == pytest.approx(1.2e-9)


def test_grid_includes_exact_baseline():
    linear = diag.grid(1.0, 3.0, 1.7, 3)
    logarithmic = diag.grid(1.0, 100.0, 3.0, 3, log=True)
    assert 1.7 in linear
    assert 3.0 in logarithmic


def test_run_preserves_iv_power_and_detects_joint_improvement(monkeypatch):
    freq = np.geomspace(1_000.0, 200_000.0, 121)
    target = np.ones_like(freq)
    target[(freq >= 5_000.0) & (freq <= 15_000.0)] = 1.25
    target[(freq >= 40_000.0) & (freq <= 100_000.0)] = 0.8
    iv_power = 1.0e-9

    candidate = {
        "thermal_link_model": "stycast_node",
        "T_c": 0.240,
        "T_bath": 0.215,
        "n": 2.0,
        "R": 0.02,
        "G_tes-bath": diag.opt.g_tes_bath_from_joule_power(
            iv_power, 0.240, 0.215, 2.0
        ),
        "G_tes-stycast": 1.0e-9,
    }

    def fake_point(trial):
        power = (
            trial["G_tes-bath"]
            * trial["T_c"]
            * (1.0 - (trial["T_bath"] / trial["T_c"]) ** trial["n"])
            / trial["n"]
        )
        return {
            "valid": True,
            "stable": True,
            "reason": None,
            "joule_power_W": power,
            "current_A": np.sqrt(power / trial["R"]),
        }

    def fake_spectrum(trial, frequency):
        model = np.ones_like(frequency)
        if trial["T_c"] > 0.240:
            model[(frequency >= 5_000.0) & (frequency <= 15_000.0)] = 1.15
            model[(frequency >= 40_000.0) & (frequency <= 100_000.0)] = 0.90
        return model, 1.0

    monkeypatch.setattr(diag.opt, "tes_operating_point", fake_point)
    monkeypatch.setattr(diag.opt, "deterministic_simulated_spectrum", fake_spectrum)

    result = diag.run(
        candidate,
        iv_power,
        target,
        freq,
        args(),
        (0.240, 0.250),
        (2.0,),
        (1.0e-9,),
    )

    assert result["diagnostic_only"] is True
    assert result["all_stable_points_preserve_iv_power"] is True
    assert result[
        "can_move_5_15k_deficit_and_40_100k_excess_in_correct_direction"
    ] is True
    assert result["can_improve_5_15k_and_40_100k_absolute_error"] is True
    assert result[
        "can_improve_5_15k_and_40_100k_without_worsening_100_200k"
    ] is True
    assert max(
        abs(row["iv_power_relative_error"])
        for row in result["rows"]
        if row.get("stable")
    ) < 1e-12
    assert candidate["T_c"] == pytest.approx(0.240)
