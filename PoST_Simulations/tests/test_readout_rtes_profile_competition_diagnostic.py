from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_rtes_profile_competition_diagnostic as diag  # noqa: E402


def test_tracked_detector_state_snapshot_parses():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    path = diag.resolve_config_path(
        config["detector_state_snapshot"],
        diag.DEFAULT_CONFIG,
    )
    snapshot = diag.load_detector_state_snapshot(
        path
    )
    assert snapshot["inherited_R_TES_Ohm"] > 0.0
    assert snapshot["best_shape_score"] > 0.0
    assert snapshot["repeat_local_shape_score"] > 0.0
    assert snapshot["best_family_name"] == (
        "transition_local_plus_bath"
    )


def test_default_profile_contains_positive_unique_R1_anchor():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    ratios = [
        float(value)
        for value in config["R_TES_profile"]["ratios"]
    ]
    assert all(value > 0.0 for value in ratios)
    assert len(ratios) == len(set(ratios))
    assert 1.0 in ratios
    assert min(ratios) == pytest.approx(0.75)
    assert max(ratios) == pytest.approx(1.50)


def test_default_nuisance_families_do_not_include_R():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    for family in config["nuisance_families"]:
        assert "R" not in family["parameters"]


def test_profile_point_holds_R_as_outer_coordinate(monkeypatch):
    inherited = {
        "R": 0.020,
        "alpha": 100.0,
        "beta": 2.0,
        "C_tes": 2.0e-13,
        "L": 1.0e-8,
        "T_bath": 0.215,
    }
    seen = {}

    monkeypatch.setattr(
        diag.opt,
        "tes_operating_point",
        lambda candidate: {
            "valid": True,
            "stable": True,
            "current_A": 1.0,
            "joule_power_W": 2.0,
        },
    )
    monkeypatch.setattr(
        diag,
        "full_model",
        lambda candidate, frequency, transfer, scale_hz, white_asd: (
            np.ones_like(np.asarray(frequency, dtype=float))
        ),
    )
    monkeypatch.setattr(
        diag.holdout,
        "model_metrics",
        lambda *args, **kwargs: {
            "shape_score": 1.0,
            "residual_metrics": {
                "rms_residual_dB": 0.1
            },
            "bands": {},
        },
    )

    def fake_fit(**kwargs):
        seen["R"] = kwargs["baseline_candidate"]["R"]
        seen["parameters"] = tuple(
            kwargs["parameter_names"]
        )
        trial = dict(kwargs["baseline_candidate"])
        return {
            "name": kwargs["name"],
            "parameters_varied": list(
                kwargs["parameter_names"]
            ),
            "shape_score": 0.5,
            "residual_metrics": {
                "rms_residual_dB": 0.2
            },
            "bands": {},
            "boundary_hits": {},
            "operating_point": {},
            "optimizer": {},
            "_candidate_full": trial,
            "_model": np.ones(3),
        }

    monkeypatch.setattr(
        diag.competition,
        "fit_nuisance_family",
        fake_fit,
    )

    result = diag.profile_point(
        ratio=1.10,
        inherited_candidate=inherited,
        families=[
            {
                "name": "transition_sensitivity",
                "parameters": ["alpha", "beta"],
            }
        ],
        reference_transfer={},
        full_frequency=np.array([1.0, 2.0, 3.0]),
        fit_frequency=np.array([1.0, 2.0]),
        fit_target=np.ones(2),
        fit_args=SimpleNamespace(),
        hold_frequency=np.array([3.0]),
        hold_target=np.ones(1),
        hold_args=SimpleNamespace(),
        full_target=np.ones(3),
        full_args=SimpleNamespace(),
        hold_mask=np.array([False, False, True]),
        scale_hz=40000.0,
        white_asd=6.0e-11,
        physical_bounds={
            "alpha": (5.0, 200.0),
            "beta": (0.0, 12.0),
        },
        optimizer_cfg={
            "DE_maxiter": 1,
            "least_squares_max_nfev": 1,
            "instability_penalty": 1e6,
        },
        seed=1,
        warm_candidates=[],
        repeat_local_score=0.1,
        repeat_local_rms=0.1,
        screen={
            "near_repeat_local_score_ratio": 5.0,
            "near_repeat_local_rms_delta_dB": 0.05,
        },
    )

    assert seen["R"] == pytest.approx(0.022)
    assert "R" not in seen["parameters"]
    assert result["R_ratio_to_inherited"] == pytest.approx(
        1.10
    )
    assert result["R_TES_Ohm"] == pytest.approx(
        0.022
    )


def test_default_bounds_match_detector_competition_schema():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    inherited = {
        "alpha": 200.0,
        "beta": 5.0,
        "C_tes": 2.0e-13,
        "L": 1.0e-8,
        "T_bath": 0.215,
    }
    bounds = diag.competition.nuisance_bounds(
        inherited,
        {"bounds": config["bounds"]},
    )
    assert bounds["alpha"][0] == pytest.approx(10.0)
    assert bounds["alpha"][1] == pytest.approx(200.0)
    assert bounds["beta"] == pytest.approx((0.0, 12.0))
    assert bounds["T_bath"] == pytest.approx(
        (0.213, 0.217)
    )


def test_fit_and_holdout_regions_are_disjoint():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    frequency = np.geomspace(
        1000.0,
        200000.0,
        601,
    )
    fit = config["fit_region_Hz"]
    hold = config["holdout_region_Hz"]
    fit_mask = (
        (frequency >= float(fit["min"]))
        & (
            frequency
            < float(fit["max_exclusive"])
        )
    )
    hold_mask = (
        (frequency >= float(hold["min"]))
        & (frequency <= float(hold["max"]))
    )
    assert not np.any(fit_mask & hold_mask)


def test_profile_point_searches_nuisance_when_inherited_state_is_unstable(monkeypatch):
    inherited = {
        "R": 0.020,
        "alpha": 100.0,
        "beta": 2.0,
        "C_tes": 2.0e-13,
        "L": 1.0e-8,
        "T_bath": 0.215,
    }
    calls = {"fit": 0}

    monkeypatch.setattr(
        diag.opt,
        "tes_operating_point",
        lambda candidate: {
            "valid": True,
            "stable": False,
            "reason": "positive_or_zero_eigenvalue",
            "current_A": 1.0,
            "joule_power_W": 2.0,
        },
    )

    def fake_fit(**kwargs):
        calls["fit"] += 1
        assert kwargs["allow_unstable_baseline"] is True
        trial = dict(kwargs["baseline_candidate"])
        trial["alpha"] = 50.0
        return {
            "name": kwargs["name"],
            "parameters_varied": list(kwargs["parameter_names"]),
            "shape_score": 0.2,
            "residual_metrics": {"rms_residual_dB": 0.2},
            "bands": {},
            "boundary_hits": {},
            "operating_point": {
                "valid": True,
                "stable": True,
            },
            "optimizer": {},
            "_candidate_full": trial,
            "_model": np.ones(3),
        }

    monkeypatch.setattr(
        diag.competition,
        "fit_nuisance_family",
        fake_fit,
    )
    monkeypatch.setattr(
        diag,
        "full_model",
        lambda candidate, frequency, transfer, scale_hz, white_asd: (
            np.ones_like(np.asarray(frequency, dtype=float))
        ),
    )
    monkeypatch.setattr(
        diag.holdout,
        "model_metrics",
        lambda *args, **kwargs: {
            "shape_score": 0.2,
            "residual_metrics": {"rms_residual_dB": 0.2},
            "bands": {},
        },
    )

    result = diag.profile_point(
        ratio=0.90,
        inherited_candidate=inherited,
        families=[
            {
                "name": "transition_sensitivity",
                "parameters": ["alpha", "beta"],
            }
        ],
        reference_transfer={},
        full_frequency=np.array([1.0, 2.0, 3.0]),
        fit_frequency=np.array([1.0, 2.0]),
        fit_target=np.ones(2),
        fit_args=SimpleNamespace(),
        hold_frequency=np.array([3.0]),
        hold_target=np.ones(1),
        hold_args=SimpleNamespace(),
        full_target=np.ones(3),
        full_args=SimpleNamespace(),
        hold_mask=np.array([False, False, True]),
        scale_hz=40000.0,
        white_asd=6.0e-11,
        physical_bounds={
            "alpha": (5.0, 200.0),
            "beta": (0.0, 12.0),
        },
        optimizer_cfg={
            "DE_maxiter": 1,
            "least_squares_max_nfev": 1,
            "instability_penalty": 1e6,
        },
        seed=1,
        warm_candidates=[],
        repeat_local_score=0.1,
        repeat_local_rms=0.1,
        screen={
            "near_repeat_local_score_ratio": 5.0,
            "near_repeat_local_rms_delta_dB": 0.2,
        },
    )

    assert calls["fit"] == 1
    assert result["status"] == "evaluated"
    assert result["inherited_nuisance_state_was_stable"] is False
    assert result["fixed_reference_readout_inherited_nuisance_metrics"] is None


def test_profile_point_reports_no_stable_solution_only_after_family_search(monkeypatch):
    inherited = {
        "R": 0.020,
        "alpha": 100.0,
        "beta": 2.0,
        "C_tes": 2.0e-13,
        "L": 1.0e-8,
        "T_bath": 0.215,
    }
    monkeypatch.setattr(
        diag.opt,
        "tes_operating_point",
        lambda candidate: {
            "valid": True,
            "stable": False,
            "reason": "positive_or_zero_eigenvalue",
            "current_A": 1.0,
            "joule_power_W": 2.0,
        },
    )
    monkeypatch.setattr(
        diag.competition,
        "fit_nuisance_family",
        lambda **kwargs: (_ for _ in ()).throw(
            RuntimeError("no stable candidate found")
        ),
    )

    result = diag.profile_point(
        ratio=0.85,
        inherited_candidate=inherited,
        families=[
            {
                "name": "transition_sensitivity",
                "parameters": ["alpha", "beta"],
            }
        ],
        reference_transfer={},
        full_frequency=np.array([1.0, 2.0, 3.0]),
        fit_frequency=np.array([1.0, 2.0]),
        fit_target=np.ones(2),
        fit_args=SimpleNamespace(),
        hold_frequency=np.array([3.0]),
        hold_target=np.ones(1),
        hold_args=SimpleNamespace(),
        full_target=np.ones(3),
        full_args=SimpleNamespace(),
        hold_mask=np.array([False, False, True]),
        scale_hz=40000.0,
        white_asd=6.0e-11,
        physical_bounds={
            "alpha": (5.0, 200.0),
            "beta": (0.0, 12.0),
        },
        optimizer_cfg={
            "DE_maxiter": 1,
            "least_squares_max_nfev": 1,
            "instability_penalty": 1e6,
        },
        seed=1,
        warm_candidates=[],
        repeat_local_score=0.1,
        repeat_local_rms=0.1,
        screen={
            "near_repeat_local_score_ratio": 5.0,
            "near_repeat_local_rms_delta_dB": 0.2,
        },
    )

    assert result["status"] == "no_stable_solution_at_fixed_R"
    assert result["best_family"] is None
    assert result["family_fits"][0]["status"] == "fit_error"
