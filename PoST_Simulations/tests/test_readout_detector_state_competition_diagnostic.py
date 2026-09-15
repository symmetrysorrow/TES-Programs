from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_detector_state_competition_diagnostic as diag  # noqa: E402


def baseline_candidate():
    return {
        "alpha": 200.0,
        "beta": 5.5,
        "C_tes": 2.0e-13,
        "L": 1.1e-8,
        "T_bath": 0.2157,
    }


def test_encode_decode_candidate_roundtrip():
    candidate = baseline_candidate()
    names = (
        "alpha",
        "beta",
        "C_tes",
        "L",
        "T_bath",
    )
    vector = diag.encode_candidate(
        names,
        candidate,
    )
    decoded = diag.decode_candidate(
        vector,
        names,
        candidate,
    )
    for name in names:
        assert decoded[name] == pytest.approx(
            candidate[name]
        )


def test_vector_bounds_log_only_for_ctes_and_l():
    candidate = baseline_candidate()
    physical = {
        "alpha": (10.0, 200.0),
        "beta": (0.0, 12.0),
        "C_tes": (1.0e-13, 1.0e-11),
        "L": (1.0e-10, 1.2e-8),
        "T_bath": (0.2137, 0.2177),
    }
    names = (
        "alpha",
        "beta",
        "C_tes",
        "L",
        "T_bath",
    )
    bounds = diag.vector_bounds(
        names,
        physical,
    )
    assert bounds[0] == pytest.approx(
        physical["alpha"]
    )
    assert bounds[1] == pytest.approx(
        physical["beta"]
    )
    assert bounds[2] == pytest.approx(
        tuple(np.log10(physical["C_tes"]))
    )
    assert bounds[3] == pytest.approx(
        tuple(np.log10(physical["L"]))
    )
    assert bounds[4] == pytest.approx(
        physical["T_bath"]
    )


def test_fixed_white_context_overrides_candidate_dependent_white(monkeypatch):
    monkeypatch.setattr(
        diag.base,
        "intrinsic_context",
        lambda candidate, frequency: {
            "candidate": dict(candidate),
            "frequency_Hz": np.asarray(
                frequency,
                dtype=float,
            ),
            "post_filter_white_asd_A_rtHz": 1.0,
        },
    )
    result = diag.fixed_white_context(
        baseline_candidate(),
        np.asarray([1000.0, 2000.0]),
        6.0e-11,
    )
    assert result[
        "post_filter_white_asd_A_rtHz"
    ] == pytest.approx(6.0e-11)


def test_parameter_boundary_hits_identifies_limits():
    physical = {
        "alpha": (10.0, 200.0),
        "beta": (0.0, 12.0),
    }
    candidate = {
        "alpha": 200.0,
        "beta": 3.0,
    }
    result = diag.parameter_boundary_hits(
        candidate,
        ("alpha", "beta"),
        physical,
    )
    assert result["alpha"]["at_upper"] is True
    assert result["alpha"]["at_lower"] is False
    assert result["beta"]["at_upper"] is False
    assert result["beta"]["at_lower"] is False


def test_default_config_uses_nested_limited_detector_families():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    families = config["nuisance_families"]
    assert [
        row["parameters"]
        for row in families
    ] == [
        ["alpha", "beta"],
        ["alpha", "beta", "C_tes", "L"],
        [
            "alpha",
            "beta",
            "C_tes",
            "L",
            "T_bath",
        ],
    ]
    assert config["fixed_readout"][
        "effective_numerator_order"
    ] == 2
    assert config["fixed_readout"][
        "reference_scale_Hz"
    ] == 40000.0


def test_tracked_snapshot_parses_reference_and_repeat_order2():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    path = diag.resolve_config_path(
        config["order2_day_snapshot"],
        diag.DEFAULT_CONFIG,
    )
    snapshot = diag.load_day_snapshot(path)
    assert snapshot["reference_transfer"][
        "pole_Hz"
    ] > 0.0
    assert snapshot["repeat_transfer"][
        "pole_Hz"
    ] > 0.0
    assert snapshot[
        "repeat_white_asd_A_rtHz"
    ] > 0.0
    assert snapshot[
        "repeat_local_expected_score"
    ] > 0.0


def test_nuisance_bounds_follow_production_constants():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    candidate = baseline_candidate()
    bounds = diag.nuisance_bounds(
        candidate,
        config,
    )
    assert bounds["alpha"][0] == pytest.approx(
        candidate["alpha"] * 0.05
    )
    assert bounds["alpha"][1] == pytest.approx(
        200.0
    )
    assert bounds["beta"] == pytest.approx(
        (0.0, 12.0)
    )
    assert bounds["C_tes"] == pytest.approx(
        (
            diag.opt.C_TES_FIT_MIN_J_PER_K,
            diag.opt.C_TES_FIT_MAX_J_PER_K,
        )
    )
    assert bounds["L"] == pytest.approx(
        (
            diag.opt.L_FIT_MIN_H,
            diag.opt.L_FIT_MAX_H,
        )
    )
    assert bounds["T_bath"] == pytest.approx(
        (
            candidate["T_bath"] - 0.002,
            candidate["T_bath"] + 0.002,
        )
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


def test_fit_nuisance_family_can_start_from_unstable_baseline(monkeypatch):
    baseline = {
        "alpha": 100.0,
        "beta": 2.0,
        "C_tes": 2.0e-13,
        "L": 1.0e-8,
        "T_bath": 0.215,
    }
    stable_trial = dict(baseline)
    stable_trial["alpha"] = 80.0

    def fake_model(candidate, frequency, fixed_transfer, scale_hz, white_asd):
        if candidate["alpha"] >= 99.0:
            return None, {
                "valid": True,
                "stable": False,
                "current_A": 1.0,
                "joule_power_W": 2.0,
            }
        return np.ones_like(np.asarray(frequency, dtype=float)), {
            "valid": True,
            "stable": True,
            "current_A": 1.0,
            "joule_power_W": 2.0,
        }

    monkeypatch.setattr(diag, "model_for_candidate", fake_model)
    monkeypatch.setattr(
        diag.opt,
        "fit_score",
        lambda model, target, frequency, args: float(
            np.sum((np.asarray(model) - 1.0) ** 2)
        ),
    )
    monkeypatch.setattr(
        diag.opt,
        "weighted_residual_vector",
        lambda model, target, frequency, args: (
            np.asarray(model) - 1.0
        ),
    )
    monkeypatch.setattr(
        diag,
        "differential_evolution",
        lambda objective, bounds, **kwargs: SimpleNamespace(
            x=np.asarray([80.0]),
            success=True,
            nfev=1,
        ),
    )
    monkeypatch.setattr(
        diag,
        "least_squares",
        lambda fun, x0, **kwargs: SimpleNamespace(
            x=np.asarray(x0, dtype=float),
            success=True,
            nfev=1,
        ),
    )
    monkeypatch.setattr(
        diag,
        "metrics",
        lambda model, target, frequency, args: {
            "shape_score": 0.0,
            "residual_metrics": {"rms_residual_dB": 0.0},
            "bands": {},
        },
    )
    monkeypatch.setattr(
        diag.opt,
        "tes_operating_point",
        lambda candidate: {
            "valid": True,
            "stable": candidate["alpha"] < 99.0,
            "current_A": 1.0,
            "joule_power_W": 2.0,
        },
    )

    result = diag.fit_nuisance_family(
        name="transition_sensitivity",
        parameter_names=["alpha"],
        baseline_candidate=baseline,
        fixed_transfer={},
        frequency=np.asarray([1.0]),
        target=np.asarray([1.0]),
        args=SimpleNamespace(),
        scale_hz=40000.0,
        white_asd=6.0e-11,
        physical_bounds={"alpha": (10.0, 200.0)},
        optimizer_cfg={
            "instability_penalty": 1e6,
            "DE_maxiter": 1,
            "least_squares_max_nfev": 1,
        },
        seed=1,
        allow_unstable_baseline=True,
    )

    assert result["baseline_was_stable"] is False
    assert result["allow_unstable_baseline"] is True
    assert result["candidate"]["alpha"] == pytest.approx(80.0)
