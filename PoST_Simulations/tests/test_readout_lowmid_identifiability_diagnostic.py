from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_lowmid_identifiability_diagnostic as diag  # noqa: E402


def test_decode_free_merges_fixed_and_free_parameters():
    fixed = {"zero_Hz": 70000.0, "zero_Q": 1.0}
    free_names = ("pole_Hz", "pole_Q", "leadlag_zero_Hz")
    vector = np.log10(np.array([12000.0, 0.6, 5000.0]))
    result = diag.decode_free(vector, free_names, fixed)
    assert result["zero_Hz"] == pytest.approx(70000.0)
    assert result["zero_Q"] == pytest.approx(1.0)
    assert result["pole_Hz"] == pytest.approx(12000.0)
    assert result["pole_Q"] == pytest.approx(0.6)
    assert result["leadlag_zero_Hz"] == pytest.approx(5000.0)


def test_summarize_profile_detects_broad_nonidentifiability():
    rows = []
    for value in (30000.0, 60000.0, 120000.0, 240000.0):
        rows.append(
            {
                "fixed_parameters": {"zero_Hz": value},
                "shape_score": 1.0,
                "near_free": True,
            }
        )
    result = diag.summarize_profile("zero_Hz", rows)
    assert result["n_near_free_points"] == 4
    assert result["near_free_span_ratio"] == pytest.approx(8.0)
    assert result["broad_nonidentifiability_over_tested_range"] is True


def test_annotate_profile_requires_score_and_rms_screen():
    rows = [
        {
            "shape_score": 2.0,
            "residual_metrics": {"rms_residual_dB": 0.12},
        },
        {
            "shape_score": 20.0,
            "residual_metrics": {"rms_residual_dB": 0.12},
        },
        {
            "shape_score": 2.0,
            "residual_metrics": {"rms_residual_dB": 0.25},
        },
    ]
    result = diag.annotate_profile(
        rows,
        free_score=1.0,
        free_rms=0.10,
        screen={
            "near_free_rms_delta_dB": 0.05,
            "near_free_score_ratio": 5.0,
        },
    )
    assert result[0]["near_free"] is True
    assert result[1]["near_free"] is False
    assert result[2]["near_free"] is False


def test_fit_with_fixed_parameters_preserves_fixed_value(monkeypatch):
    target_parameters = {
        "pole_Hz": 12000.0,
        "pole_Q": 0.6,
        "zero_Hz": 70000.0,
        "zero_Q": 1.0,
        "leadlag_zero_Hz": 5000.0,
    }

    monkeypatch.setattr(
        diag.profile,
        "pre_analysis_model",
        lambda context, parameters, fixed_pole: np.asarray(
            [
                parameters["pole_Hz"] / 12000.0,
                parameters["pole_Q"] / 0.6,
                parameters["zero_Hz"] / 70000.0,
                parameters["zero_Q"],
                parameters["leadlag_zero_Hz"] / 5000.0,
            ]
        ),
    )
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
            x=np.asarray(
                [
                    np.mean(bound)
                    for bound in bounds
                ],
                dtype=float,
            ),
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
        diag.base,
        "residual_db_metrics",
        lambda *args, **kwargs: {
            "rms_residual_dB": 0.0,
            "max_abs_residual_dB": 0.0,
        },
    )
    monkeypatch.setattr(
        diag.base,
        "band_summary",
        lambda *args, **kwargs: {},
    )

    result = diag.fit_with_fixed_parameters(
        {},
        np.ones(5),
        np.arange(5.0),
        SimpleNamespace(),
        fixed_parameters={"zero_Hz": 70000.0},
        shared_parameters=target_parameters,
        free_reference_parameters=target_parameters,
        center_min_hz=1000.0,
        center_max_hz=300000.0,
        general_q_min=0.1,
        q_max=20.0,
        seed=1,
        de_maxiter=1,
    )
    assert result["parameters"]["zero_Hz"] == pytest.approx(70000.0)
    assert result["shape_score"] == pytest.approx(0.0, abs=1e-15)
    assert result["best_candidate_source"] in {
        "shared_exact",
        "free_reference_exact",
    }


def test_default_git_tracked_inputs_parse():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    manifest_path = diag.resolve_config_path(
        config["manifest"],
        diag.DEFAULT_CONFIG,
    )
    reference_path = diag.resolve_config_path(
        config["reference_transfer"],
        diag.DEFAULT_CONFIG,
    )
    snapshot_path = diag.resolve_config_path(
        config["lowmid_free_snapshot"],
        diag.DEFAULT_CONFIG,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))

    assert manifest["cases"]
    assert reference["best_profile_row"]["parameters"]["zero_Hz"] > 0.0
    assert snapshot["transfer"]["parameters"]["leadlag_zero_Hz"] > 0.0
    assert "zero_Hz" in config["fixed_profiles"]
    assert "leadlag_zero_Hz" in config["fixed_profiles"]
    assert "pole_Hz" in config["fixed_profiles"]


def test_reduced_topology_specs_only_use_known_parameters():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    for spec in config["reduced_topologies"]:
        assert set(spec["fixed_parameters"]) <= set(
            diag.PARAMETER_NAMES
        )
