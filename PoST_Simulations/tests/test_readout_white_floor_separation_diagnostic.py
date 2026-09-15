from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_white_floor_separation_diagnostic as diag  # noqa: E402


def test_white_scaled_context_changes_only_white_floor():
    context = {
        "post_filter_white_asd_A_rtHz": 4.0,
        "frequency_Hz": np.array([1.0, 2.0]),
        "marker": "keep",
    }
    result = diag.white_scaled_context(context, 2.5)
    assert result["post_filter_white_asd_A_rtHz"] == pytest.approx(10.0)
    assert result["marker"] == "keep"
    np.testing.assert_array_equal(
        result["frequency_Hz"],
        context["frequency_Hz"],
    )
    assert context["post_filter_white_asd_A_rtHz"] == 4.0


def test_white_scaled_context_accepts_exact_zero():
    context = {"post_filter_white_asd_A_rtHz": 3.0}
    result = diag.white_scaled_context(context, 0.0)
    assert result["post_filter_white_asd_A_rtHz"] == 0.0


def test_profile_white_floor_finds_known_scale(monkeypatch):
    monkeypatch.setattr(
        diag,
        "shared_model_for_scale",
        lambda context, scale, shared_transfer: np.asarray(
            [float(scale)]
        ),
    )
    monkeypatch.setattr(
        diag.opt,
        "fit_score",
        lambda model, target, frequency, fit_args: float(
            (float(model[0]) - 2.0) ** 2 + 0.5
        ),
    )

    result = diag.profile_white_floor(
        {"post_filter_white_asd_A_rtHz": 5.0},
        np.asarray([1.0]),
        np.asarray([1000.0]),
        SimpleNamespace(),
        {
            "parameters": {},
            "fixed_leadlag_pole_Hz": None,
        },
        {
            "min_positive_scale": 0.01,
            "max_scale": 30.0,
            "grid_points": 21,
            "include_exact_zero": True,
            "include_baseline_scale_one": True,
        },
    )
    assert result["best"]["white_scale"] == pytest.approx(
        2.0, rel=1e-4
    )
    assert result["best"]["white_asd_A_rtHz"] == pytest.approx(
        10.0, rel=1e-4
    )
    assert any(row["white_scale"] == 0.0 for row in result["grid"])
    assert any(row["white_scale"] == 1.0 for row in result["grid"])


def test_load_previous_local_transfer_preserves_infinity(tmp_path):
    path = tmp_path / "local.json"
    path.write_text(
        json.dumps(
            {
                "transfer": {
                    "fixed_leadlag_pole_Hz": None,
                    "leadlag_pole_is_infinite": True,
                    "parameters": {
                        "pole_Hz": 12000.0,
                        "pole_Q": 0.5,
                        "zero_Hz": 65000.0,
                        "zero_Q": 1.1,
                        "leadlag_zero_Hz": 4000.0,
                    },
                }
            }
        ),
        encoding="utf-8",
    )
    result = diag.load_previous_local_transfer(path)
    assert result["fixed_leadlag_pole_Hz"] is None
    assert result["parameters"]["leadlag_zero_Hz"] == pytest.approx(
        4000.0
    )


def test_parameter_drift_zero_for_shared_parameters():
    parameters = {
        "pole_Hz": 13000.0,
        "pole_Q": 0.6,
        "zero_Hz": 68000.0,
        "zero_Q": 1.0,
        "leadlag_zero_Hz": 5000.0,
    }
    result = diag.parameter_drift(parameters, dict(parameters))
    assert result["rms_log10_parameter_ratio"] == pytest.approx(0.0)
    assert result["max_abs_log10_parameter_ratio"] == pytest.approx(
        0.0
    )


def test_default_git_tracked_inputs_parse():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    manifest_path = diag.resolve_config_path(
        config["manifest"],
        diag.DEFAULT_CONFIG,
    )
    reference_transfer_path = diag.resolve_config_path(
        config["reference_transfer"],
        diag.DEFAULT_CONFIG,
    )
    previous_local_path = diag.resolve_config_path(
        config["previous_repeat_local_transfer"],
        diag.DEFAULT_CONFIG,
    )

    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    reference = json.loads(
        reference_transfer_path.read_text(encoding="utf-8")
    )
    previous = json.loads(
        previous_local_path.read_text(encoding="utf-8")
    )

    assert manifest["cases"]
    assert reference["best_profile_row"]["parameters"]["pole_Hz"] > 0.0
    assert previous["transfer"]["parameters"]["pole_Hz"] > 0.0


def test_repeat_case_exists_and_is_repeat_validation():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    manifest_path = diag.resolve_config_path(
        config["manifest"],
        diag.DEFAULT_CONFIG,
    )
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    cases = diag.shared.normalize_manifest(
        manifest,
        manifest_path,
    )
    repeat = diag.case_by_label(
        cases,
        config["repeat_case_label"],
    )
    assert repeat["role"] == "repeat_validation"
