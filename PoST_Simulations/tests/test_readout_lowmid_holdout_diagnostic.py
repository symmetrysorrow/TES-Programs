from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_lowmid_holdout_diagnostic as diag  # noqa: E402


def test_subset_context_slices_frequency_shaped_arrays_only():
    context = {
        "frequency_Hz": np.array([1.0, 2.0, 3.0]),
        "alias_frequency_Hz": np.array([9.0, 8.0, 7.0]),
        "main_intrinsic_asd": np.array([1.0, 2.0, 3.0]),
        "same_bin": np.array([False, True, False]),
        "scalar": 5.0,
        "candidate": {"x": 1},
    }
    mask = np.array([True, False, True])
    result = diag.subset_context(context, mask)
    np.testing.assert_array_equal(
        result["frequency_Hz"],
        np.array([1.0, 3.0]),
    )
    np.testing.assert_array_equal(
        result["alias_frequency_Hz"],
        np.array([9.0, 7.0]),
    )
    np.testing.assert_array_equal(
        result["same_bin"],
        np.array([False, False]),
    )
    assert result["scalar"] == 5.0
    assert result["candidate"] == {"x": 1}


def test_band_args_changes_only_band_fields():
    original = SimpleNamespace(
        fit_min_hz=1000.0,
        fit_max_hz=200000.0,
        fit_points=601,
        robust_delta_dex=0.3,
        fit_weight_start_hz=40000.0,
        high_frequency_weight=4.0,
        absolute_asd_weight=0.0,
    )
    result = diag.band_args(
        original,
        1000.0,
        40000.0,
        300,
    )
    assert result.fit_min_hz == 1000.0
    assert result.fit_max_hz == 40000.0
    assert result.fit_points == 300
    assert result.robust_delta_dex == 0.3
    assert result.high_frequency_weight == 4.0


def test_drift_toward_shared_detects_parameter_return():
    shared = {
        "pole_Hz": 10.0,
        "pole_Q": 1.0,
        "zero_Hz": 100.0,
        "zero_Q": 1.0,
        "leadlag_zero_Hz": 5.0,
    }
    full = {
        "pole_Hz": 8.0,
        "pole_Q": 0.7,
        "zero_Hz": 150.0,
        "zero_Q": 0.7,
        "leadlag_zero_Hz": 3.0,
    }
    low = {
        "pole_Hz": 9.0,
        "pole_Q": 0.9,
        "zero_Hz": 110.0,
        "zero_Q": 0.9,
        "leadlag_zero_Hz": 4.0,
    }
    result = diag.drift_toward_shared(shared, full, low)
    assert result["n_parameters_closer_to_shared"] == 5
    assert result["fraction_parameters_closer_to_shared"] == pytest.approx(
        1.0
    )
    assert result["parameters"]["zero_Hz"][
        "low_mid_is_closer_to_shared"
    ] is True


def test_load_full_band_snapshot_preserves_white_and_infinite_pole(tmp_path):
    path = tmp_path / "snapshot.json"
    path.write_text(
        json.dumps(
            {
                "provenance": {
                    "stage2_full_band_local_shape_score": 0.0015
                },
                "white_floor": {
                    "best_white_scale": 1.3,
                    "best_white_asd_A_rtHz": 6.0e-11,
                },
                "full_band_transfer_fit": {
                    "fixed_leadlag_pole_Hz": None,
                    "leadlag_pole_is_infinite": True,
                    "parameters": {
                        "pole_Hz": 12000.0,
                        "pole_Q": 0.6,
                        "zero_Hz": 75000.0,
                        "zero_Q": 0.9,
                        "leadlag_zero_Hz": 4000.0,
                    },
                },
            }
        ),
        encoding="utf-8",
    )
    result = diag.load_full_band_snapshot(path)
    assert result["white_scale"] == pytest.approx(1.3)
    assert result["white_asd_A_rtHz"] == pytest.approx(6.0e-11)
    assert result["fixed_leadlag_pole_Hz"] is None
    assert result["parameters"]["zero_Hz"] == pytest.approx(75000.0)


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
        config["white_profiled_full_band_snapshot"],
        diag.DEFAULT_CONFIG,
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    reference = json.loads(reference_path.read_text(encoding="utf-8"))
    snapshot = json.loads(snapshot_path.read_text(encoding="utf-8"))

    assert manifest["cases"]
    assert reference["best_profile_row"]["parameters"]["pole_Hz"] > 0.0
    assert snapshot["white_floor"]["best_white_scale"] > 0.0


def test_repeat_case_exists_and_holdout_is_disjoint():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    manifest_path = diag.resolve_config_path(
        config["manifest"],
        diag.DEFAULT_CONFIG,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = diag.shared.normalize_manifest(manifest, manifest_path)
    repeat = diag.case_by_label(
        cases,
        config["repeat_case_label"],
    )
    assert repeat["role"] == "repeat_validation"

    frequency = np.geomspace(1000.0, 200000.0, 601)
    fit_cfg = config["low_mid_fit"]
    hold_cfg = config["high_frequency_holdout"]
    fit_mask = (
        (frequency >= float(fit_cfg["min_Hz"]))
        & (frequency < float(fit_cfg["max_Hz"]))
    )
    hold_mask = (
        (frequency >= float(hold_cfg["min_Hz"]))
        & (frequency <= float(hold_cfg["max_Hz"]))
    )
    assert not np.any(fit_mask & hold_mask)
    assert np.count_nonzero(fit_mask) > 0
    assert np.count_nonzero(hold_mask) > 0
