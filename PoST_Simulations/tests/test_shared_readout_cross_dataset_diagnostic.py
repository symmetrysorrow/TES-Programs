from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import shared_readout_cross_dataset_diagnostic as diag  # noqa: E402


def reference_profile():
    return {
        "best_profile_row": {
            "shape_score": 0.0017,
            "leadlag_pole_is_infinite": True,
            "parameters": {
                "pole_Hz": 13_342.0,
                "pole_Q": 0.622,
                "zero_Hz": 68_502.0,
                "zero_Q": 1.021,
                "leadlag_zero_Hz": 4_940.0,
                "leadlag_pole_Hz": None,
                "leadlag_pole_is_infinite": True,
            },
        }
    }


def test_load_shared_transfer_preserves_exact_infinite_pole():
    shared = diag.load_shared_transfer(reference_profile())
    assert shared["leadlag_pole_is_infinite"] is True
    assert shared["fixed_leadlag_pole_Hz"] is None
    assert shared["parameters"]["leadlag_zero_Hz"] == pytest.approx(
        4_940.0
    )
    assert shared["parameters"]["pole_Hz"] == pytest.approx(
        13_342.0
    )


def test_load_shared_transfer_accepts_finite_pole():
    value = reference_profile()
    row = value["best_profile_row"]
    row["leadlag_pole_is_infinite"] = False
    row["parameters"]["leadlag_pole_is_infinite"] = False
    row["parameters"]["leadlag_pole_Hz"] = 400_000.0
    shared = diag.load_shared_transfer(value)
    assert shared["leadlag_pole_is_infinite"] is False
    assert shared["fixed_leadlag_pole_Hz"] == pytest.approx(
        400_000.0
    )


def test_normalize_manifest_resolves_relative_paths(tmp_path):
    manifest_path = tmp_path / "manifest.json"
    manifest = {
        "cases": [
            {
                "label": "reference",
                "role": "reference",
                "summary": "a/summary.json",
                "comparison_summary": "b/comparison_summary.json",
            },
            {
                "label": "validation",
                "summary": "c/summary.json",
                "comparison_summary": "d/comparison_summary.json",
                "experiment_path": "G:/data/example",
            },
        ]
    }
    rows = diag.normalize_manifest(manifest, manifest_path)
    assert rows[0]["role"] == "reference"
    assert rows[1]["role"] == "validation"
    assert rows[0]["summary"] == (
        tmp_path / "a/summary.json"
    ).resolve()
    assert rows[0]["comparison_summary"] == (
        tmp_path / "b/comparison_summary.json"
    ).resolve()
    assert str(rows[1]["experiment_path"]).replace("\\", "/") == (
        "G:/data/example"
    )


def test_normalize_manifest_rejects_duplicate_labels(tmp_path):
    manifest = {
        "cases": [
            {
                "label": "same",
                "summary": "a.json",
                "comparison_summary": "b.json",
            },
            {
                "label": "same",
                "summary": "c.json",
                "comparison_summary": "d.json",
            },
        ]
    }
    with pytest.raises(ValueError):
        diag.normalize_manifest(manifest, tmp_path / "manifest.json")


def test_parameter_drift_zero_for_identical_transfer():
    shared = diag.load_shared_transfer(reference_profile())[
        "parameters"
    ]
    drift = diag.parameter_drift(shared, dict(shared))
    assert drift["rms_log10_parameter_ratio"] == pytest.approx(0.0)
    assert drift["max_abs_log10_parameter_ratio"] == pytest.approx(
        0.0
    )


def test_evaluation_flags_require_tail_for_strict_improvement():
    baseline = {
        "5000-15000_Hz": {"mean_model_over_target_dB": -4.0},
        "40000-100000_Hz": {"mean_model_over_target_dB": 4.0},
        "100000-200000_Hz": {"mean_model_over_target_dB": 0.5},
    }
    shared = {
        "5000-15000_Hz": {"mean_model_over_target_dB": -0.2},
        "40000-100000_Hz": {"mean_model_over_target_dB": 0.3},
        "100000-200000_Hz": {"mean_model_over_target_dB": -0.8},
    }
    flags = diag.evaluation_flags(
        0.08,
        0.01,
        baseline,
        shared,
        local_score=0.009,
        shared_local_tolerance=1.25,
    )
    assert flags["shared_improves_shape_score"] is True
    assert flags["shared_improves_5_15k"] is True
    assert flags["shared_improves_40_100k"] is True
    assert flags[
        "shared_does_not_worsen_100_200k_mean"
    ] is False
    assert flags["shared_strict_broad_improvement"] is False
    assert flags["shared_over_local_best_score_ratio"] == pytest.approx(
        0.01 / 0.009
    )
    assert flags["shared_within_local_score_tolerance"] is True


def test_aggregate_requires_independent_validation_case():
    reference_row = {
        "role": "reference",
        "shared_transfer": {"score_ratio_to_baseline": 0.02},
        "flags": {
            "shared_improves_shape_score": True,
            "shared_improves_5_15k": True,
            "shared_improves_40_100k": True,
            "shared_over_local_best_score_ratio": 1.0,
            "shared_within_local_score_tolerance": True,
        },
    }
    result = diag.aggregate_results([reference_row], 1.25)
    assert result["n_validation_cases"] == 0
    assert result["interpretation_flags"][
        "supports_shared_readout_transfer_across_validation_cases"
    ] is False
    assert result["interpretation_flags"][
        "validation_dataset_required_for_cross_validation_claim"
    ] is True


def test_aggregate_supports_shared_transfer_when_validation_passes():
    rows = [
        {
            "role": "reference",
            "shared_transfer": {"score_ratio_to_baseline": 0.02},
            "flags": {
                "shared_improves_shape_score": True,
                "shared_improves_5_15k": True,
                "shared_improves_40_100k": True,
                "shared_over_local_best_score_ratio": 1.01,
                "shared_within_local_score_tolerance": True,
            },
        },
        {
            "role": "validation",
            "shared_transfer": {"score_ratio_to_baseline": 0.2},
            "flags": {
                "shared_improves_shape_score": True,
                "shared_improves_5_15k": True,
                "shared_improves_40_100k": True,
                "shared_over_local_best_score_ratio": 1.10,
                "shared_within_local_score_tolerance": True,
            },
        },
    ]
    result = diag.aggregate_results(rows, 1.25)
    assert result["n_validation_cases"] == 1
    assert result["all_validation_cases_shared_improve_score"] is True
    assert result[
        "all_validation_cases_improve_mid_and_high"
    ] is True
    assert result[
        "all_validation_cases_within_local_score_tolerance"
    ] is True
    assert result["interpretation_flags"][
        "supports_shared_readout_transfer_across_validation_cases"
    ] is True


def test_normalized_model_from_components():
    frequency = np.array([1_000.0, 2_000.0])
    components = {
        "main_asd": np.array([3.0, 4.0]),
        "alias_asd": np.array([4.0, 3.0]),
        "white_asd": np.array([0.0, 0.0]),
    }
    result = diag.normalized_model_from_components(
        frequency,
        components,
    )
    np.testing.assert_allclose(result, np.ones_like(result))
