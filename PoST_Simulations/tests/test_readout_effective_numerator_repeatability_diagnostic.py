from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_effective_numerator_repeatability_diagnostic as diag  # noqa: E402


def test_comparison_with_indices_does_not_mutate_source():
    comparison = {
        "acquisition": {
            "accepted_record_indices": [0, 1, 2, 3]
        }
    }
    result = diag.comparison_with_indices(
        comparison,
        [1, 3],
    )
    assert result["acquisition"][
        "accepted_record_indices"
    ] == [1, 3]
    assert comparison["acquisition"][
        "accepted_record_indices"
    ] == [0, 1, 2, 3]


def test_shape_difference_db_zero_for_identical_shapes():
    values = np.array([1.0, 2.0, 3.0])
    mask = np.array([True, True, False])
    result = diag.shape_difference_db(
        values,
        values,
        mask,
    )
    assert result["rms_dB"] == pytest.approx(0.0)
    assert result["mean_dB"] == pytest.approx(0.0)
    assert result["max_abs_dB"] == pytest.approx(0.0)


def test_parameter_repeatability_detects_day_shift():
    reference = {
        "pole_Hz": 13000.0,
        "pole_Q": 0.70,
        "c2": 100.0,
        "c4": 30.0,
    }
    repeat_full = {
        "pole_Hz": 12000.0,
        "pole_Q": 0.60,
        "c2": 80.0,
        "c4": 20.0,
    }
    blocks = [
        {
            "pole_Hz": 11900.0 + 50.0 * i,
            "pole_Q": 0.59 + 0.005 * i,
            "c2": 79.0 + i,
            "c4": 19.5 + 0.25 * i,
        }
        for i in range(5)
    ]
    result = diag.parameter_repeatability(
        reference,
        repeat_full,
        blocks,
        sigma_threshold=2.0,
    )
    assert result[
        "n_reference_outside_repeat_block_range"
    ] == 4
    assert result[
        "n_reference_shift_exceeds_sigma_threshold"
    ] == 4
    for name in diag.CANONICAL_PARAMETER_NAMES:
        assert result["parameters"][name][
            "reference_outside_repeat_block_range"
        ] is True
        assert result["parameters"][name][
            "reference_shift_over_repeat_block_log10_std"
        ] > 2.0


def test_parameter_repeatability_handles_zero_coefficient_without_log():
    reference = {
        "pole_Hz": 12000.0,
        "pole_Q": 0.6,
        "c2": 80.0,
        "c4": 0.0,
    }
    repeat_full = dict(reference)
    blocks = [dict(reference) for _ in range(5)]
    result = diag.parameter_repeatability(
        reference,
        repeat_full,
        blocks,
        sigma_threshold=2.0,
    )
    assert result["parameters"]["c4"][
        "reference_shift_over_repeat_block_log10_std"
    ] is None
    assert result["parameters"]["c4"][
        "reference_shift_exceeds_sigma_threshold"
    ] is False


def test_default_config_is_order2_and_manifest_cases_exist():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    assert config["effective_numerator_order"] == 2
    assert config["reference_scale_Hz"] == 40000.0

    manifest_path = diag.ident.resolve_config_path(
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
    reference = diag.ident.case_by_label(
        cases,
        config["reference_case_label"],
    )
    repeat = diag.ident.case_by_label(
        cases,
        config["repeat_case_label"],
    )
    assert reference["role"] == "reference"
    assert repeat["role"] == "repeat_validation"


def test_tracked_lowmid_snapshot_supplies_repeat_white_scale():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    snapshot_path = diag.ident.resolve_config_path(
        config["lowmid_free_snapshot"],
        diag.DEFAULT_CONFIG,
    )
    snapshot = diag.ident.load_free_snapshot(
        snapshot_path
    )
    assert snapshot["white_scale"] > 1.0
    assert snapshot["white_asd_A_rtHz"] > 0.0


def test_block_size_mode_uses_reference_record_count():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(encoding="utf-8")
    )
    assert config["block_size"]["mode"] == (
        "reference_accepted_record_count"
    )
    blocks, remainder = diag.day_drift.make_blocks(
        list(range(1731)),
        345,
    )
    assert len(blocks) == 5
    assert all(len(block) == 345 for block in blocks)
    assert remainder == 6
