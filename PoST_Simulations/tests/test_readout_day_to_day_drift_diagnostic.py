from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_day_to_day_drift_diagnostic as diag  # noqa: E402


def test_make_blocks_uses_only_full_nonoverlapping_blocks():
    blocks, remainder = diag.make_blocks(list(range(1037)), 345)
    assert len(blocks) == 3
    assert all(len(block) == 345 for block in blocks)
    assert blocks[0][0] == 0
    assert blocks[1][0] == 345
    assert blocks[2][-1] == 1034
    assert remainder == 2


def test_normalized_transfer_is_unity_at_reference():
    frequency = np.array([1000.0, 5000.0, 10000.0])
    parameters = {
        "pole_Hz": 13000.0,
        "pole_Q": 0.6,
        "zero_Hz": 68000.0,
        "zero_Q": 1.0,
        "leadlag_zero_Hz": 5000.0,
    }
    values = diag.normalized_transfer(
        frequency,
        parameters,
        None,
        1000.0,
    )
    assert values[0] == pytest.approx(1.0)


def test_band_stats_zero_when_empirical_matches_transfer():
    frequency = np.linspace(1000.0, 200000.0, 100)
    empirical = np.sin(np.linspace(0.0, 1.0, 100))
    rows = diag.band_stats(
        frequency,
        empirical,
        empirical.copy(),
        [
            {
                "name": "all",
                "min": 1000.0,
                "max": 200000.0,
            }
        ],
    )
    assert rows["all"]["empirical_minus_transfer_mean_dB"] == pytest.approx(
        0.0
    )
    assert rows["all"]["empirical_minus_transfer_rms_dB"] == pytest.approx(
        0.0
    )


def test_block_envelope_and_day_comparison_identify_large_day_drift():
    frequency = np.array(
        [1000.0, 5000.0, 10000.0, 20000.0]
    )
    repeat = np.ones(4)
    blocks = [
        10.0 ** (np.array([0.0, 0.1, -0.1, 0.05]) / 20.0),
        10.0 ** (np.array([0.0, -0.1, 0.1, -0.05]) / 20.0),
        np.ones(4),
    ]
    fit_mask = np.ones(4, dtype=bool)
    anchors = frequency.tolist()
    envelope = diag.block_envelope(
        frequency,
        repeat,
        blocks,
        anchors,
        fit_mask,
    )
    empirical_day_db = np.array([0.0, 1.0, -1.0, 0.8])
    comparison = diag.compare_day_to_block(
        frequency,
        empirical_day_db,
        envelope["anchor_envelope"],
        anchors,
        fit_mask,
        envelope["fit_band_rms_max_dB"],
    )
    assert comparison[
        "day_fit_band_rms_exceeds_all_blocks"
    ] is True
    assert comparison[
        "n_anchors_day_exceeds_all_blocks"
    ] == 3


def test_load_repeat_transfer_preserves_infinite_pole(tmp_path):
    path = tmp_path / "repeat.json"
    path.write_text(
        json.dumps(
            {
                "transfer": {
                    "fixed_leadlag_pole_Hz": None,
                    "leadlag_pole_is_infinite": True,
                    "parameters": {
                        "pole_Hz": 12760.0,
                        "pole_Q": 0.54,
                        "zero_Hz": 64643.0,
                        "zero_Q": 1.06,
                        "leadlag_zero_Hz": 3750.0,
                    },
                }
            }
        ),
        encoding="utf-8",
    )
    result = diag.load_repeat_transfer(path)
    assert result["fixed_leadlag_pole_Hz"] is None
    assert result["parameters"]["leadlag_zero_Hz"] == pytest.approx(
        3750.0
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
    repeat_transfer_path = diag.resolve_config_path(
        config["repeat_local_transfer"],
        diag.DEFAULT_CONFIG,
    )

    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    reference = json.loads(
        reference_transfer_path.read_text(encoding="utf-8")
    )
    repeat = json.loads(
        repeat_transfer_path.read_text(encoding="utf-8")
    )

    assert manifest["cases"]
    assert reference["best_profile_row"]["parameters"]["pole_Hz"] > 0.0
    assert repeat["transfer"]["parameters"]["pole_Hz"] > 0.0


def test_config_labels_exist_in_manifest():
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
    reference = diag.case_by_label(
        cases,
        config["reference_case_label"],
    )
    repeat = diag.case_by_label(
        cases,
        config["repeat_case_label"],
    )
    assert reference["role"] == "reference"
    assert repeat["role"] == "repeat_validation"
