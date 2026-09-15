from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import readout_residual_dof_competition_diagnostic as diag  # noqa: E402
from subScript import readout_effective_numerator_diagnostic as effective  # noqa: E402


def test_canonical_to_effective_preserves_polynomial():
    readout = {
        "pole_Hz": 12000.0,
        "pole_Q": 0.7,
        "c2": 84.0,
        "c4": 24.0,
    }
    frequency = np.geomspace(
        1000.0,
        500000.0,
        500,
    )
    scale = 40000.0
    mapped = diag.canonical_to_effective(readout)
    direct = diag.canonical_polynomial(
        frequency,
        readout,
        scale,
    )
    factored = effective.effective_numerator_squared(
        frequency,
        2,
        mapped["latent"],
        scale,
    )
    np.testing.assert_allclose(
        direct,
        factored,
        rtol=2e-13,
        atol=2e-13,
    )


def test_canonical_mapping_is_positive_for_nonnegative_coefficients():
    frequency = np.geomspace(
        1.0,
        1.0e7,
        500,
    )
    for c2, c4 in (
        (0.0, 0.0),
        (42.0, 34.0),
        (84.0, 24.0),
        (300.0, 150.0),
    ):
        readout = {
            "pole_Hz": 12000.0,
            "pole_Q": 0.7,
            "c2": c2,
            "c4": c4,
        }
        values = diag.canonical_polynomial(
            frequency,
            readout,
            40000.0,
        )
        assert np.all(values > 0.0)
        mapped = diag.canonical_to_effective(
            readout
        )
        assert mapped["latent"]["v"] >= 0.0


def test_joint_encode_decode_preserves_fixed_R():
    baseline_detector = {
        "R": 0.022,
        "alpha": 180.0,
        "beta": 0.1,
        "C_tes": 3.0e-13,
        "L": 1.2e-8,
        "T_bath": 0.216,
    }
    reference_readout = {
        "pole_Hz": 12300.0,
        "pole_Q": 0.79,
        "c2": 42.0,
        "c4": 34.0,
    }
    detector_names = (
        "alpha",
        "beta",
        "C_tes",
        "L",
        "T_bath",
    )
    readout_names = (
        "pole_Hz",
        "c2",
    )
    vector = diag.encode_joint(
        detector_names,
        baseline_detector,
        readout_names,
        reference_readout,
    )
    detector, readout = diag.decode_joint(
        vector,
        detector_names,
        baseline_detector,
        readout_names,
        reference_readout,
    )
    assert detector["R"] == pytest.approx(
        baseline_detector["R"]
    )
    for name in detector_names:
        assert detector[name] == pytest.approx(
            baseline_detector[name]
        )
    for name in diag.READOUT_PARAMETER_NAMES:
        assert readout[name] == pytest.approx(
            reference_readout[name]
        )


def test_fixed_reference_family_has_no_readout_coordinates():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    family = next(
        row
        for row in config[
            "residual_readout_families"
        ]
        if row["name"] == "fixed_reference_readout"
    )
    assert family["readout_parameters"] == []


def test_main_ladder_is_nested():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    by_name = {
        row["name"]: set(
            row["readout_parameters"]
        )
        for row in config[
            "residual_readout_families"
        ]
    }
    previous = set()
    for name in diag.MAIN_NESTED_LADDER:
        current = by_name[name]
        assert previous <= current
        previous = current

    assert by_name[
        "fixed_reference_readout"
    ] == set()
    assert by_name[
        "pole_frequency_only"
    ] == {"pole_Hz"}
    assert by_name[
        "pole_section"
    ] == {"pole_Hz", "pole_Q"}
    assert by_name[
        "pole_section_plus_c2"
    ] == {"pole_Hz", "pole_Q", "c2"}
    assert by_name[
        "full_order2"
    ] == {
        "pole_Hz",
        "pole_Q",
        "c2",
        "c4",
    }


def test_branch_families_probe_c2_without_forcing_full_pole_section():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    by_name = {
        row["name"]: set(
            row["readout_parameters"]
        )
        for row in config[
            "residual_readout_families"
        ]
    }
    assert by_name["c2_only"] == {"c2"}
    assert by_name[
        "pole_frequency_plus_c2"
    ] == {"pole_Hz", "c2"}


def test_tracked_baseline_snapshot_parses():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    path = diag.resolve_config_path(
        config["residual_baseline_snapshot"],
        diag.DEFAULT_CONFIG,
    )
    snapshot = diag.load_baseline_snapshot(
        path
    )
    assert snapshot["R_ratio_to_inherited"] == pytest.approx(
        1.25
    )
    assert snapshot["R_TES_Ohm"] > 0.0
    assert snapshot["baseline_shape_score"] > 0.0
    assert snapshot["repeat_local_shape_score"] > 0.0
    assert snapshot["white_asd_A_rtHz"] > 0.0
    assert snapshot["reference_readout"]["c2"] > 0.0
    assert snapshot["local_readout"]["c2"] > 0.0


def test_readout_bounds_include_reference_and_local_values():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    path = diag.resolve_config_path(
        config["residual_baseline_snapshot"],
        diag.DEFAULT_CONFIG,
    )
    snapshot = diag.load_baseline_snapshot(
        path
    )
    for name in diag.READOUT_PARAMETER_NAMES:
        lower, upper = [
            float(value)
            for value in config[
                "readout_bounds"
            ][name]
        ]
        for source in (
            snapshot["reference_readout"],
            snapshot["local_readout"],
        ):
            assert lower <= source[name] <= upper


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


def test_detector_nuisance_never_contains_R():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    names = config[
        "detector_nuisance"
    ]["parameters"]
    assert "R" not in names
    assert names == [
        "alpha",
        "beta",
        "C_tes",
        "L",
        "T_bath",
    ]


def test_detector_bounds_are_git_tracked():
    config = json.loads(
        diag.DEFAULT_CONFIG.read_text(
            encoding="utf-8"
        )
    )
    bounds = config["detector_bounds"]
    assert bounds["alpha"]["minimum_fraction_of_baseline"] == pytest.approx(
        0.05
    )
    assert bounds["alpha"]["maximum"] == pytest.approx(200.0)
    assert bounds["beta"]["min"] == pytest.approx(0.0)
    assert bounds["beta"]["max"] == pytest.approx(12.0)
    assert bounds["T_bath"]["half_width_K"] == pytest.approx(0.002)
