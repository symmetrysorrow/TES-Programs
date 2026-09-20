from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_effective_Ctes_joint_c2_diagnostic as diag


def test_absolute_bounds_to_material_multipliers():
    low, high = diag.absolute_bounds_to_material_multipliers(
        5.0e-13,
        2.0e-11,
        1.0e-12,
    )
    assert low == pytest.approx(0.5)
    assert high == pytest.approx(20.0)


def test_absolute_bounds_reject_invalid():
    with pytest.raises(ValueError):
        diag.absolute_bounds_to_material_multipliers(
            2.0e-11,
            5.0e-13,
            1.0e-12,
        )


def test_log_position_endpoints_and_middle():
    assert diag.log_position(1.0, 1.0, 100.0) == pytest.approx(0.0)
    assert diag.log_position(100.0, 1.0, 100.0) == pytest.approx(1.0)
    assert diag.log_position(10.0, 1.0, 100.0) == pytest.approx(0.5)


def _result(*, material=True, c2_upper=False, ctes_boundary=False):
    hit = {
        "at_lower": bool(ctes_boundary),
        "at_upper": False,
    }
    return {
        "nested_test": {
            "shared_c2_material_improvement": bool(material),
        },
        "c2_zero": {
            "solution": {
                "shared_physical_boundary_hits": {
                    "C_tes": dict(hit),
                },
            },
        },
        "c2_free": {
            "solution": {
                "shared_c2_boundary_hit": {
                    "at_lower": False,
                    "at_upper": bool(c2_upper),
                },
                "shared_physical_boundary_hits": {
                    "C_tes": dict(hit),
                },
            },
        },
    }


def test_effective_classification_c2_not_material():
    result = diag.effective_classification(_result(material=False))
    assert result == (
        "broad_effective_Ctes_profile_removes_material_c2_need"
    )


def test_effective_classification_material_interior():
    result = diag.effective_classification(_result())
    assert result == (
        "shared_c2_material_after_broad_effective_Ctes_profile"
    )


def test_effective_classification_c2_upper_bound():
    result = diag.effective_classification(
        _result(c2_upper=True)
    )
    assert result == (
        "shared_c2_material_but_c2_upper_bound_limited_with_"
        "broad_effective_Ctes_profile"
    )


def test_effective_classification_ctes_boundary():
    result = diag.effective_classification(
        _result(ctes_boundary=True)
    )
    assert result == (
        "shared_c2_material_but_effective_Ctes_profile_"
        "boundary_active"
    )
