from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import line_robust_continuum_fit_diagnostic as diag


def _repair_config():
    return {
        "detect_min_Hz": 40000,
        "detect_max_Hz": 200000,
        "baseline_width_Hz": 2000,
        "minimum_excess_dB": 3.0,
        "minimum_prominence_dB": 2.0,
        "minimum_spacing_Hz": 250,
        "max_detected_peaks": 64,
        "max_line_width_Hz": 500,
        "minimum_repair_half_width_Hz": 75,
        "width_expansion_factor": 1.5,
    }


def test_repair_narrow_line_preserves_broad_hump():
    frequency = np.arange(0.0, 200005.0, 5.0)
    asd = np.ones_like(frequency)

    # Broad smooth structure that should remain part of the continuum.
    asd *= 1.0 + 2.0 * np.exp(
        -0.5 * ((frequency - 115000.0) / 2500.0) ** 2
    )

    # Native-bin-like narrow line.
    index = int(np.argmin(np.abs(frequency - 100000.0)))
    asd[index] *= 10.0

    result = diag.repair_narrow_lines(
        frequency,
        asd,
        _repair_config(),
    )
    repaired = result["continuum_normalized"]
    raw = result["raw_normalized"]

    assert any(
        abs(row["frequency_Hz"] - 100000.0) <= 5.0
        for row in result["accepted_lines"]
    )
    assert repaired[index] < raw[index] / 3.0

    hump = int(np.argmin(np.abs(frequency - 115000.0)))
    assert np.isclose(
        repaired[hump],
        raw[hump],
        rtol=1e-12,
        atol=0.0,
    )


def test_choose_recommended_keeps_simpler_family_without_material_gain():
    pole = {
        "name": "pole_section_plus_c2",
        "shape_score": 1.0,
        "residual_metrics": {"rms_residual_dB": 0.50},
        "readout_boundary_hits": {},
    }
    full = {
        "name": "full_order2",
        "shape_score": 0.90,
        "residual_metrics": {"rms_residual_dB": 0.47},
        "readout_boundary_hits": {
            "c4": {
                "at_lower": False,
                "at_upper": False,
            }
        },
    }
    screen = {
        "full_order2_material_score_ratio": 0.8,
        "full_order2_min_rms_improvement_dB": 0.05,
        "reject_if_added_c4_at_bound": True,
    }
    result = diag.choose_recommended(pole, full, screen)
    assert result["selected_family"] == "pole_section_plus_c2"
    assert not result["full_order2_material_improvement"]


def test_choose_recommended_accepts_material_unbounded_c4_gain():
    pole = {
        "name": "pole_section_plus_c2",
        "shape_score": 1.0,
        "residual_metrics": {"rms_residual_dB": 0.50},
        "readout_boundary_hits": {},
    }
    full = {
        "name": "full_order2",
        "shape_score": 0.70,
        "residual_metrics": {"rms_residual_dB": 0.40},
        "readout_boundary_hits": {
            "c4": {
                "at_lower": False,
                "at_upper": False,
            }
        },
    }
    screen = {
        "full_order2_material_score_ratio": 0.8,
        "full_order2_min_rms_improvement_dB": 0.05,
        "reject_if_added_c4_at_bound": True,
    }
    result = diag.choose_recommended(pole, full, screen)
    assert result["selected_family"] == "full_order2"
    assert result["full_order2_material_improvement"]


def test_choose_recommended_rejects_c4_at_bound():
    pole = {
        "name": "pole_section_plus_c2",
        "shape_score": 1.0,
        "residual_metrics": {"rms_residual_dB": 0.50},
        "readout_boundary_hits": {},
    }
    full = {
        "name": "full_order2",
        "shape_score": 0.60,
        "residual_metrics": {"rms_residual_dB": 0.35},
        "readout_boundary_hits": {
            "c4": {
                "at_lower": False,
                "at_upper": True,
            }
        },
    }
    screen = {
        "full_order2_material_score_ratio": 0.8,
        "full_order2_min_rms_improvement_dB": 0.05,
        "reject_if_added_c4_at_bound": True,
    }
    result = diag.choose_recommended(pole, full, screen)
    assert result["selected_family"] == "pole_section_plus_c2"
    assert result["full_order2_c4_at_bound"]
