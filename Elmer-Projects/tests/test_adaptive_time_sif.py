from __future__ import annotations

import pytest

from scripts.support.build_cases import (
    _adaptive_time_lines,
    _validate_adaptive_time_consistency,
)


def test_sif_schedule_and_events_are_explicit_and_independent() -> None:
    lines = _adaptive_time_lines(
        {
            "adaptive_time": {
                "start": "0[s]",
                "end": "1[ms]",
                "requested_output_times": {
                    "mode": "uniform",
                    "start": "0[s]",
                    "end": "1[ms]",
                    "count": 50,
                },
                "dt_initial": "1[us]",
                "dt_min": "1[ns]",
                "dt_max": "100[us]",
                "physical_event_times": ["250[us]", "2[ms]"],
            }
        },
        {},
    )
    assert "Adaptive Output Decoupling = Logical True" in [line.strip() for line in lines]
    assert any(line.startswith("  Requested Output Times(50) = Real") for line in lines)
    assert any(line.startswith("  Physical Event Times(1) = Real") for line in lines)
    assert "Adaptive Initial Timestep = Real 1e-06" in [line.strip() for line in lines]


def test_adaptive_events_are_clipped_to_integration_interval() -> None:
    lines = _adaptive_time_lines(
        {
            "adaptive_time": {
                "start": "1[ms]",
                "end": "2[ms]",
                "requested_output_times": ["1[ms]", "2[ms]"],
                "physical_event_times": ["0.5[ms]", "1.5[ms]", "2.5[ms]"],
            }
        },
        {},
    )
    assert any(line.startswith("  Physical Event Times(1) = Real") for line in lines)
    assert "Physical Event Times(1) = Real 0.0015" in [line.strip() for line in lines]


def test_adaptive_case_interval_must_match_stages_and_restart() -> None:
    spec = {
        "adaptive_time": {"start": "10[ms]", "end": "12[ms]"},
        "restart_time": "10[ms]",
        "timesteps": [["1[ms]", 2]],
    }
    _validate_adaptive_time_consistency("ok", spec, {})

    with pytest.raises(ValueError, match="restart_time"):
        _validate_adaptive_time_consistency(
            "bad_restart",
            {**spec, "restart_time": "9[ms]"},
            {},
        )
    with pytest.raises(ValueError, match="timestep stages"):
        _validate_adaptive_time_consistency(
            "bad_interval",
            {**spec, "timesteps": [["1[ms]", 1]]},
            {},
        )
