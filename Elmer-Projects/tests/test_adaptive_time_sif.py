from __future__ import annotations

from scripts.support.build_cases import _adaptive_time_lines


def test_sif_schedule_and_events_are_explicit_and_independent() -> None:
    lines = _adaptive_time_lines(
        {
            "adaptive_time": {
                "requested_output_times": {
                    "mode": "uniform",
                    "start": "0[s]",
                    "end": "1[ms]",
                    "count": 50,
                },
                "dt_initial": "1[us]",
                "dt_min": "1[ns]",
                "dt_max": "100[us]",
                "physical_event_times": ["250[us]"],
            }
        },
        {},
    )
    assert "Adaptive Output Decoupling = Logical True" in [line.strip() for line in lines]
    assert any(line.startswith("  Requested Output Times(50) = Real") for line in lines)
    assert any(line.startswith("  Physical Event Times(1) = Real") for line in lines)
    assert "Adaptive Initial Timestep = Real 1e-06" in [line.strip() for line in lines]
