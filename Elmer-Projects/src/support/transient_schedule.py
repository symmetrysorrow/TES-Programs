"""Time grids used by project TOML transient windows."""

from __future__ import annotations

import math


def _ratio_limited(end_s: float, *, cap_scale: float = 1.0, coarse_after_s: float = 0.0) -> list[list[object]]:
    caps = [(100e-6, 2.5e-6), (200e-6, 5e-6), (500e-6, 10e-6), (1e-3, 25e-6),
            (2e-3, 50e-6), (5e-3, 100e-6), (15e-3, 250e-6), (35e-3, 500e-6),
            (float("inf"), 1e-3)]
    steps: list[float] = []
    t, dt = 0.0, 1e-9
    while t < end_s - 1e-15:
        cap = next(value for limit, value in caps if t < limit)
        if t >= coarse_after_s:
            cap *= cap_scale
        dt = min(dt, cap, end_s - t)
        steps.append(dt)
        t += dt
        dt = min(2.0 * dt, cap)
    schedule: list[list[object]] = [["18[us]", 1], ["1[us]", 2]]
    for step in steps:
        token = f"{step:.12e}[s]"
        if schedule[-1][0] == token:
            schedule[-1][1] = int(schedule[-1][1]) + 1
        else:
            schedule.append([token, 1])
    return schedule


def _fixed_post_pulse(end_s: float, max_step_s: float) -> list[list[object]]:
    pulse_steps: list[float] = []
    t, dt = 0.0, 1e-9
    while dt <= max_step_s and t < end_s - 1e-15:
        dt = min(dt, end_s - t)
        pulse_steps.append(dt)
        t += dt
        dt *= 2.0
    if not pulse_steps or pulse_steps[-1] < 2e-9:
        raise ValueError("fixed post-pulse schedule must resolve the pulse ramp")
    remaining = end_s - t
    count = max(1, math.ceil(remaining / max_step_s))
    step = remaining / count
    schedule: list[list[object]] = [["18[us]", 1], ["1[us]", 2]]
    for value in pulse_steps:
        token = f"{value:.12e}[s]"
        if schedule[-1][0] == token:
            schedule[-1][1] = int(schedule[-1][1]) + 1
        else:
            schedule.append([token, 1])
    schedule.append([f"{step:.12e}[s]", count])
    return schedule


def _fixed_segments(segments: list[tuple[float, float]]) -> list[list[object]]:
    pulse_steps: list[float] = []
    t, dt = 0.0, 1e-9
    while dt <= segments[0][1] and t < segments[0][0] - 1e-15:
        dt = min(dt, segments[0][0] - t)
        pulse_steps.append(dt)
        t += dt
        dt *= 2.0
    schedule: list[list[object]] = [["18[us]", 1], ["1[us]", 2]]
    for value in pulse_steps:
        token = f"{value:.12e}[s]"
        if schedule[-1][0] == token:
            schedule[-1][1] = int(schedule[-1][1]) + 1
        else:
            schedule.append([token, 1])
    for end_s, max_step_s in segments:
        remaining = end_s - t
        if remaining <= 0:
            continue
        count = max(1, math.ceil(remaining / max_step_s))
        token = f"{remaining / count:.12e}[s]"
        if schedule[-1][0] == token:
            schedule[-1][1] = int(schedule[-1][1]) + count
        else:
            schedule.append([token, count])
        t = end_s
    return schedule


def project_schedule(window: str) -> list[list[object]]:
    """Return the configured project window's validated timestep schedule."""
    if window == "75ms_bdf2h15":
        return _ratio_limited(0.055, cap_scale=1.5, coarse_after_s=500e-6)
    if window == "1ms_bdf2h15":
        return _ratio_limited(1e-3, cap_scale=1.5, coarse_after_s=500e-6)
    if window == "1ms_bdf2fixed2p5":
        return _fixed_post_pulse(1e-3, 2.5e-6)
    if window == "1ms_bdf2fixed1p25":
        return _fixed_post_pulse(1e-3, 1.25e-6)
    if window == "1ms_bdf2fixed2p5_5":
        return _fixed_segments([(500e-6, 2.5e-6), (1e-3, 5e-6)])
    if window == "singlepixel_hybrid_1p48ms":
        return [["18[us]", 1], ["1[us]", 2], ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9],
                ["1[us]", 9], ["0.625[us]", 144], ["5.26134395146352e-05[s]", 1],
                ["1.37771551611708e-04[s]", 1], ["1.89614008873657e-04[s]", 1],
                ["2.5e-04[s]", 1], ["2.5e-04[s]", 1]]
    raise ValueError(f"unknown project time window: {window}")
