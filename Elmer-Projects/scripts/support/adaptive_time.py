"""Stage 11 time/output control primitives.

The native Elmer driver owns the FEM solve and its persistent linear-solver
objects.  This module owns the *policy* that sits above that solve: requested
sample times, physical events, variable-step BDF2 coefficients, weighted
error control, and counters.  Keeping these concerns separate prevents an
output schedule from accidentally becoming a solver schedule.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Iterable, Sequence


_TIME_EPS = 1.0e-12
_POST_REJECT_MAX_GROWTH = 1.2


def _finite_times(values: Iterable[float], *, name: str) -> tuple[float, ...]:
    result = tuple(float(value) for value in values)
    if any(not math.isfinite(value) for value in result):
        raise ValueError(f"{name} must contain only finite times")
    if any(b <= a for a, b in zip(result, result[1:])):
        raise ValueError(f"{name} must be strictly increasing")
    return result


@dataclass(frozen=True)
class OutputSchedule:
    """Absolute physical times at which observables are requested."""

    times: tuple[float, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "times", _finite_times(self.times, name="output times"))

    @classmethod
    def explicit(cls, times: Iterable[float]) -> "OutputSchedule":
        return cls(tuple(times))

    @classmethod
    def uniform(cls, start: float, end: float, count: int) -> "OutputSchedule":
        if count < 1:
            raise ValueError("uniform output count must be positive")
        if end < start:
            raise ValueError("output end must not precede output start")
        if end == start and count != 1:
            raise ValueError("a zero-length uniform schedule has exactly one sample")
        if count == 1:
            values = (float(start),)
        else:
            step = (end - start) / (count - 1)
            values = tuple(start + i * step for i in range(count))
            values = (*values[:-1], float(end))
        return cls(values)

    @classmethod
    def interval(cls, start: float, end: float, interval: float) -> "OutputSchedule":
        if interval <= 0.0 or not math.isfinite(interval):
            raise ValueError("output interval must be finite and positive")
        if end < start:
            raise ValueError("output end must not precede output start")
        values = []
        current = float(start)
        while current < end - _TIME_EPS * max(1.0, abs(end)):
            values.append(current)
            current += interval
        values.append(float(end))
        return cls(tuple(values))

    @classmethod
    def logarithmic(
        cls, start: float, end: float, count: int, *, origin: float = 0.0
    ) -> "OutputSchedule":
        """Return a log-spaced schedule for positive ``t-origin``."""
        if count < 2:
            raise ValueError("logarithmic output count must be at least two")
        a = start - origin
        b = end - origin
        if a <= 0.0 or b <= a:
            raise ValueError("logarithmic schedule requires origin < start < end")
        ratio = (b / a) ** (1.0 / (count - 1))
        values = tuple(origin + a * ratio**i for i in range(count))
        return cls((*values[:-1], float(end)))

    def between(
        self, start: float, end: float, *, include_start: bool = True
    ) -> tuple[float, ...]:
        """Return requested times in ``[start, end]`` including endpoints.

        Set ``include_start=False`` when this span follows another accepted
        span whose endpoint is the same physical time.  This keeps the
        endpoint available to the first span while preventing duplicate
        accounting across adjacent spans.
        """
        scale = max(1.0, abs(start), abs(end))
        lower = (
            (lambda t: t >= start - _TIME_EPS * scale)
            if include_start
            else (lambda t: t > start + _TIME_EPS * scale)
        )
        return tuple(
            t
            for t in self.times
            if lower(t) and t <= end + _TIME_EPS * scale
        )


@dataclass(frozen=True)
class AdaptiveConfig:
    """Conservative defaults suitable for TES thermal transients."""

    dt_initial: float
    dt_min: float
    dt_max: float
    relative_tolerance: float = 1.0e-4
    absolute_tolerance: float = 1.0e-8
    r_min: float = 0.5
    r_max: float = 2.0
    max_growth: float = 1.5
    max_shrink: float = 0.5
    max_rejected: int = 12

    def __post_init__(self) -> None:
        positive = {
            "dt_initial": self.dt_initial,
            "dt_min": self.dt_min,
            "dt_max": self.dt_max,
            "relative_tolerance": self.relative_tolerance,
            "absolute_tolerance": self.absolute_tolerance,
            "r_min": self.r_min,
            "r_max": self.r_max,
            "max_growth": self.max_growth,
            "max_shrink": self.max_shrink,
        }
        if any(not math.isfinite(value) or value <= 0.0 for value in positive.values()):
            raise ValueError("adaptive time configuration values must be finite and positive")
        if not self.r_min <= 1.0 <= self.r_max:
            raise ValueError("step-ratio bounds must contain one")
        if self.max_shrink > 1.0 or self.max_growth < 1.0:
            raise ValueError("max_shrink must be <= 1 and max_growth must be >= 1")
        if not self.dt_min <= self.dt_initial <= self.dt_max:
            raise ValueError("dt_min <= dt_initial <= dt_max is required")
        if self.max_rejected < 0:
            raise ValueError("max_rejected must be non-negative")


def bdf1_coefficients(h_new: float) -> tuple[float, float]:
    if h_new <= 0.0:
        raise ValueError("new timestep must be positive")
    return 1.0 / h_new, -1.0 / h_new


def variable_bdf2_coefficients(h_new: float, h_old: float) -> tuple[float, float, float]:
    """Return derivative coefficients for ``u[n+1], u[n], u[n-1]``.

    With ``r = h_new / h_old`` this is
    ``((1+2r)/(1+r), -(1+r), r²/(1+r)) / h_new``.
    """
    if h_new <= 0.0 or h_old <= 0.0:
        raise ValueError("BDF2 step sizes must be positive")
    r = h_new / h_old
    return (
        (1.0 + 2.0 * r) / ((1.0 + r) * h_new),
        -(1.0 + r) / h_new,
        r * r / ((1.0 + r) * h_new),
    )


def dense_linear(u0: Sequence[float], u1: Sequence[float], theta: float) -> tuple[float, ...]:
    """Interpolate an accepted state without invoking a FEM solve."""
    if len(u0) != len(u1):
        raise ValueError("interpolation states must have equal length")
    if not 0.0 <= theta <= 1.0:
        raise ValueError("interpolation fraction must lie in [0, 1]")
    return tuple((1.0 - theta) * a + theta * b for a, b in zip(u0, u1))


def weighted_error(
    coarse: Sequence[float], fine: Sequence[float], *, atol: float, rtol: float
) -> float:
    """Infinity-norm embedded error using absolute/relative tolerances."""
    if len(coarse) != len(fine):
        raise ValueError("error-estimator states must have equal length")
    if atol <= 0.0 or rtol <= 0.0:
        raise ValueError("error tolerances must be positive")
    return max(
        (abs(a - b) / (atol + rtol * max(abs(a), abs(b))))
        for a, b in zip(coarse, fine)
    ) if coarse else 0.0


@dataclass
class StepCounters:
    accepted_internal_steps: int = 0
    rejected_internal_steps: int = 0
    requested_outputs: int = 0
    interpolation_only_outputs: int = 0
    event_forced_steps: int = 0
    bdf1_steps: int = 0
    bdf2_steps: int = 0
    matrix_refreshes: int = 0
    amg_setups: int = 0
    linear_solves: int = 0
    recoveries: int = 0


@dataclass
class AdaptiveController:
    """Small policy object used by the native driver and unit tests."""

    config: AdaptiveConfig
    dt: float | None = None
    previous_dt: float | None = None
    bdf_order: int = 1
    rejected_in_row: int = 0
    growth_cooldown_remaining: int = 0
    counters: StepCounters = field(default_factory=StepCounters)

    def __post_init__(self) -> None:
        if self.dt is None:
            self.dt = self.config.dt_initial
        self.dt = min(max(self.dt, self.config.dt_min), self.config.dt_max)

    def propose(self, time: float, end: float, *, next_event: float | None = None) -> tuple[float, bool]:
        """Return ``(dt, event_forced)``; output times are never considered."""
        if end < time:
            raise ValueError("integration end must not precede current time")
        dt = min(self.dt or self.config.dt_initial, end - time)
        forced = False
        if self.previous_dt is not None:
            dt = min(dt, self.previous_dt * self.config.r_max)
            dt = max(dt, min(end - time, self.previous_dt * self.config.r_min))
        if next_event is not None and time < next_event < time + dt + _TIME_EPS:
            dt = next_event - time
            forced = True
        if dt <= 0.0:
            raise ValueError("proposed timestep is not positive")
        return dt, forced

    def accept(
        self,
        dt: float,
        error: float,
        *,
        event_forced: bool = False,
        nonlinear_difficulty: float = 0.0,
        discontinuity: bool = False,
    ) -> None:
        if dt <= 0.0:
            raise ValueError("accepted timestep must be positive")
        had_history = self.previous_dt is not None
        self.previous_dt = dt
        self.counters.accepted_internal_steps += 1
        self.counters.event_forced_steps += int(event_forced)
        self.bdf_order = 1 if event_forced or discontinuity or not had_history else 2
        forced_floor = dt <= self.config.dt_min * (1.0 + 1.0e-10) and error > 1.0
        self.rejected_in_row = 0
        if forced_floor:
            # A floor accept is a safety valve, not evidence that the error
            # target was met.  Hold the floor until an ordinary accept gives
            # the controller a reliable recovery point.
            self.dt = self.config.dt_min
            self.growth_cooldown_remaining = 1
        elif self.growth_cooldown_remaining > 0:
            # After reject -> accept, avoid immediately replaying the old
            # max-growth proposal.  One conservative recovery step is enough
            # to preserve normal behavior once the local transient settles.
            self.dt = min(self.config.dt_max, dt * min(self.config.max_growth, _POST_REJECT_MAX_GROWTH))
            self.growth_cooldown_remaining -= 1
        else:
            difficulty_factor = 0.75 if nonlinear_difficulty > 1.0 else 1.0
            if error <= 0.1:
                factor = self.config.max_growth
            elif error >= 1.0:
                factor = self.config.max_shrink
            else:
                factor = min(self.config.max_growth, max(self.config.max_shrink, 0.9 * error ** -0.5))
            self.dt = min(self.config.dt_max, max(self.config.dt_min, dt * factor * difficulty_factor))
        if self.bdf_order == 1:
            self.counters.bdf1_steps += 1
        else:
            self.counters.bdf2_steps += 1

    def reject(self, dt: float) -> float:
        if dt <= 0.0:
            raise ValueError("rejected timestep must be positive")
        self.counters.rejected_internal_steps += 1
        self.rejected_in_row += 1
        if self.rejected_in_row > self.config.max_rejected:
            raise RuntimeError("maximum rejected adaptive timesteps exceeded")
        self.dt = max(self.config.dt_min, dt * self.config.max_shrink)
        self.growth_cooldown_remaining = 1
        self.bdf_order = 1
        return self.dt


def add_output_counters(
    counters: StepCounters,
    schedule: OutputSchedule,
    start: float,
    end: float,
    *,
    include_start: bool = True,
) -> None:
    """Count outputs for one accepted span without double-counting endpoints."""
    outputs = schedule.between(start, end, include_start=include_start)
    counters.requested_outputs += len(outputs)
    if end > start:
        counters.interpolation_only_outputs += sum(
            start < t < end for t in outputs
        )
