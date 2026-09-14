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


def _event_time_tolerance(event: float, endpoint: float) -> float:
    """Allow a couple of ulps when testing a decimal epsilon boundary."""
    return _TIME_EPS + 2.0 * math.ulp(max(1.0, abs(event), abs(endpoint)))


def is_event_landing(time: float, dt: float, event: float | None) -> bool:
    """Return whether the trial endpoint lands on a physical event.

    Landing is independent of whether the proposal had to be clipped.  In
    particular, a proposal whose endpoint is already within ``_TIME_EPS`` of
    the event is a landing even when no clipping is required.
    """
    if event is None or time >= event or dt <= 0.0:
        return False
    endpoint = time + dt
    return abs(endpoint - event) <= _event_time_tolerance(event, endpoint)


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
    event_landing_steps: int = 0
    bdf1_steps: int = 0
    bdf2_steps: int = 0
    matrix_refreshes: int = 0
    amg_setups: int = 0
    linear_solves: int = 0
    recoveries: int = 0


@dataclass(frozen=True)
class TrialLogEntry:
    """Stage 11C-0 observability record for one propose/accept-or-reject cycle.

    Populated by :meth:`AdaptiveController.accept`/``reject`` as a pure
    side-channel: nothing here feeds back into a control decision, so adding
    or reading fields cannot change accepted trajectories.
    """

    trial_id: int
    time: float | None
    proposed_dt: float
    final_dt: float
    previous_accepted_dt: float | None
    step_ratio: float | None
    bdf_order: int
    event_clipped: bool
    event_landing: bool
    force_bdf1_active: bool
    rejection_streak_before: int
    cooldown_before: int
    floor_accept: bool
    estimator_raw_error: float | None
    estimator_control_error: float | None
    accepted: bool
    proposed_next_dt: float | None
    applied_growth_factor: float | None
    accepted_history_length_before: int
    bdf2_reentry: bool
    bdf2_reentry_pending_before: bool
    linear_iterations: int | None
    linear_residual: float | None


@dataclass(frozen=True)
class TrialMetadata:
    """Immutable proposal facts consumed by one accepted or rejected trial."""

    proposed_dt: float
    final_dt: float
    event_clipped: bool
    event_landing: bool
    event_endpoint_snapped: bool
    event_forced: bool
    bdf_order: int

    def __iter__(self):
        """Preserve the historical ``dt, event_forced = propose(...)`` API."""
        yield self.final_dt
        yield self.event_forced

    def __getitem__(self, index: int):
        return (self.final_dt, self.event_forced)[index]


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
    last_event_clipped: bool = field(default=False, init=False)
    last_event_landing: bool = field(default=False, init=False)
    last_event_endpoint_snapped: bool = field(default=False, init=False)
    force_bdf1_steps_remaining: int = 0
    bdf2_reentry_pending: bool = False
    log: list[TrialLogEntry] = field(default_factory=list, repr=False, compare=False)
    _last_proposed_time: float | None = field(default=None, init=False, repr=False)
    _next_trial_id: int = field(default=1, init=False, repr=False)

    def __post_init__(self) -> None:
        if self.dt is None:
            self.dt = self.config.dt_initial
        self.dt = min(max(self.dt, self.config.dt_min), self.config.dt_max)

    def propose(self, time: float, end: float, *, next_event: float | None = None) -> TrialMetadata:
        """Return immutable trial facts; output times are never considered.

        ``event_forced`` retains its historical meaning: only event clipping
        sets it.  A numerical endpoint snap is reported separately.
        """
        if end < time:
            raise ValueError("integration end must not precede current time")
        dt = min(self.dt or self.config.dt_initial, end - time)
        forced = False
        self.last_event_clipped = False
        self.last_event_landing = False
        self.last_event_endpoint_snapped = False
        if self.previous_dt is not None:
            dt = min(dt, self.previous_dt * self.config.r_max)
            # A rejected retry owns the shrunken controller dt.  The
            # accepted-step ratio lower bound must not expand it before event
            # clipping, otherwise a 0.5 ns retry can replay the old 1 ns
            # event-landing step.
            if self.rejected_in_row == 0:
                dt = max(dt, min(end - time, self.previous_dt * self.config.r_min))
        proposed_dt = dt
        if next_event is not None and time < next_event < time + dt - _TIME_EPS:
            dt = next_event - time
            forced = True
        if dt <= 0.0:
            raise ValueError("proposed timestep is not positive")
        self.last_event_clipped = forced
        self.last_event_landing = is_event_landing(time, dt, next_event)
        snapped = False
        if self.last_event_landing and next_event is not None:
            snapped = (time + dt) != next_event
            dt = next_event - time
            self.last_event_endpoint_snapped = snapped
        trial_bdf_order = (
            1
            if self.force_bdf1_steps_remaining > 0
            or self.rejected_in_row > 0
            or self.previous_dt is None
            or self.last_event_landing
            else 2
        )
        self._last_proposed_time = time
        return TrialMetadata(
            proposed_dt=proposed_dt,
            final_dt=dt,
            event_clipped=forced,
            event_landing=self.last_event_landing,
            event_endpoint_snapped=snapped,
            event_forced=forced,
            bdf_order=trial_bdf_order,
        )

    def accept(
        self,
        dt: float | TrialMetadata,
        error: float,
        *,
        event_forced: bool = False,
        event_landing: bool | None = None,
        nonlinear_difficulty: float = 0.0,
        discontinuity: bool = False,
        raw_error: float | None = None,
        linear_iterations: int | None = None,
        linear_residual: float | None = None,
    ) -> None:
        trial = dt if isinstance(dt, TrialMetadata) else None
        proposed_dt = trial.proposed_dt if trial is not None else dt
        event_clipped = trial.event_clipped if trial is not None else False
        if trial is not None:
            if event_landing is None:
                event_landing = trial.event_landing
            event_forced = trial.event_forced
            dt = trial.final_dt
        if dt <= 0.0:
            raise ValueError("accepted timestep must be positive")
        # --- Stage 11C-0 observability: snapshot pre-mutation state only. ---
        previous_dt_before = self.previous_dt
        rejected_in_row_before = self.rejected_in_row
        cooldown_before = self.growth_cooldown_remaining
        history_length_before = self.counters.accepted_internal_steps
        force_bdf1_active_before = self.force_bdf1_steps_remaining > 0
        previous_bdf_order = self.bdf_order
        reentry_pending_before = self.bdf2_reentry_pending
        had_history = self.previous_dt is not None
        if event_landing is None:
            event_landing = False
        self.previous_dt = dt
        self.counters.accepted_internal_steps += 1
        self.counters.event_forced_steps += int(event_forced)
        self.counters.event_landing_steps += int(event_landing)
        held_bdf1 = self.force_bdf1_steps_remaining > 0
        is_event_boundary = event_forced or event_landing or discontinuity
        if trial is not None:
            # Trust the order the trial actually used -- propose()'s
            # trial_bdf_order already accounts for rejected_in_row (a
            # post-reject retry) as well as force_bdf1/no-history/landing.
            # Recomputing independently here (as below) previously missed
            # the rejected_in_row case: a non-event post-reject BDF1 retry
            # would be mislabeled order=2 the instant it was accepted, which
            # silently defeated Stage 11C-2 Candidate A's same-dt re-entry
            # check. Native avoids this by counting from AdaptiveStepOrder,
            # the value actually used to build the trial's stencil, instead
            # of a second guess made at accept time -- mirrored here.
            self.bdf_order = trial.bdf_order
        else:
            self.bdf_order = 1 if is_event_boundary or held_bdf1 or not had_history else 2
        if held_bdf1:
            self.force_bdf1_steps_remaining -= 1
        if is_event_boundary:
            self.force_bdf1_steps_remaining = 1
        forced_floor = dt <= self.config.dt_min * (1.0 + 1.0e-10) and error > 1.0
        self.rejected_in_row = 0
        # Stage 11C-2 Candidate A: same-dt BDF2 re-entry.  A BDF2 reject sets
        # bdf2_reentry_pending (see reject()).  While it is pending, the
        # accepted BDF1 recovery step must not also grow dt -- otherwise the
        # very next trial both flips back to BDF2 AND tries a larger,
        # unvalidated dt in the same step, which is the reject/floor/re-entry
        # cycle this candidate exists to break.  An event boundary always
        # takes priority and clears any stale pending state, since the old
        # recovery attempt no longer applies once the physics has moved on.
        if is_event_boundary:
            self.bdf2_reentry_pending = False
        elif reentry_pending_before and self.bdf_order == 2:
            # The same-dt BDF2 validation trial was just accepted: re-entry
            # succeeded.  Only now does ordinary growth policy resume.
            self.bdf2_reentry_pending = False
        if reentry_pending_before and self.bdf_order == 1 and not is_event_boundary:
            # Hold dt exactly where it is; the next trial must test BDF2 at
            # this same dt, not a grown one.
            self.dt = dt
        elif forced_floor:
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
        bdf2_reentry = self.bdf_order == 2 and previous_bdf_order == 1 and had_history
        self.log.append(
            TrialLogEntry(
                trial_id=self._next_trial_id,
                time=self._last_proposed_time,
                proposed_dt=proposed_dt,
                final_dt=dt,
                previous_accepted_dt=previous_dt_before,
                step_ratio=(dt / previous_dt_before) if previous_dt_before else None,
                bdf_order=self.bdf_order,
                event_clipped=event_clipped,
                event_landing=event_landing,
                force_bdf1_active=force_bdf1_active_before,
                rejection_streak_before=rejected_in_row_before,
                cooldown_before=cooldown_before,
                floor_accept=forced_floor,
                estimator_raw_error=raw_error if raw_error is not None else error,
                estimator_control_error=error,
                bdf2_reentry_pending_before=reentry_pending_before,
                accepted=True,
                proposed_next_dt=self.dt,
                applied_growth_factor=(self.dt / dt) if dt else None,
                accepted_history_length_before=history_length_before,
                bdf2_reentry=bdf2_reentry,
                linear_iterations=linear_iterations,
                linear_residual=linear_residual,
            )
        )
        self._next_trial_id += 1

    def reject(
        self,
        dt: float | TrialMetadata,
        *,
        raw_error: float | None = None,
        linear_iterations: int | None = None,
        linear_residual: float | None = None,
        error: float | None = None,
    ) -> float:
        trial = dt if isinstance(dt, TrialMetadata) else None
        proposed_dt = trial.proposed_dt if trial is not None else dt
        event_clipped = trial.event_clipped if trial is not None else False
        event_landing = trial.event_landing if trial is not None else False
        trial_bdf_order = trial.bdf_order if trial is not None else self.bdf_order
        if trial is not None:
            dt = trial.final_dt
        if dt <= 0.0:
            raise ValueError("rejected timestep must be positive")
        # --- Stage 11C-0 observability: snapshot pre-mutation state only. ---
        previous_dt_before = self.previous_dt
        rejected_in_row_before = self.rejected_in_row
        cooldown_before = self.growth_cooldown_remaining
        history_length_before = self.counters.accepted_internal_steps
        force_bdf1_active_before = self.force_bdf1_steps_remaining > 0
        reentry_pending_before = self.bdf2_reentry_pending
        self.counters.rejected_internal_steps += 1
        self.rejected_in_row += 1
        if self.rejected_in_row > self.config.max_rejected:
            raise RuntimeError("maximum rejected adaptive timesteps exceeded")
        self.dt = max(self.config.dt_min, dt * self.config.max_shrink)
        self.growth_cooldown_remaining = 1
        self.bdf_order = 1
        if trial_bdf_order == 2:
            # Stage 11C-2 Candidate A: a rejected BDF2 trial arms same-dt
            # re-entry validation (see accept()). Stays armed across repeated
            # BDF2 rejects during the same recovery attempt.
            self.bdf2_reentry_pending = True
        self.last_event_clipped = False
        self.last_event_landing = False
        self.last_event_endpoint_snapped = False
        self.log.append(
            TrialLogEntry(
                trial_id=self._next_trial_id,
                time=self._last_proposed_time,
                proposed_dt=proposed_dt,
                final_dt=dt,
                previous_accepted_dt=previous_dt_before,
                step_ratio=(dt / previous_dt_before) if previous_dt_before else None,
                bdf_order=trial_bdf_order,
                event_clipped=event_clipped,
                event_landing=event_landing,
                force_bdf1_active=force_bdf1_active_before,
                rejection_streak_before=rejected_in_row_before,
                cooldown_before=cooldown_before,
                floor_accept=False,
                estimator_raw_error=raw_error if raw_error is not None else error,
                estimator_control_error=error,
                bdf2_reentry_pending_before=reentry_pending_before,
                accepted=False,
                proposed_next_dt=self.dt,
                applied_growth_factor=(self.dt / dt) if dt else None,
                accepted_history_length_before=history_length_before,
                bdf2_reentry=False,
                linear_iterations=linear_iterations,
                linear_residual=linear_residual,
            )
        )
        self._next_trial_id += 1
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
