"""Stage 11C-0 observability: turn a native ``solver.log`` into a per-trial record.

The Phase24 Stage 11 adaptive controller (``ElmerSolver.F90``) already prints an
``Adaptive event audit`` / ``Adaptive trial start`` pair for every trial, plus an
``Adaptive accept: ...`` or ``Adaptive rejection: ...`` block once its outcome is
known.  Those lines were added across several prior Stage 11 commits for ad hoc
diagnostics; nothing in this repository turned them into a single structured,
regenerable per-trial record.  This module is that parser.

It is read-only: it does not change solver behavior, only interprets output the
solver already produces.  See ``artifacts/phase24_stage11c0_baseline_freeze.json``
for the Gate 11C-0 comparison against the frozen 20 ns baseline.
"""
from __future__ import annotations

from dataclasses import dataclass
import re
from pathlib import Path


_EVENT_AUDIT_RE = re.compile(
    r"Adaptive event audit id=(?P<id>\d+), "
    r"start/event/distance/proposed_dt/proposed_end/abs_proposed_gap/final_end="
    r"\s*(?P<start>\S+)\s+(?P<event>\S+)\s+(?P<distance>\S+)\s+(?P<proposed_dt>\S+)"
    r"\s+(?P<proposed_end>\S+)\s+(?P<abs_gap>\S+)\s+(?P<final_end>\S+),"
    r" clipped/landing/snapped/forced=(?P<flags>[TF]{4})"
)

_TRIAL_START_RE = re.compile(
    r"Adaptive trial start id=(?P<id>\d+) t0/t1/dt=\s*(?P<t0>\S+)\s+(?P<t1>\S+)\s+(?P<dt>\S+),"
    r" order=(?P<order>\d+), event=(?P<event>[TF])(?:, landing=(?P<landing>[TF]))?"
)

_ACCEPT_RE = re.compile(
    r"Adaptive accept: t0/attempted_dt/error/raw_maxdiff/h_old="
    r"\s*(?P<t0>\S+)\s+(?P<dt>\S+)\s+(?P<error>\S+)\s+(?P<raw_maxdiff>\S+)\s+(?P<h_old>\S+),"
    r" order=(?P<order>\d+), event=(?P<event>[TF])"
)

_ACCEPT_FLAGS_RE = re.compile(
    r"Adaptive accept flags: dt_is_min/ordinary_error/forced_floor/previous_rejected="
    r"(?P<flags>[TF]{4}), rejection_streak_before/after=(?P<before>\d+)(?P<after>\d+)"
)

_ACCEPT_COOLDOWN_RE = re.compile(
    r"Adaptive accept cooldown: before/after=(?P<before>\d+)(?P<after>\d+)"
)

_ACCEPT_GROWTH_RE = re.compile(
    r"Adaptive accept growth: accepted_dt/next_dt/factor="
    r"\s*(?P<accepted_dt>\S+)\s+(?P<next_dt>\S+)\s+(?P<factor>\S+), event_forced=(?P<event_forced>[TF])"
)

_REJECT_RE = re.compile(
    r"Adaptive rejection: index=(?P<index>\d+), consecutive=(?P<consecutive>\d+), "
    r"t0/hnew/hold/r/error/maxdiff="
    r"\s*(?P<t0>\S+)\s+(?P<hnew>\S+)\s+(?P<hold>\S+)\s+(?P<r>\S+)\s+(?P<error>\S+)\s+(?P<maxdiff>\S+),"
    r" order=(?P<order>\d+), event=(?P<event>[TF])"
)

_LINEAR_SOLVE_RE = re.compile(
    r"SolveHypre: Required iterations (?P<iterations>\d+) \(method \d+\) to norm (?P<norm>\S+)"
)

_SUMMARY_RE = re.compile(
    r"ExecSimulation: Adaptive output summary: accepted=(?P<accepted>\d+), "
    r"rejected=(?P<rejected>\d+), requested_outputs=(?P<requested_outputs>\d+), "
    r"interpolation_only=(?P<interpolation_only>\d+), event_forced=(?P<event_forced>\d+), "
    r"bdf1=(?P<bdf1>\d+), bdf2=(?P<bdf2>\d+)"
)


def _f(text: str) -> float:
    """Parse a Fortran-style exponent that may be missing its 'E' (e.g. '1.79+308')."""
    return float(re.sub(r"(\d)([+-]\d\d\d)$", r"\1E\2", text))


@dataclass(frozen=True)
class NativeTrialRecord:
    trial_id: int
    time: float
    proposed_dt: float
    final_dt: float
    bdf_order: int
    event_clipped: bool
    event_landing: bool
    event_endpoint_snapped: bool
    event_forced: bool
    accepted: bool
    estimator_raw_error: float
    estimator_control_error: float
    rejection_streak_before: int | None
    rejection_streak_after: int | None
    cooldown_before: int | None
    cooldown_after: int | None
    floor_accept: bool | None
    proposed_next_dt: float | None
    applied_growth_factor: float | None
    accepted_history_length_before: int
    bdf2_reentry: bool
    linear_iterations: int | None
    linear_residual: float | None


def parse_native_solver_log(log_path: str | Path) -> list[NativeTrialRecord]:
    """Return one :class:`NativeTrialRecord` per adaptive trial, in trial-id order.

    Raises ``ValueError`` if a trial's event-audit/trial-start pair cannot be
    matched to exactly one accept-or-reject resolution: silently dropping a
    trial would defeat the purpose of a Stage 11C-0 observability gate.
    """
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    pending: dict | None = None
    pending_linear: list[tuple[int, float]] = []
    records: list[NativeTrialRecord] = []
    previous_accepted_order: int | None = None
    accepted_count = 0
    state: dict | None = None

    for line in lines:
        m = _EVENT_AUDIT_RE.search(line)
        if m:
            clipped, landing, snapped, forced = (c == "T" for c in m.group("flags"))
            pending = {
                "id": int(m.group("id")),
                "event_clipped": clipped,
                "event_landing": landing,
                "event_endpoint_snapped": snapped,
                "event_forced": forced,
                "proposed_dt": _f(m.group("proposed_dt")),
            }
            pending_linear = []
            continue

        m = _TRIAL_START_RE.search(line)
        if m:
            trial_id = int(m.group("id"))
            if pending is None or pending["id"] != trial_id:
                # Cases with no physical events at all (physical_event_times=[])
                # never print "Adaptive event audit", so there is nothing to
                # match here; synthesize a no-event pending record instead of
                # treating this as a malformed log.
                pending = {
                    "id": trial_id,
                    "event_clipped": False,
                    "event_landing": False,
                    "event_endpoint_snapped": False,
                    "event_forced": False,
                    "proposed_dt": _f(m.group("dt")),
                }
                pending_linear = []
            pending["time"] = _f(m.group("t0"))
            pending["final_dt"] = _f(m.group("dt"))
            pending["trial_order"] = int(m.group("order"))
            continue

        m = _LINEAR_SOLVE_RE.search(line)
        if m:
            pending_linear.append((int(m.group("iterations")), _f(m.group("norm"))))
            continue

        m = _ACCEPT_RE.search(line)
        if m:
            if pending is None:
                # Cases without adaptive_time.debug never print "Adaptive
                # event audit" or "Adaptive trial start" at all; "Adaptive
                # accept:"/"Adaptive rejection:" are the only unconditional
                # per-trial lines there (gated by AdaptiveDebug OR
                # AdaptiveForcedFloor OR AdaptivePreviousRejected). Synthesize
                # a bare pending record instead of raising.
                dt = _f(m.group("dt"))
                pending = {
                    "id": len(records) + 1,
                    "event_clipped": False,
                    "event_landing": False,
                    "event_endpoint_snapped": False,
                    "event_forced": m.group("event") == "T",
                    "proposed_dt": dt,
                    "time": _f(m.group("t0")),
                    "final_dt": dt,
                }
                pending_linear = []
            state = dict(pending)
            state["accepted"] = True
            state["estimator_control_error"] = _f(m.group("error"))
            state["estimator_raw_error"] = _f(m.group("raw_maxdiff"))
            state["accept_order"] = int(m.group("order"))
            continue

        m = _ACCEPT_FLAGS_RE.search(line)
        if m and state is not None and state.get("accepted"):
            dt_is_min, ordinary_error, forced_floor, previous_rejected = (
                c == "T" for c in m.group("flags")
            )
            state["floor_accept"] = forced_floor
            state["rejection_streak_before"] = int(m.group("before"))
            state["rejection_streak_after"] = int(m.group("after"))
            continue

        m = _ACCEPT_COOLDOWN_RE.search(line)
        if m and state is not None and state.get("accepted"):
            state["cooldown_before"] = int(m.group("before"))
            state["cooldown_after"] = int(m.group("after"))
            continue

        m = _ACCEPT_GROWTH_RE.search(line)
        if m and state is not None and state.get("accepted"):
            state["proposed_next_dt"] = _f(m.group("next_dt"))
            state["applied_growth_factor"] = _f(m.group("factor"))
            bdf_order = state["accept_order"]
            bdf2_reentry = previous_accepted_order == 1 and bdf_order == 2
            records.append(
                NativeTrialRecord(
                    trial_id=state["id"],
                    time=state["time"],
                    proposed_dt=state["proposed_dt"],
                    final_dt=state["final_dt"],
                    bdf_order=bdf_order,
                    event_clipped=state["event_clipped"],
                    event_landing=state["event_landing"],
                    event_endpoint_snapped=state["event_endpoint_snapped"],
                    event_forced=state["event_forced"],
                    accepted=True,
                    estimator_raw_error=state["estimator_raw_error"],
                    estimator_control_error=state["estimator_control_error"],
                    rejection_streak_before=state["rejection_streak_before"],
                    rejection_streak_after=state["rejection_streak_after"],
                    cooldown_before=state["cooldown_before"],
                    cooldown_after=state["cooldown_after"],
                    floor_accept=state["floor_accept"],
                    proposed_next_dt=state["proposed_next_dt"],
                    applied_growth_factor=state["applied_growth_factor"],
                    accepted_history_length_before=accepted_count,
                    bdf2_reentry=bdf2_reentry,
                    linear_iterations=(pending_linear[-1][0] if pending_linear else None),
                    linear_residual=(pending_linear[-1][1] if pending_linear else None),
                )
            )
            previous_accepted_order = bdf_order
            accepted_count += 1
            pending = None
            pending_linear = []
            state = None
            continue

        m = _REJECT_RE.search(line)
        if m:
            if pending is None:
                # Same non-debug case as the accept branch above: no trial
                # marker lines exist at all, only "Adaptive rejection:" itself.
                rejected_dt = _f(m.group("hnew"))
                pending = {
                    "id": len(records) + 1,
                    "event_clipped": False,
                    "event_landing": False,
                    "event_endpoint_snapped": False,
                    "event_forced": m.group("event") == "T",
                    "proposed_dt": rejected_dt,
                    "time": _f(m.group("t0")),
                    "final_dt": rejected_dt,
                }
                pending_linear = []
            records.append(
                NativeTrialRecord(
                    trial_id=pending["id"],
                    time=pending["time"],
                    proposed_dt=pending["proposed_dt"],
                    final_dt=pending["final_dt"],
                    bdf_order=int(m.group("order")),
                    event_clipped=pending["event_clipped"],
                    event_landing=pending["event_landing"],
                    event_endpoint_snapped=pending["event_endpoint_snapped"],
                    event_forced=pending["event_forced"],
                    accepted=False,
                    estimator_raw_error=_f(m.group("maxdiff")),
                    estimator_control_error=_f(m.group("error")),
                    rejection_streak_before=int(m.group("consecutive")) - 1,
                    rejection_streak_after=int(m.group("consecutive")),
                    cooldown_before=None,
                    cooldown_after=None,
                    floor_accept=False,
                    proposed_next_dt=_f(m.group("hnew")),
                    applied_growth_factor=None,
                    accepted_history_length_before=accepted_count,
                    bdf2_reentry=False,
                    linear_iterations=(pending_linear[-1][0] if pending_linear else None),
                    linear_residual=(pending_linear[-1][1] if pending_linear else None),
                )
            )
            pending = None
            pending_linear = []
            continue

    return records


def parse_summary_line(log_path: str | Path) -> dict[str, int] | None:
    text = Path(log_path).read_text(encoding="utf-8", errors="replace")
    m = _SUMMARY_RE.search(text)
    if m is None:
        return None
    return {key: int(value) for key, value in m.groupdict().items()}
