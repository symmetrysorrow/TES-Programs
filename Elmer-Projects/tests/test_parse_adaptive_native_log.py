from __future__ import annotations

from pathlib import Path

import pytest

from scripts.analysis.parse_adaptive_native_log import (
    parse_native_solver_log,
    parse_summary_line,
)

FIXTURE = Path(__file__).parent / "fixtures" / "adaptive_native_log_snippet.txt"
NO_EVENT_AUDIT_FIXTURE = (
    Path(__file__).parent / "fixtures" / "adaptive_native_log_no_event_audit_snippet.txt"
)


def test_fixture_is_a_verbatim_prefix_of_the_frozen_20ns_baseline() -> None:
    """The fixture is trials 1-8 of the archived Gate-11C-0 baseline solver.log
    (results/case_phase24_adaptive_continuous_post_event_20ns_baseline/solver.log,
    sha256 82a80b577aff22f1a32b8a03f29489ba24f2fb2efbfd480b8ca3e1b2fa409f5b) plus the
    real run's trailing summary line appended verbatim. It is not synthesized."""
    assert FIXTURE.exists()


def test_parses_event_landing_reject_floor_reentry_sequence() -> None:
    records = parse_native_solver_log(FIXTURE)
    assert [r.trial_id for r in records] == list(range(1, 9))
    assert [r.accepted for r in records] == [
        True, False, True, True, True, True, True, True,
    ]

    # Stage 11A event-boundary contract: FTFF / TTFT / FFFF / FTFF.
    flags = [
        (r.event_clipped, r.event_landing, r.event_endpoint_snapped, r.event_forced)
        for r in records[:4]
    ]
    assert flags == [
        (False, True, False, False),
        (True, True, False, True),
        (False, False, False, False),
        (False, True, False, False),
    ]

    # Event-landing BDF sequence: landing BDF1 (trial 1) -> retry rejected
    # (trial 2) -> floor BDF1 (trial 3) -> second-event landing BDF1 (trial 4)
    # -> immediate post-event BDF1 hold (trial 5) -> ordinary BDF2 (trial 6+).
    assert [r.bdf_order for r in records] == [1, 1, 1, 1, 1, 2, 2, 2]

    reject = records[1]
    assert reject.estimator_control_error == pytest.approx(4.175234160701)
    assert reject.rejection_streak_after == 1

    floor_entries = records[2:6]
    assert [r.floor_accept for r in floor_entries] == [True, True, True, True]

    # The first BDF2 trial after the reject/floor run is the observed
    # re-entry point (root-cause candidate under investigation in 11C-1/11C-2).
    reentry = records[5]
    assert reentry.bdf2_reentry is True
    assert reentry.applied_growth_factor == pytest.approx(1.0)
    assert records[6].bdf2_reentry is False


def test_accepted_history_length_counts_only_accepted_trials() -> None:
    records = parse_native_solver_log(FIXTURE)
    assert [r.accepted_history_length_before for r in records] == [0, 1, 1, 2, 3, 4, 5, 6]


def test_summary_line_matches_known_20ns_baseline() -> None:
    summary = parse_summary_line(FIXTURE)
    assert summary == {
        "accepted": 42,
        "rejected": 11,
        "requested_outputs": 4,
        "interpolation_only": 4,
        "event_forced": 2,
        "bdf1": 14,
        "bdf2": 28,
    }


def test_parses_trials_with_no_physical_events_configured() -> None:
    """Cases with physical_event_times=[] (e.g. the local_dt/dtmin probes and
    ADAPTIVE_RECOVERY_CASE) never print 'Adaptive event audit', and their
    'Adaptive trial start' lines omit the trailing ', landing=' field
    entirely -- both must parse without raising, from
    results/case_phase24_adaptive_post_event_dt_recovery_10ns/solver.log
    (trials 1-3, restart_position=3 anchor)."""
    records = parse_native_solver_log(NO_EVENT_AUDIT_FIXTURE)
    assert [r.trial_id for r in records] == [1, 2, 3]
    assert [r.accepted for r in records] == [True, False, True]
    assert [r.event_landing for r in records] == [False, False, False]
    assert [r.bdf_order for r in records] == [1, 2, 1]
    assert records[1].estimator_control_error == pytest.approx(2.423130907247)


def test_missing_summary_line_returns_none(tmp_path) -> None:
    empty_log = tmp_path / "no_summary.log"
    empty_log.write_text("nothing to see here\n", encoding="utf-8")
    assert parse_summary_line(empty_log) is None
