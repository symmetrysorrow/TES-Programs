from __future__ import annotations

import pytest

from scripts.support.adaptive_time import (
    AdaptiveConfig,
    AdaptiveController,
    OutputSchedule,
    add_output_counters,
    bdf1_coefficients,
    dense_linear,
    is_event_landing,
    TrialLogEntry,
    TrialMetadata,
    variable_bdf2_coefficients,
    weighted_error,
)


def test_variable_bdf2_reduces_to_constant_step_formula() -> None:
    a0, a1, a2 = variable_bdf2_coefficients(2.0, 2.0)
    assert (a0, a1, a2) == pytest.approx((3.0 / 4.0, -1.0, 1.0 / 4.0))


def test_variable_bdf2_uses_new_step_ratio() -> None:
    h_new, h_old = 0.25, 1.0
    a0, a1, a2 = variable_bdf2_coefficients(h_new, h_old)
    r = h_new / h_old
    assert (a0, a1, a2) == pytest.approx(
        ((1.0 + 2.0 * r) / ((1.0 + r) * h_new), -(1.0 + r) / h_new, r * r / ((1.0 + r) * h_new))
    )
    assert bdf1_coefficients(h_new) == pytest.approx((4.0, -4.0))


def test_schedules_are_independent_of_internal_steps() -> None:
    sparse = OutputSchedule.uniform(0.0, 1.0, 5)
    dense = OutputSchedule.uniform(0.0, 1.0, 5000)
    assert sparse.times[0] == dense.times[0] == 0.0
    assert sparse.times[-1] == dense.times[-1] == 1.0
    assert len(dense.between(0.2, 0.4)) == 1000


def test_events_clip_but_outputs_do_not() -> None:
    config = AdaptiveConfig(0.1, 0.01, 0.5)
    controller = AdaptiveController(config)
    trial = controller.propose(0.0, 1.0, next_event=0.25)
    assert isinstance(trial, TrialMetadata)
    assert (trial.final_dt, trial.event_forced) == pytest.approx((0.1, False))
    trial = controller.propose(0.2, 1.0, next_event=0.25)
    assert (trial.final_dt, trial.event_forced) == pytest.approx((0.05, True))
    assert controller.last_event_clipped is True
    assert controller.last_event_landing is True


def test_exact_event_landing_without_clipping_resets_to_bdf1() -> None:
    t0 = 0.0200200005
    event = 0.020020001
    dt = 0.5e-9
    controller = AdaptiveController(AdaptiveConfig(dt, 0.5e-9, 1.0e-4))
    controller.dt = dt

    trial = controller.propose(t0, event, next_event=event)

    assert trial.final_dt == pytest.approx(dt)
    assert trial.event_forced is False
    assert controller.last_event_clipped is False
    assert controller.last_event_landing is True
    controller.accept(trial, 25.0)
    assert controller.bdf_order == 1
    assert trial.final_dt + t0 == pytest.approx(event)


@pytest.mark.parametrize("offset", (-1.1, -1.0, -0.5, 0.0, 0.5, 1.0, 1.1))
def test_event_epsilon_contract(offset: float) -> None:
    event = 1.0
    t0 = 0.0
    epsilon = 1.0e-12
    endpoint = event + offset * epsilon
    controller = AdaptiveController(AdaptiveConfig(1.0, 0.1, 2.0))
    trial = controller.propose(t0, 2.0, next_event=event)
    # Drive the final endpoint through the public controller using a fresh
    # controller whose proposal is the requested boundary case.
    controller.dt = endpoint
    trial = controller.propose(t0, 2.0, next_event=event)
    within = abs(offset) <= 1.0
    clipped = offset > 1.0
    assert trial.event_clipped is clipped
    assert trial.event_forced is clipped
    assert trial.event_landing is (within or clipped)
    assert trial.event_endpoint_snapped is (within and offset != 0.0)
    assert trial.final_dt <= event
    if within or clipped:
        assert trial.final_dt == pytest.approx(event)
    else:
        assert trial.event_endpoint_snapped is False


def test_propose_forces_bdf1_for_a_landing_trial_from_ordinary_history() -> None:
    """Regression for a Stage 11C-0 finding: propose()'s trial_bdf_order used to
    check only prior force/reject/no-history state, never whether THIS trial
    itself lands on the event -- so a landing trial reached from ordinary
    history (no active force-BDF1 or reject streak) was mislabeled order=2,
    unlike native's proactive 'IF (AdaptiveEventLanding) AdaptiveForceBDF1 =
    .TRUE.' before it picks AdaptiveStepOrder. accept()'s own recomputation
    already had this condition; propose() must match it."""
    controller = AdaptiveController(AdaptiveConfig(1.0e-9, 0.5e-9, 1.0e-4))
    controller.accept(1.0e-9, 0.05)  # ordinary history, no force/reject state active
    trial = controller.propose(1.0e-9, 5.0e-9, next_event=2.0e-9)
    assert trial.event_landing is True
    assert trial.bdf_order == 1


def test_event_landing_holds_bdf1_for_immediate_post_event_trial() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.5e-9, 0.5e-9, 1.0e-4))
    trial = controller.propose(0.5e-9, 1.0e-9, next_event=1.0e-9)
    controller.accept(trial, 0.1)
    assert controller.bdf_order == 1
    post = controller.propose(1.0e-9, 2.0e-9)
    assert post.bdf_order == 1
    controller.accept(post, 0.1)
    following = controller.propose(1.0e-9 + post.final_dt, 2.0e-9)
    assert following.bdf_order == 2


def test_event_landing_requires_endpoint_at_event() -> None:
    t0 = 0.0200200005
    event = 0.020020001
    assert is_event_landing(t0, 0.5e-9, event) is True
    assert is_event_landing(t0, 0.25e-9, event) is False


def test_step_ratio_bounds_and_event_restart_bdf1() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 1.0, r_min=0.5, r_max=1.5))
    controller.accept(0.1, 0.01)
    assert controller.propose(0.1, 1.0)[0] == pytest.approx(0.15)
    controller.accept(0.15, 0.01)
    assert controller.bdf_order == 2
    controller.accept(0.1, 0.01, event_forced=True)
    assert controller.bdf_order == 1


def test_rejection_does_not_advance_accepted_count() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 0.5))
    controller.reject(0.1)
    assert controller.counters.accepted_internal_steps == 0
    assert controller.counters.rejected_internal_steps == 1
    controller.accept(0.05, 0.01)
    assert controller.counters.accepted_internal_steps == 1
    assert controller.counters.bdf1_steps == 1


def test_rejected_event_retry_keeps_controller_shrink_below_ratio_floor() -> None:
    controller = AdaptiveController(
        AdaptiveConfig(1.0e-9, 0.5e-9, 1.0e-4, r_min=0.5, r_max=2.0)
    )
    controller.accept(6.8125e-6, 0.01)
    controller.dt = 1.0e-9
    controller.reject(1.0e-9)

    retry_dt, event_forced = controller.propose(
        0.02002, 0.020020001, next_event=0.020020001
    )

    assert retry_dt == pytest.approx(0.5e-9)
    assert event_forced is False


def test_rejection_budget_is_consecutive_not_cumulative() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 0.5, max_rejected=2))

    # Two reject/accept cycles must not consume the next rejection streak.
    controller.reject(0.1)
    controller.accept(0.05, 0.01)
    controller.reject(0.1)
    controller.accept(0.05, 0.01)

    controller.reject(0.1)
    controller.reject(0.05)
    with pytest.raises(RuntimeError, match="maximum rejected"):
        controller.reject(0.025)

    assert controller.counters.rejected_internal_steps == 5
    assert controller.rejected_in_row == 3


def test_reject_accept_does_not_immediately_use_max_growth() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 1.0))
    controller.reject(0.02)
    controller.accept(0.01, 0.2)
    assert controller.dt == pytest.approx(0.012)
    assert controller.growth_cooldown_remaining == 0


def test_reject_accept_accept_resumes_growth_after_cooldown() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 1.0))
    controller.reject(0.02)
    controller.accept(0.01, 0.2)
    controller.accept(0.012, 0.01)
    assert controller.dt == pytest.approx(0.018)


def test_normal_accept_keeps_max_growth_behavior() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 1.0))
    controller.accept(0.1, 0.01)
    assert controller.dt == pytest.approx(0.15)


def test_event_restart_preserves_reject_cooldown() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 1.0))
    controller.reject(0.02)
    controller.accept(0.01, 0.2, event_forced=True)
    assert controller.bdf_order == 1
    assert controller.dt == pytest.approx(0.012)


def test_dt_min_forced_accept_stays_at_floor_without_oscillation() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 1.0))
    controller.reject(0.015)
    controller.accept(0.01, 2.0)
    assert controller.dt == pytest.approx(0.01)
    controller.accept(0.01, 2.0)
    assert controller.dt == pytest.approx(0.01)


def test_dense_output_and_weighted_error() -> None:
    assert dense_linear((0.0, 10.0), (10.0, 20.0), 0.25) == pytest.approx((2.5, 12.5))
    assert weighted_error((1.0, 10.0), (1.001, 10.0), atol=1.0e-4, rtol=1.0e-3) == pytest.approx(
        1.0e-3 / (1.0e-4 + 1.0e-3 * 1.001)
    )


def test_interpolation_counter_is_explicit() -> None:
    counters = AdaptiveController(AdaptiveConfig(0.1, 0.01, 0.5)).counters
    add_output_counters(counters, OutputSchedule.uniform(0.0, 1.0, 11), 0.0, 0.5)
    assert counters.requested_outputs == 6
    assert counters.interpolation_only_outputs == 4


def test_adjacent_output_spans_do_not_double_count_shared_endpoint() -> None:
    counters = AdaptiveController(AdaptiveConfig(0.1, 0.01, 0.5)).counters
    schedule = OutputSchedule.uniform(0.0, 1.0, 11)
    add_output_counters(counters, schedule, 0.0, 0.5)
    add_output_counters(counters, schedule, 0.5, 1.0, include_start=False)
    assert counters.requested_outputs == 11
    assert counters.interpolation_only_outputs == 8


# --- Stage 11C-0: observability must not perturb the accepted trajectory. ---


def _drive_reject_floor_reentry_sequence(controller: AdaptiveController) -> None:
    """Replay the reject -> floor-accept -> BDF2-reentry shape seen in the
    Stage 11 continuous post-event trace (phase24_stage11_continuous_post_event_bdf2_validation.json)."""
    controller.reject(1.0e-9)
    controller.accept(0.5e-9, 1.75)  # forced-floor accept (dt==dt_min, error>1)
    controller.accept(0.5e-9, 25.7)  # forced-floor accept again; re-arms cooldown
    controller.accept(0.5e-9, 0.4)  # ordinary accept: cooldown-limited growth fires here


def test_observability_log_does_not_change_accepted_trajectory() -> None:
    quiet = AdaptiveController(AdaptiveConfig(1.0e-9, 0.5e-9, 1.0e-4))
    watched = AdaptiveController(AdaptiveConfig(1.0e-9, 0.5e-9, 1.0e-4))
    _drive_reject_floor_reentry_sequence(quiet)
    _drive_reject_floor_reentry_sequence(watched)
    # Reading .log after every call must not be observable in the numbers.
    assert len(watched.log) == 4
    assert quiet.dt == pytest.approx(watched.dt)
    assert quiet.bdf_order == watched.bdf_order
    assert quiet.growth_cooldown_remaining == watched.growth_cooldown_remaining
    assert quiet.rejected_in_row == watched.rejected_in_row


def test_log_records_reject_then_floor_accept_then_reentry_growth() -> None:
    controller = AdaptiveController(AdaptiveConfig(1.0e-9, 0.5e-9, 1.0e-4))
    controller.accept(1.0e-9, 0.05)  # ordinary history-establishing accept, BDF1
    controller.reject(1.0e-9)
    controller.accept(0.5e-9, 1.75)  # forced-floor accept after the reject retry
    controller.accept(0.5e-9, 0.4)  # first ordinary accept: grows dt via cooldown

    entries = controller.log
    assert len(entries) == 4
    reject_entry = entries[1]
    assert reject_entry.accepted is False
    assert reject_entry.rejection_streak_before == 0

    floor_entry = entries[2]
    assert floor_entry.floor_accept is True
    assert floor_entry.rejection_streak_before == 1
    assert floor_entry.cooldown_before == 1
    # KNOWN PARITY GAP (found while writing this Stage 11C-0 observability
    # test, not introduced by it): AdaptiveController.accept's own bdf_order
    # recomputation checks event_forced/event_landing/discontinuity/held_bdf1
    # but -- unlike propose()'s trial_bdf_order -- does not check
    # rejected_in_row. A non-event post-reject retry is therefore proposed
    # and solved at BDF1 (propose() forces it via rejected_in_row > 0), yet
    # this Python reference labels the accepted step order 2 here, unlike
    # native's ElmerSolver.F90 which increments its BDF1/BDF2 counters from
    # the AdaptiveStepOrder actually used to build the trial's stencil,
    # not from an independent post-hoc recomputation. Left unresolved by
    # deliberate Stage 11C-0 scope (observability only, no control-logic
    # changes); tracked for Stage 11C-2's BDF2 re-entry state-machine work.
    assert floor_entry.bdf_order == 2

    reentry_entry = entries[3]
    assert reentry_entry.bdf_order == 2
    assert reentry_entry.bdf2_reentry is False
    assert reentry_entry.floor_accept is False
    assert reentry_entry.applied_growth_factor == pytest.approx(1.2)
    assert reentry_entry.accepted_history_length_before == 2


def test_log_captures_estimator_and_linear_diagnostics_when_supplied() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.1, 0.01, 0.5))
    controller.accept(
        0.1,
        0.2,
        raw_error=0.6,
        linear_iterations=12,
        linear_residual=4.5e-7,
    )
    entry = controller.log[-1]
    assert entry.estimator_control_error == pytest.approx(0.2)
    assert entry.estimator_raw_error == pytest.approx(0.6)
    assert entry.linear_iterations == 12
    assert entry.linear_residual == pytest.approx(4.5e-7)


def test_log_entry_from_trial_metadata_records_event_flags() -> None:
    controller = AdaptiveController(AdaptiveConfig(0.5e-9, 0.5e-9, 1.0e-4))
    trial = controller.propose(0.5e-9, 1.0e-9, next_event=1.0e-9)
    controller.accept(trial, 0.1)
    entry = controller.log[-1]
    assert entry.event_landing is True
    assert entry.time == pytest.approx(0.5e-9)
    assert isinstance(entry, TrialLogEntry)
