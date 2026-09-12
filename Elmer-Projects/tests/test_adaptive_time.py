from __future__ import annotations

import pytest

from scripts.support.adaptive_time import (
    AdaptiveConfig,
    AdaptiveController,
    OutputSchedule,
    add_output_counters,
    bdf1_coefficients,
    dense_linear,
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
    dt, forced = controller.propose(0.0, 1.0, next_event=0.25)
    assert (dt, forced) == pytest.approx((0.1, False))
    dt, forced = controller.propose(0.2, 1.0, next_event=0.25)
    assert (dt, forced) == pytest.approx((0.05, True))


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
