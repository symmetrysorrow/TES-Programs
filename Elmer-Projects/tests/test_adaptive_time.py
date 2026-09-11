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
