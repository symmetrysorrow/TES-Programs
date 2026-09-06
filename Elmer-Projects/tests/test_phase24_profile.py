from __future__ import annotations

import json

import pytest

from scripts.analysis.phase24_profile import (
    ProfileError,
    classify,
    make_profile,
    parse_phase24_log,
    validate_gate,
)


def test_nested_wall_events_are_exclusive() -> None:
    profile = make_profile(
        {
            "events": [
                {"component": "element traversal", "start_s": 0.0, "end_s": 10.0},
                {"component": "local element calculation", "start_s": 2.0, "end_s": 5.0},
            ]
        },
        case="smoke",
        total_wall=12.0,
    )
    assert profile["components_seconds"]["element_traversal"] == pytest.approx(7.0)
    assert profile["components_seconds"]["local_element_calculation"] == pytest.approx(3.0)
    assert profile["unclassified_seconds"] == pytest.approx(2.0)
    assert profile["accounting_ok"] is True


def test_cpu_only_measurement_is_rejected() -> None:
    with pytest.raises(ProfileError, match="CPU-only"):
        make_profile([{"component": "hypre solve", "cpu_seconds": 1.0}], case="x", total_wall=1.0)


def test_gate_requires_all_correctness_observables() -> None:
    baseline = {"case": "base"}
    candidate = {
        "case": "candidate",
        "matrix_sparsity_identical": True,
        "matrix_values_identical": True,
        "rhs_identical": True,
        "temperature_field_identical": True,
        "tes_observables_identical": False,
        "solver_converged": True,
        "numeric_parity": {"rhs_linf": 1.0e-12},
    }
    result = validate_gate(baseline, candidate, 1.0e-10)
    assert result["passed"] is False
    assert result["checks"]["tes_observables_identical"] is False


def test_classification_is_not_allowed_to_overrule_correctness() -> None:
    profile = make_profile([], case="x", total_wall=1.0)
    gate = {"passed": False}
    result = classify(profile, profile, gate)
    assert result["classification"] == "REJECTED_CORRECTNESS_GATE"


def test_phase24_log_events_are_read_as_wall_seconds(tmp_path) -> None:
    path = tmp_path / "solver.log"
    path.write_text(
        "PHASE24_EVENT: component=global_sparse_insertion wall_seconds= 1.2500E+00\n"
        "PHASE24_EVENT: component=udf_circuit wall_seconds= 2.5E-01\n",
        encoding="utf-8",
    )
    events = parse_phase24_log(path)
    assert events[0]["component"] == "global_sparse_insertion"
    assert events[0]["wall_seconds"] == pytest.approx(1.25)
