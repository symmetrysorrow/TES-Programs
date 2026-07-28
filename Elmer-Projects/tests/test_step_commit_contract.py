from __future__ import annotations

from pathlib import Path

import pytest

from scripts.support.build_cases import solver1_block


def commit(
    temperature: float,
    previous_current: float,
    step_dt: float,
    *,
    bias_current: float = 0.000715,
    shunt_resistance: float = 0.0039,
    r0: float = 0.015527,
    rmin: float = 1.0e-6,
    alpha: float = 256.46,
    beta: float = 5.03,
    i0: float = 0.000143537344932311,
    t0: float = 0.16857,
    inductance: float = 1.23e-8,
) -> tuple[float, float, float]:
    """Python transcription of the opt-in Fortran step-boundary equation."""
    a = r0 * (1.0 + alpha * (temperature - t0) / t0 - beta)
    b = r0 * beta / i0
    c = shunt_resistance + a + inductance / step_dt
    discriminant = c * c + 4.0 * b * (
        bias_current * shunt_resistance
        + inductance * previous_current / step_dt
    )
    current = max(
        min((discriminant**0.5 - c) / (2.0 * b), bias_current),
        0.0,
    )
    resistance = a + b * abs(current)
    if resistance < rmin:
        resistance = rmin
        current = (
            bias_current * shunt_resistance
            + inductance * previous_current / step_dt
        ) / (shunt_resistance + resistance + inductance / step_dt)
    return current, resistance, current * current * resistance


def test_commit_backward_euler_snapshot() -> None:
    assert commit(0.168, 0.000164, 1.0e-5) == pytest.approx(
        (
            0.0001638579017319079,
            0.013118860557804382,
            3.5223369156206987e-10,
        ),
        rel=1.0e-14,
    )


def test_commit_rmin_path() -> None:
    assert commit(0.01, 0.000164, 1.0e-5) == pytest.approx(
        (0.0005827752874683297, 1.0e-6, 3.3962703568379435e-13),
        rel=1.0e-14,
    )


def test_commit_branch_does_not_use_sweep_history() -> None:
    source = (
        Path(__file__).resolve().parents[2]
        / "tools/elmer-hypre/src/fem/src/modules/HeatSolve.F90"
    ).read_text(encoding="utf-8")
    branch = source.split("IF (StepCommit) THEN", 1)[1].split("ELSE", 1)[0]
    assert "CircuitCallsInStep" not in branch
    assert "CircuitIterInStep" not in branch
    assert "MidTemperature" not in branch


def test_sif_step_commit_is_opt_in() -> None:
    solver = {
        "nonlinear_max_iterations": 25,
        "nonlinear_convergence_tolerance": 3.0e-7,
        "nonlinear_relaxation_factor": 1.0,
        "steady_state_convergence_tolerance": 1.0e-8,
        "linear_system": "mumps",
    }
    off = "\n".join(
        solver1_block(solver, inner_circuit=True, tes_body_id=8)
    )
    on = "\n".join(
        solver1_block(
            solver,
            inner_circuit=True,
            tes_body_id=8,
            inner_circuit_step_commit=True,
        )
    )
    assert "TES Inner Circuit Step Commit" not in off
    assert '"TES Inner Circuit Step Commit" = Logical True' in on
