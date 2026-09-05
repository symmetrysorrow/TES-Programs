"""Regression test for the common TES electrical/thermal linearization.

The C++ solver is a distributed absorber model while Python's noise solver is
the reduced model.  The absorber boundary conductance therefore differs
(`G_abs-tes` versus `G_eff`); this test compares the shared intrinsic TES
sub-block and records that reduction boundary explicitly.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SIMULATION_ROOT = ROOT / "PoST_Simulations"
sys.path.insert(0, str(SIMULATION_ROOT))

from PoST_Simulation import tes_linearized_time_matrix, tes_operating_point  # noqa: E402


INPUT = ROOT / "PoST_Simulations" / "input.json"
EXE_CANDIDATES = (
    ROOT / "tes_cpp" / "build" / "posi2pulse" / "Release" / "posi2pulse.exe",
    ROOT / "tes_cpp" / "build-release" / "posi2pulse.exe",
)


def _cpp_summary():
    executable = next((path for path in EXE_CANDIDATES if path.exists()), None)
    if executable is None:
        pytest.skip("C++ posi2pulse debug executable is not built")
    completed = subprocess.run(
        [str(executable), "--dump-linearization", str(INPUT)],
        check=True,
        capture_output=True,
        text=True,
    )
    return json.loads(completed.stdout)


def test_cpp_python_intrinsic_tes_subblock_parity():
    parameters = json.loads(INPUT.read_text())
    point = tes_operating_point(parameters)
    matrix = tes_linearized_time_matrix(parameters)
    cpp = _cpp_summary()

    assert np.isclose(cpp["current_A"], point["current_A"], rtol=1e-12)
    assert np.isclose(cpp["tau_el_s"], parameters["L"] / (parameters["R_l"] + parameters["R"] * (1 + parameters["beta"])), rtol=1e-12)
    assert np.isclose(cpp["loop_gain"], parameters["alpha"] * point["current_A"] ** 2 * parameters["R"] / (parameters["G_tes-bath"] * parameters["T_c"]), rtol=1e-12)

    expected_common = np.array([matrix[0, 0], matrix[0, 1], matrix[1, 0]])
    actual_common = np.asarray(cpp["tes1_time_block"][:3])
    np.testing.assert_allclose(actual_common, expected_common, rtol=2e-12, atol=1e-7)

    intrinsic_python_thermal = -(1.0 - cpp["loop_gain"]) * parameters["G_tes-bath"] / parameters["C_tes"]
    assert np.isclose(cpp["tes_intrinsic_thermal_diag_per_s"], intrinsic_python_thermal, rtol=2e-12)
    assert np.isclose(
        cpp["tes1_time_block"][3],
        intrinsic_python_thermal - cpp["tes_boundary_rate_per_s"],
        rtol=2e-12,
    )
