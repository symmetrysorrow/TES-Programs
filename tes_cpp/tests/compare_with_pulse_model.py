"""Numerically compare C++ ``posi2pulse`` with simu_py/pulse_model.py.

Run from the repository root after building ``posi2pulse``::

    python tes_cpp/tests/compare_with_pulse_model.py

The default input is simu_py/input.json, so both implementations receive the
same physical parameters and absorber positions.  Use --input and
--executable to compare another built executable or input file.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np


TES_CPP = Path(__file__).resolve().parents[1]
WORKSPACE = TES_CPP.parent
sys.path.insert(0, str(WORKSPACE / "simu_py"))
from pulse_model import model  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=WORKSPACE / "simu_py" / "input.json")
    parser.add_argument(
        "--executable",
        type=Path,
        default=Path(os.environ.get("TES_CPP_POSI2PULSE_EXECUTABLE", TES_CPP / "build-release" / "posi2pulse.exe")),
    )
    parser.add_argument("--rtol", type=float, default=1e-8)
    # Eigen and SciPy may reconstruct the t=0 state with a few e-18 of
    # round-off error even when every meaningful signal sample agrees.
    parser.add_argument("--atol", type=float, default=1e-16)
    args = parser.parse_args()

    parameters = json.loads(args.input.read_text(encoding="utf-8"))
    positions = parameters["position"]
    expected_ch0, expected_ch1 = model(parameters)

    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "pulses.json"
        subprocess.run(
            [str(args.executable), str(args.input), str(output), "--positions", ",".join(map(str, positions))],
            check=True,
        )
        actual_document = json.loads(output.read_text(encoding="utf-8"))

    assert actual_document["input"] == parameters
    actual = actual_document["pulses"]
    assert [int(position) for position in actual] == positions
    expected_time = np.linspace(0, parameters["samples"] / parameters["rate"], int(parameters["samples"]))
    np.testing.assert_allclose(actual_document["time"], expected_time, rtol=0, atol=0)
    for index, position in enumerate(positions):
        pulse = actual[str(position)]
        np.testing.assert_allclose(pulse["ch0"], expected_ch0[index], rtol=args.rtol, atol=args.atol)
        np.testing.assert_allclose(pulse["ch1"], expected_ch1[index], rtol=args.rtol, atol=args.atol)

    print(f"PASS: {len(positions)} positions, {int(parameters['samples'])} samples; rtol={args.rtol:g}, atol={args.atol:g}")


if __name__ == "__main__":
    main()
