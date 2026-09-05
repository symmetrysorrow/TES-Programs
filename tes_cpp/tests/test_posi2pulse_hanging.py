"""Smoke tests for the optional hanging TES state in the C++ pulse solver."""

import json
import os
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[2]


def _executable():
    configured = os.environ.get("TES_CPP_POSI2PULSE_EXECUTABLE")
    candidates = [
        Path(configured) if configured else None,
        ROOT / "tes_cpp" / "build" / "posi2pulse" / "Release" / "posi2pulse.exe",
        ROOT / "tes_cpp" / "build-release" / "posi2pulse.exe",
        ROOT / "tes_cpp" / "build-msvc" / "Release" / "posi2pulse.exe",
        ROOT / "tes_cpp" / "build-ninja-msvc" / "posi2pulse.exe",
    ]
    for candidate in candidates:
        if candidate is not None and candidate.is_file():
            return candidate
    return None


@pytest.mark.parametrize("model", ["none", "hanging"])
def test_position_one_and_ch1_index_for_single_absorber(tmp_path, model):
    executable = _executable()
    if executable is None:
        pytest.skip("posi2pulse executable is not built")

    source = ROOT / "PoST_Simulations" / "input.json"
    parameters = json.loads(source.read_text(encoding="utf-8"))
    parameters.update({"n_abs": 1, "samples": 16, "rate": 100000.0})
    if model == "hanging":
        parameters.update(
            {
                "tes_internal_model": "hanging",
                "C_tes_hanging": 2.0e-13,
                "G_tes-hanging": 2.0e-7,
            }
        )
    else:
        parameters["tes_internal_model"] = "none"

    input_path = tmp_path / f"{model}.json"
    output_path = tmp_path / f"{model}-pulses.json"
    input_path.write_text(json.dumps(parameters), encoding="utf-8")
    subprocess.run(
        [str(executable), str(input_path), str(output_path), "--positions", "1"],
        check=True,
        capture_output=True,
        text=True,
    )

    document = json.loads(output_path.read_text(encoding="utf-8"))
    pulse = document["pulses"]["1"]
    assert len(pulse["ch0"]) == 16
    assert len(pulse["ch1"]) == 16
    assert any(value != 0.0 for value in pulse["ch0"])
    assert any(value != 0.0 for value in pulse["ch1"])
