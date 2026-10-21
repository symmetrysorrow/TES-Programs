"""Migration guards: signed reaction-difference converter must be unreachable."""

from __future__ import annotations

import ast
import os
from pathlib import Path
import shutil
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[3]
NATIVE = ROOT / "tes_cpp" / "dump2event"


def test_legacy_implementation_not_in_native_build():
    cmake = (NATIVE / "CMakeLists.txt").read_text(encoding="utf-8")
    assert "src/converter.cpp" not in cmake
    assert not (NATIVE / "src/converter.cpp").exists()
    assert not (NATIVE / "include/tes_cpp/converter.hpp").exists()
    assert (NATIVE / "src/physical_deposit.cpp").exists()
    assert (NATIVE / "src/physical_dump_reconstruct.cpp").exists()


def test_only_one_dump2event_action_and_it_uses_pb_v2_ui():
    tree = ast.parse((ROOT / "PoST_Simulations/PoST_Simulation.py").read_text(encoding="utf-8"))
    function = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "Dump2Event")
    calls = [n.func.id for n in ast.walk(function) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)]
    assert "run_pb_egs5_ui" in calls
    assert "Dump2EventPbEgs5" not in calls
    source = (ROOT / "PoST_Simulations/PoST_Simulation.py").read_text(encoding="utf-8")
    assert 'value="dump_to_event_pb_egs5"' not in source
    assert '"--normalize-deposits"' not in source


@pytest.mark.parametrize("args", [
    ["--input-energy", "2"],
    ["--normalize-deposits"],
    ["--full-energy-only"],
    ["--save-all"],
    ["--summary-only"],
    ["--physics", "legacy"],
])
def test_old_native_flags_are_rejected(args, tmp_path):
    executable = os.environ.get("TES_CPP_DUMP2EVENT_EXECUTABLE") or shutil.which(
        "dump2event.exe" if os.name == "nt" else "dump2event"
    )
    if not executable:
        pytest.skip("native dump2event executable unavailable")
    dump = tmp_path / "dumpall.dat"
    dump.write_text("test", encoding="utf-8")
    output = tmp_path / "output.h5"
    result = subprocess.run([executable, str(dump), str(output), *args],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert not output.exists()


def test_native_without_contract_never_produces_file(tmp_path):
    executable = os.environ.get("TES_CPP_DUMP2EVENT_EXECUTABLE") or shutil.which(
        "dump2event.exe" if os.name == "nt" else "dump2event"
    )
    if not executable:
        pytest.skip("native dump2event executable unavailable")
    dump = tmp_path / "dumpall.dat"
    dump.write_text("test", encoding="utf-8")
    output = tmp_path / "output.h5"
    result = subprocess.run([executable, str(dump), str(output)],
                            capture_output=True, text=True)
    assert result.returncode != 0
    assert "physics-config" in result.stderr
    assert not output.exists()
