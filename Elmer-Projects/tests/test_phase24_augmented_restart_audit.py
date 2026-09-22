from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "support" / "audit_phase24_restart_state.py"
SPEC = importlib.util.spec_from_file_location("phase24_restart_state_audit", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
audit = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)

from phase24_restart_state_audit import (
    materialize,
    read_stage_csv,
    resolve_runtime_state_path,
)


HEADER = (
    "stage,physical_time,timestep,dt,bdf_order,state_file,state_exists,state_open_success," 
    "state_parse_success,state_loaded,current,previous_current,resistance,power,average_temperature," 
    "sweep_temperature,omega,omega_cap,prev_residual,last_circuit_dt,last_circuit_step," 
    "circuit_iter_in_step,circuit_calls_in_step,circuit_initialized,circuit_power_initialized,source"
)


def write_audit(path: Path, loaded: str = "T") -> None:
    rows = [
        "load_after,1.02,1,1,0,work/state.state,T,T,T," + loaded + ",1e-4,1e-4,1e-2,3e-10,0,0,.5,.5,0,-1,-2147483647,0,0,F,F,checkpoint_or_default",
        "circuit_init_after,1.02,1,1,0,work/state.state,T,T,T," + loaded + ",1e-4,1e-4,1e-2,3e-10,1.68,1.68,.5,.5,0,-1,1,1,1,T,T,loaded_or_default_then_init",
        "pre_first_assembly,1.02,1,1,0,work/state.state,T,T,T," + loaded + ",1e-4,1e-4,1e-2,3e-10,1.68,1.68,.5,.5,0,-1,1,1,1,T,T,preserved_before_first_assembly",
    ]
    path.write_text("\n".join((HEADER, *rows)) + "\n", encoding="utf-8")


def test_checkpoint_path_uses_solver_cwd() -> None:
    root = Path("D:/repo")
    assert resolve_runtime_state_path("work/meshes/gate3.state", root) == root / "work/meshes/gate3.state"
    assert resolve_runtime_state_path("../../work/meshes/gate3.state", root) == Path("D:/work/meshes/gate3.state")


def test_state_audit_detects_load_failure_and_restoration(tmp_path: Path) -> None:
    audit = tmp_path / "audit.csv"
    write_audit(audit)
    saved = tmp_path / "saved.state"
    actual = tmp_path / "work" / "state.state"
    actual.parent.mkdir()
    saved.write_text("1.68 1e-4 1e-2 3e-10 1e-4\n", encoding="utf-8")
    actual.write_text(saved.read_text(encoding="utf-8"), encoding="utf-8")
    report = materialize(audit, saved, tmp_path / "out", runtime_cwd=tmp_path)
    assert report["all_required_fields_loaded"] is True
    assert report["first_mismatch"] is None
    assert report["path_audit"]["diagnostic_snapshot_isolated"] is True
    assert report["stages"]["load_after"]["current_A"] == report["saved_checkpoint"]["current_A"]

    failed = tmp_path / "failed.csv"
    write_audit(failed, loaded="F")
    failed_report = materialize(failed, saved, tmp_path / "failed_out", runtime_cwd=tmp_path)
    assert failed_report["all_required_fields_loaded"] is False


def test_original_checkpoint_hash_is_protected(tmp_path: Path) -> None:
    audit = tmp_path / "audit.csv"
    write_audit(audit)
    saved = tmp_path / "saved.state"
    before = tmp_path / "before.state"
    actual = tmp_path / "work" / "state.state"
    actual.parent.mkdir()
    content = "1.68 1e-4 1e-2 3e-10 1e-4\n"
    saved.write_text(content, encoding="utf-8")
    before.write_text(content, encoding="utf-8")
    actual.write_text(content, encoding="utf-8")
    report = materialize(audit, saved, tmp_path / "out", original_before=before, runtime_cwd=tmp_path)
    assert report["path_audit"]["original_checkpoint_overwritten"] is False
