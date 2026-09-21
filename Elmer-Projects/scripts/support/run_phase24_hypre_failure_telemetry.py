"""Re-run the 1e-8 Phase24 case and collect native HYPRE failure telemetry.

This diagnostic assumes the native Elmer/HYPRE binary has been rebuilt from
the ``SolveHypre.c`` telemetry change.  It intentionally reuses the existing
1-us ``1e-8`` project, disables the unrelated epoch-44 capture selector, and
records the first native failure without changing the numerical solver policy.

Run from the repository root::

    python scripts/support/run_phase24_hypre_failure_telemetry.py

The script returns non-zero when Elmer fails; that is expected for the
diagnostic if the current ``1e-8`` tolerance remains unattainable.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CASE = "case_phase24_short_hypre_tol1e8_1us"
PROJECT = ROOT / "artifacts/phase24_hypre_tol1e8_diagnostic_1us/phase24_hypre_tol1e8_diagnostic.json"
RESULT_DIR = ROOT / "results" / CASE
LAUNCHER_LOG = RESULT_DIR / "failure_telemetry_launcher.log"
SOLVER_LOG = RESULT_DIR / "solver.log"
ARTIFACT_DIR = ROOT / "artifacts/phase24_hypre_failure_telemetry_1us"


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def command_for(solver: Path, runtime_bin: Path, toolchain_bin: Path) -> list[str]:
    return [
        sys.executable,
        str(ROOT / "run.py"), CASE,
        "--project", str(PROJECT), "--skip-sync",
        "--elmer-solver", str(solver),
        "--runtime-bin", str(runtime_bin),
        "--toolchain-bin", str(toolchain_bin),
        "--mpi-procs", "1",
    ]


def parse_failure_telemetry(path: Path) -> dict:
    if not path.is_file():
        return {"solver_log_exists": False, "telemetry_records": []}
    text = path.read_text(encoding="utf-8", errors="replace")
    records = []
    pattern = re.compile(r"PHASE24_HYPRE_SOLVE_FAILURE\s+([^\n\r]+)")
    for match in pattern.finditer(text):
        record = {}
        for item in match.group(1).split():
            if "=" not in item:
                continue
            key, value = item.split("=", 1)
            try:
                record[key] = int(value)
            except ValueError:
                try:
                    record[key] = float(value.replace("D", "E").replace("d", "e"))
                except ValueError:
                    record[key] = value
        records.append(record)
    return {
        "solver_log_exists": True,
        "all_done": "ALL DONE" in text,
        "stop_1": text.count("STOP 1"),
        "telemetry_records": records,
        "telemetry_missing": "PHASE24_HYPRE_SOLVE_FAILURE" not in text,
        "selector_warning_present": "NATIVE_XVEC_CAPTURE_CANDIDATE_SEEN_BUT_SELECTOR_MISMATCH" in text,
        "fatal_message_present": "HYPRE failed to satisfy the production linear tolerance" in text,
    }


def write_summary(payload: dict) -> None:
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    (ARTIFACT_DIR / "summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    records = payload["audit"].get("telemetry_records", [])
    lines = [
        "# Phase24 HYPRE failure telemetry diagnostic",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        "- Case: `case_phase24_short_hypre_tol1e8_1us`",
        "- Numerical policy: unchanged; only failure telemetry and capture suppression are diagnostic changes",
        f"- Exit code: `{payload['run'].get('exit_code')}`",
        f"- ALL DONE: `{payload['audit'].get('all_done')}`",
        f"- Telemetry records: `{len(records)}`",
        f"- Telemetry missing: `{payload['audit'].get('telemetry_missing')}`",
        "",
        "## Native failure records",
        "",
    ]
    if records:
        lines.append("```json")
        lines.append(json.dumps(records, indent=2, ensure_ascii=False))
        lines.append("```")
    else:
        lines.append("No `PHASE24_HYPRE_SOLVE_FAILURE` record was found. Verify that the rebuilt binary is the one being executed.")
    lines.extend([
        "",
        f"Solver log: `{SOLVER_LOG}`",
        f"Launcher log: `{LAUNCHER_LOG}`",
        "",
    ])
    (ARTIFACT_DIR / "summary.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\bin\ElmerSolver_mpi.exe"))
    parser.add_argument("--runtime-bin", type=Path, default=Path(r"D:\Github\TES-Programs\tools\elmer-hypre\install-stage11\lib\elmersolver"))
    parser.add_argument("--toolchain-bin", type=Path, default=Path(r"C:\msys64\ucrt64\bin"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    if not PROJECT.is_file():
        raise SystemExit(f"diagnostic project not found: {PROJECT}")
    solver = absolute(args.solver)
    runtime_bin = absolute(args.runtime_bin)
    toolchain_bin = absolute(args.toolchain_bin)
    if not args.dry_run:
        for path, label in ((solver, "solver"), (runtime_bin, "runtime DLL directory"), (toolchain_bin, "toolchain DLL directory")):
            if not path.exists():
                raise SystemExit(f"{label} not found: {path}")

    command = command_for(solver, runtime_bin, toolchain_bin)
    print("[failure-telemetry] command: " + " ".join(command))
    if args.dry_run:
        run = {"command": command, "exit_code": None, "log": str(LAUNCHER_LOG)}
        audit = {}
    else:
        ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
        started = datetime.now(timezone.utc).isoformat()
        start = time.monotonic()
        env = os.environ.copy()
        env["PHASE24_DISABLE_NATIVE_CAPTURE"] = "1"
        with LAUNCHER_LOG.open("w", encoding="utf-8") as handle:
            process = subprocess.run(command, cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
        run = {
            "command": command,
            "started": started,
            "finished": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic() - start,
            "exit_code": process.returncode,
            "log": str(LAUNCHER_LOG),
        }
        audit = parse_failure_telemetry(SOLVER_LOG)

    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(PROJECT),
        "solver_log": str(SOLVER_LOG),
        "run": run,
        "audit": audit,
    }
    write_summary(payload)
    print(f"summary={ARTIFACT_DIR / 'summary.md'}")
    return 1 if run.get("exit_code") not in (0, None) else 0


if __name__ == "__main__":
    raise SystemExit(main())
