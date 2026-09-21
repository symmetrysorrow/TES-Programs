"""Collect the native FlexGMRES residual progression for the 1e-8 case.

The existing ``1e-8 / max-iterations=10000`` case is reused.  The rebuilt
solver is run with ``PHASE24_HYPRE_RESIDUAL_TRACE=1`` so HYPRE print level 3
emits its native Krylov trace.  This is a diagnostic-only run and does not
change the production SIF policy.
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
CASE = "case_phase24_short_hypre_tol1e8_max10000_1us"
PROJECT = ROOT / "artifacts/phase24_hypre_max_iterations_diagnostic_1us/phase24_hypre_max_iterations_diagnostic.json"
SIF = ROOT / "generated/cases" / f"{CASE}.sif"
RESULT_DIR = ROOT / "results" / CASE
LAUNCHER_LOG = RESULT_DIR / "residual_trace_launcher.log"
SOLVER_LOG = RESULT_DIR / "solver.log"
ARTIFACT_DIR = ROOT / "artifacts/phase24_hypre_residual_trace_1us"


def absolute(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def command_for(solver: Path, runtime_bin: Path, toolchain_bin: Path) -> list[str]:
    return [
        sys.executable, str(ROOT / "run.py"), CASE,
        "--project", str(PROJECT), "--skip-sync",
        "--elmer-solver", str(solver),
        "--runtime-bin", str(runtime_bin),
        "--toolchain-bin", str(toolchain_bin),
        "--mpi-procs", "1",
    ]


def parse_trace(path: Path) -> dict:
    if not path.is_file():
        return {"solver_log_exists": False, "trace_lines": [], "failure_records": []}
    text = path.read_text(encoding="utf-8", errors="replace")
    iteration_pattern = re.compile(
        r"^\s*(\d+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s+([-+0-9.EeDd]+)\s*$"
    )
    trace_lines = []
    native_iterations = []
    for line in text.splitlines():
        match = iteration_pattern.match(line)
        if match:
            values = [float(item.replace("D", "E").replace("d", "e")) for item in match.groups()]
            native_iterations.append({
                "iteration": int(values[0]),
                "residual": values[1],
                "relative_step": values[2],
                "reported_residual": values[3],
            })
            trace_lines.append(line)
        elif (
            "PHASE24_HYPRE_RESIDUAL_TRACE" in line
            or "FlexGMRES" in line
            or "flexgmres" in line
            or re.search(r"(?:iteration|residual|norm)\s*[=:]", line, re.IGNORECASE)
        ):
            trace_lines.append(line)
    failure_records = []
    for match in re.finditer(r"PHASE24_HYPRE_SOLVE_FAILURE\s+([^\n\r]+)", text):
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
        failure_records.append(record)
    return {
        "solver_log_exists": True,
        "all_done": "ALL DONE" in text,
        "stop_1": text.count("STOP 1"),
        "trace_enabled_marker": "PHASE24_HYPRE_RESIDUAL_TRACE_ENABLED" in text,
        "trace_line_count": len(trace_lines),
        "native_iteration_count": len(native_iterations),
        "native_iterations_head": native_iterations[:10],
        "native_iterations_tail": native_iterations[-10:],
        "native_residual_min": min((item["residual"] for item in native_iterations), default=None),
        "native_residual_final": native_iterations[-1]["residual"] if native_iterations else None,
        "trace_head": trace_lines[:20],
        "trace_tail": trace_lines[-20:],
        "failure_records": failure_records,
        "trace_missing": "PHASE24_HYPRE_RESIDUAL_TRACE_ENABLED" not in text,
    }


def write_summary(payload: dict) -> None:
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    (ARTIFACT_DIR / "summary.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    audit = payload["audit"]
    lines = [
        "# Phase24 HYPRE native residual trace",
        "",
        f"Generated: `{payload['generated_utc']}`",
        "",
        "- Case: `1e-8`, `Linear System Max Iterations=10000`",
        "- Diagnostic environment: `PHASE24_HYPRE_RESIDUAL_TRACE=1`",
        f"- Exit code: `{payload['run'].get('exit_code')}`",
        f"- ALL DONE: `{audit.get('all_done')}`",
        f"- Trace enabled marker: `{audit.get('trace_enabled_marker')}`",
        f"- Candidate trace lines: `{audit.get('trace_line_count')}`",
        f"- Native iteration records: `{audit.get('native_iteration_count')}`",
        f"- Native minimum residual: `{audit.get('native_residual_min')}`",
        f"- Native final residual: `{audit.get('native_residual_final')}`",
        "",
        "## Interpretation",
        "",
        "The complete native trace is in `solver.log`. The head and tail below are only a compact index. Inspect whether residuals decrease steadily, stagnate, or become unstable before the 10000-iteration stop.",
        "",
        "### Trace head",
        "",
        "```text",
        *audit.get("trace_head", []),
        "```",
        "",
        "### Trace tail",
        "",
        "```text",
        *audit.get("trace_tail", []),
        "```",
        "",
        f"Solver log: `{SOLVER_LOG}`",
        f"Launcher log: `{LAUNCHER_LOG}`",
        "",
    ]
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
    if not args.dry_run and (not SIF.is_file() or "Linear System Max Iterations = 10000" not in SIF.read_text(encoding="utf-8", errors="replace")):
        raise SystemExit(f"expected 10000-iteration SIF not found: {SIF}")
    solver, runtime_bin, toolchain_bin = map(absolute, (args.solver, args.runtime_bin, args.toolchain_bin))
    if not args.dry_run:
        for path, label in ((solver, "solver"), (runtime_bin, "runtime DLL directory"), (toolchain_bin, "toolchain DLL directory")):
            if not path.exists():
                raise SystemExit(f"{label} not found: {path}")

    command = command_for(solver, runtime_bin, toolchain_bin)
    print("[residual-trace] command: " + " ".join(command))
    if args.dry_run:
        run = {"command": command, "exit_code": None, "log": str(LAUNCHER_LOG)}
        audit = {}
    else:
        RESULT_DIR.mkdir(parents=True, exist_ok=True)
        started = datetime.now(timezone.utc).isoformat()
        start = time.monotonic()
        env = os.environ.copy()
        env["PHASE24_DISABLE_NATIVE_CAPTURE"] = "1"
        env["PHASE24_HYPRE_RESIDUAL_TRACE"] = "1"
        with LAUNCHER_LOG.open("w", encoding="utf-8") as handle:
            process = subprocess.run(command, cwd=ROOT, env=env, stdout=handle, stderr=subprocess.STDOUT)
        run = {
            "command": command, "started": started,
            "finished": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": time.monotonic() - start,
            "exit_code": process.returncode, "log": str(LAUNCHER_LOG),
        }
        audit = parse_trace(SOLVER_LOG)

    payload = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "project": str(PROJECT), "sif": str(SIF),
        "run": run, "audit": audit,
        "solver_log": str(SOLVER_LOG),
    }
    write_summary(payload)
    print(f"summary={ARTIFACT_DIR / 'summary.md'}")
    return 1 if run.get("exit_code") not in (0, None) else 0


if __name__ == "__main__":
    raise SystemExit(main())
