"""Run one bounded epoch-44 capture from the existing Phase24 restart.

This runner creates a temporary SIF with the existing solver settings plus
the already-supported Phase24 diagnostic keyword.  It stops the Elmer process
as soon as the native epoch-44 A/b/x capture is complete, so it never runs the
remaining 31 us interval.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path


CAPTURE_FILES = (
    "native_exact_A_rank0000.csrbin",
    "native_hypre_indices_rank0000.dat",
    "native_hypre_b_rank0000.dat",
    "native_hypre_x_before_rank0000.dat",
    "native_hypre_x_after_rank0000.dat",
    "native_hypre_metadata_rank0000.txt",
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-sif", type=Path, default=Path("generated/cases/case_phase24_adaptive_output_smoke.sif"))
    parser.add_argument("--elmer-solver", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--timeout-seconds", type=float, default=1800.0)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    source_sif = (root / args.source_sif).resolve() if not args.source_sif.is_absolute() else args.source_sif.resolve()
    solver = args.elmer_solver.resolve()
    output_dir = (root / args.output_dir).resolve() if not args.output_dir.is_absolute() else args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    capture_sif = output_dir / "epoch44_capture.sif"
    run_log = output_dir / "epoch44_capture_run.log"
    prefix = (output_dir / "epoch44_state").relative_to(root).as_posix()
    source = source_sif.read_text(encoding="utf-8")
    marker = "  Phase24 HYPRE Reuse = Logical True\n"
    injection = (
        marker
        + "  Phase24 HYPRE Diagnostic = Logical True\n"
        + f'  "Phase24 HYPRE Diagnostic Prefix" = String "{prefix}"\n'
    )
    if marker not in source:
        raise RuntimeError("Phase24 HYPRE Reuse marker was not found in source SIF")
    capture_sif.write_text(source.replace(marker, injection, 1), encoding="utf-8")

    env = dict(__import__("os").environ)
    env["PHASE24_HYPRE_CAPTURE_DIR"] = output_dir.relative_to(root).as_posix()
    start = time.monotonic()
    with run_log.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            [str(solver), str(capture_sif)],
            cwd=root,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        captured = False
        try:
            while True:
                elapsed = time.monotonic() - start
                if all((output_dir / name).exists() for name in CAPTURE_FILES):
                    captured = True
                    process.terminate()
                    break
                if process.poll() is not None:
                    break
                if elapsed > args.timeout_seconds:
                    process.kill()
                    break
                time.sleep(1.0)
            try:
                process.wait(timeout=30.0)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=30.0)
        finally:
            if process.poll() is None:
                process.kill()
                process.wait(timeout=30.0)
    elapsed = time.monotonic() - start
    result = {
        "capture_complete": captured,
        "runtime_seconds": elapsed,
        "solver_returncode": process.returncode,
        "capture_dir": str(output_dir),
        "capture_sif": str(capture_sif),
        "run_log": str(run_log),
        "missing_files": [name for name in CAPTURE_FILES if not (output_dir / name).exists()],
    }
    (output_dir / "epoch44_runner_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if captured else 1


if __name__ == "__main__":
    raise SystemExit(main())
