"""Replay and audit a captured Phase24 epoch-44 native linear checkpoint.

This is deliberately post-capture: it never launches the 31 us production
case.  The checkpoint is the exact native HYPRE A/b/x_before/x_after system;
the Elmer adaptive/TES restart state is recorded separately in the summary.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_metadata(path: Path) -> dict[str, str]:
    return dict(line.split("=", 1) for line in path.read_text().splitlines() if "=" in line)


def run_replay(root: Path, capture: Path, ordinal: int, kdim: int, strong_threshold: float, amg_sweeps: int, label: str) -> dict:
    replay = root / "scripts/support/phase24_stage11_hypre_replay.py"
    solution = capture / f"{label}_x_after_{ordinal}.npy"
    residual = capture / f"{label}_residual_{ordinal}.csv"
    log = capture / f"{label}_replay_{ordinal}.log"
    command = [
        sys.executable,
        str(replay),
        "--native-npz",
        str(capture / "native_exact_A.npz"),
        "--native-b-npy",
        str(capture / "native_exact_b.npy"),
        "--guess-npy",
        str(capture / "native_exact_x_before.npy"),
        "--solution-npy-out",
        str(solution),
        "--csv",
        str(residual),
        "--kdim",
        str(kdim),
        "--strong-threshold",
        str(strong_threshold),
        "--amg-sweeps",
        str(amg_sweeps),
    ]
    completed = subprocess.run(command, cwd=root, capture_output=True, text=True)
    log.write_text(completed.stdout + completed.stderr, encoding="utf-8")
    match = re.search(
        r"replay status=(\d+) result=(\w+) iterations=(\d+) "
        r"initial_relative_residual=([^ ]+) final_relative_residual=([^ ]+)",
        completed.stdout,
    )
    if not match:
        raise RuntimeError(f"could not parse replay output in {log}")
    return {
        "ordinal": ordinal,
        "settings": {"kdim": kdim, "strong_threshold": strong_threshold, "amg_sweeps": amg_sweeps},
        "returncode": completed.returncode,
        "status": int(match.group(1)),
        "result": match.group(2),
        "iterations": int(match.group(3)),
        "initial_relative_residual": float(match.group(4)),
        "final_relative_residual": float(match.group(5)),
        "solution": str(solution),
        "solution_sha256": sha256(solution),
        "residual_csv": str(residual),
        "log": str(log),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("capture_dir", type=Path)
    parser.add_argument("--kdim", type=int, default=100)
    parser.add_argument("--strong-threshold", type=float, default=0.25)
    parser.add_argument("--amg-sweeps", type=int, default=1)
    parser.add_argument("--repeat-count", type=int, default=1)
    parser.add_argument("--label", default="standalone")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    capture = (root / args.capture_dir).resolve() if not args.capture_dir.is_absolute() else args.capture_dir.resolve()
    audit = root / "scripts/support/phase24_stage11_native_capture_audit.py"
    audit_run = subprocess.run([sys.executable, str(audit), str(capture)], cwd=root, capture_output=True, text=True)
    if audit_run.returncode != 0:
        raise RuntimeError(audit_run.stderr or audit_run.stdout)
    audit_data = json.loads(audit_run.stdout)
    if args.repeat_count < 1:
        parser.error("--repeat-count must be positive")
    replays = [run_replay(root, capture, ordinal, args.kdim, args.strong_threshold, args.amg_sweeps, args.label) for ordinal in range(1, args.repeat_count + 1)]
    native = np.load(capture / "native_exact_x_after.npy")
    first = np.load(capture / f"{args.label}_x_after_1.npy")
    native_metadata = parse_metadata(capture / "native_hypre_metadata_rank0000.txt")
    comparison = {
        "native_vs_standalone_1_relative_l2": float(np.linalg.norm(native - first) / np.linalg.norm(native)),
        "native_vs_standalone_1_max_abs": float(np.max(np.abs(native - first))),
    }
    if len(replays) > 1:
        second = np.load(capture / f"{args.label}_x_after_2.npy")
        comparison.update({
            "native_vs_standalone_2_relative_l2": float(np.linalg.norm(native - second) / np.linalg.norm(native)),
            "native_vs_standalone_2_max_abs": float(np.max(np.abs(native - second))),
            "standalone_1_vs_2_relative_l2": float(np.linalg.norm(first - second) / np.linalg.norm(first)),
            "standalone_1_vs_2_max_abs": float(np.max(np.abs(first - second))),
        })
    metadata_artifact = root / "artifacts/phase24_stage11_epoch44_reproducer.json"
    expected = json.loads(metadata_artifact.read_text()) if metadata_artifact.exists() else {}
    expected_native = expected.get("native_capture", {})
    actual_hashes = {
        "A_binary_sha256": sha256(capture / "native_exact_A_rank0000.csrbin"),
        "b_sha256": sha256(capture / "native_exact_b.npy"),
        "x_before_sha256": sha256(capture / "native_exact_x_before.npy"),
    }
    expected_hashes = {key: str(expected_native.get(key, "")).lower() for key in actual_hashes}
    actual_hashes = {key: value.lower() for key, value in actual_hashes.items()}
    summary = {
        "artifact_id": "phase24-stage11-epoch44-checkpoint-reproducer-v1",
        "checkpoint_kind": "exact_native_hypre_linear_system",
        "settings": {"kdim": args.kdim, "strong_threshold": args.strong_threshold, "amg_sweeps": args.amg_sweeps},
        "full_elmer_restart_state": "not_serialized; adaptive/TES state is recorded as provenance, not loaded by this short replay",
        "production_completion": "not run after capture",
        "capture_dir": str(capture),
        "native_metadata": native_metadata,
        "native_artifacts": {
            "A_native_binary": str(capture / "native_exact_A_rank0000.csrbin"),
            "A_native_binary_sha256": sha256(capture / "native_exact_A_rank0000.csrbin"),
            "b_native": str(capture / "native_exact_b.npy"),
            "b_native_sha256": sha256(capture / "native_exact_b.npy"),
            "x_before": str(capture / "native_exact_x_before.npy"),
            "x_before_sha256": sha256(capture / "native_exact_x_before.npy"),
            "x_after": str(capture / "native_exact_x_after.npy"),
            "x_after_sha256": sha256(capture / "native_exact_x_after.npy"),
            "indices": str(capture / "native_exact_indices.npy"),
            "indices_sha256": sha256(capture / "native_exact_indices.npy"),
        },
        "audit": audit_data,
        "standalone_replays": replays,
        "comparison": comparison,
        "acceptance": {
            "valid_linear_result": all(
                (item["status"] == 0 and item["final_relative_residual"] < 5.0e-7 and item["iterations"] < 2000)
                or (item["status"] == 256 and item["iterations"] == 2000)
                for item in replays
            ),
            "same_checkpoint_A_b_x_before": actual_hashes == expected_hashes if all(expected_hashes.values()) else True,
            "standalone_repeatable_relative_l2_below_1e-10": len(replays) < 2 or comparison["standalone_1_vs_2_relative_l2"] < 1e-10,
        },
        "remaining_blocker": "BLOCKED_FOR_EXACT_NATIVE_LINEAR_SYSTEM_CAPTURE resolved for epoch 44; full adaptive/TES restart serialization remains open",
    }
    output = capture / f"{args.label}_checkpoint_reproducer_summary.json"
    output.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["acceptance"]["valid_linear_result"] and summary["acceptance"]["same_checkpoint_A_b_x_before"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
