"""Automatically choose SinglePixel spatial and early-time resolutions.

The tuner varies one factor at a time, keeps every previously accepted factor
fixed, and stops each spatial axis at its first accepted adjacent refinement.
Use ``--seed stycast_layers=32`` to continue from an already completed pilot.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.prep import run_singlepixel_resolution_pilot as pilot


CONFIG = ROOT / "singlepixel_resolution_optimization.json"
PILOT = Path(pilot.__file__).resolve()
OUT = ROOT / "artifacts" / "optimization" / "singlepixel_resolution_autotune.json"


def slug(value: float) -> str:
    return (f"{value:g}").replace(".", "p").replace("-", "m")


def parse_seed(values: list[str], config: dict) -> dict[str, float]:
    selected: dict[str, float] = {}
    for value in values:
        axis, sep, raw = value.partition("=")
        if not sep or axis not in config["spatial_axes"]:
            raise ValueError(f"--seed must name a spatial axis: {value!r}")
        number = float(raw)
        if number <= 0:
            raise ValueError(f"--seed must be positive: {value!r}")
        selected[axis] = number
    return selected


def context_tag(selected: dict[str, float]) -> str:
    """Return a short deterministic namespace for a chain of fixed settings.

    Mesh/case/state names incorporate this tag.  Keeping it short is required
    because the TES UDF accepts state-file paths of at most 128 characters.
    The full setting dictionary remains in the pilot JSON and manifest.
    """
    if not selected:
        return "auto_base"
    encoded = json.dumps(sorted(selected.items()), separators=(",", ":"), ensure_ascii=True).encode("ascii")
    return "auto_" + hashlib.sha256(encoded).hexdigest()[:6]


def selected_spatial_level(manifest: dict, levels: list[float]) -> tuple[float, str]:
    for index, comparison in enumerate(manifest.get("comparisons", [])):
        if comparison["accepted"]:
            return levels[index + 1], "accepted"
    return levels[-1], "unresolved_at_finest_candidate"


def selected_time_step(manifest: dict, levels: list[float]) -> tuple[float, str]:
    for index, comparison in enumerate(manifest.get("comparisons", [])):
        if comparison["accepted"]:
            return levels[index + 1], "accepted"
    return levels[0], "finest_candidate_required"


def run_axis(args: argparse.Namespace, config: dict, axis: str, selected: dict[str, float]) -> dict:
    levels = [float(value) for value in (config["time_step_candidates_us"] if axis == "time" else config["spatial_axes"][axis])]
    tag = context_tag(selected)
    command = [sys.executable, str(PILOT), "--axis", axis, "--tag", tag,
               "--end-us", f"{args.end_us:g}", "--early-step-us", f"{args.early_step_us:g}",
               "--linear-system", args.linear_system]
    for name, value in sorted(selected.items()):
        command.extend(["--fixed-setting", f"{name}={value:g}"])
    if args.tail_start_us is not None:
        command.extend(["--tail-start-us", f"{args.tail_start_us:g}"])
    if args.elmer_solver is not None:
        command.extend(["--elmer-solver", str(args.elmer_solver)])
    if args.runtime_bin is not None:
        command.extend(["--runtime-bin", str(args.runtime_bin)])
    if args.execute:
        command.append("--execute")
    print("autotune:", subprocess.list2cmdline(command))
    if not args.execute:
        return {"axis": axis, "levels": levels, "tag": tag, "status": "planned", "command": command}
    subprocess.run(command, cwd=ROOT, check=True)
    _, manifest_path = pilot.pilot_output_paths(axis, tag)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    value, status = (
        selected_time_step(manifest, levels) if axis == "time"
        else selected_spatial_level(manifest, levels)
    )
    return {"axis": axis, "levels": levels, "tag": tag, "selected": value,
            "status": status, "manifest": str(manifest_path.relative_to(ROOT)),
            "comparisons": manifest.get("comparisons", [])}


def main() -> None:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", action="append", default=[], metavar="AXIS=VALUE")
    parser.add_argument("--end-us", type=float, default=105.0)
    parser.add_argument("--early-step-us", type=float, default=0.625)
    parser.add_argument("--tail-start-us", type=float)
    parser.add_argument("--linear-system", default="mumps", choices=["direct", "iterative", "iterative_hypre_boomeramg", "mumps"])
    parser.add_argument("--elmer-solver", type=Path)
    parser.add_argument("--runtime-bin", type=Path)
    parser.add_argument("--execute", action="store_true")
    args = parser.parse_args()

    selected = parse_seed(args.seed, config)
    report = {"execute": args.execute, "seed": selected.copy(), "axes": [], "selected": selected}
    for axis in config["spatial_axes"]:
        if axis in selected:
            report["axes"].append({"axis": axis, "selected": selected[axis], "status": "seeded"})
            continue
        outcome = run_axis(args, config, axis, selected)
        report["axes"].append(outcome)
        if args.execute:
            selected[axis] = outcome["selected"]
            report["selected"] = selected.copy()
            if outcome["status"] != "accepted":
                report["status"] = "needs_more_refinement"
                break
    else:
        outcome = run_axis(args, config, "time", selected)
        report["time"] = outcome
        if args.execute:
            report["selected"]["early_step_us"] = outcome["selected"]
            report["status"] = outcome["status"]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
