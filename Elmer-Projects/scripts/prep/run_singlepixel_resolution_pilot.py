"""Prepare and optionally execute one-axis SinglePixel convergence pilots.

Examples
--------
Prepare a TES 1-vs-2-layer pilot (no solver execution)::

    python scripts/prep/run_singlepixel_resolution_pilot.py --axis tes_layers

Create meshes, run the two cases, compare them, and stop when converged::

    python scripts/prep/run_singlepixel_resolution_pilot.py --axis tes_layers --execute

Only one axis is varied in one invocation.  This is deliberate: changing mesh
and time resolution together cannot identify the source of a waveform change.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from generate_hybrid_prism_geometry import inspect_mesh_file
from scripts.support.build_cases import TES_STATE_FILE_MAX_LEN


CONFIG = ROOT / "singlepixel_resolution_optimization.json"
SOURCE = ROOT / "elmer_project_hybrid_prism_phase23_pulse_tight.json"
OUT_DIR = ROOT / "artifacts" / "optimization" / "pilots"
SOURCE_STEADY = "case_tes_steady_hybrid_prism_stack17_abs35r50_noextend_inner_1rank_tight"
SOURCE_PULSE = "case_p19_pulse_phase23_tight"
PULSE_PREFIX = [["18[us]", 1], ["1[us]", 2], ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9]]
# The staged prefix ends at 20.030001 ms, i.e. 10.001 us after the 20.020 ms pulse.
POST_PREFIX_US = 10.001
# The steady circuit solve can converge slowly near the superconducting
# transition.  It is the restart source for every pulse, so give it enough
# iterations to reach the same fixed point while retaining the shorter cap
# for the much more numerous transient steps.
PILOT_STEADY_NONLINEAR_MAX_ITERATIONS = 60
PILOT_PULSE_NONLINEAR_MAX_ITERATIONS = 25
# The production/reference projects retain their 1e-8 nonlinear tolerance.
# Resolution pilots use 1e-7 only as a QoI pre-screen: the observed steady
# iteration noise floor is about 1e-7, so requiring 1e-8 rejects otherwise
# stable pilot waveforms without improving the resolution comparison.
PILOT_NONLINEAR_CONVERGENCE_TOLERANCE = 1e-7
PILOT_IMPLEMENTATION_SUFFIX = "alg5qg1"
STEADY_TEMPERATURE_LIMITS_K = (0.16, 0.18)
STEADY_CURRENT_LIMITS_A = (50e-6, 300e-6)
STEADY_MIN_RESISTANCE_OHM = 1e-4
PREFERRED_ELMER = ROOT.parent / "tools" / "elmer-hypre" / "install-phase13-step-commit" / "bin" / "ElmerSolver.exe"
PREFERRED_RUNTIME_BIN = Path(r"C:\msys64\ucrt64\bin")
LAYER_AXIS_TO_OPTION = {
    "stycast_layers": "stycast",
    "tes_layers": "tes",
    "sinx_layers": "sinx",
    "sio2_1_layers": "sio2-1",
    "si_1_layers": "si-1",
    "si_2_layers": "si-2",
}


def parse_args() -> argparse.Namespace:
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    parser = argparse.ArgumentParser(description=__doc__)
    selection = parser.add_mutually_exclusive_group(required=True)
    selection.add_argument("--axis", choices=[*config["spatial_axes"], "time"])
    selection.add_argument("--all", action="store_true", help="run every spatial axis, then the time-step axis")
    parser.add_argument("--levels", type=float, nargs="+", help="override the first-to-last candidates")
    parser.add_argument("--end-us", type=float, default=105.0, help="post-pulse comparison end time")
    parser.add_argument("--early-step-us", type=float, default=0.625)
    parser.add_argument("--tail-start-us", type=float, help="for --axis time: retain the finest step until this post-pulse time, then test coarser steps")
    parser.add_argument(
        "--fixed-setting", action="append", default=[], metavar="AXIS=VALUE",
        help="hold a previously accepted spatial setting fixed, e.g. stycast_layers=32; repeat as needed",
    )
    parser.add_argument(
        "--tag", default="", help="safe suffix used to isolate meshes and results for a chained pilot",
    )
    default_solver = "mumps" if PREFERRED_ELMER.is_file() else "direct"
    parser.add_argument(
        "--linear-system", choices=["direct", "iterative", "iterative_hypre_boomeramg", "mumps"], default=default_solver,
        help="Elmer linear solver; MUMPS is selected when the project-local MUMPS build is available",
    )
    parser.add_argument("--elmer-solver", type=Path, default=PREFERRED_ELMER if PREFERRED_ELMER.is_file() else None,
                        help="ElmerSolver executable; defaults to the project's proven MUMPS build when present")
    parser.add_argument("--runtime-bin", type=Path, default=PREFERRED_RUNTIME_BIN if PREFERRED_RUNTIME_BIN.is_dir() else None,
                        help="runtime DLL directory paired with --elmer-solver")
    parser.add_argument("--execute", action="store_true", help="run Elmer after preparing meshes and cases")
    parser.add_argument("--force", action="store_true", help="pass --force-deps to run.py")
    return parser.parse_args()


def slug(value: float) -> str:
    return (f"{value:g}").replace(".", "p").replace("-", "m")


def pilot_output_paths(axis: str, tag: str = "") -> tuple[Path, Path]:
    """Return versioned project/manifest paths that cannot replace old pilots."""
    tag_suffix = f"_{tag}" if tag else ""
    base = f"singlepixel{tag_suffix}_{axis}_{PILOT_IMPLEMENTATION_SUFFIX}_pilot"
    return OUT_DIR / f"{base}.json", OUT_DIR / f"{base}_manifest.json"


def parse_fixed_settings(values: list[str], config: dict) -> dict[str, float]:
    """Parse stable ``axis=value`` settings passed from the auto-tuner."""
    settings: dict[str, float] = {}
    for value in values:
        axis, separator, raw = value.partition("=")
        if not separator or axis not in config["spatial_axes"]:
            raise ValueError(f"--fixed-setting must name a spatial axis: {value!r}")
        number = float(raw)
        if number <= 0:
            raise ValueError(f"--fixed-setting must be positive: {value!r}")
        if axis in LAYER_AXIS_TO_OPTION and not number.is_integer():
            raise ValueError(f"{axis} must be an integer layer count: {value!r}")
        settings[axis] = number
    return settings


def fixed_generator_options(settings: dict[str, float]) -> list[str]:
    """Translate accepted settings into geometry-generator arguments."""
    options: list[str] = []
    for axis, value in sorted(settings.items()):
        if axis in LAYER_AXIS_TO_OPTION:
            options.extend([f"--{LAYER_AXIS_TO_OPTION[axis]}-layers", str(int(value))])
        elif axis == "stack_local_size_um":
            options.extend(["--stack-local-size", f"{value}e-6"])
        elif axis == "absorber_local_size_um":
            options.extend(["--absorber-local-size", f"{value}e-6"])
        elif axis == "global_mesh_size_um":
            options.extend(["--global-mesh-size", f"{value}e-6"])
        elif axis == "stack_local_half_width_um":
            options.extend(["--stack-local-half-width", f"{value}e-6"])
        elif axis == "absorber_local_radius_um":
            options.extend(["--absorber-local-radius", f"{value}e-6"])
        else:
            raise ValueError(axis)
    return options


def read_response(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    time_s = np.asarray([float(row["time_s"]) for row in rows])
    current_uA = np.asarray([float(row["tes_current_A"]) * 1e6 for row in rows])
    pre = (time_s >= 19.5e-3) & (time_s < 20.020e-3)
    if not np.any(pre):
        raise ValueError(f"No pre-pulse baseline in {path}")
    return (time_s - 20.020e-3) * 1e6, float(np.mean(current_uA[pre])) - current_uA


def compare_series(left: Path, right: Path, end_us: float, peak_uA: float) -> dict[str, float]:
    tx, yx = read_response(left)
    ty, yy = read_response(right)
    grid = np.linspace(0.0, end_us, 1001)
    difference = np.interp(grid, tx, yx) - np.interp(grid, ty, yy)
    i = int(np.argmax(np.abs(difference)))
    return {
        "max_abs_uA": float(abs(difference[i])),
        "time_us": float(grid[i]),
        "pct_comsol_peak": float(abs(difference[i]) / peak_uA * 100.0),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def mesh_complete(mesh: str, raw_mesh: Path | None = None) -> bool:
    """Return whether a converted mesh exists and its source MSH is safe.

    When ``raw_mesh`` is supplied, a converted mesh is never reused solely on
    the presence of mesh.header: the source MSH must pass the current quality
    contract first. Inspection failures are fatal to avoid overwriting or
    silently reusing an unverifiable mesh.
    """
    if not (ROOT / "work" / "meshes" / mesh / "mesh.header").is_file():
        return False
    if raw_mesh is None:
        return True
    try:
        quality = inspect_mesh_file(raw_mesh)
    except Exception as error:
        raise RuntimeError(f"[{mesh}] cannot inspect source mesh {raw_mesh}: {error}") from error
    reasons = quality["reasons"]
    if reasons:
        detail = "; ".join(str(reason) for reason in reasons)
        raise RuntimeError(f"[{mesh}] existing mesh rejected by quality gate ({raw_mesh}): {detail}")
    return True


def raw_mesh_from_recipe(commands: list[str]) -> Path:
    """Resolve the generator's --output path from a stored mesh recipe."""
    for command in commands:
        tokens = shlex.split(command, posix=False)
        if "generate_hybrid_prism_geometry.py" not in command or "--output" not in tokens:
            continue
        output = Path(tokens[tokens.index("--output") + 1].strip('"'))
        return output if output.is_absolute() else ROOT / output
    raise ValueError("mesh recipe has no generate_hybrid_prism_geometry.py --output command")


def _iteration_rows(case: str) -> list[dict[str, str]]:
    path = ROOT / "results" / case / f"{case}_iterations.csv"
    if not path.is_file():
        raise ValueError(f"{case}: missing iteration CSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{case}: iteration CSV is empty")
    required = {"time_step", "nonlinear_iter"}
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{case}: iteration CSV missing columns: {', '.join(sorted(missing))}")
    return rows


def nonlinear_iteration_cap(case: str) -> int:
    """Return the configured pilot cap for a generated steady or pulse case."""
    return (
        PILOT_STEADY_NONLINEAR_MAX_ITERATIONS
        if case.endswith("_steady")
        else PILOT_PULSE_NONLINEAR_MAX_ITERATIONS
    )


def iteration_failure_reason(case: str, max_iterations: int | None = None) -> str | None:
    """Return an explicit reason when any timestep exhausted that case's cap.

    ``max_iterations`` is optional for callers checking a nonstandard case;
    generated pilot names otherwise select the steady or pulse cap directly.
    """
    max_iterations = nonlinear_iteration_cap(case) if max_iterations is None else max_iterations
    try:
        rows = _iteration_rows(case)
        reached = [row for row in rows if int(row["nonlinear_iter"]) >= max_iterations]
    except (ValueError, KeyError) as error:
        return str(error)
    if reached:
        steps = sorted({row["time_step"] for row in reached}, key=int)
        return (
            f"{case}: nonlinear iteration cap {max_iterations} reached "
            f"at timestep(s) {', '.join(steps)}"
        )
    return None


def steady_state_failure_reason(case: str) -> str | None:
    """Validate the final steady circuit state recorded in its iteration CSV."""
    try:
        final = _iteration_rows(case)[-1]
        temperature = float(final["tes_temperature_K"])
        current = float(final.get("tes_current_A", final.get("raw_current_A", "")))
        resistance = float(final["tes_resistance_ohm"])
    except (ValueError, KeyError) as error:
        return f"{case}: cannot read final steady TES state: {error}"
    if not STEADY_TEMPERATURE_LIMITS_K[0] <= temperature <= STEADY_TEMPERATURE_LIMITS_K[1]:
        return f"{case}: final steady TES temperature {temperature:.6g} K is outside 0.16..0.18 K"
    if not STEADY_CURRENT_LIMITS_A[0] <= current <= STEADY_CURRENT_LIMITS_A[1]:
        return f"{case}: final steady TES current {current:.6g} A is outside 50..300 uA"
    if resistance <= STEADY_MIN_RESISTANCE_OHM:
        return f"{case}: final steady TES resistance {resistance:.6g} ohm is not above 1e-4 ohm"
    return None


def case_failure_reason(case: str, project_path: Path) -> str | None:
    """Return why a pilot pulse is not safe to reuse or compare.

    Checking the project hash prevents an old 105-us pulse from being reused
    accidentally after the user changes, for example, --end-us to 225 us.
    """
    result_dir = ROOT / "results" / case
    manifest_path = result_dir / "manifest.json"
    series = result_dir / f"{case}_series.csv"
    if not manifest_path.is_file() or not series.is_file():
        return f"{case}: missing manifest or pulse series"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return f"{case}: invalid manifest JSON"
    if manifest.get("exit_code") != 0:
        return f"{case}: solver exit code was {manifest.get('exit_code')!r}"
    if manifest.get("inputs_sha256", {}).get(project_path.name) != sha256(project_path):
        return f"{case}: result was made from a different pilot project"
    if reason := iteration_failure_reason(case):
        return reason
    steady_case = manifest.get("restart_from") or (case[:-6] + "steady" if case.endswith("_pulse") else None)
    if not steady_case:
        return f"{case}: no linked steady case to validate"
    if reason := iteration_failure_reason(steady_case):
        return reason
    return steady_state_failure_reason(steady_case)


def case_complete(case: str, project_path: Path) -> bool:
    return case_failure_reason(case, project_path) is None


def resolve_elmergrid() -> str:
    """Return an executable ElmerGrid command independent of shell PATH lookup."""
    executable = shutil.which("ElmerGrid") or shutil.which("ElmerGrid.exe")
    if executable is None:
        raise RuntimeError(
            "ElmerGrid executable was not found. Add its bin directory to PATH, "
            "or start this command from an Elmer environment."
        )
    return executable


def recipe_command(command: list[str]) -> str:
    """Serialize an argv list for the Windows ``shell=True`` recipe runner."""
    return subprocess.list2cmdline(command)


def mesh_command(mesh_file: str, mesh_name: str, options: list[str]) -> tuple[list[str], list[str]]:
    return [
        sys.executable, "generate_hybrid_prism_geometry.py", "elmer_project_comsol_timegrid.json",
        "--output", f"gmsh/{mesh_file}", "--stack-local-size", "16.6666666666667e-6",
        "--absorber-local-size", "35e-6", "--absorber-local-radius", "50e-6",
        "--disable-mesh-size-extend-from-boundary", "--mesh-algorithm", "5", *options,
    ], [resolve_elmergrid(), "14", "2", f"gmsh/{mesh_file}", "-merge", "1e-10", "-out", f"work/meshes/{mesh_name}"]


def step_count(start_us: float, end_us: float, step_us: float) -> int:
    count = round((end_us - start_us) / step_us)
    if count < 1 or not np.isclose(start_us + count * step_us, end_us, atol=0.01):
        raise ValueError(f"interval {start_us:g}..{end_us:g} us is not compatible with {step_us:g} us")
    return count


def build_project(axis: str, levels: list[float], step_us: float, end_us: float, tail_start_us: float | None, linear_system: str, output: Path, fixed_settings: dict[str, float] | None = None, tag: str = "") -> tuple[dict, list[tuple[str, str]]]:
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    fixed_settings = fixed_settings or {}
    if axis in fixed_settings:
        raise ValueError(f"{axis} cannot be both fixed and varied")
    if tag and not re.fullmatch(r"[A-Za-z0-9_]+", tag):
        raise ValueError("--tag may contain only letters, numbers, and underscores")
    cases: list[tuple[str, str]] = []
    for level in levels:
        # Namespace every generated artifact (raw MSH, converted mesh, cases,
        # results and restart state) away from pre-quality-gate pilot runs.
        # The "pilot_" literal is dropped (redundant with the "auto_" tag and
        # this directory's own name) to keep the identifier well under the
        # TES UDF's 128-character state-file limit for the longest axis name
        # (absorber_local_size_um) chained with a full autotune tag.
        ident = f"{tag + '_' if tag else 'pilot_'}{axis}_{slug(level)}_{PILOT_IMPLEMENTATION_SUFFIX}"
        mesh = f"mesh_{ident}"
        state_file = f"work/meshes/{mesh}/{ident}.state"
        if len(state_file) > TES_STATE_FILE_MAX_LEN:
            raise ValueError(
                f"[{axis}={level}] generated TES state file would be {len(state_file)} "
                f"characters, over the UDF's {TES_STATE_FILE_MAX_LEN}-character limit: "
                f"{state_file!r}. Shorten --tag or the axis/level naming before running "
                "mesh generation."
            )
        mesh_file = f"{ident}.msh"
        options: list[str] = []
        if axis in LAYER_AXIS_TO_OPTION:
            options = [f"--{LAYER_AXIS_TO_OPTION[axis]}-layers", str(int(level))]
        elif axis == "stack_local_size_um":
            options = ["--stack-local-size", f"{level}e-6"]
        elif axis == "absorber_local_size_um":
            options = ["--absorber-local-size", f"{level}e-6"]
        elif axis == "global_mesh_size_um":
            options = ["--global-mesh-size", f"{level}e-6"]
        elif axis == "stack_local_half_width_um":
            options = ["--stack-local-half-width", f"{level}e-6"]
        elif axis == "absorber_local_radius_um":
            options = ["--absorber-local-radius", f"{level}e-6"]
        elif axis == "time":
            options = ["--stycast-layers", "32"]
        else:
            raise ValueError(axis)
        options = [*fixed_generator_options(fixed_settings), *options]
        gmsh, grid = mesh_command(mesh_file, mesh, options)
        project["meshes"][mesh] = {
            "geometry": "single_pixel", "dir": mesh,
            "recipe": {
                "generator": "generate_hybrid_prism_geometry.py",
                "commands": [recipe_command(gmsh), recipe_command(grid)],
            },
            "notes": f"Automatic one-axis convergence pilot: {axis}={level}; fixed={fixed_settings}.",
        }
        steady_name = f"case_{ident}_steady"
        steady = copy.deepcopy(project["cases"][SOURCE_STEADY])
        steady.update({"mesh": mesh, "state_file": f"work/meshes/{mesh}/{ident}.state",
                       "output_file_path": f"../work/meshes/{mesh}/{steady_name}.result",
                       "series_file": f"{steady_name}_series.csv", "iteration_series_file": f"{steady_name}_iterations.csv"})
        steady["solver"] = {
            **steady.get("solver", {}),
            "linear_system": linear_system,
            "nonlinear_max_iterations": PILOT_STEADY_NONLINEAR_MAX_ITERATIONS,
            "nonlinear_convergence_tolerance": PILOT_NONLINEAR_CONVERGENCE_TOLERANCE,
        }
        project["cases"][steady_name] = steady
        pulse_name = f"case_{ident}_pulse"
        case_step_us = level if axis == "time" else step_us
        if axis == "time" and tail_start_us is not None and level != levels[0]:
            timesteps = [*PULSE_PREFIX, [f"{levels[0]:g}[us]", step_count(POST_PREFIX_US, tail_start_us, levels[0])],
                         [f"{case_step_us:g}[us]", step_count(tail_start_us, end_us, case_step_us)]]
        else:
            timesteps = [*PULSE_PREFIX, [f"{case_step_us:g}[us]", step_count(POST_PREFIX_US, end_us, case_step_us)]]
        pulse = copy.deepcopy(project["cases"][SOURCE_PULSE])
        pulse.update({"mesh": mesh, "restart_from": steady_name,
                      "restart_file_path": f"../work/meshes/{mesh}/{steady_name}.result",
                      "state_file": f"work/meshes/{mesh}/{ident}.state",
                      "series_file": f"{pulse_name}_series.csv", "iteration_series_file": f"{pulse_name}_iterations.csv",
                      "timesteps": timesteps,
                      "output_intervals": [1] * len(timesteps),
                      "comparison_time_grid": {"mode": "automatic one-axis pilot", "axis": axis, "value": level, "end_after_pulse": f"{end_us:g}[us]"}})
        pulse["solver"] = {
            **pulse.get("solver", {}),
            "linear_system": linear_system,
            "nonlinear_max_iterations": PILOT_PULSE_NONLINEAR_MAX_ITERATIONS,
            "nonlinear_convergence_tolerance": PILOT_NONLINEAR_CONVERGENCE_TOLERANCE,
        }
        project["cases"][pulse_name] = pulse
        cases.append((pulse_name, mesh))
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    return project, cases


def main() -> None:
    args = parse_args()
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    if args.all:
        if args.fixed_setting or args.tag:
            raise ValueError("--all does not support --fixed-setting/--tag; use run_singlepixel_resolution_autotune.py")
        if args.levels:
            raise ValueError("--levels is only valid with one --axis")
        for axis in [*config["spatial_axes"], "time"]:
            command = [sys.executable, str(Path(__file__).resolve()), "--axis", axis,
                       "--end-us", f"{args.end_us:g}", "--early-step-us", f"{args.early_step_us:g}"]
            command.extend(["--linear-system", args.linear_system])
            if args.elmer_solver is not None:
                command.extend(["--elmer-solver", str(args.elmer_solver)])
            if args.runtime_bin is not None:
                command.extend(["--runtime-bin", str(args.runtime_bin)])
            if args.tail_start_us is not None:
                command.extend(["--tail-start-us", f"{args.tail_start_us:g}"])
            if args.execute:
                command.append("--execute")
            if args.force:
                command.append("--force")
            subprocess.run(command, cwd=ROOT, check=True)
        return
    levels = args.levels or (config["time_step_candidates_us"] if args.axis == "time" else config["spatial_axes"][args.axis])
    if len(levels) < 2:
        raise ValueError("At least two candidate levels are required")
    levels = [float(value) for value in levels]
    # For time, each level is itself the post-pulse step; for spatial axes it is fixed at early-step-us.
    step_us = args.early_step_us if args.axis != "time" else levels[0]
    fixed_settings = parse_fixed_settings(args.fixed_setting, config)
    project_path, manifest_path = pilot_output_paths(args.axis, args.tag)
    _, cases = build_project(args.axis, levels, step_us, args.end_us, args.tail_start_us, args.linear_system, project_path, fixed_settings, args.tag)
    manifest = {"axis": args.axis, "levels": levels, "tail_start_us": args.tail_start_us, "fixed_settings": fixed_settings, "tag": args.tag or None, "linear_system": args.linear_system, "elmer_solver": str(args.elmer_solver) if args.elmer_solver else "run.py default", "runtime_bin": str(args.runtime_bin) if args.runtime_bin else None, "project": str(project_path.relative_to(ROOT)), "execute": args.execute, "cases": [case for case, _ in cases]}

    if args.execute:
        # Meshes are generated before sync_elmer_parameters.py, because SIF generation
        # needs the converted mesh to compute the discrete absorber pulse normalisation.
        for _, mesh in cases:
            recipe = json.loads(project_path.read_text(encoding="utf-8"))["meshes"][mesh]["recipe"]["commands"]
            raw_mesh = raw_mesh_from_recipe(recipe)
            if mesh_complete(mesh, raw_mesh):
                print(f"[{mesh}] converted mesh exists - skipping generation")
                continue
            for command in recipe:
                subprocess.run(command, cwd=ROOT, shell=True, check=True)
        # Execute only as far as the first accepted adjacent pair.  Meshes are
        # cheap compared with the transient solves; the unused later meshes
        # make the next refinement immediately runnable if it is required.
        previous: tuple[str, str] | None = None
        comparisons = []
        peak = float(json.loads((ROOT / "artifacts/comparison/stycast_z16_resolution/metrics.json").read_text(encoding="utf-8"))["comsol_peak_uA"])
        for current in cases:
            pulse, _ = current
            if case_complete(pulse, project_path):
                print(f"[{pulse}] successful matching result exists - skipping pulse")
            else:
                command = [sys.executable, "run.py", pulse, "--project", str(project_path)]
                if args.elmer_solver is not None:
                    command.extend(["--elmer-solver", str(args.elmer_solver)])
                if args.runtime_bin is not None:
                    command.extend(["--runtime-bin", str(args.runtime_bin)])
                if args.force:
                    command.append("--force-deps")
                subprocess.run(command, cwd=ROOT, check=True)
            if reason := case_failure_reason(pulse, project_path):
                raise RuntimeError(f"[{pulse}] pilot result rejected; comparison skipped: {reason}")
            if previous is not None:
                left, _ = previous
                result = compare_series(ROOT / "results" / left / f"{left}_series.csv", ROOT / "results" / pulse / f"{pulse}_series.csv", args.end_us, peak)
                result.update({"left": left, "right": pulse, "accepted": result["pct_comsol_peak"] <= config["waveform_change_limit_pct_of_peak"]})
                comparisons.append(result)
                if result["accepted"]:
                    break
            previous = current
        manifest["comparisons"] = comparisons
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(project_path)
    print(manifest_path)


if __name__ == "__main__":
    main()
