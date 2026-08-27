"""Prepare, run, and score PoST dual-TES multivariate searches.

Unlike ``single_pixel_search.py``, this runner compares the PoST geometry
directly with paired experimental CH0/CH1 voltage pulses.  It keeps the
single-pixel scan as a model-sensitivity study and uses its selected shared
material parameters only as bounded priors for the PoST fit.

Typical use from ``Elmer-Projects``::

    python scripts/search/post_multivariate_search.py prepare-references
    python scripts/search/post_multivariate_search.py prepare
    python scripts/search/post_multivariate_search.py dry-run mv_000_<hash>
    python scripts/search/post_multivariate_search.py run-all --limit 3
    python scripts/search/post_multivariate_search.py run-all --skip-scored

The experimental reference builder uses the selected-event peak sharing
between CH0 and CH1 to form five ordered average pulses.  The corresponding
Elmer cases cover left edge, left middle, centre, and their mirrored partners.
Both simulated TES currents receive the same experimental Bessel low-pass
before scoring.  A target pair is normalized by the sum of its two channel
peaks, preserving left/right energy sharing while removing the unknown global
current-to-voltage gain.  One common bounded time shift is fitted across every
position and both channels.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy import signal as scipy_signal


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "post_multivariate_search_config.json"


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    factors: dict[str, float]


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def resolve_from_root(value: str | Path) -> Path:
    expanded = os.path.expandvars(str(value))
    path = Path(expanded)
    return path if path.is_absolute() else ROOT / path


def resolve_from_experiment(config: dict[str, Any], value: str | Path) -> Path:
    path = Path(os.path.expandvars(str(value)))
    if path.is_absolute():
        return path
    return resolve_from_root(config["experiment"]["root"]) / path


def nested_get(data: dict[str, Any], path: Iterable[str]) -> Any:
    node: Any = data
    for key in path:
        node = node[key]
    return node


def nested_set(data: dict[str, Any], path: Iterable[str], value: Any) -> None:
    keys = list(path)
    node: Any = data
    for key in keys[:-1]:
        node = node[key]
    node[keys[-1]] = value


def scaled_expression(expression: str, factor: float) -> str:
    if factor <= 0:
        raise ValueError(f"multiplicative factor must be positive, got {factor}")
    return f"({expression})*{factor:.12g}"


def short_hash(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:8]


def latin_hypercube_log_factors(
    config: dict[str, Any],
    *,
    sample_count: int | None = None,
    seed: int | None = None,
) -> list[dict[str, float]]:
    settings = config["multivariate_search"]
    variable_names = list(settings["variables"])
    count = int(settings["samples"] if sample_count is None else sample_count)
    random_seed = int(settings["seed"] if seed is None else seed)
    if count <= 0:
        raise ValueError("multivariate sample count must be positive")
    if not variable_names:
        raise ValueError("multivariate search requires at least one variable")

    rng = np.random.default_rng(random_seed)
    unit = np.empty((count, len(variable_names)), dtype=float)
    for column in range(len(variable_names)):
        strata = (np.arange(count, dtype=float) + rng.random(count)) / count
        rng.shuffle(strata)
        unit[:, column] = strata

    rows: list[dict[str, float]] = []
    for sample in unit:
        factors: dict[str, float] = {}
        for index, variable_name in enumerate(variable_names):
            spec = config["variables"][variable_name]
            lower, upper = (float(value) for value in spec["search_range"])
            if not 0 < lower < upper:
                raise ValueError(
                    f"invalid positive search_range for {variable_name}: {lower}, {upper}"
                )
            factors[variable_name] = math.exp(
                math.log(lower)
                + sample[index] * (math.log(upper) - math.log(lower))
            )
        rows.append(factors)
    return rows


def candidates_from_config(
    config: dict[str, Any],
    *,
    sample_count: int | None = None,
    seed: int | None = None,
) -> list[Candidate]:
    return [
        Candidate(f"mv_{index:03d}_{short_hash(factors)}", factors)
        for index, factors in enumerate(
            latin_hypercube_log_factors(
                config, sample_count=sample_count, seed=seed
            )
        )
    ]


def output_root(config: dict[str, Any]) -> Path:
    return resolve_from_root(config["output_dir"])


def candidate_dir(config: dict[str, Any], candidate_id: str) -> Path:
    return output_root(config) / "candidates" / candidate_id


def reference_path(config: dict[str, Any], target: dict[str, Any]) -> Path:
    reference_root = config.get("reference_output_dir", config["output_dir"])
    return resolve_from_root(reference_root) / "references" / f"{target['name']}.csv"


def side_series_file(series_file: str, side: str) -> str:
    suffix = "_series.csv"
    if series_file.endswith(suffix):
        return f"{series_file[:-len(suffix)]}_{side}{suffix}"
    stem, dot, extension = series_file.rpartition(".")
    return f"{stem}_{side}.{extension}" if dot else f"{series_file}_{side}"


def case_stem(candidate: Candidate, trial_index: int) -> str:
    return f"postsearch_{candidate.candidate_id}_{trial_index:02d}"


def apply_case_overrides(case: dict[str, Any], overrides: dict[str, Any]) -> None:
    for key, value in overrides.items():
        case[key] = copy.deepcopy(value)


def mutate_project(
    base_project: dict[str, Any],
    config: dict[str, Any],
    candidate: Candidate,
    g0_factor: float,
    trial_index: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    project = copy.deepcopy(base_project)
    for variable_name, factor in candidate.factors.items():
        spec = config["variables"][variable_name]
        path = list(spec["path"])
        original = str(nested_get(project, path))
        nested_set(project, path, scaled_expression(original, factor))

    original_g0 = str(project["parameter_expressions"]["G0"])
    project["parameter_expressions"]["G0"] = scaled_expression(
        original_g0, g0_factor
    )

    source_cases = base_project["cases"]
    stem = case_stem(candidate, trial_index)
    steady_name = f"{stem}_steady"
    steady_source = str(config["cases"]["steady"])
    if steady_source not in source_cases:
        raise KeyError(f"steady source case does not exist: {steady_source}")
    steady_case = copy.deepcopy(source_cases[steady_source])
    steady_case["series_file"] = f"{steady_name}_series.csv"
    steady_case["output_result"] = True
    steady_case["vtu"] = False
    steady_case.pop("output_file_path", None)
    apply_case_overrides(steady_case, config["cases"].get("steady_overrides", {}))

    mesh_name = str(steady_case["mesh"])
    mesh_dir = str(project["meshes"][mesh_name]["dir"])
    steady_result_path = f"../work/meshes/{mesh_dir}/{steady_name}.result"
    steady_case.setdefault("output_file_path", steady_result_path)

    pulse_cases: dict[str, dict[str, str]] = {}
    cases: dict[str, Any] = {steady_name: steady_case}
    # Two targets sharing a source_case (e.g. a mirrored "swap_simulation_channels"
    # pair) describe the same physical Elmer run — solve it once and let both
    # targets reference the same case/series files instead of re-running Elmer.
    source_case_to_pulse_name: dict[str, str] = {}
    for target in config["targets"]:
        target_name = str(target["name"])
        source_name = str(target["source_case"])
        if source_name not in source_cases:
            raise KeyError(f"pulse source case does not exist: {source_name}")
        if source_name in source_case_to_pulse_name:
            pulse_name = source_case_to_pulse_name[source_name]
        else:
            pulse_name = f"{stem}_{target_name}"
            pulse_case = copy.deepcopy(source_cases[source_name])
            pulse_case["restart_from"] = steady_name
            pulse_case.pop("restart_file_base", None)
            pulse_case.pop("restart_file_path", None)
            pulse_case["series_file"] = f"{pulse_name}_series.csv"
            pulse_case["vtu"] = False
            pulse_overrides = copy.deepcopy(config["cases"].get("pulse_overrides", {}))
            solver_overrides = pulse_overrides.pop("solver_overrides", {})
            apply_case_overrides(pulse_case, pulse_overrides)
            if solver_overrides:
                pulse_case["solver"] = copy.deepcopy(pulse_case.get("solver", {}))
                pulse_case["solver"].update(copy.deepcopy(solver_overrides))
            pulse_case.setdefault("restart_file_path", steady_result_path)
            cases[pulse_name] = pulse_case
            source_case_to_pulse_name[source_name] = pulse_name
        pulse_cases[target_name] = {
            "case": pulse_name,
            "series_file": cases[pulse_name]["series_file"],
            "source_case": source_name,
        }

    project["cases"] = cases
    metadata = {
        "candidate_id": candidate.candidate_id,
        "candidate_factors": candidate.factors,
        "g0_factor": g0_factor,
        "trial_index": trial_index,
        "steady_case": steady_name,
        "steady_series": steady_case["series_file"],
        "mesh": mesh_name,
        "mesh_dir": mesh_dir,
        "state_files": {
            side: f"work/meshes/{mesh_dir}/{steady_name}_{side}.state"
            for side in ("L", "R")
        },
        "pulse_cases": pulse_cases,
        "project_hash": short_hash(
            {"factors": candidate.factors, "g0_factor": g0_factor}
        ),
    }
    return project, metadata


def trial_paths(
    config: dict[str, Any], candidate: Candidate, trial_index: int
) -> tuple[Path, Path]:
    directory = candidate_dir(config, candidate.candidate_id) / f"trial_{trial_index:02d}"
    return directory / "project.json", directory / "metadata.json"


def write_trial(
    base_project: dict[str, Any],
    config: dict[str, Any],
    candidate: Candidate,
    g0_factor: float,
    trial_index: int,
) -> tuple[Path, dict[str, Any]]:
    project, metadata = mutate_project(
        base_project, config, candidate, g0_factor, trial_index
    )
    project_path, metadata_path = trial_paths(config, candidate, trial_index)
    write_json(project_path, project)
    write_json(metadata_path, metadata)
    return project_path, metadata


def run_command(
    project_path: Path,
    case_name: str,
    *,
    dry_run: bool = False,
    force_deps: bool = False,
) -> None:
    command = [
        sys.executable,
        str(ROOT / "run.py"),
        case_name,
        "--project",
        str(project_path),
        "--mpi-procs",
        "1",
    ]
    if dry_run:
        command.append("--dry-run")
    if force_deps:
        command.append("--force-deps")
    print("$ " + subprocess.list2cmdline(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def read_state(path: Path) -> dict[str, float]:
    values = np.fromstring(path.read_text(encoding="utf-8"), sep=" ")
    if len(values) != 5 or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid TES state file {path}: expected five finite values")
    return {
        "temperature_K": float(values[0]),
        "current_A": float(values[1]),
        "resistance_ohm": float(values[2]),
        "power_W": float(values[3]),
        "previous_current_A": float(values[4]),
    }


def target_temperature(project: dict[str, Any]) -> float:
    sys.path.insert(0, str(ROOT))
    from scripts.support.reconcile_project import reconcile_project

    return float(reconcile_project(project)["parameters"]["T_0"])


def evaluate_steady_trial(
    base_project: dict[str, Any],
    config: dict[str, Any],
    candidate: Candidate,
    g0_factor: float,
    trial_index: int,
    *,
    dry_run: bool = False,
) -> tuple[float, Path, dict[str, Any], dict[str, dict[str, float]] | None]:
    project_path, metadata = write_trial(
        base_project, config, candidate, g0_factor, trial_index
    )
    run_command(project_path, metadata["steady_case"], dry_run=dry_run)
    if dry_run:
        return math.nan, project_path, metadata, None

    states = {
        side: read_state(ROOT / relative_path)
        for side, relative_path in metadata["state_files"].items()
    }
    target = target_temperature(load_json(project_path))
    temperatures = np.asarray(
        [states["L"]["temperature_K"], states["R"]["temperature_K"]]
    )
    metric = str(config["g0_calibration"].get("temperature_metric", "mean"))
    if metric == "mean":
        fitted_temperature = float(np.mean(temperatures))
    elif metric == "max":
        fitted_temperature = float(np.max(temperatures))
    else:
        raise ValueError("g0_calibration.temperature_metric must be mean or max")
    error_K = fitted_temperature - target
    write_json(
        project_path.parent / "steady_result.json",
        {
            **metadata,
            "target_temperature_K": target,
            "fitted_temperature_K": fitted_temperature,
            "temperature_error_mK": error_K * 1e3,
            "states": states,
        },
    )
    return error_K, project_path, metadata, states


def calibrate_g0(
    base_project: dict[str, Any],
    config: dict[str, Any],
    candidate: Candidate,
) -> tuple[Path, dict[str, Any]]:
    settings = config["g0_calibration"]
    tolerance_K = float(settings["temperature_tolerance_mK"]) * 1e-3
    initial = float(settings["initial_factor"])
    lower = float(settings["lower_factor"])
    upper = float(settings["upper_factor"])
    max_iterations = int(settings["max_iterations"])

    def accept(
        error: float, project_path: Path, metadata: dict[str, Any]
    ) -> tuple[Path, dict[str, Any]]:
        write_json(
            candidate_dir(config, candidate.candidate_id) / "calibration.json",
            {
                "candidate_id": candidate.candidate_id,
                "g0_factor": metadata["g0_factor"],
                "temperature_error_mK": error * 1e3,
                "project": str(project_path.relative_to(ROOT)),
                "steady_case": metadata["steady_case"],
            },
        )
        return project_path, metadata

    trial_index = 0
    points: list[tuple[float, float, Path, dict[str, Any]]] = []
    for factor in (initial, lower, upper):
        if any(math.isclose(factor, existing[0]) for existing in points):
            continue
        error, project_path, metadata, _ = evaluate_steady_trial(
            base_project, config, candidate, factor, trial_index
        )
        trial_index += 1
        points.append((factor, error, project_path, metadata))
        if abs(error) <= tolerance_K:
            return accept(error, project_path, metadata)

    points.sort(key=lambda item: item[0])
    bracket: tuple[
        tuple[float, float, Path, dict[str, Any]],
        tuple[float, float, Path, dict[str, Any]],
    ] | None = None
    for left, right in zip(points, points[1:]):
        if left[1] * right[1] <= 0:
            bracket = (left, right)
            break
    if bracket is None:
        details = ", ".join(
            f"factor={factor:.6g} error={error*1e3:+.6f} mK"
            for factor, error, _, _ in points
        )
        raise RuntimeError(f"G0 calibration did not bracket T0: {details}")

    left, right = bracket
    best = min(points, key=lambda item: abs(item[1]))
    while trial_index < max_iterations:
        lower_factor, lower_error = left[0], left[1]
        upper_factor, upper_error = right[0], right[1]
        log_lower = math.log(lower_factor)
        log_upper = math.log(upper_factor)
        denominator = upper_error - lower_error
        if denominator == 0:
            log_middle = 0.5 * (log_lower + log_upper)
        else:
            log_middle = log_lower - lower_error * (
                log_upper - log_lower
            ) / denominator
            guard = 0.12 * (log_upper - log_lower)
            log_middle = min(max(log_middle, log_lower + guard), log_upper - guard)
        middle_factor = math.exp(log_middle)
        error, project_path, metadata, _ = evaluate_steady_trial(
            base_project, config, candidate, middle_factor, trial_index
        )
        trial_index += 1
        middle = (middle_factor, error, project_path, metadata)
        if abs(error) < abs(best[1]):
            best = middle
        if abs(error) <= tolerance_K:
            return accept(error, project_path, metadata)
        if left[1] * error <= 0:
            right = middle
        else:
            left = middle

    if abs(best[1]) <= tolerance_K:
        return accept(best[1], best[2], best[3])
    raise RuntimeError(
        f"G0 calibration failed to reach {tolerance_K*1e3:.4f} mK; "
        f"best error was {abs(best[1])*1e3:.6f} mK"
    )


def read_peak_summary(path: Path) -> dict[int, float]:
    peaks: dict[int, float] = {}
    with path.open("r", newline="", encoding="utf-8-sig") as file:
        for row in csv.DictReader(file):
            key = int(float(row["key"]))
            peak = float(row["Peak"])
            if math.isfinite(peak) and peak > 0:
                peaks[key] = peak
    return peaks


def load_selected_keys(path: Path) -> list[int]:
    values = np.atleast_1d(np.loadtxt(path, comments="#"))
    return [int(round(float(value))) for value in values]


def select_reference_keys(
    selected_keys: list[int],
    left_peaks: dict[int, float],
    right_peaks: dict[int, float],
    targets: list[dict[str, Any]],
    *,
    left_gain: float = 1.0,
    right_gain: float = 1.0,
) -> dict[str, dict[str, Any]]:
    if left_gain <= 0 or right_gain <= 0:
        raise ValueError("channel peak gains must be positive")
    records: list[tuple[int, float, float]] = []
    for key in selected_keys:
        if key not in left_peaks or key not in right_peaks:
            continue
        corrected_left = left_peaks[key] / left_gain
        corrected_right = right_peaks[key] / right_gain
        total = corrected_left + corrected_right
        if total > 0:
            records.append((key, corrected_left, corrected_left / total))
    if not records:
        raise RuntimeError("no selected events have finite positive peaks in both channels")

    fractions = np.asarray([fraction for _, _, fraction in records])
    selected: dict[str, dict[str, Any]] = {}
    for target in targets:
        selector = target["selector"]

        if "extremum" in selector:
            # Fitting_LSE.py-style selection: rank by the raw CH0 (left) peak
            # alone, not by the left/right fraction, and take the single
            # most extreme event (count defaults to 1).
            extremum = str(selector["extremum"])
            count = int(selector.get("count", 1))
            if count <= 0:
                raise ValueError(f"target {target['name']}: selector.count must be positive")
            if extremum == "max":
                ranked = sorted(records, key=lambda item: item[1], reverse=True)[:count]
            elif extremum == "min":
                ranked = sorted(records, key=lambda item: item[1])[:count]
            else:
                raise ValueError(f"target {target['name']}: selector.extremum must be max or min")
            selected[str(target["name"])] = {
                "selection": f"left_peak_{extremum}",
                "keys": [key for key, _, _ in ranked],
                "mean_left_peak": float(np.mean([peak for _, peak, _ in ranked])),
                "mean_left_fraction": float(np.mean([fraction for _, _, fraction in ranked])),
                "min_left_fraction": float(min(fraction for _, _, fraction in ranked)),
                "max_left_fraction": float(max(fraction for _, _, fraction in ranked)),
            }
            continue

        count = int(selector["count"])
        if count <= 0:
            raise ValueError(f"target {target['name']}: selector.count must be positive")
        if "fraction_quantile" in selector:
            quantile = float(selector["fraction_quantile"])
            if not 0 <= quantile <= 1:
                raise ValueError(
                    f"target {target['name']}: fraction_quantile must be in [0, 1]"
                )
            desired = float(np.quantile(fractions, quantile))
        elif "fraction" in selector:
            desired = float(selector["fraction"])
        else:
            raise ValueError(
                f"target {target['name']}: selector needs 'extremum', 'fraction_quantile', or 'fraction'"
            )
        ranked = sorted(records, key=lambda item: abs(item[2] - desired))[:count]
        selected[str(target["name"])] = {
            "selection": "fraction_quantile",
            "desired_left_fraction": desired,
            "keys": [key for key, _, _ in ranked],
            "mean_left_fraction": float(np.mean([fraction for _, _, fraction in ranked])),
            "min_left_fraction": float(min(fraction for _, _, fraction in ranked)),
            "max_left_fraction": float(max(fraction for _, _, fraction in ranked)),
        }
    return selected


def load_binary_pulse(path: Path, *, skip_bytes: int = 4) -> np.ndarray:
    with path.open("rb") as file:
        file.seek(skip_bytes)
        values = np.frombuffer(file.read(), dtype="float64").copy()
    if len(values) == 0 or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid pulse data: {path}")
    return values


def bessel_lowpass(values: np.ndarray, rate_hz: float, cutoff_hz: float, order: int) -> np.ndarray:
    if not 0 < cutoff_hz < 0.5 * rate_hz:
        raise ValueError("Bessel cutoff must be between zero and Nyquist")
    sos = scipy_signal.bessel(
        order, cutoff_hz / (0.5 * rate_hz), btype="low", output="sos", norm="phase"
    )
    return scipy_signal.sosfiltfilt(sos, values)


def experimental_response(
    raw: np.ndarray,
    *,
    rate_hz: float,
    cutoff_hz: float,
    filter_order: int,
    presample: int,
    direction: str,
    channel_gain: float = 1.0,
) -> np.ndarray:
    if channel_gain <= 0:
        raise ValueError("experimental channel_gain must be positive")
    filtered = bessel_lowpass(raw, rate_hz, cutoff_hz, filter_order)
    baseline = float(np.mean(filtered[:presample]))
    if direction == "rise":
        return (filtered - baseline) / channel_gain
    if direction == "drop":
        return (baseline - filtered) / channel_gain
    raise ValueError("experiment.response_direction must be rise or drop")


def prepare_references(config: dict[str, Any]) -> dict[str, Path]:
    experiment = config["experiment"]
    pulse_config = load_json(resolve_from_experiment(config, experiment["pulse_config"]))
    rate_hz = float(pulse_config["Readout"]["Rate"])
    sample_count = int(pulse_config["Readout"]["Sample"])
    presample = int(pulse_config["Readout"]["PreSample"])
    cutoff_hz = float(
        experiment.get(
            "cutoff_hz", pulse_config["Analysis"]["CutoffFrequency"]
        )
    )
    filter_order = int(experiment.get("filter_order", 2))
    direction = str(experiment.get("response_direction", "rise"))
    skip_bytes = int(experiment.get("binary_skip_bytes", 4))

    selected_keys = load_selected_keys(
        resolve_from_experiment(config, experiment["selected_keys"])
    )
    left_peaks = read_peak_summary(
        resolve_from_experiment(config, experiment["left_summary"])
    )
    right_peaks = read_peak_summary(
        resolve_from_experiment(config, experiment["right_summary"])
    )
    common_keys = [
        key for key in selected_keys if key in left_peaks and key in right_peaks
    ]
    equalization = str(
        experiment.get("channel_gain_equalization", "selected_median_peak")
    )
    if equalization == "selected_median_peak":
        left_gain = float(np.median([left_peaks[key] for key in common_keys]))
        right_gain = float(np.median([right_peaks[key] for key in common_keys]))
    elif equalization == "none":
        left_gain = right_gain = 1.0
    else:
        raise ValueError(
            "experiment.channel_gain_equalization must be "
            "selected_median_peak or none"
        )
    selections = select_reference_keys(
        selected_keys,
        left_peaks,
        right_peaks,
        config["targets"],
        left_gain=left_gain,
        right_gain=right_gain,
    )

    time_ms = (np.arange(sample_count, dtype=float) - presample) / rate_hz * 1e3
    window_start, window_end = (
        float(value) for value in experiment.get("reference_window_ms", [-20.0, 100.0])
    )
    window_mask = (time_ms >= window_start) & (time_ms <= window_end)
    outputs: dict[str, Path] = {}
    manifest_targets: dict[str, Any] = {}
    for target in config["targets"]:
        name = str(target["name"])
        keys = selections[name]["keys"]
        left_sum = np.zeros(sample_count, dtype=float)
        right_sum = np.zeros(sample_count, dtype=float)
        for key in keys:
            left_path = resolve_from_experiment(
                config, str(experiment["left_pulse_pattern"]).format(key=key)
            )
            right_path = resolve_from_experiment(
                config, str(experiment["right_pulse_pattern"]).format(key=key)
            )
            left_raw = load_binary_pulse(left_path, skip_bytes=skip_bytes)
            right_raw = load_binary_pulse(right_path, skip_bytes=skip_bytes)
            if len(left_raw) != sample_count or len(right_raw) != sample_count:
                raise ValueError(
                    f"pulse {key} length mismatch: "
                    f"left={len(left_raw)}, right={len(right_raw)}, expected={sample_count}"
                )
            left_sum += experimental_response(
                left_raw,
                rate_hz=rate_hz,
                cutoff_hz=cutoff_hz,
                filter_order=filter_order,
                presample=presample,
                direction=direction,
                channel_gain=left_gain,
            )
            right_sum += experimental_response(
                right_raw,
                rate_hz=rate_hz,
                cutoff_hz=cutoff_hz,
                filter_order=filter_order,
                presample=presample,
                direction=direction,
                channel_gain=right_gain,
            )
        left_average = left_sum / len(keys)
        right_average = right_sum / len(keys)
        path = reference_path(config, target)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file)
            writer.writerow(["time_ms", "left_signal", "right_signal"])
            for t, left, right in zip(
                time_ms[window_mask],
                left_average[window_mask],
                right_average[window_mask],
            ):
                writer.writerow([f"{t:.12g}", f"{left:.12g}", f"{right:.12g}"])
        outputs[name] = path
        manifest_targets[name] = {
            **selections[name],
            "reference": str(path.relative_to(ROOT)),
            "source_case": target["source_case"],
            "swap_simulation_channels": bool(
                target.get("swap_simulation_channels", False)
            ),
        }

    manifest_path = output_root(config) / "reference_manifest.json"
    write_json(
        manifest_path,
        {
            "experiment_root": str(resolve_from_root(experiment["root"])),
            "rate_hz": rate_hz,
            "sample_count": sample_count,
            "presample": presample,
            "cutoff_hz": cutoff_hz,
            "filter_order": filter_order,
            "channel_gain_equalization": equalization,
            "channel_peak_gains": {"left": left_gain, "right": right_gain},
            "targets": manifest_targets,
        },
    )
    outputs["manifest"] = manifest_path
    return outputs


def load_reference_pair(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return (
        np.asarray(data["time_ms"], dtype=float),
        np.asarray(data["left_signal"], dtype=float),
        np.asarray(data["right_signal"], dtype=float),
    )


def load_simulation_current(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return (
        np.asarray(data["time_s"], dtype=float) * 1e3,
        np.asarray(data["tes_current_A"], dtype=float),
    )


def baseline_correct_current_drop(
    time_ms: np.ndarray, current_A: np.ndarray, pulse_start_ms: float
) -> np.ndarray:
    mask = time_ms < pulse_start_ms
    if np.count_nonzero(mask) < 2:
        raise ValueError("simulation needs at least two pre-pulse samples")
    baseline = float(np.median(current_A[mask]))
    return baseline - current_A


def waveform_regions(
    grid_ms: np.ndarray,
    envelope: np.ndarray,
    peak_time_ms: float,
    score: dict[str, Any],
) -> np.ndarray:
    half_width = float(score.get("peak_half_width_ms", 1.0))
    tail_threshold = float(score.get("tail_threshold", 0.2))
    regions = np.full(len(grid_ms), "decay", dtype="U8")
    regions[grid_ms < peak_time_ms - half_width] = "rise"
    regions[np.abs(grid_ms - peak_time_ms) <= half_width] = "peak"
    regions[
        (grid_ms > peak_time_ms + half_width) & (envelope < tail_threshold)
    ] = "tail"
    return regions


def regional_pair_error(
    left_residual: np.ndarray,
    right_residual: np.ndarray | None,
    regions: np.ndarray,
    weights: dict[str, float],
) -> tuple[float, dict[str, float]]:
    region_rmse: dict[str, float] = {}
    weighted = 0.0
    total_weight = 0.0
    for name, raw_weight in weights.items():
        weight = float(raw_weight)
        mask = regions == name
        if weight <= 0 or np.count_nonzero(mask) == 0:
            continue
        series = [left_residual[mask]]
        if right_residual is not None:
            series.append(right_residual[mask])
        combined = np.concatenate(series)
        rmse = float(np.sqrt(np.mean(combined**2)))
        region_rmse[name] = rmse
        weighted += weight * rmse
        total_weight += weight
    if total_weight <= 0:
        raise ValueError("at least one waveform region needs a positive weight")
    return weighted / total_weight, region_rmse


def pair_scale(left: np.ndarray, right: np.ndarray) -> tuple[float, float, float]:
    left_peak = float(np.max(left))
    right_peak = float(np.max(right))
    scale = left_peak + right_peak
    if not math.isfinite(scale) or scale <= 0:
        raise ValueError("paired waveform has no positive peak sum")
    return scale, left_peak, right_peak


def signal_peak(signal: np.ndarray, *, label: str) -> float:
    peak = float(np.max(signal))
    if not math.isfinite(peak) or peak <= 0:
        raise ValueError(f"{label} has no positive peak")
    return peak


def comparison_mode(score: dict[str, Any]) -> str:
    mode = str(score.get("comparison_mode", "paired"))
    if mode not in {"paired", "ch0_extrema"}:
        raise ValueError(
            "score.comparison_mode must be 'paired' or 'ch0_extrema', "
            f"got {mode!r}"
        )
    return mode


def selected_channel(score: dict[str, Any]) -> str:
    channel = str(score.get("comparison_channel", "left"))
    if channel not in {"left", "right"}:
        raise ValueError(
            "score.comparison_channel must be 'left' or 'right', "
            f"got {channel!r}"
        )
    return channel


def prepare_target_score_data(
    config: dict[str, Any],
    target: dict[str, Any],
    pulse_metadata: dict[str, str],
) -> dict[str, Any]:
    reference = reference_path(config, target)
    if not reference.exists():
        raise FileNotFoundError(
            f"reference is missing: {reference}; run prepare-references first"
        )
    ref_time, ref_left, ref_right = load_reference_pair(reference)

    results_dir = ROOT / "results" / pulse_metadata["case"]
    left_series = results_dir / side_series_file(pulse_metadata["series_file"], "L")
    right_series = results_dir / side_series_file(pulse_metadata["series_file"], "R")
    sim_time_left, sim_left_current = load_simulation_current(left_series)
    sim_time_right, sim_right_current = load_simulation_current(right_series)
    pulse_start_ms = float(config["pulse_start_ms"])
    sim_left = baseline_correct_current_drop(
        sim_time_left, sim_left_current, pulse_start_ms
    )
    sim_right = baseline_correct_current_drop(
        sim_time_right, sim_right_current, pulse_start_ms
    )
    sim_time_left = sim_time_left - pulse_start_ms
    sim_time_right = sim_time_right - pulse_start_ms
    if target.get("swap_simulation_channels", False):
        sim_time_left, sim_time_right = sim_time_right, sim_time_left
        sim_left, sim_right = sim_right, sim_left

    score = config["score"]
    rate_hz = float(score["readout_rate_hz"])
    cutoff_hz = float(score["readout_cutoff_hz"])
    filter_order = int(score.get("filter_order", 2))
    window_start, window_end = (
        float(value) for value in score["comparison_window_ms"]
    )
    padding_ms = float(score.get("filter_padding_ms", 10.0))
    max_shift_ms = float(score["max_shift_ms"])
    step_ms = 1e3 / rate_hz
    filter_grid = np.arange(
        window_start - padding_ms - max_shift_ms,
        window_end + padding_ms + max_shift_ms + 0.5 * step_ms,
        step_ms,
    )
    if filter_grid[0] < max(sim_time_left.min(), sim_time_right.min()):
        raise ValueError(
            f"target {target['name']}: simulation lacks pre-pulse padding for readout filtering"
        )
    if filter_grid[-1] > min(sim_time_left.max(), sim_time_right.max()):
        raise ValueError(
            f"target {target['name']}: simulation ends before the requested comparison/filter window"
        )
    sim_left_uniform = np.interp(filter_grid, sim_time_left, sim_left)
    sim_right_uniform = np.interp(filter_grid, sim_time_right, sim_right)
    sim_left_filtered = bessel_lowpass(
        sim_left_uniform, rate_hz, cutoff_hz, filter_order
    )
    sim_right_filtered = bessel_lowpass(
        sim_right_uniform, rate_hz, cutoff_hz, filter_order
    )

    comparison_grid = np.arange(
        window_start, window_end + 0.5 * step_ms, step_ms
    )
    ref_left_grid = np.interp(comparison_grid, ref_time, ref_left)
    ref_right_grid = np.interp(comparison_grid, ref_time, ref_right)
    mode = comparison_mode(score)
    ref_scale, ref_left_peak, ref_right_peak = pair_scale(
        ref_left_grid, ref_right_grid
    )
    sim_scale, sim_left_peak, sim_right_peak = pair_scale(
        sim_left_filtered, sim_right_filtered
    )
    ref_left_norm = ref_left_grid / ref_scale
    ref_right_norm = ref_right_grid / ref_scale
    sim_left_norm = sim_left_filtered / sim_scale
    sim_right_norm = sim_right_filtered / sim_scale
    envelope = ref_left_norm + ref_right_norm
    comparison_channel = selected_channel(score)
    if comparison_channel == "left":
        reference_selected_raw = ref_left_grid
        simulation_selected_raw = sim_left_filtered
    else:
        reference_selected_raw = ref_right_grid
        simulation_selected_raw = sim_right_filtered
    reference_selected_peak = signal_peak(
        reference_selected_raw,
        label=f"target {target['name']} reference {comparison_channel}",
    )
    simulation_selected_peak = signal_peak(
        simulation_selected_raw,
        label=f"target {target['name']} simulation {comparison_channel}",
    )
    reference_selected = reference_selected_raw / reference_selected_peak
    simulation_selected = simulation_selected_raw / simulation_selected_peak
    if mode == "ch0_extrema":
        envelope = reference_selected
    peak_time = float(comparison_grid[int(np.argmax(envelope))])
    regions = waveform_regions(comparison_grid, envelope, peak_time, score)
    return {
        "name": str(target["name"]),
        "weight": float(target.get("weight", 1.0)),
        "comparison_mode": mode,
        "comparison_channel": comparison_channel,
        "comparison_grid": comparison_grid,
        "filter_grid": filter_grid,
        "reference_left": ref_left_norm,
        "reference_right": ref_right_norm,
        "simulation_left": sim_left_norm,
        "simulation_right": sim_right_norm,
        "regions": regions,
        "reference_peaks": {"left": ref_left_peak, "right": ref_right_peak},
        "simulation_peaks": {"left": sim_left_peak, "right": sim_right_peak},
        "reference_left_fraction": ref_left_peak / ref_scale,
        "simulation_left_fraction": sim_left_peak / sim_scale,
        "reference_selected": reference_selected,
        "simulation_selected": simulation_selected,
        "reference_selected_peak": reference_selected_peak,
        "simulation_selected_peak": simulation_selected_peak,
        "reference_file": str(reference),
        "left_series": str(left_series),
        "right_series": str(right_series),
    }


def score_target_shift(
    data: dict[str, Any], shift_ms: float, score: dict[str, Any]
) -> dict[str, Any]:
    grid = data["comparison_grid"]
    query = grid - shift_ms
    sim_left = np.interp(query, data["filter_grid"], data["simulation_left"])
    sim_right = np.interp(query, data["filter_grid"], data["simulation_right"])
    weights = {name: float(weight) for name, weight in score["region_weights"].items()}
    if data.get("comparison_mode", "paired") == "ch0_extrema":
        simulation_selected = np.interp(
            query, data["filter_grid"], data["simulation_selected"]
        )
        selected_residual = simulation_selected - data["reference_selected"]
        waveform, region_rmse = regional_pair_error(
            selected_residual,
            None,
            data["regions"],
            weights,
        )
        return {
            "waveform_objective": waveform,
            "region_rmse": region_rmse,
            "selected_rmse": float(np.sqrt(np.mean(selected_residual**2))),
            "simulation_selected_aligned": simulation_selected,
        }

    left_residual = sim_left - data["reference_left"]
    right_residual = sim_right - data["reference_right"]
    waveform, region_rmse = regional_pair_error(
        left_residual,
        right_residual,
        data["regions"],
        weights,
    )
    peak_share_error = (
        data["simulation_left_fraction"] - data["reference_left_fraction"]
    ) ** 2
    return {
        "waveform_objective": waveform,
        "peak_share_error": peak_share_error,
        "region_rmse": region_rmse,
        "left_rmse": float(np.sqrt(np.mean(left_residual**2))),
        "right_rmse": float(np.sqrt(np.mean(right_residual**2))),
        "simulation_left_aligned": sim_left,
        "simulation_right_aligned": sim_right,
    }


def write_aligned_pair(path: Path, data: dict[str, Any], metrics: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "time_after_pulse_ms",
                "reference_left_normalized",
                "simulation_left_normalized",
                "reference_right_normalized",
                "simulation_right_normalized",
                "left_residual",
                "right_residual",
                "region",
            ]
        )
        for index, time_ms in enumerate(data["comparison_grid"]):
            ref_left = data["reference_left"][index]
            ref_right = data["reference_right"][index]
            sim_left = metrics["simulation_left_aligned"][index]
            sim_right = metrics["simulation_right_aligned"][index]
            writer.writerow(
                [
                    f"{time_ms:.12g}",
                    f"{ref_left:.12g}",
                    f"{sim_left:.12g}",
                    f"{ref_right:.12g}",
                    f"{sim_right:.12g}",
                    f"{sim_left-ref_left:.12g}",
                    f"{sim_right-ref_right:.12g}",
                    data["regions"][index],
                ]
            )


def write_aligned_selected(
    path: Path, data: dict[str, Any], metrics: dict[str, Any]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "time_after_pulse_ms",
                "comparison_channel",
                "reference_normalized",
                "simulation_normalized",
                "residual",
                "region",
            ]
        )
        for index, time_ms in enumerate(data["comparison_grid"]):
            reference = data["reference_selected"][index]
            simulation = metrics["simulation_selected_aligned"][index]
            writer.writerow(
                [
                    f"{time_ms:.12g}",
                    data["comparison_channel"],
                    f"{reference:.12g}",
                    f"{simulation:.12g}",
                    f"{simulation-reference:.12g}",
                    data["regions"][index],
                ]
            )


def relative_peak_ratio_metrics(
    target_data: list[dict[str, Any]], score: dict[str, Any]
) -> dict[str, float]:
    names = score.get("peak_ratio_targets", ["high", "low"])
    if not isinstance(names, list) or len(names) != 2:
        raise ValueError("score.peak_ratio_targets must name exactly two targets")
    by_name = {str(data["name"]): data for data in target_data}
    numerator, denominator = (str(name) for name in names)
    if numerator not in by_name or denominator not in by_name:
        raise ValueError(
            "score.peak_ratio_targets must refer to configured targets, "
            f"got {names!r}"
        )
    reference_ratio = (
        by_name[numerator]["reference_selected_peak"]
        / by_name[denominator]["reference_selected_peak"]
    )
    simulation_ratio = (
        by_name[numerator]["simulation_selected_peak"]
        / by_name[denominator]["simulation_selected_peak"]
    )
    log_error = math.log(simulation_ratio / reference_ratio)
    return {
        "reference_peak_ratio": float(reference_ratio),
        "simulation_peak_ratio": float(simulation_ratio),
        "peak_ratio_error": float(log_error**2),
    }


def score_filename(config: dict[str, Any]) -> str:
    mode = comparison_mode(config["score"])
    return "score.json" if mode == "paired" else f"score_{mode}.json"


def leaderboard_filename(config: dict[str, Any]) -> str:
    mode = comparison_mode(config["score"])
    return "leaderboard.csv" if mode == "paired" else f"leaderboard_{mode}.csv"


def score_candidate_outputs(
    config: dict[str, Any],
    candidate: Candidate,
    metadata: dict[str, Any],
) -> dict[str, Any]:
    target_data = [
        prepare_target_score_data(
            config, target, metadata["pulse_cases"][str(target["name"])]
        )
        for target in config["targets"]
    ]
    score = config["score"]
    max_shift = float(score["max_shift_ms"])
    shift_step = float(score["shift_step_ms"])
    shift_sigma = float(score["shift_sigma_ms"])
    shift_penalty_weight = float(score["shift_penalty_weight"])
    mode = comparison_mode(score)
    peak_share_weight = float(score.get("peak_share_weight", 0.0))
    peak_ratio = (
        relative_peak_ratio_metrics(target_data, score)
        if mode == "ch0_extrema"
        else None
    )
    peak_ratio_weight = float(score.get("peak_ratio_weight", 0.0))
    shifts = np.arange(-max_shift, max_shift + 0.5 * shift_step, shift_step)
    best: dict[str, Any] | None = None
    no_shift: dict[str, Any] | None = None
    for shift in shifts:
        target_metrics = [score_target_shift(data, float(shift), score) for data in target_data]
        weight_sum = sum(data["weight"] for data in target_data)
        waveform = sum(
            data["weight"] * metrics["waveform_objective"]
            for data, metrics in zip(target_data, target_metrics)
        ) / weight_sum
        peak_share = (
            sum(
                data["weight"] * metrics["peak_share_error"]
                for data, metrics in zip(target_data, target_metrics)
            ) / weight_sum
            if mode == "paired"
            else 0.0
        )
        shift_penalty = shift_penalty_weight * (float(shift) / shift_sigma) ** 2
        peak_ratio_objective = (
            peak_ratio_weight * peak_ratio["peak_ratio_error"]
            if peak_ratio is not None
            else 0.0
        )
        aligned = (
            waveform
            + peak_share_weight * peak_share
            + peak_ratio_objective
            + shift_penalty
        )
        record = {
            "shift_ms": float(shift),
            "waveform_objective": waveform,
            "peak_share_objective": peak_share,
            "peak_ratio_objective": peak_ratio_objective,
            "shift_penalty": shift_penalty,
            "aligned_objective": aligned,
            "target_metrics": target_metrics,
        }
        if abs(shift) < 0.5 * shift_step:
            no_shift = record
        if best is None or aligned < best["aligned_objective"]:
            best = record
    if best is None or no_shift is None:
        raise RuntimeError("failed to evaluate the requested common time shifts")

    prior = 0.0
    for variable_name, factor in candidate.factors.items():
        sigma = float(config["variables"][variable_name]["prior_sigma_log"])
        prior += (math.log(factor) / sigma) ** 2
    prior_weight = float(score["prior_weight"])
    objective = float(best["aligned_objective"]) + prior_weight * prior

    targets: dict[str, Any] = {}
    for data, metrics in zip(target_data, best["target_metrics"]):
        if mode == "ch0_extrema":
            trace_path = (
                candidate_dir(config, candidate.candidate_id)
                / f"aligned_{mode}_{data['name']}.csv"
            )
            write_aligned_selected(trace_path, data, metrics)
            targets[data["name"]] = {
                "weight": data["weight"],
                "comparison_channel": data["comparison_channel"],
                "waveform_objective": metrics["waveform_objective"],
                "region_rmse": metrics["region_rmse"],
                "selected_rmse": metrics["selected_rmse"],
                "reference_peak": data["reference_selected_peak"],
                "simulation_peak": data["simulation_selected_peak"],
                "reference_file": data["reference_file"],
                "left_series": data["left_series"],
                "right_series": data["right_series"],
                "aligned_trace": str(trace_path),
            }
        else:
            trace_path = candidate_dir(config, candidate.candidate_id) / f"aligned_{data['name']}.csv"
            write_aligned_pair(trace_path, data, metrics)
            targets[data["name"]] = {
                "weight": data["weight"],
                "waveform_objective": metrics["waveform_objective"],
                "peak_share_error": metrics["peak_share_error"],
                "region_rmse": metrics["region_rmse"],
                "left_rmse": metrics["left_rmse"],
                "right_rmse": metrics["right_rmse"],
                "reference_left_fraction": data["reference_left_fraction"],
                "simulation_left_fraction": data["simulation_left_fraction"],
                "reference_peaks": data["reference_peaks"],
                "simulation_peaks": data["simulation_peaks"],
                "reference_file": data["reference_file"],
                "left_series": data["left_series"],
                "right_series": data["right_series"],
                "aligned_trace": str(trace_path),
            }
    return {
        "candidate_id": candidate.candidate_id,
        "factors": candidate.factors,
        "comparison_mode": mode,
        "comparison_channel": target_data[0]["comparison_channel"],
        "best_shift_ms": best["shift_ms"],
        "waveform_objective": best["waveform_objective"],
        "peak_share_objective": best["peak_share_objective"],
        "peak_ratio_objective": best["peak_ratio_objective"],
        "shift_penalty": best["shift_penalty"],
        "aligned_objective": best["aligned_objective"],
        "objective_no_shift": no_shift["aligned_objective"],
        "prior_penalty": prior,
        "objective": objective,
        **(peak_ratio or {}),
        "targets": targets,
    }


def prepare(
    config: dict[str, Any],
    *,
    sample_count: int | None = None,
    seed: int | None = None,
) -> list[Candidate]:
    candidates = candidates_from_config(
        config, sample_count=sample_count, seed=seed
    )
    base_project = load_json(resolve_from_root(config["base_project"]))
    root = output_root(config)
    root.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for candidate in candidates:
        project_path, metadata = write_trial(
            base_project,
            config,
            candidate,
            float(config["g0_calibration"]["initial_factor"]),
            0,
        )
        rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "factors": json.dumps(candidate.factors, sort_keys=True),
                "project": str(project_path.relative_to(ROOT)),
                "steady_case": metadata["steady_case"],
                "pulse_cases": json.dumps(metadata["pulse_cases"], sort_keys=True),
            }
        )
    table_path = root / "candidates.csv"
    with table_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    write_json(
        root / "search_manifest.json",
        {
            "config": str(DEFAULT_CONFIG.relative_to(ROOT)),
            "base_project": config["base_project"],
            "candidate_count": len(candidates),
            "seed": int(
                config["multivariate_search"]["seed"] if seed is None else seed
            ),
            "variables": list(config["multivariate_search"]["variables"]),
            "targets": [target["name"] for target in config["targets"]],
            "candidates": [candidate.candidate_id for candidate in candidates],
        },
    )
    return candidates


def update_leaderboard(config: dict[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for path in output_root(config).glob(f"candidates/*/{score_filename(config)}"):
        score = load_json(path)
        rows.append(
            {
                "candidate_id": score["candidate_id"],
                "comparison_mode": score.get("comparison_mode", "paired"),
                "objective": score["objective"],
                "waveform_objective": score["waveform_objective"],
                "peak_share_objective": score["peak_share_objective"],
                "peak_ratio_objective": score.get("peak_ratio_objective", 0.0),
                "best_shift_ms": score["best_shift_ms"],
                "g0_factor": score.get("g0_factor", ""),
                "factors": json.dumps(score["factors"], sort_keys=True),
            }
        )
    if not rows:
        return
    rows.sort(key=lambda row: float(row["objective"]))
    with (output_root(config) / leaderboard_filename(config)).open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def run_candidate(config: dict[str, Any], candidate: Candidate) -> dict[str, Any]:
    base_project = load_json(resolve_from_root(config["base_project"]))
    if bool(config["g0_calibration"].get("enabled", True)):
        project_path, metadata = calibrate_g0(base_project, config, candidate)
    else:
        project_path, metadata = write_trial(
            base_project,
            config,
            candidate,
            float(config["g0_calibration"]["initial_factor"]),
            0,
        )
        run_command(project_path, metadata["steady_case"])

    run_pulse_cases: set[str] = set()
    for target in config["targets"]:
        pulse = metadata["pulse_cases"][str(target["name"])]
        # run.py always (re-)executes its final requested case even if the
        # result already exists, so targets sharing a deduped case (see
        # mutate_project) must only be run once here, not once per target.
        if pulse["case"] in run_pulse_cases:
            continue
        run_pulse_cases.add(pulse["case"])
        run_command(project_path, pulse["case"])
    score = score_candidate_outputs(config, candidate, metadata)
    score.update(
        {
            "project": str(project_path.relative_to(ROOT)),
            "steady_case": metadata["steady_case"],
            "pulse_cases": metadata["pulse_cases"],
            "g0_factor": metadata["g0_factor"],
        }
    )
    write_json(candidate_dir(config, candidate.candidate_id) / score_filename(config), score)
    update_leaderboard(config)
    return score


def rescore_candidate(config: dict[str, Any], candidate: Candidate) -> dict[str, Any]:
    """Score an already-run candidate without executing Elmer again."""
    legacy_path = candidate_dir(config, candidate.candidate_id) / "score.json"
    if not legacy_path.exists():
        raise FileNotFoundError(
            f"candidate {candidate.candidate_id} has no legacy score metadata at "
            f"{legacy_path}; run the candidate before rescoring it"
        )
    existing = load_json(legacy_path)
    required = ("project", "steady_case", "pulse_cases", "g0_factor")
    missing = [key for key in required if key not in existing]
    if missing:
        raise ValueError(
            f"candidate {candidate.candidate_id} cannot be rescored; "
            f"legacy score lacks {missing}"
        )
    metadata = {key: existing[key] for key in required}
    score = score_candidate_outputs(config, candidate, metadata)
    score.update(metadata)
    write_json(candidate_dir(config, candidate.candidate_id) / score_filename(config), score)
    update_leaderboard(config)
    return score


def print_candidate_table(candidates: list[Candidate]) -> None:
    print(f"prepared {len(candidates)} PoST multivariate candidates")
    for candidate in candidates:
        print(f"  {candidate.candidate_id:20s} factors={candidate.factors}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("prepare-references")
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--samples", type=int)
    prepare_parser.add_argument("--seed", type=int)
    dry_parser = subparsers.add_parser("dry-run")
    dry_parser.add_argument("candidate")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("candidate")
    rescore_parser = subparsers.add_parser("rescore")
    rescore_parser.add_argument("candidate")
    run_all_parser = subparsers.add_parser("run-all")
    run_all_parser.add_argument("--limit", type=int)
    run_all_parser.add_argument("--samples", type=int)
    run_all_parser.add_argument("--seed", type=int)
    run_all_parser.add_argument("--skip-scored", action="store_true")
    rescore_all_parser = subparsers.add_parser("rescore-all")
    rescore_all_parser.add_argument("--skip-scored", action="store_true")
    args = parser.parse_args()

    config_path = args.config if args.config.is_absolute() else ROOT / args.config
    config = load_json(config_path)
    if args.command == "prepare-references":
        paths = prepare_references(config)
        print(json.dumps({key: str(path) for key, path in paths.items()}, indent=2))
        return 0

    sample_count = getattr(args, "samples", None)
    seed = getattr(args, "seed", None)
    candidates = candidates_from_config(
        config, sample_count=sample_count, seed=seed
    )
    by_id = {candidate.candidate_id: candidate for candidate in candidates}
    if args.command == "prepare":
        print_candidate_table(prepare(config, sample_count=sample_count, seed=seed))
        return 0
    if args.command == "dry-run":
        if args.candidate not in by_id:
            raise SystemExit(f"unknown candidate: {args.candidate}")
        base_project = load_json(resolve_from_root(config["base_project"]))
        project_path, metadata = write_trial(
            base_project,
            config,
            by_id[args.candidate],
            float(config["g0_calibration"]["initial_factor"]),
            0,
        )
        first_target = str(config["targets"][0]["name"])
        run_command(
            project_path,
            metadata["pulse_cases"][first_target]["case"],
            dry_run=True,
        )
        return 0
    if args.command == "run":
        if args.candidate not in by_id:
            raise SystemExit(f"unknown candidate: {args.candidate}")
        print(json.dumps(run_candidate(config, by_id[args.candidate]), indent=2))
        return 0
    if args.command == "rescore":
        if args.candidate not in by_id:
            raise SystemExit(f"unknown candidate: {args.candidate}")
        print(json.dumps(rescore_candidate(config, by_id[args.candidate]), indent=2))
        return 0
    if args.command == "run-all":
        selected = candidates
        if args.skip_scored:
            selected = [
                candidate
                for candidate in selected
                if not (candidate_dir(config, candidate.candidate_id) / score_filename(config)).exists()
            ]
        if args.limit is not None:
            selected = selected[: args.limit]
        for index, candidate in enumerate(selected, start=1):
            print(f"[{index}/{len(selected)}] {candidate.candidate_id}", flush=True)
            run_candidate(config, candidate)
        return 0
    if args.command == "rescore-all":
        selected = [
            candidate
            for candidate in candidates
            if (candidate_dir(config, candidate.candidate_id) / "score.json").exists()
        ]
        if args.skip_scored:
            selected = [
                candidate
                for candidate in selected
                if not (candidate_dir(config, candidate.candidate_id) / score_filename(config)).exists()
            ]
        for index, candidate in enumerate(selected, start=1):
            print(f"[{index}/{len(selected)}] {candidate.candidate_id}", flush=True)
            rescore_candidate(config, candidate)
        return 0
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
