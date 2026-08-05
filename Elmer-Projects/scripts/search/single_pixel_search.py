"""Prepare, run, score, and analyse single-pixel parameter searches.

The first stage is an one-factor-at-a-time sensitivity scan on the coarse
single-pixel mesh.  Every candidate gets its own steady and pulse case names,
so results are cacheable and cannot collide with the main project cases.

Typical use from ``Elmer-Projects``::

    python scripts/search/single_pixel_search.py prepare
    python scripts/search/single_pixel_search.py score-existing
    python scripts/search/single_pixel_search.py dry-run baseline
    python scripts/search/single_pixel_search.py run baseline
    python scripts/search/single_pixel_search.py run-all --limit 3
    python scripts/search/single_pixel_search.py analyze-sensitivity

Positive material parameters are changed as multiplicative factors in their
expression fields.  For steady-sensitive candidates, G0 is calibrated by a
bracketed bisection on the converged circuit state before the pulse is run.
The pulse objective uses baseline correction, independent peak normalization,
phase-balanced waveform errors, one bounded time shift, and weak log-space
priors.  Experimental multivariate fitting is intentionally handled by
``post_multivariate_search.py`` on the PoST dual-TES geometry.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "single_pixel_search_config.json"


@dataclass(frozen=True)
class Candidate:
    candidate_id: str
    factors: dict[str, float]
    steady_sensitive: bool


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def resolve_from_root(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


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


def factor_tag(factor: float) -> str:
    percent = round((factor - 1.0) * 100)
    return f"p{percent}" if percent >= 0 else f"m{abs(percent)}"


def candidates_from_config(config: dict[str, Any]) -> list[Candidate]:
    candidates = [Candidate("baseline", {}, False)]
    for variable_name, spec in config["variables"].items():
        for factor in spec["factors"]:
            candidates.append(
                Candidate(
                    f"{variable_name}_{factor_tag(float(factor))}",
                    {variable_name: float(factor)},
                    bool(spec.get("steady_sensitive", False)),
                )
            )
    return candidates


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
            log_factor = math.log(lower) + sample[index] * (
                math.log(upper) - math.log(lower)
            )
            factors[variable_name] = math.exp(log_factor)
        rows.append(factors)
    return rows


def multivariate_candidates_from_config(
    config: dict[str, Any],
    *,
    sample_count: int | None = None,
    seed: int | None = None,
) -> list[Candidate]:
    candidates: list[Candidate] = []
    for index, factors in enumerate(
        latin_hypercube_log_factors(config, sample_count=sample_count, seed=seed)
    ):
        steady_sensitive = any(
            bool(config["variables"][name].get("steady_sensitive", False))
            for name in factors
        )
        candidate_id = f"mv_{index:03d}_{short_hash(factors)}"
        candidates.append(Candidate(candidate_id, factors, steady_sensitive))
    return candidates


def candidate_map(config: dict[str, Any]) -> dict[str, Candidate]:
    return {candidate.candidate_id: candidate for candidate in candidates_from_config(config)}


def short_hash(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()[:8]


def case_names(candidate: Candidate, trial_index: int) -> tuple[str, str]:
    stem = f"spsearch_{candidate.candidate_id}_{trial_index:02d}"
    return f"{stem}_steady", f"{stem}_pulse"


def mutate_project(
    base_project: dict[str, Any],
    config: dict[str, Any],
    candidate: Candidate,
    g0_factor: float,
    trial_index: int,
) -> tuple[dict[str, Any], dict[str, Any]]:
    project = copy.deepcopy(base_project)
    for variable_name, factor in candidate.factors.items():
        path = config["variables"][variable_name]["path"]
        original = str(nested_get(project, path))
        nested_set(project, path, scaled_expression(original, factor))

    original_g0 = str(project["parameter_expressions"]["G0"])
    project["parameter_expressions"]["G0"] = scaled_expression(
        original_g0, g0_factor
    )

    steady_name, pulse_name = case_names(candidate, trial_index)
    mesh_name = str(config["mesh"])
    state_file = f"work/meshes/{project['meshes'][mesh_name]['dir']}/{steady_name}.state"
    steady_series = f"{steady_name}_series.csv"
    pulse_series = f"{pulse_name}_series.csv"
    timesteps = config["screen_timesteps"]
    output_intervals = [1] * len(timesteps)

    steady_case = {
        "template": "steady",
        "mesh": mesh_name,
        "heat_source": "circuit_implicit",
        "series_file": steady_series,
        "state_file": state_file,
        "initial_temperature": "T_0",
        "output_result": True,
        "output_file_path": (
            f"../work/meshes/{project['meshes'][mesh_name]['dir']}/"
            f"{steady_name}.result"
        ),
        "vtu": False,
        "steady_state_max_iterations": 1,
        "output_intervals": 1,
        "solver": {
            "nonlinear_max_iterations": 120,
            "nonlinear_convergence_tolerance": 1e-8,
            "nonlinear_relaxation_factor": 1.0,
            "steady_state_convergence_tolerance": 1e-8,
        },
    }
    pulse_case = {
        "template": "pulse",
        "mesh": mesh_name,
        "restart_from": steady_name,
        "restart_file_path": (
            f"../work/meshes/{project['meshes'][mesh_name]['dir']}/"
            f"{steady_name}.result"
        ),
        "restart_time": 0.0,
        "series_file": pulse_series,
        "state_file": state_file,
        "initial_temperature": "T_0",
        "lumped_mass": True,
        "vtu": False,
        "timesteps": timesteps,
        "output_intervals": output_intervals,
        "pulse": {
            "energy": "1332[keV]",
            "start": config["pulse_start"],
            "duration": "1[ns]",
            "sigma": "50[um]",
            "center": "auto",
        },
        "steady_state_max_iterations": 1,
        "solver_comment": "short coarse-mesh parameter-screening pulse",
        "solver": {
            "nonlinear_max_iterations": 25,
            "nonlinear_convergence_tolerance": 3e-7,
            "nonlinear_relaxation_factor": 1.0,
            "steady_state_convergence_tolerance": 1e-9,
        },
    }
    project["cases"] = {steady_name: steady_case, pulse_name: pulse_case}
    metadata = {
        "candidate_id": candidate.candidate_id,
        "candidate_factors": candidate.factors,
        "steady_sensitive": candidate.steady_sensitive,
        "g0_factor": g0_factor,
        "trial_index": trial_index,
        "steady_case": steady_name,
        "pulse_case": pulse_name,
        "state_file": state_file,
        "pulse_series": pulse_series,
        "project_hash": short_hash(
            {"factors": candidate.factors, "g0_factor": g0_factor}
        ),
    }
    return project, metadata


def output_root(config: dict[str, Any]) -> Path:
    return resolve_from_root(config["output_dir"])


def candidate_dir(config: dict[str, Any], candidate_id: str) -> Path:
    return output_root(config) / "candidates" / candidate_id


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


def prepare_candidates(
    config: dict[str, Any],
    candidates: list[Candidate],
    *,
    manifest_name: str,
    table_name: str,
) -> list[Candidate]:
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
                "steady_sensitive": candidate.steady_sensitive,
                "project": str(project_path.relative_to(ROOT)),
                "steady_case": metadata["steady_case"],
                "pulse_case": metadata["pulse_case"],
            }
        )
    with (root / table_name).open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    write_json(
        root / manifest_name,
        {
            "config": str(DEFAULT_CONFIG.relative_to(ROOT)),
            "base_project": config["base_project"],
            "candidate_count": len(candidates),
            "candidates": [row["candidate_id"] for row in rows],
        },
    )
    return candidates


def prepare(config: dict[str, Any]) -> list[Candidate]:
    return prepare_candidates(
        config,
        candidates_from_config(config),
        manifest_name="search_manifest.json",
        table_name="candidates.csv",
    )


def prepare_multivariate(
    config: dict[str, Any],
    *,
    sample_count: int | None = None,
    seed: int | None = None,
) -> list[Candidate]:
    candidates = multivariate_candidates_from_config(
        config, sample_count=sample_count, seed=seed
    )
    prepared = prepare_candidates(
        config,
        candidates,
        manifest_name="multivariate_manifest.json",
        table_name="multivariate_candidates.csv",
    )
    manifest_path = output_root(config) / "multivariate_manifest.json"
    manifest = load_json(manifest_path)
    manifest["sample_count"] = len(candidates)
    manifest["seed"] = int(
        config["multivariate_search"]["seed"] if seed is None else seed
    )
    manifest["variables"] = list(config["multivariate_search"]["variables"])
    write_json(manifest_path, manifest)
    return prepared


def run_command(
    project_path: Path,
    case_name: str,
    *,
    dry_run: bool = False,
    force_deps: bool = False,
    mpi_procs: int = 1,
) -> None:
    command = [
        sys.executable,
        str(ROOT / "run.py"),
        case_name,
        "--project",
        str(project_path),
        "--mpi-procs",
        str(mpi_procs),
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
) -> tuple[float, Path, dict[str, Any], dict[str, float] | None]:
    project_path, metadata = write_trial(
        base_project, config, candidate, g0_factor, trial_index
    )
    run_command(project_path, metadata["steady_case"], dry_run=dry_run)
    if dry_run:
        return math.nan, project_path, metadata, None
    state_path = ROOT / metadata["state_file"]
    state = read_state(state_path)
    project = load_json(project_path)
    error_K = state["temperature_K"] - target_temperature(project)
    trial_record = {
        **metadata,
        **state,
        "temperature_error_mK": error_K * 1e3,
    }
    write_json(project_path.parent / "steady_result.json", trial_record)
    return error_K, project_path, metadata, state


def calibrate_g0(
    base_project: dict[str, Any],
    config: dict[str, Any],
    candidate: Candidate,
    initial_factor: float | None = None,
) -> tuple[Path, dict[str, Any]]:
    settings = config["g0_calibration"]
    tolerance_K = float(settings["temperature_tolerance_mK"]) * 1e-3
    initial = (
        float(settings["initial_factor"])
        if initial_factor is None
        else float(initial_factor)
    )
    lower = float(settings["lower_factor"])
    upper = float(settings["upper_factor"])
    max_iterations = int(settings["max_iterations"])

    def accept(
        accepted_error: float,
        accepted_project: Path,
        accepted_metadata: dict[str, Any],
    ) -> tuple[Path, dict[str, Any]]:
        write_json(
            candidate_dir(config, candidate.candidate_id) / "calibration.json",
            {
                "candidate_id": candidate.candidate_id,
                "g0_factor": accepted_metadata["g0_factor"],
                "temperature_error_mK": accepted_error * 1e3,
                "project": str(accepted_project.relative_to(ROOT)),
                "steady_case": accepted_metadata["steady_case"],
            },
        )
        return accepted_project, accepted_metadata

    error, project_path, metadata, _ = evaluate_steady_trial(
        base_project, config, candidate, initial, 0
    )
    if abs(error) <= tolerance_K:
        return accept(error, project_path, metadata)

    trial_index = 1
    low_error, low_project, low_metadata, _ = evaluate_steady_trial(
        base_project, config, candidate, lower, trial_index
    )
    trial_index += 1
    if low_error == 0:
        return accept(low_error, low_project, low_metadata)

    # Prefer the initial point as one bracket endpoint.  This avoids an
    # unnecessary upper-bound run in the common case where the nominal coarse
    # mesh is too cold and reducing G0 crosses T0.
    if low_error * error <= 0:
        lower_factor, lower_error = lower, low_error
        lower_project, lower_metadata = low_project, low_metadata
        upper_factor, upper_error = initial, error
        upper_project, upper_metadata = project_path, metadata
    else:
        high_error, high_project, high_metadata, _ = evaluate_steady_trial(
            base_project, config, candidate, upper, trial_index
        )
        trial_index += 1
        if high_error == 0:
            return accept(high_error, high_project, high_metadata)
        if low_error * high_error > 0:
            raise RuntimeError(
                "G0 calibration did not bracket T0: "
                f"factor={lower} error={low_error*1e3:+.6f} mK, "
                f"factor={initial} error={error*1e3:+.6f} mK, "
                f"factor={upper} error={high_error*1e3:+.6f} mK"
            )
        lower_factor, lower_error = lower, low_error
        lower_project, lower_metadata = low_project, low_metadata
        upper_factor, upper_error = upper, high_error
        upper_project, upper_metadata = high_project, high_metadata

    best = min(
        [
            (abs(lower_error), lower_error, lower_project, lower_metadata),
            (abs(upper_error), upper_error, upper_project, upper_metadata),
        ],
        key=lambda item: item[0],
    )
    while trial_index < max_iterations:
        log_lower = math.log(lower_factor)
        log_upper = math.log(upper_factor)
        denominator = upper_error - lower_error
        if denominator == 0:
            log_middle = 0.5 * (log_lower + log_upper)
        else:
            log_middle = log_lower - lower_error * (
                log_upper - log_lower
            ) / denominator
            # Keep the regula-falsi proposal away from an endpoint.  This
            # preserves a shrinking bracket when the response is curved.
            guard = 0.12 * (log_upper - log_lower)
            log_middle = min(
                max(log_middle, log_lower + guard), log_upper - guard
            )
        middle = math.exp(log_middle)
        middle_error, middle_project, middle_metadata, _ = evaluate_steady_trial(
            base_project, config, candidate, middle, trial_index
        )
        trial_index += 1
        if abs(middle_error) < best[0]:
            best = (
                abs(middle_error),
                middle_error,
                middle_project,
                middle_metadata,
            )
        if abs(middle_error) <= tolerance_K:
            return accept(middle_error, middle_project, middle_metadata)
        if lower_error * middle_error <= 0:
            upper_factor, upper_error = middle, middle_error
            upper_project, upper_metadata = middle_project, middle_metadata
        else:
            lower_factor, lower_error = middle, middle_error
            lower_project, lower_metadata = middle_project, middle_metadata
    if best[0] <= tolerance_K:
        return accept(best[1], best[2], best[3])
    raise RuntimeError(
        f"G0 calibration failed to reach {tolerance_K*1e3:.4f} mK; "
        f"best error was {best[0]*1e3:.6f} mK"
    )


def pulse_metrics(
    time_ms: np.ndarray,
    signal: np.ndarray,
    pulse_start_ms: float,
    baseline_window: tuple[float, float] | None = None,
    *,
    response_direction: str = "drop",
    baseline_statistic: str = "median",
) -> dict[str, float | np.ndarray]:
    if baseline_window is None:
        baseline_mask = time_ms < pulse_start_ms
    else:
        baseline_mask = (time_ms >= baseline_window[0]) & (time_ms <= baseline_window[1])
    if np.count_nonzero(baseline_mask) < 2:
        raise ValueError("at least two pre-pulse samples are required for baseline correction")
    if baseline_statistic == "median":
        baseline = float(np.median(signal[baseline_mask]))
    elif baseline_statistic == "mean":
        baseline = float(np.mean(signal[baseline_mask]))
    else:
        raise ValueError(
            f"unknown baseline_statistic {baseline_statistic!r}; use 'median' or 'mean'"
        )

    post = np.flatnonzero(time_ms >= pulse_start_ms)
    if len(post) == 0:
        raise ValueError("series has no post-pulse samples")

    rising_response = signal - baseline
    dropping_response = baseline - signal
    if response_direction == "rise":
        response = rising_response
    elif response_direction == "drop":
        response = dropping_response
    elif response_direction == "auto":
        response = (
            rising_response
            if np.max(rising_response[post]) >= np.max(dropping_response[post])
            else dropping_response
        )
        response_direction = "rise" if response is rising_response else "drop"
    else:
        raise ValueError(
            f"unknown response_direction {response_direction!r}; use rise, drop, or auto"
        )

    peak_index = int(post[np.argmax(response[post])])
    peak = float(response[peak_index])
    if peak <= 0:
        raise ValueError("pulse response has no positive peak after polarity correction")
    normalized = response / peak

    def crossing(level: float) -> float:
        y = normalized[: peak_index + 1]
        indices = np.flatnonzero((y[:-1] < level) & (y[1:] >= level))
        if len(indices) == 0:
            return math.nan
        index = int(indices[-1])
        fraction = (level - y[index]) / (y[index + 1] - y[index])
        return float(time_ms[index] + fraction * (time_ms[index + 1] - time_ms[index]))

    t10 = crossing(0.1)
    t90 = crossing(0.9)
    return {
        "baseline_signal": baseline,
        "response_signal": response,
        "normalized": normalized,
        "peak_signal": peak,
        "response_direction": response_direction,
        "peak_time_ms": float(time_ms[peak_index]),
        "peak_delay_ms": float(time_ms[peak_index] - pulse_start_ms),
        "rise_time_10_90_ms": float(t90 - t10),
    }


def load_reference_series(config: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, str]:
    spec = config.get("reference_data", {})
    delimiter = spec.get("delimiter")
    reference = np.loadtxt(
        resolve_from_root(config["reference"]),
        comments=str(spec.get("comments", "%")),
        delimiter=delimiter,
        encoding="utf-8",
    )
    time_column = int(spec.get("time_column", 0))
    signal_column = int(spec.get("signal_column", 4))
    time_ms = reference[:, time_column] * float(spec.get("time_scale_to_ms", 1.0))
    signal = reference[:, signal_column] * float(spec.get("signal_scale", 1.0))
    return time_ms, signal, str(spec.get("response_direction", "drop"))


def load_simulation_series(config: dict[str, Any], path: Path) -> tuple[np.ndarray, np.ndarray, str]:
    spec = config.get("simulation_data", {})
    simulation = np.genfromtxt(path, delimiter=",", names=True)
    time_column = str(spec.get("time_column", "time_s"))
    signal_column = str(spec.get("signal_column", "tes_current_A"))
    time_ms = simulation[time_column] * float(spec.get("time_scale_to_ms", 1e3))
    signal = simulation[signal_column] * float(spec.get("signal_scale", 1e6))
    return time_ms, signal, str(spec.get("response_direction", "drop"))


def waveform_regions(
    grid_ms: np.ndarray,
    reference_normalized: np.ndarray,
    reference_peak_delay_ms: float,
    score_settings: dict[str, Any],
) -> np.ndarray:
    half_width = float(score_settings.get("peak_half_width_ms", 0.03))
    tail_threshold = float(score_settings.get("tail_threshold", 0.2))
    regions = np.full(len(grid_ms), "decay", dtype="U8")
    regions[grid_ms < reference_peak_delay_ms - half_width] = "rise"
    peak_mask = np.abs(grid_ms - reference_peak_delay_ms) <= half_width
    regions[peak_mask] = "peak"
    tail_mask = (
        (grid_ms > reference_peak_delay_ms + half_width)
        & (reference_normalized < tail_threshold)
    )
    regions[tail_mask] = "tail"
    return regions


def regional_waveform_error(
    residual: np.ndarray,
    regions: np.ndarray,
    weights: dict[str, float],
) -> tuple[float, dict[str, float]]:
    region_rmse: dict[str, float] = {}
    weighted_sum = 0.0
    weight_sum = 0.0
    for region_name, weight_value in weights.items():
        weight = float(weight_value)
        mask = regions == region_name
        if weight <= 0 or np.count_nonzero(mask) == 0:
            continue
        rmse = float(np.sqrt(np.mean(residual[mask] ** 2)))
        region_rmse[region_name] = rmse
        weighted_sum += weight * rmse
        weight_sum += weight
    if weight_sum <= 0:
        raise ValueError("at least one waveform region must have a positive weight")
    return weighted_sum / weight_sum, region_rmse


def write_aligned_trace(
    path: Path,
    time_ms: np.ndarray,
    reference: np.ndarray,
    simulation: np.ndarray,
    regions: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "time_after_pulse_ms",
                "reference_normalized",
                "simulation_normalized",
                "residual",
                "region",
            ]
        )
        for time_value, ref_value, sim_value, region in zip(
            time_ms, reference, simulation, regions
        ):
            writer.writerow(
                [
                    f"{time_value:.12g}",
                    f"{ref_value:.12g}",
                    f"{sim_value:.12g}",
                    f"{sim_value-ref_value:.12g}",
                    region,
                ]
            )


def score_series(
    config: dict[str, Any],
    series_path: Path,
    candidate: Candidate,
    simulation_pulse_start_ms: float,
    *,
    trace_path: Path | None = None,
) -> dict[str, Any]:
    ref_time, ref_signal, ref_direction = load_reference_series(config)
    sim_time, sim_signal, sim_direction = load_simulation_series(config, series_path)
    reference_pulse_start = float(config["reference_pulse_start_ms"])
    baseline_window = tuple(float(value) for value in config["baseline_window_ms"])
    score_settings = config["score"]
    baseline_statistic = str(score_settings.get("baseline_statistic", "median"))
    ref_metrics = pulse_metrics(
        ref_time,
        ref_signal,
        reference_pulse_start,
        baseline_window,
        response_direction=ref_direction,
        baseline_statistic=baseline_statistic,
    )
    sim_metrics = pulse_metrics(
        sim_time,
        sim_signal,
        simulation_pulse_start_ms,
        response_direction=sim_direction,
        baseline_statistic=baseline_statistic,
    )

    ref_relative = ref_time - reference_pulse_start
    sim_relative = sim_time - simulation_pulse_start_ms
    window_start, window_end = (
        float(value) for value in config["comparison_window_ms"]
    )
    window_end = min(window_end, float(sim_relative.max()), float(ref_relative.max()))
    grid = np.linspace(
        window_start, window_end, int(score_settings["grid_points"])
    )
    ref_norm = np.interp(grid, ref_relative, ref_metrics["normalized"])
    regions = waveform_regions(
        grid,
        ref_norm,
        float(ref_metrics["peak_delay_ms"]),
        score_settings,
    )
    region_weights = {
        name: float(weight)
        for name, weight in score_settings.get(
            "region_weights", {"rise": 1.0, "peak": 1.0, "decay": 1.0, "tail": 1.0}
        ).items()
    }
    max_shift = float(score_settings["max_shift_ms"])
    shift_step = float(score_settings["shift_step_ms"])
    shift_sigma = float(score_settings.get("shift_sigma_ms", max_shift))
    shift_weight = float(score_settings.get("shift_penalty_weight", 0.0))
    if shift_sigma <= 0:
        raise ValueError("shift_sigma_ms must be positive")
    shifts = np.arange(-max_shift, max_shift + 0.5 * shift_step, shift_step)
    best: dict[str, Any] | None = None
    no_shift_rmse = math.nan
    no_shift_waveform_objective = math.nan
    no_shift_region_rmse: dict[str, float] = {}
    for shift in shifts:
        # Positive shift means the simulated waveform is delayed.
        shifted_query = grid - shift
        valid = (shifted_query >= sim_relative.min()) & (
            shifted_query <= sim_relative.max()
        )
        if np.count_nonzero(valid) < 0.8 * len(grid):
            continue
        sim_norm = np.interp(
            shifted_query[valid], sim_relative, sim_metrics["normalized"]
        )
        residual = sim_norm - ref_norm[valid]
        rmse = float(np.sqrt(np.mean(residual**2)))
        mae = float(np.mean(np.abs(residual)))
        waveform_objective, region_rmse = regional_waveform_error(
            residual, regions[valid], region_weights
        )
        shift_penalty = shift_weight * (float(shift) / shift_sigma) ** 2
        aligned_objective = waveform_objective + shift_penalty
        if abs(shift) < 0.5 * shift_step:
            no_shift_rmse = rmse
            no_shift_waveform_objective = waveform_objective
            no_shift_region_rmse = region_rmse
        if best is None or aligned_objective < best["aligned_objective"]:
            best = {
                "aligned_objective": aligned_objective,
                "waveform_objective": waveform_objective,
                "shift_penalty": shift_penalty,
                "shift_ms": float(shift),
                "rmse": rmse,
                "mae": mae,
                "region_rmse": region_rmse,
                "valid": valid,
                "simulation_normalized": sim_norm,
            }
    if best is None:
        raise ValueError("no valid time shift for the requested comparison window")

    prior_scale = float(score_settings["prior_scale"])
    prior = 0.0
    for variable_name, factor in candidate.factors.items():
        sigma_log = float(
            config["variables"][variable_name].get(
                "prior_sigma_log", math.log(prior_scale)
            )
        )
        if sigma_log <= 0:
            raise ValueError(f"prior_sigma_log must be positive for {variable_name}")
        prior += (math.log(factor) / sigma_log) ** 2
    prior_weight = float(score_settings["prior_weight"])
    objective = float(best["aligned_objective"]) + prior_weight * prior

    if trace_path is not None:
        valid = best["valid"]
        write_aligned_trace(
            trace_path,
            grid[valid],
            ref_norm[valid],
            best["simulation_normalized"],
            regions[valid],
        )
    return {
        "candidate_id": candidate.candidate_id,
        "factors": candidate.factors,
        "series": str(series_path),
        "aligned_trace": str(trace_path) if trace_path is not None else None,
        "comparison_window_ms": [window_start, window_end],
        "best_shift_ms": best["shift_ms"],
        "shift_penalty": best["shift_penalty"],
        "waveform_objective": best["waveform_objective"],
        "region_rmse": best["region_rmse"],
        "normalized_rmse": best["rmse"],
        "normalized_mae": best["mae"],
        "normalized_rmse_no_shift": no_shift_rmse,
        "waveform_objective_no_shift": no_shift_waveform_objective,
        "region_rmse_no_shift": no_shift_region_rmse,
        "prior_penalty": prior,
        "objective": objective,
        "reference": {
            key: value
            for key, value in ref_metrics.items()
            if key not in {"response_signal", "normalized"}
        },
        "simulation": {
            key: value
            for key, value in sim_metrics.items()
            if key not in {"response_signal", "normalized"}
        },
    }


def run_candidate(config: dict[str, Any], candidate: Candidate) -> dict[str, Any]:
    base_project = load_json(resolve_from_root(config["base_project"]))
    if config["g0_calibration"]["enabled"]:
        baseline_calibration_path = (
            candidate_dir(config, "baseline") / "calibration.json"
        )
        baseline_g0_factor: float | None = None
        if baseline_calibration_path.exists():
            baseline_g0_factor = float(
                load_json(baseline_calibration_path)["g0_factor"]
            )

        if (
            candidate.candidate_id != "baseline"
            and not candidate.steady_sensitive
            and baseline_g0_factor is not None
        ):
            # Heat-capacity changes do not alter the mathematical steady
            # equilibrium.  Reuse the calibrated coarse-mesh G0, but still run
            # and verify this candidate's own steady dependency before pulse.
            error, project_path, metadata, _ = evaluate_steady_trial(
                base_project,
                config,
                candidate,
                baseline_g0_factor,
                0,
            )
            tolerance_K = (
                float(config["g0_calibration"]["temperature_tolerance_mK"])
                * 1e-3
            )
            if abs(error) > tolerance_K:
                project_path, metadata = calibrate_g0(
                    base_project,
                    config,
                    candidate,
                    initial_factor=baseline_g0_factor,
                )
            else:
                write_json(
                    candidate_dir(config, candidate.candidate_id)
                    / "calibration.json",
                    {
                        "candidate_id": candidate.candidate_id,
                        "g0_factor": baseline_g0_factor,
                        "temperature_error_mK": error * 1e3,
                        "project": str(project_path.relative_to(ROOT)),
                        "steady_case": metadata["steady_case"],
                        "reused_baseline_g0": True,
                    },
                )
        else:
            project_path, metadata = calibrate_g0(
                base_project,
                config,
                candidate,
                initial_factor=baseline_g0_factor,
            )
    else:
        project_path, metadata = write_trial(
            base_project,
            config,
            candidate,
            float(config["g0_calibration"]["initial_factor"]),
            0,
        )
        run_command(project_path, metadata["steady_case"])

    # The calibrated steady result already exists, so run.py will reuse it.
    run_command(project_path, metadata["pulse_case"])
    series_path = ROOT / "results" / metadata["pulse_case"] / metadata["pulse_series"]
    trace_path = candidate_dir(config, candidate.candidate_id) / "aligned_waveform.csv"
    score = score_series(
        config,
        series_path,
        candidate,
        float(config["pulse_start_ms"]),
        trace_path=trace_path,
    )
    score.update(
        {
            "project": str(project_path.relative_to(ROOT)),
            "steady_case": metadata["steady_case"],
            "pulse_case": metadata["pulse_case"],
            "g0_factor": metadata["g0_factor"],
        }
    )
    write_json(candidate_dir(config, candidate.candidate_id) / "score.json", score)
    update_leaderboard(config)
    return score


def update_leaderboard(config: dict[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for score_path in output_root(config).glob("candidates/*/score.json"):
        score = load_json(score_path)
        rows.append(
            {
                "candidate_id": score["candidate_id"],
                "objective": score["objective"],
                "waveform_objective": score.get("waveform_objective", ""),
                "normalized_rmse": score["normalized_rmse"],
                "normalized_rmse_no_shift": score["normalized_rmse_no_shift"],
                "best_shift_ms": score["best_shift_ms"],
                "g0_factor": score.get("g0_factor", ""),
                "factors": json.dumps(score["factors"], sort_keys=True),
            }
        )
    rows.sort(key=lambda row: float(row["objective"]))
    if not rows:
        return
    with (output_root(config) / "leaderboard.csv").open(
        "w", newline="", encoding="utf-8"
    ) as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def score_existing(config: dict[str, Any]) -> dict[str, Any]:
    existing = (
        ROOT
        / "results"
        / "raw"
        / "legacy-root-output"
        / "tes_mpi_comsol_grid_series.csv"
    )
    score = score_series(
        config,
        existing,
        Candidate("existing_baseline", {}, False),
        float(config["reference_pulse_start_ms"]),
        trace_path=output_root(config) / "existing_baseline_aligned_waveform.csv",
    )
    write_json(output_root(config) / "existing_baseline_score.json", score)
    return score


def read_aligned_trace(path: Path) -> dict[str, np.ndarray]:
    rows: list[dict[str, str]] = []
    with path.open("r", newline="", encoding="utf-8") as file:
        rows.extend(csv.DictReader(file))
    if not rows:
        raise ValueError(f"aligned trace is empty: {path}")
    return {
        "time_ms": np.array([float(row["time_after_pulse_ms"]) for row in rows]),
        "simulation": np.array([float(row["simulation_normalized"]) for row in rows]),
        "reference": np.array([float(row["reference_normalized"]) for row in rows]),
        "region": np.array([row["region"] for row in rows], dtype="U8"),
    }


def score_trace_path(score: dict[str, Any]) -> Path:
    value = score.get("aligned_trace")
    if not value:
        raise ValueError(
            f"score for {score.get('candidate_id', 'unknown')} has no aligned trace; rerun scoring"
        )
    path = Path(str(value))
    return path if path.is_absolute() else ROOT / path


def analyze_sensitivity(config: dict[str, Any]) -> dict[str, Path]:
    root = output_root(config)
    summaries: list[dict[str, Any]] = []
    sensitivities: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    for variable_name, spec in config["variables"].items():
        factors = sorted(float(value) for value in spec["factors"])
        if len(factors) < 2:
            continue
        low_factor, high_factor = factors[0], factors[-1]
        low_id = f"{variable_name}_{factor_tag(low_factor)}"
        high_id = f"{variable_name}_{factor_tag(high_factor)}"
        low_score_path = candidate_dir(config, low_id) / "score.json"
        high_score_path = candidate_dir(config, high_id) / "score.json"
        if not low_score_path.exists() or not high_score_path.exists():
            continue

        low_score = load_json(low_score_path)
        high_score = load_json(high_score_path)
        low_trace = read_aligned_trace(score_trace_path(low_score))
        high_trace = read_aligned_trace(score_trace_path(high_score))
        time_start = max(float(low_trace["time_ms"].min()), float(high_trace["time_ms"].min()))
        time_end = min(float(low_trace["time_ms"].max()), float(high_trace["time_ms"].max()))
        common_time = low_trace["time_ms"]
        common_time = common_time[(common_time >= time_start) & (common_time <= time_end)]
        high_simulation = np.interp(
            common_time, high_trace["time_ms"], high_trace["simulation"]
        )
        low_simulation = np.interp(
            common_time, low_trace["time_ms"], low_trace["simulation"]
        )
        denominator = math.log(high_factor) - math.log(low_factor)
        sensitivity = (high_simulation - low_simulation) / denominator
        sensitivities[variable_name] = (common_time, sensitivity)

        region_labels = np.array(
            [
                low_trace["region"][
                    int(np.argmin(np.abs(low_trace["time_ms"] - time_value)))
                ]
                for time_value in common_time
            ],
            dtype="U8",
        )
        region_norms: dict[str, float] = {}
        for region_name in ("rise", "peak", "decay", "tail"):
            mask = region_labels == region_name
            region_norms[f"{region_name}_sensitivity_rms"] = (
                float(np.sqrt(np.mean(sensitivity[mask] ** 2)))
                if np.count_nonzero(mask)
                else math.nan
            )
        summaries.append(
            {
                "variable": variable_name,
                "low_factor": low_factor,
                "high_factor": high_factor,
                "waveform_sensitivity_rms": float(np.sqrt(np.mean(sensitivity**2))),
                "objective_log_slope": (
                    float(high_score["objective"]) - float(low_score["objective"])
                )
                / denominator,
                "low_objective": float(low_score["objective"]),
                "high_objective": float(high_score["objective"]),
                **region_norms,
            }
        )

    if not summaries:
        raise RuntimeError(
            "no complete low/high OFAT score pairs were found; run the first-stage candidates first"
        )

    summary_path = root / "sensitivity_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=list(summaries[0].keys()))
        writer.writeheader()
        writer.writerows(
            sorted(summaries, key=lambda row: row["waveform_sensitivity_rms"], reverse=True)
        )

    names = sorted(sensitivities)
    correlation_path = root / "sensitivity_correlation.csv"
    with correlation_path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["variable", *names])
        for left_name in names:
            left_time, left_values = sensitivities[left_name]
            row: list[Any] = [left_name]
            for right_name in names:
                right_time, right_values = sensitivities[right_name]
                start = max(float(left_time.min()), float(right_time.min()))
                end = min(float(left_time.max()), float(right_time.max()))
                grid = left_time[(left_time >= start) & (left_time <= end)]
                left_interp = np.interp(grid, left_time, left_values)
                right_interp = np.interp(grid, right_time, right_values)
                if np.std(left_interp) == 0 or np.std(right_interp) == 0:
                    correlation = math.nan
                else:
                    correlation = float(np.corrcoef(left_interp, right_interp)[0, 1])
                row.append(correlation)
            writer.writerow(row)
    return {"summary": summary_path, "correlation": correlation_path}


def print_candidate_table(candidates: list[Candidate]) -> None:
    print(f"prepared {len(candidates)} candidates")
    for candidate in candidates:
        print(
            f"  {candidate.candidate_id:20s} "
            f"steady_sensitive={str(candidate.steady_sensitive):5s} "
            f"factors={candidate.factors}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("prepare")
    prepare_multi_parser = subparsers.add_parser("prepare-multivariate")
    prepare_multi_parser.add_argument("--samples", type=int)
    prepare_multi_parser.add_argument("--seed", type=int)
    subparsers.add_parser("score-existing")
    subparsers.add_parser("analyze-sensitivity")
    dry_parser = subparsers.add_parser("dry-run")
    dry_parser.add_argument("candidate")
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("candidate")
    run_all_parser = subparsers.add_parser("run-all")
    run_all_parser.add_argument("--limit", type=int)
    run_all_parser.add_argument(
        "--skip-baseline",
        action="store_true",
        help="run only the one-factor candidates, preserving an existing baseline calibration",
    )
    run_all_parser.add_argument(
        "--skip-scored",
        action="store_true",
        help="skip candidates that already have score.json",
    )
    run_multi_parser = subparsers.add_parser("run-multivariate")
    run_multi_parser.add_argument("--limit", type=int)
    run_multi_parser.add_argument("--samples", type=int)
    run_multi_parser.add_argument("--seed", type=int)
    args = parser.parse_args()

    config_path = args.config if args.config.is_absolute() else ROOT / args.config
    config = load_json(config_path)
    candidates = candidates_from_config(config)
    by_id = {candidate.candidate_id: candidate for candidate in candidates}

    if args.command == "prepare":
        print_candidate_table(prepare(config))
        return 0
    if args.command == "prepare-multivariate":
        raise SystemExit(
            "single-pixel multivariate fitting is disabled: use "
            "scripts/search/post_multivariate_search.py prepare"
        )
    if args.command == "score-existing":
        score = score_existing(config)
        print(json.dumps(score, indent=2, ensure_ascii=False))
        return 0
    if args.command == "analyze-sensitivity":
        paths = analyze_sensitivity(config)
        print(json.dumps({key: str(value) for key, value in paths.items()}, indent=2))
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
        run_command(project_path, metadata["pulse_case"], dry_run=True)
        return 0
    if args.command == "run":
        if args.candidate not in by_id:
            raise SystemExit(f"unknown candidate: {args.candidate}")
        score = run_candidate(config, by_id[args.candidate])
        print(json.dumps(score, indent=2, ensure_ascii=False))
        return 0
    if args.command == "run-all":
        available = candidates[1:] if args.skip_baseline else candidates
        if args.skip_scored:
            available = [
                candidate
                for candidate in available
                if not (candidate_dir(config, candidate.candidate_id) / "score.json").exists()
            ]
        selected = available[: args.limit] if args.limit is not None else available
        for index, candidate in enumerate(selected, start=1):
            print(f"[{index}/{len(selected)}] {candidate.candidate_id}", flush=True)
            run_candidate(config, candidate)
        return 0
    if args.command == "run-multivariate":
        raise SystemExit(
            "single-pixel multivariate fitting is disabled: use "
            "scripts/search/post_multivariate_search.py run-all"
        )
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
