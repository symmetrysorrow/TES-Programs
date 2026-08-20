"""Standalone PoST optimizer: differential evolution -> Powell.

For each trial, this script writes H:/hata2025/new/input.json, runs
PoST_Simulation.py, then compares the newly generated
noise_total-bessel100k.dat with tagawa CH0_noise/modelnoise.txt from 1--30 kHz.
The .dat file must include the 100 kHz hardware Bessel and a 10 kHz analysis
Bessel, so each candidate fixes input.json["cutoff"] to 10000.

By default input.json is restored after optimization.  Add --apply-final to
keep the best parameters and run PoST one final time.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, minimize


# ---------- Paths ----------
POST_SCRIPT = Path(r"D:\Github\TES-Programs\PoST_Simulations\PoST_Simulation.py")
INPUT_PATH = Path(r"H:\hata2025\new\input.json")
# This is the original, constraint-defining input.  Do not use a previously
# optimized input.json as the reference, or R's 50--200% range will compound.
DEFAULT_REFERENCE_INPUT_PATH = Path(r"H:\hata2025\new\input.before_noise_optimization.20260818-011257.json")
NOISE_DAT_PATH = Path(r"H:\hata2025\new\noise_total-bessel100k.dat")
EXPERIMENT_ROOT = Path(
    r"G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_"
    r"difftrig5e-5_rate500k_samples100k_gain5_day2"
)
MODELNOISE_PATH = EXPERIMENT_ROOT / "CH0_noise" / "modelnoise.txt"
PULSE_CONFIG_PATH = EXPERIMENT_ROOT / "PulseConfig.json"

# ---------- Comparison settings ----------
# 1--10 kHz already agrees reasonably well; optimize the remaining mismatch.
FIT_MIN_HZ = 1000.0
FIT_MAX_HZ = 30000.0
SIM_ANALYSIS_CUTOFF_HZ = 10000  # modelnoise.txt is filtered at this cutoff.
OPTIMIZATION_SAMPLES = 4096


@dataclass(frozen=True)
class Bound:
    lower: float
    upper: float
    logarithmic: bool = True


def arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("--de-maxiter", type=int, default=8)
    parser.add_argument("--de-popsize", type=int, default=3)
    parser.add_argument("--powell-maxfev", type=int, default=250)
    parser.add_argument("--skip-de", action="store_true")
    parser.add_argument("--apply-final", action="store_true")
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--reference-input", type=Path, default=DEFAULT_REFERENCE_INPUT_PATH,
                        help="Original JSON that defines all parameter bounds and fixed values.")
    return parser.parse_args()


def load_json(path: Path):
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def write_json_atomically(path: Path, data: dict):
    temporary = path.with_suffix(path.suffix + ".optimize-tmp")
    with temporary.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=4)
        f.write("\n")
    temporary.replace(path)


def normalize_at(freq, asd, reference_hz=1000.0):
    value = asd[np.argmin(np.abs(freq - reference_hz))]
    if not np.isfinite(value) or value <= 0:
        raise ValueError(f"Invalid ASD at {reference_hz} Hz")
    return asd / value


def target_spectrum():
    pulse_config = load_json(PULSE_CONFIG_PATH)
    exp_asd = np.loadtxt(MODELNOISE_PATH)
    exp_rate = float(pulse_config["Readout"]["Rate"])
    exp_freq = np.arange(len(exp_asd)) * (exp_rate / 2.0) / len(exp_asd)
    fit_freq = np.geomspace(FIT_MIN_HZ, FIT_MAX_HZ, 121)
    target = np.interp(fit_freq, exp_freq, normalize_at(exp_freq, exp_asd))
    return fit_freq, target


def run_post(timeout_seconds: int, output_dir: Path, noise_path: Path):
    before = noise_path.stat().st_mtime_ns if noise_path.exists() else -1
    result = subprocess.run(
        [
            sys.executable,
            str(POST_SCRIPT),
            "--noise-only",
            "--output",
            str(output_dir),
        ],
        cwd=POST_SCRIPT.parent,
        capture_output=True,
        text=True,
        timeout=timeout_seconds,
    )
    if result.returncode != 0:
        raise RuntimeError("PoST failed:\n" + result.stderr[-4000:])
    if not noise_path.exists() or noise_path.stat().st_mtime_ns <= before:
        raise RuntimeError(
            "noise_total-bessel100k.dat was not updated. "
            "PoST_Simulation.py must end with MakeNoise() followed by _ShowNoiseSpectrum()."
        )


def simulated_spectrum(candidate: dict, fit_freq: np.ndarray, noise_path: Path):
    sim_asd = np.loadtxt(noise_path)
    sim_rate = float(candidate["rate"])
    sim_freq = np.arange(len(sim_asd)) * (sim_rate / 2.0) / len(sim_asd)
    return np.interp(fit_freq, sim_freq, normalize_at(sim_freq, sim_asd))


def parameter_bounds(original: dict):
    """Only parameters the user has allowed to vary are listed here."""
    return {
        "C_abs": Bound(original["C_abs"] * 1e-3, original["C_abs"] * 1e3),
        "C_tes": Bound(original["C_tes"] * 1e-3, original["C_tes"] * 1e3),
        "G_abs-abs": Bound(original["G_abs-abs"] * 1e-3, original["G_abs-abs"] * 1e3),
        "G_abs-tes": Bound(original["G_abs-tes"] * 1e-3, original["G_abs-tes"] * 1e3),
        "G_tes-bath": Bound(original["G_tes-bath"] * 1e-3, original["G_tes-bath"] * 1e3),
        "R": Bound(original["R"] * 0.5, original["R"] * 2.0),
        "T_bath": Bound(1e-6, original["T_c"]),
        "alpha": Bound(50.0, 1000.0),
        "beta": Bound(0.0, 5.0, logarithmic=False),
        # L is now intentionally variable.  Tighten this upper bound if known.
        "L": Bound(original["L"], 1e-5),
    }


def vector_bounds(bounds: dict):
    keys = list(bounds)
    transformed = []
    for key in keys:
        bound = bounds[key]
        if bound.logarithmic:
            transformed.append((np.log10(bound.lower), np.log10(bound.upper)))
        else:
            transformed.append((bound.lower, bound.upper))
    return keys, transformed


def encode(candidate: dict, keys, bounds):
    return np.asarray([
        np.log10(candidate[key]) if bounds[key].logarithmic else candidate[key]
        for key in keys
    ], dtype=float)


def decode(
    vector,
    original: dict,
    keys,
    bounds,
    sample_rate: float,
    post_filter_white_asd: float,
):
    candidate = original.copy()
    for key, value in zip(keys, vector):
        candidate[key] = float(10.0 ** value if bounds[key].logarithmic else value)
    # Match the experimental modelnoise.txt analysis filter; this is not fitted.
    candidate["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    # Filter coefficients and the frequency grid must use the acquisition rate
    # of modelnoise.txt, even when an old optimization reference used 300 kHz.
    candidate["rate"] = sample_rate
    candidate["post_filter_white_asd_A_rtHz"] = post_filter_white_asd
    candidate.pop("readout_white_asd_A_rtHz", None)
    return candidate


def main():
    args = arguments()
    for path in (POST_SCRIPT, INPUT_PATH, MODELNOISE_PATH, PULSE_CONFIG_PATH, args.reference_input):
        if not path.is_file():
            raise FileNotFoundError(f"Missing required file: {path}")

    original = load_json(INPUT_PATH)  # state to restore when --apply-final is absent
    reference = load_json(args.reference_input)
    pulse_config = load_json(PULSE_CONFIG_PATH)
    experimental_rate = float(pulse_config["Readout"]["Rate"])
    post_filter_white_asd = float(
        original.get(
            "post_filter_white_asd_A_rtHz",
            original.get("readout_white_asd_A_rtHz", 0.0),
        )
    )
    backup = INPUT_PATH.with_name(f"input.before_scipy_optimization.{time.strftime('%Y%m%d-%H%M%S')}.json")
    shutil.copy2(INPUT_PATH, backup)
    print("backup:", backup)
    work_dir = INPUT_PATH.parent / ".noise_optimization_work"
    work_dir.mkdir(exist_ok=True)
    work_input_path = work_dir / "input.json"
    work_noise_path = work_dir / NOISE_DAT_PATH.name

    bounds = parameter_bounds(reference)
    keys, scipy_bounds = vector_bounds(bounds)
    fit_freq, target = target_spectrum()
    initial = original.copy()
    initial["cutoff"] = SIM_ANALYSIS_CUTOFF_HZ
    initial["rate"] = experimental_rate
    initial["post_filter_white_asd_A_rtHz"] = post_filter_white_asd
    initial.pop("readout_white_asd_A_rtHz", None)
    initial["samples"] = OPTIMIZATION_SAMPLES
    initial_x = encode(initial, keys, bounds)

    cache = {}
    evaluation_count = 0
    best_score = np.inf
    best_candidate = initial.copy()

    def objective(vector):
        nonlocal evaluation_count, best_score, best_candidate
        cache_key = tuple(np.round(vector, 12))
        if cache_key in cache:
            return cache[cache_key]
        candidate = decode(
            vector,
            original,
            keys,
            bounds,
            experimental_rate,
            post_filter_white_asd,
        )
        candidate["samples"] = OPTIMIZATION_SAMPLES
        evaluation_count += 1
        try:
            write_json_atomically(work_input_path, candidate)
            run_post(args.timeout, work_dir, work_noise_path)
            model = simulated_spectrum(candidate, fit_freq, work_noise_path)
            score = float(np.mean((np.log10(model) - np.log10(target)) ** 2))
            print(f"evaluation {evaluation_count:4d}: {score:.6g}")
        except Exception as error:
            print(f"evaluation {evaluation_count:4d}: rejected ({error})")
            score = 1e12
        cache[cache_key] = score
        if score < best_score:
            best_score = score
            best_candidate = candidate.copy()
            print("  new best score:", best_score)
        return score

    try:
        print("Initial evaluation")
        objective(initial_x)

        if args.skip_de:
            starting_x = initial_x
        else:
            print("\nStage 1/2: differential_evolution")
            global_result = differential_evolution(
                objective,
                scipy_bounds,
                maxiter=args.de_maxiter,
                popsize=args.de_popsize,
                seed=args.seed,
                workers=1,       # input.json and output directory are shared.
                updating="immediate",
                polish=False,
            )
            starting_x = global_result.x
            print("DE result:", global_result.fun)

        print("\nStage 2/2: Powell")
        local_result = minimize(
            objective,
            starting_x,
            method="Powell",
            bounds=scipy_bounds,
            options={"maxfev": args.powell_maxfev, "disp": True},
        )
        objective(local_result.x)  # ensures the returned point is cached/tracked

        print("\nPoST evaluations:", evaluation_count)
        print("Best score:", best_score)
        print("Best fitted parameters:")
        print(json.dumps({key: best_candidate[key] for key in keys}, indent=2))

        if args.apply_final:
            final_candidate = best_candidate.copy()
            final_candidate["samples"] = original["samples"]
            write_json_atomically(INPUT_PATH, final_candidate)
            run_post(args.timeout, INPUT_PATH.parent, NOISE_DAT_PATH)
            print("Best input.json applied and final .dat regenerated.")
        else:
            print("Original input.json restored; use --apply-final to keep the best candidate.")
    except Exception:
        write_json_atomically(INPUT_PATH, original)
        print("Original input.json restored after error.")
        raise


if __name__ == "__main__":
    main()
