"""Profile the line-robust continuum fit against the c2 upper bound.

This diagnostic keeps the line-repaired 1--200 kHz continuum target, detector
nuisance family, R_TES, white-floor profiling procedure, optimizer settings,
and readout parameterization fixed.  Only the allowed upper bound of c2 is
changed across a tracked grid (default: 300, 600, 1000).

The purpose is to distinguish:
  * a merely too-tight c2 bound (solution settles at a finite interior value),
  * from c2 acting as a surrogate for missing readout physics (solution keeps
    following the expanded upper bound).

No production input is modified.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "c2_bound_profile_continuum_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "c2_bound_profile_continuum_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "c2_bound_profile_continuum_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import line_robust_continuum_fit_diagnostic as continuum  # noqa: E402


def _float(value):
    return float(value)


def _at_upper(value, upper, *, relative_tolerance=1.0e-4):
    scale = max(abs(float(upper)), 1.0)
    return bool(abs(float(value) - float(upper)) <= scale * relative_tolerance)


def _profile_row(c2_upper, result):
    recommended = result["recommended"]
    selected_family = recommended["family"]
    fit_row = result["fits"][selected_family]
    c2_value = float(recommended["readout"]["c2"])

    return {
        "c2_upper_bound": float(c2_upper),
        "selected_family": selected_family,
        "c2_value": c2_value,
        "c2_fraction_of_upper_bound": float(c2_value / float(c2_upper)),
        "c2_at_upper_bound": _at_upper(c2_value, c2_upper),
        "pole_Hz": float(recommended["readout"]["pole_Hz"]),
        "pole_Q": float(recommended["readout"]["pole_Q"]),
        "c4": float(recommended["readout"]["c4"]),
        "white_asd_A_rtHz": float(recommended["white_asd_A_rtHz"]),
        "continuum_rms_dB": float(recommended["continuum_rms_dB"]),
        "continuum_1_40k_rms_dB": float(
            recommended["continuum_1_40k_rms_dB"]
        ),
        "continuum_40_200k_rms_dB": float(
            recommended["continuum_40_200k_rms_dB"]
        ),
        "shape_score": float(fit_row["shape_score"]),
        "detector_candidate": {
            key: float(value)
            for key, value in recommended["detector_candidate"].items()
        },
        "selection": result["selection"],
    }


def classify_profile(rows, config):
    if len(rows) < 2:
        return {
            "classification": "insufficient_grid",
            "reason": "Need at least two c2 upper bounds.",
        }

    last = rows[-1]
    previous = rows[-2]
    all_at_bound = all(row["c2_at_upper_bound"] for row in rows)
    last_at_bound = bool(last["c2_at_upper_bound"])

    settle_margin = float(
        config["interpretation"]["interior_fraction_threshold"]
    )
    last_is_interior = bool(
        last["c2_fraction_of_upper_bound"] <= settle_margin
    )

    rms_gain = float(
        rows[0]["continuum_rms_dB"] - last["continuum_rms_dB"]
    )
    last_step_gain = float(
        previous["continuum_rms_dB"] - last["continuum_rms_dB"]
    )
    material_gain = float(
        config["interpretation"]["material_total_rms_improvement_dB"]
    )
    negligible_step = float(
        config["interpretation"]["negligible_last_step_rms_improvement_dB"]
    )

    if last_is_interior:
        classification = "finite_interior_c2_supported"
        reason = (
            "The solution no longer follows the expanded c2 ceiling; "
            "the previous c2<=300 constraint was likely too tight."
        )
    elif all_at_bound and last_step_gain > negligible_step:
        classification = "c2_tracks_bound_missing_shape_freedom_candidate"
        reason = (
            "c2 remains ceiling-limited and the newest expansion still "
            "improves the continuum fit, consistent with c2 absorbing "
            "missing smooth readout shape freedom."
        )
    elif last_at_bound and last_step_gain <= negligible_step:
        classification = "c2_bound_limited_but_fit_saturated"
        reason = (
            "c2 remains at the ceiling, but the newest expansion gives "
            "negligible RMS improvement; more c2 alone is not buying much."
        )
    else:
        classification = "mixed_or_transitioning_profile"
        reason = (
            "The c2 profile is neither clearly interior nor monotonically "
            "ceiling-limited across the tested grid."
        )

    return {
        "classification": classification,
        "reason": reason,
        "all_tested_points_at_upper_bound": bool(all_at_bound),
        "highest_bound_point_at_upper_bound": bool(last_at_bound),
        "highest_bound_fraction": float(
            last["c2_fraction_of_upper_bound"]
        ),
        "total_rms_improvement_first_to_last_dB": rms_gain,
        "last_step_rms_improvement_dB": last_step_gain,
        "material_total_improvement": bool(rms_gain >= material_gain),
        "thresholds": config["interpretation"],
    }


def run(config, config_path: Path, experiment_path_override=None):
    base_continuum_config_path = continuum.resolve_config_path(
        config["base_continuum_config"],
        config_path,
    )
    base_config = json.loads(
        base_continuum_config_path.read_text(encoding="utf-8")
    )

    base_residual_config_path = continuum.resolve_config_path(
        base_config["base_residual_config"],
        base_continuum_config_path,
    )
    residual_config = json.loads(
        base_residual_config_path.read_text(encoding="utf-8")
    )

    grid = [float(value) for value in config["c2_upper_bounds"]]
    if len(grid) < 2 or any(value <= 0.0 for value in grid):
        raise ValueError("c2_upper_bounds must contain at least two positive values")
    if grid != sorted(set(grid)):
        raise ValueError("c2_upper_bounds must be strictly increasing and unique")

    rows = []
    full_results = {}
    for upper in grid:
        working_residual = json.loads(json.dumps(residual_config))
        working_residual["readout_bounds"]["c2"] = [0.0, float(upper)]

        # Write a temporary tracked-equivalent config next to the output so the
        # existing continuum fitter can consume it without changing production
        # config files.
        temp_dir = DEFAULT_OUTPUT.parent
        temp_dir.mkdir(parents=True, exist_ok=True)
        temp_residual_path = (
            temp_dir / f"_c2_profile_residual_bound_{int(upper)}.json"
        )
        temp_continuum_path = (
            temp_dir / f"_c2_profile_continuum_bound_{int(upper)}.json"
        )
        temp_residual_path.write_text(
            json.dumps(working_residual, indent=2) + "\n",
            encoding="utf-8",
        )

        working_continuum = json.loads(json.dumps(base_config))
        working_continuum["base_residual_config"] = str(
            temp_residual_path.resolve()
        )
        temp_continuum_path.write_text(
            json.dumps(working_continuum, indent=2) + "\n",
            encoding="utf-8",
        )

        try:
            result = continuum.run(
                working_continuum,
                temp_continuum_path,
                experiment_path_override=experiment_path_override,
            )
        finally:
            for path in (temp_residual_path, temp_continuum_path):
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass

        row = _profile_row(upper, result)
        rows.append(row)
        full_results[str(int(upper))] = result

    interpretation = classify_profile(rows, config)

    return {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "base_continuum_config": str(base_continuum_config_path),
        "base_residual_config": str(base_residual_config_path),
        "only_profiled_change": "readout_bounds.c2 upper bound",
        "c2_upper_bounds": grid,
        "profile": rows,
        "interpretation": interpretation,
        "guardrail": config["guardrail"],
        "_plot": {
            "bounds": [row["c2_upper_bound"] for row in rows],
            "c2": [row["c2_value"] for row in rows],
            "rms": [row["continuum_rms_dB"] for row in rows],
            "rms_low": [row["continuum_1_40k_rms_dB"] for row in rows],
            "rms_high": [row["continuum_40_200k_rms_dB"] for row in rows],
            "pole": [row["pole_Hz"] for row in rows],
            "q": [row["pole_Q"] for row in rows],
        },
    }


def make_plot(result, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    bounds = np.asarray(p["bounds"], dtype=float)
    c2 = np.asarray(p["c2"], dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(9.5, 8.5), sharex=True)
    top, middle, bottom = axes

    top.plot(bounds, c2, marker="o", label="fitted c2")
    top.plot(bounds, bounds, linestyle="--", label="c2 upper bound")
    top.set_ylabel("c2")
    top.set_title("Does c2 settle, or follow the bound?")
    top.legend(frameon=False)
    top.grid(True, alpha=0.25)

    middle.plot(bounds, p["rms"], marker="o", label="1-200 kHz")
    middle.plot(bounds, p["rms_low"], marker="o", label="1-40 kHz")
    middle.plot(bounds, p["rms_high"], marker="o", label="40-200 kHz")
    middle.set_ylabel("Continuum RMS [dB]")
    middle.legend(frameon=False)
    middle.grid(True, alpha=0.25)

    bottom.plot(bounds, np.asarray(p["pole"]) / 1000.0, marker="o", label="pole [kHz]")
    bottom.plot(bounds, p["q"], marker="o", label="Q")
    bottom.set_xlabel("c2 upper bound")
    bottom.set_ylabel("Readout coordinates")
    bottom.legend(frameon=False)
    bottom.grid(True, alpha=0.25)

    fig.suptitle(
        "Line-robust continuum c2-bound profile | "
        + result["interpretation"]["classification"],
        fontsize=10,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result):
    return {
        key: value
        for key, value in result.items()
        if key != "_plot"
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--experiment-path", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(
        config,
        args.config,
        experiment_path_override=args.experiment_path,
    )
    make_plot(result, args.figure, show=args.show)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(cleaned_result(result), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "profile": result["profile"],
                "interpretation": result["interpretation"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
