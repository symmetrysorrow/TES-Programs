"""Test the documented Magnicon XXF-1 10 kHz output LPF against the continuum fit.

The XXF-1 manual documents a Connector Box anti-alias filter that is second-order
Bessel with nominal 10 kHz cutoff and +/-2.5% tolerance.  This diagnostic asks a
minimal physical question: can that known filter replace the previously free
~12 kHz effective pole without materially degrading the line-repaired 1--200 kHz
continuum fit?

The measurement-time LPF state is not known from the noise files.  Therefore a
successful fit is only consistency evidence for the LPF-ON hypothesis, not proof
that the switch was ON during the acquisition.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "magnicon_xxf1_lpf_continuum_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_continuum_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_continuum_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import line_robust_continuum_fit_diagnostic as continuum  # noqa: E402
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402
from subScript import readout_effective_numerator_diagnostic as effective  # noqa: E402
from subScript import readout_lowmid_holdout_diagnostic as holdout  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


def second_order_bessel_canonical(cutoff_hz: float, norm: str) -> dict:
    """Return the canonical pole/Q coordinates for an analog 2nd-order Bessel.

    ``scipy.signal.bessel`` returns H(s)=B(s)/A(s).  After normalizing A(s) by
    its constant term, the denominator is written as

        (s/w0)^2 + s/(Q*w0) + 1,

    which is exactly the pole denominator used by the effective readout model.
    """
    cutoff = float(cutoff_hz)
    if cutoff <= 0.0:
        raise ValueError("cutoff_hz must be positive")
    if norm not in {"mag", "phase", "delay"}:
        raise ValueError("unsupported Bessel normalization")

    b, a = signal.bessel(
        2,
        2.0 * np.pi * cutoff,
        btype="low",
        analog=True,
        output="ba",
        norm=norm,
    )
    b = np.asarray(b, dtype=float)
    a = np.asarray(a, dtype=float)
    if a.size != 3:
        raise RuntimeError("second-order Bessel denominator is not quadratic")

    # Normalize to a2=1 so the constant denominator term is unity.
    a0, a1, a2 = (float(value) for value in a)
    if a0 <= 0.0 or a1 <= 0.0 or a2 <= 0.0:
        raise RuntimeError("unexpected non-positive Bessel denominator")
    omega0 = float(np.sqrt(a2 / a0))
    q = float(np.sqrt(a0 * a2) / a1)
    pole_hz = float(omega0 / (2.0 * np.pi))
    dc_gain = float(b[-1] / a2)
    return {
        "cutoff_Hz": cutoff,
        "norm": norm,
        "pole_Hz": pole_hz,
        "pole_Q": q,
        "pole_over_cutoff": float(pole_hz / cutoff),
        "dc_gain": dc_gain,
        "denominator_coefficients": a.tolist(),
        "numerator_coefficients": b.tolist(),
    }


def equivalent_cutoff_from_pole(pole_hz: float, canonical: dict) -> float:
    ratio = float(canonical["pole_over_cutoff"])
    if ratio <= 0.0:
        raise ValueError("invalid pole/cutoff ratio")
    return float(pole_hz) / ratio


def _problem(config: dict, config_path: Path, experiment_path_override=None):
    continuum_config_path = continuum.resolve_config_path(
        config["base_continuum_config"],
        config_path,
    )
    continuum_config = json.loads(
        continuum_config_path.read_text(encoding="utf-8")
    )
    residual_config_path = continuum.resolve_config_path(
        continuum_config["base_residual_config"],
        continuum_config_path,
    )
    residual_config = json.loads(
        residual_config_path.read_text(encoding="utf-8")
    )

    manifest_path = residual.resolve_config_path(
        residual_config["manifest"],
        residual_config_path,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = shared.normalize_manifest(manifest, manifest_path)
    repeat_case = ident.case_by_label(
        cases,
        residual_config["repeat_case_label"],
    )
    summary = json.loads(repeat_case["summary"].read_text(encoding="utf-8"))
    comparison, experiment_path, comparison_source = ident.comparison_for_case(
        repeat_case
    )
    if experiment_path_override is not None:
        experiment_path = Path(experiment_path_override)

    snapshot_path = residual.resolve_config_path(
        residual_config["residual_baseline_snapshot"],
        residual_config_path,
    )
    snapshot = residual.load_baseline_snapshot(snapshot_path)

    inherited_candidate = dict(summary["best_case_parameters"])
    baseline_detector = dict(inherited_candidate)
    baseline_detector["R"] = float(snapshot["R_TES_Ohm"])
    baseline_detector.update(snapshot["detector_candidate"])

    full_args = base.fit_args(summary)
    fit_min = float(continuum_config["fit_region_Hz"]["min"])
    fit_max = float(continuum_config["fit_region_Hz"]["max"])
    frequency = np.geomspace(fit_min, fit_max, int(full_args.fit_points))
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    repair = continuum.repair_narrow_lines(
        target_context["full_frequency_Hz"],
        target_context["pre_analysis_asd"],
        continuum_config["line_repair"],
    )
    raw_target = np.interp(
        frequency,
        target_context["full_frequency_Hz"],
        repair["raw_normalized"],
    )
    target = np.interp(
        frequency,
        target_context["full_frequency_Hz"],
        repair["continuum_normalized"],
    )
    fit_args = competition.band_args(
        full_args,
        fit_min,
        fit_max,
        len(frequency),
    )

    detector_names = tuple(
        residual_config["detector_nuisance"]["parameters"]
    )
    detector_bounds = competition.nuisance_bounds(
        inherited_candidate,
        {"bounds": residual_config["detector_bounds"]},
    )
    readout_bounds = {
        name: tuple(float(value) for value in residual_config["readout_bounds"][name])
        for name in residual.READOUT_PARAMETER_NAMES
    }
    readout_bounds["c2"] = (
        0.0,
        float(config["readout_profile"]["c2_upper_bound"]),
    )

    return {
        "continuum_config_path": continuum_config_path,
        "residual_config_path": residual_config_path,
        "manifest_path": manifest_path,
        "snapshot_path": snapshot_path,
        "continuum_config": continuum_config,
        "residual_config": residual_config,
        "snapshot": snapshot,
        "repeat_case": repeat_case,
        "comparison_source": comparison_source,
        "experiment_path": experiment_path,
        "accepted_records": int(target_context["accepted_records"]),
        "baseline_detector": baseline_detector,
        "detector_names": detector_names,
        "detector_bounds": detector_bounds,
        "readout_bounds": readout_bounds,
        "frequency": frequency,
        "target": target,
        "raw_target": raw_target,
        "full_args": full_args,
        "fit_args": fit_args,
        "scale_hz": float(
            residual_config["readout_parameterization"]["reference_scale_Hz"]
        ),
        "tracked_white": float(snapshot["white_asd_A_rtHz"]),
        "optimizer_cfg": residual_config["optimizer"],
        "repair": repair,
    }


def _fit_variant(
    *,
    problem: dict,
    name: str,
    readout_parameters,
    reference_readout: dict,
    local_readout: dict | None,
    readout_bounds: dict,
    seed_offset: int,
    warm_solutions=(),
):
    row = continuum.fit_with_profiled_white(
        family={
            "name": name,
            "readout_parameters": list(readout_parameters),
        },
        baseline_detector=problem["baseline_detector"],
        reference_readout=reference_readout,
        local_readout=(
            reference_readout if local_readout is None else local_readout
        ),
        detector_names=problem["detector_names"],
        detector_bounds=problem["detector_bounds"],
        readout_bounds=readout_bounds,
        frequency=problem["frequency"],
        target=problem["target"],
        fit_args=problem["fit_args"],
        scale_hz=problem["scale_hz"],
        baseline_white_asd=problem["tracked_white"],
        optimizer_cfg=problem["optimizer_cfg"],
        profile_cfg=problem["continuum_config"]["white_profile"],
        seed=int(problem["optimizer_cfg"]["seed"]) + int(seed_offset),
        warm_solutions=warm_solutions,
    )
    white = float(row["profiled_white_asd_A_rtHz"])
    model, point = residual.model_for_candidate(
        row["_detector_full"],
        row["_readout_full"],
        problem["frequency"],
        problem["scale_hz"],
        white,
    )
    if model is None:
        raise RuntimeError(f"model reconstruction failed for {name}: {point}")
    row["_model_full"] = model
    row["continuum_metrics_full"] = holdout.model_metrics(
        model,
        problem["target"],
        problem["frequency"],
        problem["fit_args"],
    )
    row["raw_experiment_metrics_full"] = holdout.model_metrics(
        model,
        problem["raw_target"],
        problem["frequency"],
        problem["fit_args"],
    )
    row["continuum_metrics_1_40k"] = continuum.metrics_in_region(
        model,
        problem["target"],
        problem["frequency"],
        problem["full_args"],
        1000.0,
        40000.0,
    )
    row["continuum_metrics_40_200k"] = continuum.metrics_in_region(
        model,
        problem["target"],
        problem["frequency"],
        problem["full_args"],
        40000.0,
        200000.0,
    )
    return row


def _rms(row: dict, key="continuum_metrics_full") -> float:
    return float(row[key]["residual_metrics"]["rms_residual_dB"])


def _comparison(row: dict, reference: dict) -> dict:
    ref_score = float(reference["shape_score"])
    return {
        "shape_score_ratio_to_free_pole": float(row["shape_score"] / ref_score),
        "continuum_rms_delta_to_free_pole_dB": float(
            _rms(row) - _rms(reference)
        ),
        "continuum_1_40k_rms_delta_to_free_pole_dB": float(
            _rms(row, "continuum_metrics_1_40k")
            - _rms(reference, "continuum_metrics_1_40k")
        ),
        "continuum_40_200k_rms_delta_to_free_pole_dB": float(
            _rms(row, "continuum_metrics_40_200k")
            - _rms(reference, "continuum_metrics_40_200k")
        ),
    }


def _clean(row: dict) -> dict:
    return {
        key: value
        for key, value in row.items()
        if not key.startswith("_")
    }


def _direct_pass(row: dict, reference: dict, screen: dict) -> bool:
    comp = _comparison(row, reference)
    return bool(
        comp["continuum_rms_delta_to_free_pole_dB"]
        <= float(screen["known_filter_max_rms_degradation_dB"])
        and comp["shape_score_ratio_to_free_pole"]
        <= float(screen["known_filter_max_score_ratio_to_free_pole"])
    )


def run(config: dict, config_path: Path, experiment_path_override=None):
    problem = _problem(config, config_path, experiment_path_override)
    snapshot = problem["snapshot"]
    bounds = dict(problem["readout_bounds"])

    free = _fit_variant(
        problem=problem,
        name="free_pole_plus_c2",
        readout_parameters=("pole_Hz", "pole_Q", "c2"),
        reference_readout=dict(snapshot["reference_readout"]),
        local_readout=dict(snapshot["local_readout"]),
        readout_bounds=bounds,
        seed_offset=6100,
    )

    nominal = float(config["magnicon_filter"]["nominal_cutoff_Hz"])
    tolerance = float(config["magnicon_filter"]["cutoff_tolerance_fraction"])
    screen = config["comparison_screen"]
    legacy_c4 = float(snapshot["reference_readout"]["c4"])

    variants = {}
    canonical = {}
    for index, norm in enumerate(config["magnicon_filter"]["normalization_candidates"]):
        coordinates = second_order_bessel_canonical(nominal, norm)
        canonical[norm] = coordinates
        reference = {
            "pole_Hz": float(coordinates["pole_Hz"]),
            "pole_Q": float(coordinates["pole_Q"]),
            "c2": float(snapshot["reference_readout"]["c2"]),
            "c4": legacy_c4,
        }

        fixed_name = f"magnicon_{norm}_fixed_plus_c2_legacy_c4"
        fixed = _fit_variant(
            problem=problem,
            name=fixed_name,
            readout_parameters=("c2",),
            reference_readout=reference,
            local_readout=reference,
            readout_bounds=bounds,
            seed_offset=6200 + 200 * index,
            warm_solutions=[
                ("free_pole_solution", free["_detector_full"], free["_readout_full"])
            ],
        )
        fixed["magnicon_interpretation"] = {
            "filter_norm": norm,
            "nominal_cutoff_Hz": nominal,
            "equivalent_fixed_pole_Hz": float(coordinates["pole_Hz"]),
            "equivalent_fixed_pole_Q": float(coordinates["pole_Q"]),
            "comparison_to_free_pole": _comparison(fixed, free),
        }
        variants[fixed_name] = fixed

        if bool(config["readout_profile"].get("also_test_c4_zero", True)):
            zero_reference = dict(reference)
            zero_reference["c4"] = 0.0
            zero_name = f"magnicon_{norm}_fixed_plus_c2_c4_zero"
            zero = _fit_variant(
                problem=problem,
                name=zero_name,
                readout_parameters=("c2",),
                reference_readout=zero_reference,
                local_readout=zero_reference,
                readout_bounds=bounds,
                seed_offset=6250 + 200 * index,
                warm_solutions=[
                    ("magnicon_fixed_legacy_c4", fixed["_detector_full"], fixed["_readout_full"])
                ],
            )
            zero["magnicon_interpretation"] = {
                "filter_norm": norm,
                "nominal_cutoff_Hz": nominal,
                "c4_forced_to_zero": True,
                "comparison_to_free_pole": _comparison(zero, free),
                "rms_delta_to_same_filter_legacy_c4_dB": float(
                    _rms(zero) - _rms(fixed)
                ),
            }
            variants[zero_name] = zero

        if bool(config["readout_profile"].get("also_test_filter_only", True)):
            only_reference = dict(reference)
            only_reference["c2"] = 0.0
            only_reference["c4"] = 0.0
            only_name = f"magnicon_{norm}_fixed_filter_only"
            only = _fit_variant(
                problem=problem,
                name=only_name,
                readout_parameters=(),
                reference_readout=only_reference,
                local_readout=only_reference,
                readout_bounds=bounds,
                seed_offset=6300 + 200 * index,
                warm_solutions=[
                    ("magnicon_fixed_legacy_c4", fixed["_detector_full"], fixed["_readout_full"])
                ],
            )
            only["magnicon_interpretation"] = {
                "filter_norm": norm,
                "nominal_cutoff_Hz": nominal,
                "numerator_forced_to_unity": True,
                "comparison_to_free_pole": _comparison(only, free),
            }
            variants[only_name] = only

        if bool(
            config["readout_profile"].get(
                "also_profile_cutoff_within_manual_tolerance", True
            )
        ):
            low_cutoff = nominal * (1.0 - tolerance)
            high_cutoff = nominal * (1.0 + tolerance)
            low_coordinates = second_order_bessel_canonical(low_cutoff, norm)
            high_coordinates = second_order_bessel_canonical(high_cutoff, norm)
            tolerance_bounds = dict(bounds)
            tolerance_bounds["pole_Hz"] = (
                float(low_coordinates["pole_Hz"]),
                float(high_coordinates["pole_Hz"]),
            )
            tol_name = f"magnicon_{norm}_tolerance_plus_c2_legacy_c4"
            tol_row = _fit_variant(
                problem=problem,
                name=tol_name,
                readout_parameters=("pole_Hz", "c2"),
                reference_readout=reference,
                local_readout=reference,
                readout_bounds=tolerance_bounds,
                seed_offset=6350 + 200 * index,
                warm_solutions=[
                    ("magnicon_nominal", fixed["_detector_full"], fixed["_readout_full"])
                ],
            )
            fitted_cutoff = equivalent_cutoff_from_pole(
                tol_row["readout"]["pole_Hz"],
                coordinates,
            )
            tol_row["magnicon_interpretation"] = {
                "filter_norm": norm,
                "manual_cutoff_range_Hz": [low_cutoff, high_cutoff],
                "fitted_cutoff_Hz": float(fitted_cutoff),
                "fitted_cutoff_fraction_from_nominal": float(
                    fitted_cutoff / nominal - 1.0
                ),
                "comparison_to_free_pole": _comparison(tol_row, free),
                "rms_improvement_over_nominal_same_norm_dB": float(
                    _rms(fixed) - _rms(tol_row)
                ),
            }
            variants[tol_name] = tol_row

    primary = str(config["magnicon_filter"]["primary_normalization"])
    primary_fixed = variants[f"magnicon_{primary}_fixed_plus_c2_legacy_c4"]
    primary_zero = variants.get(f"magnicon_{primary}_fixed_plus_c2_c4_zero")
    primary_only = variants.get(f"magnicon_{primary}_fixed_filter_only")
    primary_tol = variants.get(
        f"magnicon_{primary}_tolerance_plus_c2_legacy_c4"
    )

    direct_pass_by_norm = {
        norm: _direct_pass(
            variants[f"magnicon_{norm}_fixed_plus_c2_legacy_c4"],
            free,
            screen,
        )
        for norm in canonical
    }
    if direct_pass_by_norm.get(primary, False):
        classification = "documented_lpf_can_replace_free_pole_numerically"
    elif any(direct_pass_by_norm.values()):
        classification = "replacement_is_bessel_normalization_sensitive"
    else:
        classification = "documented_lpf_does_not_replace_free_pole_at_current_screen"

    c4_zero_close = None
    if primary_zero is not None:
        c4_zero_close = bool(
            _rms(primary_zero) - _rms(primary_fixed)
            <= float(screen["c4_zero_max_rms_degradation_dB"])
        )
    filter_only_close = None
    if primary_only is not None:
        filter_only_close = bool(
            _rms(primary_only) - _rms(free)
            <= float(screen["filter_only_max_rms_degradation_dB"])
        )
    tolerance_material = None
    if primary_tol is not None:
        tolerance_material = bool(
            _rms(primary_fixed) - _rms(primary_tol)
            >= float(screen["cutoff_tolerance_material_rms_improvement_dB"])
        )

    rows = {"free_pole_plus_c2": free, **variants}
    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "manual_reference": config["manual_reference"],
        "measurement_time_filter_state_known": False,
        "tested_hypothesis": "Magnicon XXF-1 Connector Box 10 kHz LPF was ON",
        "repeat_case": {
            "label": problem["repeat_case"]["label"],
            "comparison_source": problem["comparison_source"],
            "experiment_path": str(problem["experiment_path"]),
            "accepted_records": problem["accepted_records"],
        },
        "line_repair": {
            "accepted_narrow_lines": problem["repair"]["accepted_lines"],
            "rejected_broad_features": problem["repair"]["rejected_broad_features"],
            "native_fft_bin_Hz": float(problem["repair"]["native_bin_Hz"]),
        },
        "magnicon_filter": {
            "manual_nominal_cutoff_Hz": nominal,
            "manual_tolerance_fraction": tolerance,
            "normalization_guardrail": config["magnicon_filter"][
                "normalization_guardrail"
            ],
            "canonical_coordinates": canonical,
        },
        "free_pole_reference": {
            "readout": {
                key: float(value)
                for key, value in free["readout"].items()
            },
            "white_asd_A_rtHz": float(free["profiled_white_asd_A_rtHz"]),
            "continuum_rms_dB": _rms(free),
            "continuum_1_40k_rms_dB": _rms(free, "continuum_metrics_1_40k"),
            "continuum_40_200k_rms_dB": _rms(
                free, "continuum_metrics_40_200k"
            ),
            "shape_score": float(free["shape_score"]),
        },
        "fits": {name: _clean(row) for name, row in rows.items()},
        "comparisons": {
            name: _comparison(row, free)
            for name, row in variants.items()
        },
        "interpretation": {
            "classification": classification,
            "direct_replacement_pass_by_normalization": direct_pass_by_norm,
            "primary_normalization": primary,
            "primary_c4_zero_close_to_legacy_c4": c4_zero_close,
            "primary_filter_only_close_to_free_pole": filter_only_close,
            "primary_cutoff_tolerance_material_improvement": tolerance_material,
            "screen": screen,
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_continuum_config": str(problem["continuum_config_path"]),
            "base_residual_config": str(problem["residual_config_path"]),
            "manifest": str(problem["manifest_path"]),
            "residual_baseline_snapshot": str(problem["snapshot_path"]),
        },
        "_plot": {
            "frequency_Hz": problem["frequency"].tolist(),
            "continuum_target": problem["target"].tolist(),
            "raw_target": problem["raw_target"].tolist(),
            "models": {
                name: row["_model_full"].tolist()
                for name, row in rows.items()
                if (
                    name == "free_pole_plus_c2"
                    or name
                    in {
                        f"magnicon_{primary}_fixed_plus_c2_legacy_c4",
                        f"magnicon_{primary}_fixed_plus_c2_c4_zero",
                        f"magnicon_{primary}_fixed_filter_only",
                        f"magnicon_{primary}_tolerance_plus_c2_legacy_c4",
                    }
                )
            },
        },
    }
    return result


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    frequency = np.asarray(p["frequency_Hz"], dtype=float)
    target = np.asarray(p["continuum_target"], dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    top, bottom = axes
    top.loglog(frequency, target, label="line-repaired continuum", linewidth=2)
    for name, model in p["models"].items():
        top.loglog(frequency, model, label=name)
    top.set_ylabel("Normalized ASD")
    top.legend(fontsize=8, frameon=False)
    top.grid(True, which="both", alpha=0.25)

    for name, model in p["models"].items():
        model = np.asarray(model, dtype=float)
        residual_db = 20.0 * np.log10(model / target)
        bottom.semilogx(frequency, residual_db, label=name)
    bottom.axhline(0.0, linewidth=1)
    bottom.set_xlabel("Frequency [Hz]")
    bottom.set_ylabel("Model / continuum [dB]")
    bottom.grid(True, which="both", alpha=0.25)

    fig.suptitle(
        "Magnicon XXF-1 10 kHz LPF validation | "
        + result["interpretation"]["classification"],
        fontsize=10,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def json_safe(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return [json_safe(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    return value


def cleaned_result(result: dict) -> dict:
    return json_safe(
        {
            key: value
            for key, value in result.items()
            if key != "_plot"
        }
    )


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
                "classification": result["interpretation"]["classification"],
                "direct_replacement_pass_by_normalization": result[
                    "interpretation"
                ]["direct_replacement_pass_by_normalization"],
                "free_pole_reference": result["free_pole_reference"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
