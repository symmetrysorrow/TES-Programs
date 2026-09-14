"""Audit stored modelnoise.txt provenance against a fresh canonical raw-record target.

This diagnostic keeps the detector model frozen. It reconstructs the canonical
post-analysis ASD from the exact accepted CH0 record mask, compares that fresh
shape with the stored optimizer target, separates broad and narrow target-file
differences, and evaluates one frozen detector spectrum against both targets.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import ndimage, signal

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    estimate_one_sided_asd,
)

REFERENCE_HZ = 1_000.0
SMOOTH_WIDTH_HZ_DEFAULT = 1_000.0
NARROW_EXCURSION_THRESHOLD_DB = 0.5
FIT_BANDS_HZ = (
    (1_000.0, 5_000.0),
    (5_000.0, 15_000.0),
    (15_000.0, 40_000.0),
    (40_000.0, 100_000.0),
    (100_000.0, 200_000.0),
)


def fit_args(summary):
    fit = summary["fit"]
    return SimpleNamespace(
        fit_min_hz=float(fit["min_hz"]),
        fit_max_hz=float(fit["max_hz"]),
        fit_weight_start_hz=float(fit["weight_start_hz"]),
        high_frequency_weight=float(fit["high_frequency_weight"]),
        robust_delta_dex=float(fit["robust_delta_dex"]),
        fit_points=int(fit["points"]),
        absolute_asd_weight=float(fit.get("absolute_asd_weight", 0.0)),
    )


def read_record(path: Path) -> np.ndarray:
    return np.frombuffer(path.read_bytes()[4:], dtype=np.float64).copy()


def normalize_at(frequency, values, reference_hz=REFERENCE_HZ):
    frequency = np.asarray(frequency, dtype=float)
    values = np.asarray(values, dtype=float)
    index = int(np.argmin(np.abs(frequency - float(reference_hz))))
    reference = float(values[index])
    if not np.isfinite(reference) or reference <= 0.0:
        raise ValueError("invalid normalization reference")
    return values / reference


def reconstruct_fresh_post_asd(
    comparison_summary: dict,
    experiment_path: Path,
):
    acquisition = comparison_summary["acquisition"]
    rate = float(acquisition["rate_Hz"])
    sample = int(acquisition["samples"])
    cutoff = float(acquisition["cutoff_Hz"])
    accepted_indices = [
        int(value) for value in acquisition["accepted_record_indices"]
    ]
    raw_dir = experiment_path / "CH0_noise" / "rawdata"
    paths = sorted(raw_dir.glob("CH0_*.dat"))
    if not paths:
        raise FileNotFoundError(f"no CH0 raw records found in {raw_dir}")
    if not accepted_indices:
        raise ValueError("comparison summary contains no accepted records")
    if max(accepted_indices) >= len(paths):
        raise IndexError("accepted record index exceeds raw record list")

    records = [read_record(paths[index]) for index in accepted_indices]
    fresh, count = estimate_one_sided_asd(
        records,
        sample,
        rate,
        cutoff=cutoff,
        remove_mean=True,
    )
    if count != len(accepted_indices):
        raise RuntimeError("fresh reconstruction did not retain accepted mask")
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    return frequency, fresh, count


def band_metrics(residual_db, frequency):
    residual_db = np.asarray(residual_db, dtype=float)
    frequency = np.asarray(frequency, dtype=float)
    out = {}
    for low, high in FIT_BANDS_HZ:
        mask = (frequency >= low) & (frequency <= high)
        values = residual_db[mask]
        out[f"{low:g}-{high:g}_Hz"] = {
            "mean_dB": float(np.mean(values)),
            "rms_dB": float(np.sqrt(np.mean(values**2))),
            "max_abs_dB": float(np.max(np.abs(values))),
        }
    return out


def shape_difference_analysis(
    stored,
    fresh,
    frequency,
    smooth_width_hz=SMOOTH_WIDTH_HZ_DEFAULT,
):
    stored = np.asarray(stored, dtype=float)
    fresh = np.asarray(fresh, dtype=float)
    frequency = np.asarray(frequency, dtype=float)
    if stored.shape != fresh.shape or stored.shape != frequency.shape:
        raise ValueError("stored, fresh, and frequency must have equal shapes")

    stored_norm = normalize_at(frequency, stored)
    fresh_norm = normalize_at(frequency, fresh)
    valid = (
        (frequency >= 1_000.0)
        & (frequency <= 200_000.0)
        & np.isfinite(stored_norm)
        & np.isfinite(fresh_norm)
        & (stored_norm > 0.0)
        & (fresh_norm > 0.0)
    )
    freq = frequency[valid]
    residual = 20.0 * np.log10(stored_norm[valid] / fresh_norm[valid])
    if len(freq) < 3:
        raise ValueError("not enough valid target bins")

    spacing = float(np.median(np.diff(freq)))
    kernel = max(3, int(round(float(smooth_width_hz) / spacing)))
    if kernel % 2 == 0:
        kernel += 1
    broad = ndimage.median_filter(residual, size=kernel, mode="nearest")
    narrow = residual - broad

    min_peak_distance = max(1, int(round(250.0 / spacing)))
    peaks, properties = signal.find_peaks(
        np.abs(narrow),
        height=NARROW_EXCURSION_THRESHOLD_DB,
        distance=min_peak_distance,
    )
    order = sorted(
        peaks,
        key=lambda index: abs(float(narrow[index])),
        reverse=True,
    )[:12]
    excursions = [
        {
            "frequency_Hz": float(freq[index]),
            "total_difference_dB": float(residual[index]),
            "broad_component_dB": float(broad[index]),
            "narrow_component_dB": float(narrow[index]),
        }
        for index in order
    ]

    return {
        "definition": "20*log10(stored_normalized/fresh_normalized)",
        "reference_Hz": REFERENCE_HZ,
        "smooth_median_width_Hz": float(smooth_width_hz),
        "smooth_kernel_bins": int(kernel),
        "total": {
            "rms_dB": float(np.sqrt(np.mean(residual**2))),
            "mean_dB": float(np.mean(residual)),
            "max_abs_dB": float(np.max(np.abs(residual))),
            "bands": band_metrics(residual, freq),
        },
        "broad_component": {
            "rms_dB": float(np.sqrt(np.mean(broad**2))),
            "mean_dB": float(np.mean(broad)),
            "max_abs_dB": float(np.max(np.abs(broad))),
            "bands": band_metrics(broad, freq),
        },
        "narrow_component": {
            "rms_dB": float(np.sqrt(np.mean(narrow**2))),
            "mean_dB": float(np.mean(narrow)),
            "max_abs_dB": float(np.max(np.abs(narrow))),
            "fraction_of_bins_abs_gt_0p5dB": float(
                np.mean(np.abs(narrow) > NARROW_EXCURSION_THRESHOLD_DB)
            ),
            "excursion_count_over_0p5dB": int(
                np.count_nonzero(
                    np.abs(narrow) > NARROW_EXCURSION_THRESHOLD_DB
                )
            ),
            "top_separated_excursions": excursions,
        },
    }


def optimizer_stored_target(stored, rate, fit_frequency):
    stored = np.asarray(stored, dtype=float)
    stored_frequency = (
        np.arange(len(stored), dtype=float)
        * (float(rate) / 2.0)
        / len(stored)
    )
    return np.interp(
        fit_frequency,
        stored_frequency,
        opt.normalize_at(stored_frequency, stored),
    )


def fresh_target(fresh, fresh_frequency, fit_frequency):
    return np.interp(
        fit_frequency,
        fresh_frequency,
        opt.normalize_at(fresh_frequency, fresh),
    )


def model_target_metrics(model, target, fit_frequency, args):
    bands = opt.band_fit_diagnostics(
        model, target, fit_frequency, args
    )
    converted = {}
    for key, row in bands.items():
        mean = float(row["mean_log10_ratio"])
        converted[key] = {
            "mean_log10_model_over_target": mean,
            "mean_model_over_target_dB": float(20.0 * mean),
            "geometric_mean_model_over_target": float(10.0**mean),
            "rms_log10_ratio": float(row["rms_log10_ratio"]),
            "max_abs_log10_ratio": float(row["max_abs_log10_ratio"]),
        }
    return {
        "shape_score": float(
            opt.fit_score(model, target, fit_frequency, args)
        ),
        "bands": converted,
    }


def anchor_rows(
    frequency,
    stored,
    fresh,
):
    stored_norm = normalize_at(frequency, stored)
    fresh_norm = normalize_at(frequency, fresh)
    anchors = (
        1_000.0,
        5_000.0,
        10_000.0,
        20_000.0,
        40_000.0,
        70_000.0,
        100_000.0,
        150_000.0,
        200_000.0,
    )
    rows = []
    for anchor in anchors:
        index = int(np.argmin(np.abs(frequency - anchor)))
        ratio = float(stored_norm[index] / fresh_norm[index])
        rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "stored_normalized": float(stored_norm[index]),
                "fresh_normalized": float(fresh_norm[index]),
                "stored_over_fresh": ratio,
                "stored_over_fresh_dB": float(20.0 * np.log10(ratio)),
            }
        )
    return rows


def run(
    comparison_summary_path: Path,
    optimizer_summary_path: Path,
    experiment_path: Path | None = None,
    smooth_width_hz: float = SMOOTH_WIDTH_HZ_DEFAULT,
):
    comparison = json.loads(
        comparison_summary_path.read_text(encoding="utf-8")
    )
    optimizer_summary = json.loads(
        optimizer_summary_path.read_text(encoding="utf-8")
    )
    if experiment_path is None:
        experiment_path = Path(comparison["experiment_path"])

    frequency, fresh, accepted_count = reconstruct_fresh_post_asd(
        comparison, experiment_path
    )
    modelnoise_path = experiment_path / "CH0_noise" / "modelnoise.txt"
    if not modelnoise_path.exists():
        raise FileNotFoundError(f"stored target not found: {modelnoise_path}")
    stored = np.asarray(np.loadtxt(modelnoise_path), dtype=float)
    if stored.shape != fresh.shape:
        raise ValueError(
            "stored modelnoise and fresh canonical ASD have different lengths"
        )

    target_difference = shape_difference_analysis(
        stored,
        fresh,
        frequency,
        smooth_width_hz=smooth_width_hz,
    )

    args = fit_args(optimizer_summary)
    fit_frequency = np.geomspace(
        float(args.fit_min_hz),
        float(args.fit_max_hz),
        int(args.fit_points),
    )
    stored_target = optimizer_stored_target(
        stored,
        comparison["acquisition"]["rate_Hz"],
        fit_frequency,
    )
    fresh_fit_target = fresh_target(
        fresh,
        frequency,
        fit_frequency,
    )
    candidate = dict(optimizer_summary["best_case_parameters"])
    model, _ = opt.deterministic_simulated_spectrum(
        candidate.copy(), fit_frequency
    )
    stored_eval = model_target_metrics(
        model, stored_target, fit_frequency, args
    )
    fresh_eval = model_target_metrics(
        model, fresh_fit_target, fit_frequency, args
    )

    stored_mid = stored_eval["bands"]["5000-15000_Hz"][
        "mean_log10_model_over_target"
    ]
    stored_high = stored_eval["bands"]["40000-100000_Hz"][
        "mean_log10_model_over_target"
    ]
    fresh_mid = fresh_eval["bands"]["5000-15000_Hz"][
        "mean_log10_model_over_target"
    ]
    fresh_high = fresh_eval["bands"]["40000-100000_Hz"][
        "mean_log10_model_over_target"
    ]

    summary_score = optimizer_summary.get("best_case_score")
    score_difference = (
        float(stored_eval["shape_score"] - float(summary_score))
        if summary_score is not None
        else None
    )
    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "detector_parameters_held_fixed": True,
        "stored_modelnoise_unchanged": True,
        "comparison_summary": str(comparison_summary_path),
        "optimizer_summary": str(optimizer_summary_path),
        "experiment_path": str(experiment_path),
        "stored_modelnoise_path": str(modelnoise_path),
        "accepted_records": int(accepted_count),
        "fresh_target_semantics": (
            "same accepted CH0 record mask; mean removal; production 10 kHz "
            "second-order Bessel filtfilt; Hann; power average; one-sided ASD"
        ),
        "target_shape_difference": target_difference,
        "anchor_comparison": anchor_rows(
            frequency, stored, fresh
        ),
        "frozen_detector_evaluation": {
            "stored_target": stored_eval,
            "fresh_canonical_target": fresh_eval,
            "optimizer_summary_best_case_score": (
                float(summary_score) if summary_score is not None else None
            ),
            "stored_score_minus_summary_score": score_difference,
            "stored_score_reproduces_summary": (
                bool(abs(score_difference) <= 1e-10)
                if score_difference is not None
                else None
            ),
            "score_change_fresh_minus_stored": float(
                fresh_eval["shape_score"] - stored_eval["shape_score"]
            ),
            "score_ratio_fresh_over_stored": float(
                fresh_eval["shape_score"] / stored_eval["shape_score"]
            ),
            "stored_has_5_15_deficit_and_40_100_excess": bool(
                stored_mid < 0.0 and stored_high > 0.0
            ),
            "fresh_has_5_15_deficit_and_40_100_excess": bool(
                fresh_mid < 0.0 and fresh_high > 0.0
            ),
            "mid_mean_shift_dex_fresh_minus_stored": float(
                fresh_mid - stored_mid
            ),
            "high_mean_shift_dex_fresh_minus_stored": float(
                fresh_high - stored_high
            ),
        },
        "guardrail": (
            "A stored/fresh mismatch is target-file provenance evidence, not "
            "a detector-physics result. The stored target is never overwritten."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--comparison-summary",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--optimizer-summary",
        type=Path,
        required=True,
    )
    parser.add_argument("--experiment-path", type=Path)
    parser.add_argument(
        "--smooth-width-hz",
        type=float,
        default=SMOOTH_WIDTH_HZ_DEFAULT,
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = run(
        args.comparison_summary,
        args.optimizer_summary,
        experiment_path=args.experiment_path,
        smooth_width_hz=args.smooth_width_hz,
    )
    output = args.output or args.optimizer_summary.with_name(
        "modelnoise_provenance_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    evaluation = result["frozen_detector_evaluation"]
    print(
        json.dumps(
            {
                "output": str(output),
                "accepted_records": result["accepted_records"],
                "target_total_rms_dB": result[
                    "target_shape_difference"
                ]["total"]["rms_dB"],
                "target_broad_rms_dB": result[
                    "target_shape_difference"
                ]["broad_component"]["rms_dB"],
                "target_narrow_rms_dB": result[
                    "target_shape_difference"
                ]["narrow_component"]["rms_dB"],
                "stored_shape_score": evaluation[
                    "stored_target"
                ]["shape_score"],
                "fresh_shape_score": evaluation[
                    "fresh_canonical_target"
                ]["shape_score"],
                "stored_pattern": evaluation[
                    "stored_has_5_15_deficit_and_40_100_excess"
                ],
                "fresh_pattern": evaluation[
                    "fresh_has_5_15_deficit_and_40_100_excess"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
