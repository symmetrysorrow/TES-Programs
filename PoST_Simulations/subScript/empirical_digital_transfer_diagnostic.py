"""Measure the experimental digital-analysis transfer directly from raw records.

This diagnostic is intentionally TES-model free. It reuses the exact accepted
CH0 record mask from comparison_summary.json, reconstructs pre/post-analysis
ASDs from the same raw records, and compares their ratio with analytic and
finite-record Bessel transfer predictions.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    estimate_one_sided_asd,
    one_sided_asd_from_power,
    preprocess_noise_record,
    windowed_rfft_power,
)

FIT_BANDS_HZ = (
    (1_000.0, 5_000.0),
    (5_000.0, 15_000.0),
    (15_000.0, 40_000.0),
    (40_000.0, 100_000.0),
    (100_000.0, 200_000.0),
)
REFERENCE_HZ = 1_000.0
PROXY_DIGITAL_ORDER = 4
PROXY_DIGITAL_PASSES = 1
DEFAULT_FINITE_SEED = 20260914


def read_record(path: Path) -> np.ndarray:
    return np.frombuffer(path.read_bytes()[4:], dtype=np.float64).copy()


def bessel_magnitude(frequency_hz, rate_hz, cutoff_hz, order, passes):
    frequency = np.asarray(frequency_hz, dtype=float)
    rate = float(rate_hz)
    cutoff = float(cutoff_hz)
    order = int(order)
    passes = int(passes)
    if rate <= 0.0 or cutoff <= 0.0 or cutoff >= rate / 2.0:
        raise ValueError("invalid rate/cutoff")
    if order < 1 or passes < 1:
        raise ValueError("order and passes must be positive")
    b, a = signal.bessel(order, cutoff / (rate / 2.0), "low")
    _, response = signal.freqz(
        b,
        a,
        worN=2.0 * np.pi * frequency / rate,
    )
    return np.abs(response) ** passes


def normalize_at(frequency, values, reference_hz=REFERENCE_HZ):
    frequency = np.asarray(frequency, dtype=float)
    values = np.asarray(values, dtype=float)
    index = int(np.argmin(np.abs(frequency - float(reference_hz))))
    reference = float(values[index])
    if not np.isfinite(reference) or reference <= 0.0:
        raise ValueError("invalid normalization reference")
    return values / reference


def paired_finite_record_transfer(
    input_asd,
    sample,
    rate_hz,
    cutoff_hz,
    records,
    seed=DEFAULT_FINITE_SEED,
):
    """Simulate paired pre/post estimates from identical synthetic records."""

    sample = int(sample)
    rate = float(rate_hz)
    records = int(records)
    input_asd = np.asarray(input_asd, dtype=float)
    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    if input_asd.shape != frequency.shape:
        raise ValueError("input_asd must use the exact rFFT grid")
    if np.any(~np.isfinite(input_asd)) or np.any(input_asd < 0.0):
        raise ValueError("input_asd must be finite and non-negative")
    if records <= 0:
        raise ValueError("records must be positive")

    df = rate / sample
    window = np.hanning(sample)
    window_gain = float(np.sqrt(np.mean(window**2)))
    pre_power = np.zeros_like(frequency)
    post_power = np.zeros_like(frequency)
    rng = np.random.default_rng(int(seed))

    for _ in range(records):
        spectrum = np.zeros_like(frequency, dtype=np.complex128)
        if len(frequency) > 2:
            sigma = input_asd[1:-1] * sample * np.sqrt(df) / 2.0
            spectrum[1:-1] = (
                rng.normal(size=len(sigma)) * sigma
                + 1j * rng.normal(size=len(sigma)) * sigma
            )
        spectrum[0] = (
            rng.normal() * input_asd[0] * sample * np.sqrt(df)
        )
        if sample % 2 == 0 and len(frequency) > 1:
            spectrum[-1] = (
                rng.normal() * input_asd[-1] * sample * np.sqrt(df)
            )
        elif len(frequency) > 1:
            sigma_last = input_asd[-1] * sample * np.sqrt(df) / 2.0
            spectrum[-1] = (
                rng.normal() * sigma_last
                + 1j * rng.normal() * sigma_last
            )

        record = np.fft.irfft(spectrum, n=sample)
        pre = preprocess_noise_record(
            record, rate, cutoff=0.0, remove_mean=True
        )
        post = preprocess_noise_record(
            record, rate, cutoff=cutoff_hz, remove_mean=True
        )
        pre_power += windowed_rfft_power(pre, window)
        post_power += windowed_rfft_power(post, window)

    pre_asd = one_sided_asd_from_power(
        pre_power / records,
        sample,
        rate,
        window_gain,
    )
    post_asd = one_sided_asd_from_power(
        post_power / records,
        sample,
        rate,
        window_gain,
    )
    return pre_asd, post_asd, post_asd / np.maximum(
        pre_asd, np.finfo(float).tiny
    )


def residual_metrics(empirical, predicted, frequency, reference_hz=REFERENCE_HZ):
    empirical = np.asarray(empirical, dtype=float)
    predicted = np.asarray(predicted, dtype=float)
    frequency = np.asarray(frequency, dtype=float)
    valid = (
        (frequency >= 1_000.0)
        & (frequency <= 200_000.0)
        & np.isfinite(empirical)
        & np.isfinite(predicted)
        & (empirical > 0.0)
        & (predicted > 0.0)
    )
    absolute_db = 20.0 * np.log10(
        predicted[valid] / empirical[valid]
    )
    empirical_norm = normalize_at(
        frequency, empirical, reference_hz=reference_hz
    )
    predicted_norm = normalize_at(
        frequency, predicted, reference_hz=reference_hz
    )
    normalized_db_full = 20.0 * np.log10(
        predicted_norm / empirical_norm
    )
    normalized_db = normalized_db_full[valid]

    bands = {}
    for low, high in FIT_BANDS_HZ:
        mask = valid & (frequency >= low) & (frequency <= high)
        values = normalized_db_full[mask]
        bands[f"{low:g}-{high:g}_Hz"] = {
            "mean_residual_dB": float(np.mean(values)),
            "rms_residual_dB": float(np.sqrt(np.mean(values**2))),
            "max_abs_residual_dB": float(np.max(np.abs(values))),
        }

    return {
        "absolute_rms_residual_dB": float(
            np.sqrt(np.mean(absolute_db**2))
        ),
        "absolute_max_abs_residual_dB": float(
            np.max(np.abs(absolute_db))
        ),
        "shape_normalized_rms_residual_dB": float(
            np.sqrt(np.mean(normalized_db**2))
        ),
        "shape_normalized_max_abs_residual_dB": float(
            np.max(np.abs(normalized_db))
        ),
        "shape_normalized_mean_residual_dB": float(
            np.mean(normalized_db)
        ),
        "bands": bands,
    }


def stored_shape_metrics(stored, recomputed, frequency):
    stored = np.asarray(stored, dtype=float)
    recomputed = np.asarray(recomputed, dtype=float)
    if stored.shape != recomputed.shape:
        raise ValueError("stored modelnoise and recomputed ASD length differ")
    stored_norm = normalize_at(frequency, stored)
    recomputed_norm = normalize_at(frequency, recomputed)
    valid = (
        (frequency >= 1_000.0)
        & (frequency <= 200_000.0)
        & (stored_norm > 0.0)
        & (recomputed_norm > 0.0)
    )
    residual_db = 20.0 * np.log10(
        stored_norm[valid] / recomputed_norm[valid]
    )
    return {
        "rms_shape_difference_dB": float(
            np.sqrt(np.mean(residual_db**2))
        ),
        "max_abs_shape_difference_dB": float(
            np.max(np.abs(residual_db))
        ),
        "mean_shape_difference_dB": float(np.mean(residual_db)),
    }


def curve_sample(frequency, empirical, predictions):
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
        row = {
            "frequency_Hz": float(frequency[index]),
            "empirical_post_over_pre": float(empirical[index]),
            "empirical_transfer_dB": float(
                20.0 * np.log10(empirical[index])
            ),
        }
        for name, values in predictions.items():
            row[name] = float(values[index])
            row[f"{name}_dB"] = float(
                20.0 * np.log10(values[index])
            )
        rows.append(row)
    return rows


def run(
    comparison_summary_path: Path,
    experiment_path: Path | None = None,
    finite_records: int = 0,
    finite_seed: int = DEFAULT_FINITE_SEED,
):
    report = json.loads(
        comparison_summary_path.read_text(encoding="utf-8")
    )
    acquisition = report["acquisition"]
    rate = float(acquisition["rate_Hz"])
    sample = int(acquisition["samples"])
    cutoff = float(acquisition["cutoff_Hz"])
    accepted_indices = [
        int(value) for value in acquisition["accepted_record_indices"]
    ]
    if experiment_path is None:
        experiment_path = Path(report["experiment_path"])

    raw_dir = experiment_path / "CH0_noise" / "rawdata"
    record_paths = sorted(raw_dir.glob("CH0_*.dat"))
    if not record_paths:
        raise FileNotFoundError(f"no CH0 raw records found in {raw_dir}")
    if not accepted_indices:
        raise ValueError("comparison summary has no accepted records")
    if max(accepted_indices) >= len(record_paths):
        raise IndexError("accepted record index exceeds raw record list")

    accepted_records = [
        read_record(record_paths[index]) for index in accepted_indices
    ]
    pre_asd, pre_count = estimate_one_sided_asd(
        accepted_records,
        sample,
        rate,
        cutoff=0.0,
        remove_mean=True,
    )
    post_asd, post_count = estimate_one_sided_asd(
        accepted_records,
        sample,
        rate,
        cutoff=cutoff,
        remove_mean=True,
    )
    if pre_count != post_count or pre_count != len(accepted_indices):
        raise RuntimeError("pre/post accepted-record counts differ")

    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    empirical = post_asd / np.maximum(
        pre_asd, np.finfo(float).tiny
    )
    production_analytic = bessel_magnitude(
        frequency, rate, cutoff, order=2, passes=2
    )
    proxy_4th_single = bessel_magnitude(
        frequency,
        rate,
        cutoff,
        order=PROXY_DIGITAL_ORDER,
        passes=PROXY_DIGITAL_PASSES,
    )

    finite_count = (
        int(finite_records) if int(finite_records) > 0 else pre_count
    )
    finite_pre, finite_post, finite_transfer = (
        paired_finite_record_transfer(
            pre_asd,
            sample,
            rate,
            cutoff,
            finite_count,
            seed=finite_seed,
        )
    )

    predictions = {
        "production_analytic_order2_filtfilt": production_analytic,
        "proxy_order4_single_pass": proxy_4th_single,
        "production_finite_record": finite_transfer,
    }
    comparisons = {
        name: residual_metrics(empirical, values, frequency)
        for name, values in predictions.items()
    }
    best_name = min(
        comparisons,
        key=lambda name: comparisons[name][
            "shape_normalized_rms_residual_dB"
        ],
    )

    modelnoise_path = experiment_path / "CH0_noise" / "modelnoise.txt"
    stored_metrics = None
    if modelnoise_path.exists():
        stored_modelnoise = np.asarray(
            np.loadtxt(modelnoise_path), dtype=float
        )
        stored_metrics = stored_shape_metrics(
            stored_modelnoise, post_asd, frequency
        )

    return {
        "diagnostic_only": True,
        "tes_model_used": False,
        "detector_parameters_used": False,
        "comparison_summary": str(comparison_summary_path),
        "experiment_path": str(experiment_path),
        "raw_record_semantics": (
            "same accepted CH0 record indices as comparison_summary.json; "
            "no re-selection is performed"
        ),
        "acquisition": {
            "rate_Hz": rate,
            "samples": sample,
            "cutoff_Hz": cutoff,
            "accepted_records": pre_count,
        },
        "empirical_transfer_definition": (
            "power-averaged post-analysis ASD divided by power-averaged "
            "pre-analysis ASD from the same accepted raw records"
        ),
        "production_digital_convention": {
            "order": 2,
            "passes": 2,
            "implementation": "scipy.signal.bessel followed by filtfilt",
        },
        "proxy_from_measurement_chain_screen": {
            "order": PROXY_DIGITAL_ORDER,
            "passes": PROXY_DIGITAL_PASSES,
            "interpretation": (
                "shape proxy only; not claimed to be the actual analysis path"
            ),
        },
        "finite_record_simulation": {
            "records": finite_count,
            "seed": int(finite_seed),
            "input_spectrum": (
                "empirical pre-analysis ASD; paired synthetic records are "
                "estimated before and after the production filter"
            ),
        },
        "comparisons_to_empirical_transfer": comparisons,
        "best_shape_match": best_name,
        "stored_modelnoise_vs_recomputed_post_shape": stored_metrics,
        "stored_modelnoise_path": (
            str(modelnoise_path) if modelnoise_path.exists() else None
        ),
        "curve_sample": curve_sample(
            frequency, empirical, predictions
        ),
        "guardrail": (
            "The empirical post/pre ratio tests the digital analysis stage "
            "and finite-record estimator only. It does not identify analog "
            "hardware response or TES physics."
        ),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--comparison-summary",
        type=Path,
        required=True,
    )
    p.add_argument("--experiment-path", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--finite-records", type=int, default=0)
    p.add_argument("--finite-seed", type=int, default=DEFAULT_FINITE_SEED)
    a = p.parse_args()

    result = run(
        a.comparison_summary,
        experiment_path=a.experiment_path,
        finite_records=a.finite_records,
        finite_seed=a.finite_seed,
    )
    output = a.output or a.comparison_summary.with_name(
        "empirical_digital_transfer_diagnostic.json"
    )
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(output),
                "accepted_records": result["acquisition"][
                    "accepted_records"
                ],
                "best_shape_match": result["best_shape_match"],
                "production_analytic_shape_rms_dB": result[
                    "comparisons_to_empirical_transfer"
                ]["production_analytic_order2_filtfilt"][
                    "shape_normalized_rms_residual_dB"
                ],
                "finite_record_shape_rms_dB": result[
                    "comparisons_to_empirical_transfer"
                ]["production_finite_record"][
                    "shape_normalized_rms_residual_dB"
                ],
                "proxy_shape_rms_dB": result[
                    "comparisons_to_empirical_transfer"
                ]["proxy_order4_single_pass"][
                    "shape_normalized_rms_residual_dB"
                ],
                "stored_modelnoise_shape_rms_dB": (
                    result[
                        "stored_modelnoise_vs_recomputed_post_shape"
                    ]["rms_shape_difference_dB"]
                    if result[
                        "stored_modelnoise_vs_recomputed_post_shape"
                    ]
                    is not None
                    else None
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
