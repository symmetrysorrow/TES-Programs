"""Diagnose narrow high-frequency line noise in the raw CH0/CH1 records.

This diagnostic is intentionally model-free.  It uses the native FFT grid from
the acquisition (5 Hz for 500 kS/s and 100000 samples), keeps the production
CH0 accepted-record mask, and treats CH1 only as an auxiliary simultaneous
channel using exact event-key pairing plus finite/length/no-gross-clipping
checks.

The goal is to distinguish a smooth detector/readout noise continuum from
narrow external/interference lines.  No notch is applied and no simulation
parameter is fitted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import find_peaks, peak_widths

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_residual_dof_competition_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "high_frequency_line_noise_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "high_frequency_line_noise_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402
from subScript.pulse_contamination_common import noise_paths, record_key  # noqa: E402


def read_record(path: Path, samples: int) -> np.ndarray:
    values = np.fromfile(path, dtype=np.float64, offset=4)
    if values.size != int(samples) or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid record: {path}")
    return values


def one_sided_spectra(
    values0: np.ndarray,
    rate_hz: float,
    window: np.ndarray,
    values1: np.ndarray | None = None,
):
    x0 = np.asarray(values0, dtype=float) - float(np.mean(values0))
    f0 = np.fft.rfft(x0 * window)
    scale = float(rate_hz) * float(np.sum(window**2))
    p00 = np.abs(f0) ** 2 / scale
    p00[1:-1] *= 2.0
    if values1 is None:
        return p00, None, None

    x1 = np.asarray(values1, dtype=float) - float(np.mean(values1))
    f1 = np.fft.rfft(x1 * window)
    p11 = np.abs(f1) ** 2 / scale
    p01 = f0 * np.conjugate(f1) / scale
    p11[1:-1] *= 2.0
    p01[1:-1] *= 2.0
    return p00, p11, p01


def odd_kernel_bins(width_hz: float, bin_hz: float, minimum: int = 5) -> int:
    bins = max(int(round(float(width_hz) / float(bin_hz))), int(minimum))
    if bins % 2 == 0:
        bins += 1
    return bins


def local_baseline_db(asd: np.ndarray, kernel_bins: int) -> np.ndarray:
    db = 20.0 * np.log10(np.maximum(asd, np.finfo(float).tiny))
    return median_filter(db, size=int(kernel_bins), mode="nearest")


def detect_lines(
    frequency_hz: np.ndarray,
    asd: np.ndarray,
    *,
    baseline_width_hz: float,
    min_excess_db: float,
    min_prominence_db: float,
    min_spacing_hz: float,
    max_peaks: int,
):
    bin_hz = float(np.median(np.diff(frequency_hz)))
    kernel = odd_kernel_bins(baseline_width_hz, bin_hz)
    db = 20.0 * np.log10(np.maximum(asd, np.finfo(float).tiny))
    baseline = local_baseline_db(asd, kernel)
    excess = db - baseline
    spacing_bins = max(1, int(round(float(min_spacing_hz) / bin_hz)))
    peaks, properties = find_peaks(
        excess,
        height=float(min_excess_db),
        prominence=float(min_prominence_db),
        distance=spacing_bins,
    )
    if peaks.size:
        widths = peak_widths(excess, peaks, rel_height=0.5)[0] * bin_hz
    else:
        widths = np.asarray([], dtype=float)

    order = sorted(
        range(len(peaks)),
        key=lambda i: float(excess[peaks[i]]),
        reverse=True,
    )[: int(max_peaks)]
    rows = []
    for rank, i in enumerate(order, start=1):
        index = int(peaks[i])
        rows.append(
            {
                "rank": rank,
                "band_index": index,
                "frequency_Hz": float(frequency_hz[index]),
                "asd": float(asd[index]),
                "excess_over_local_baseline_dB": float(excess[index]),
                "prominence_dB": float(properties["prominences"][i]),
                "half_prominence_width_Hz": float(widths[i]),
            }
        )
    return rows, baseline, excess


def harmonic_relations(lines: list[dict], relative_tolerance: float = 0.01):
    relations = []
    by_frequency = sorted(lines, key=lambda row: row["frequency_Hz"])
    for i, lower in enumerate(by_frequency):
        for upper in by_frequency[i + 1 :]:
            ratio = float(upper["frequency_Hz"] / lower["frequency_Hz"])
            harmonic = int(round(ratio))
            if harmonic < 2 or harmonic > 6:
                continue
            fractional_error = abs(ratio - harmonic) / harmonic
            if fractional_error <= float(relative_tolerance):
                relations.append(
                    {
                        "fundamental_Hz": float(lower["frequency_Hz"]),
                        "candidate_harmonic_Hz": float(upper["frequency_Hz"]),
                        "integer_multiple": harmonic,
                        "frequency_ratio": ratio,
                        "fractional_error": float(fractional_error),
                    }
                )
    return relations


def auxiliary_ch1_ok(values: np.ndarray, max_range: float) -> bool:
    if not np.all(np.isfinite(values)):
        return False
    return float(np.max(values) - np.min(values)) < float(max_range)


def line_baseline_psd(
    frequency_hz: np.ndarray,
    psd: np.ndarray,
    center_hz: float,
    *,
    baseline_half_width_hz: float = 1000.0,
    exclusion_half_width_hz: float = 100.0,
):
    distance = np.abs(frequency_hz - float(center_hz))
    mask = (
        (distance <= float(baseline_half_width_hz))
        & (distance >= float(exclusion_half_width_hz))
    )
    if not np.any(mask):
        return float("nan")
    return float(np.median(psd[mask]))


def run(args):
    config = json.loads(args.config.read_text(encoding="utf-8"))
    manifest_path = ident.resolve_config_path(
        config["manifest"],
        args.config,
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = shared.normalize_manifest(manifest, manifest_path)
    repeat_case = ident.case_by_label(cases, config["repeat_case_label"])
    comparison, experiment_path, comparison_source = ident.comparison_for_case(
        repeat_case
    )
    if args.experiment_path is not None:
        experiment_path = args.experiment_path

    acquisition = comparison["acquisition"]
    rate_hz = float(acquisition["rate_Hz"])
    samples = int(acquisition["samples"])
    accepted_indices = [
        int(value) for value in acquisition["accepted_record_indices"]
    ]
    if not accepted_indices:
        raise RuntimeError("no accepted CH0 records")

    paths0 = noise_paths(experiment_path, "CH0")
    paths1 = noise_paths(experiment_path, "CH1")
    if max(accepted_indices) >= len(paths0):
        raise IndexError("CH0 accepted index exceeds available raw records")

    accepted_paths0 = [paths0[index] for index in accepted_indices]
    map1 = {record_key(path): path for path in paths1}

    frequency = np.fft.rfftfreq(samples, d=1.0 / rate_hz)
    bin_hz = float(rate_hz / samples)
    band_mask = (
        (frequency >= float(args.min_hz))
        & (frequency <= float(args.max_hz))
    )
    band_frequency = frequency[band_mask]
    if band_frequency.size < 10:
        raise RuntimeError("requested frequency band contains too few FFT bins")

    window = np.hanning(samples)
    ch0_sum = np.zeros(band_frequency.size, dtype=float)
    pair_p00 = np.zeros(band_frequency.size, dtype=float)
    pair_p11 = np.zeros(band_frequency.size, dtype=float)
    pair_p01 = np.zeros(band_frequency.size, dtype=np.complex128)
    paired_keys = []

    block_sums = []
    block_counts = []
    block_sum = np.zeros(band_frequency.size, dtype=float)
    block_count = 0

    for path0 in accepted_paths0:
        values0 = read_record(path0, samples)
        p00, _, _ = one_sided_spectra(values0, rate_hz, window)
        band_p00 = p00[band_mask]
        ch0_sum += band_p00
        block_sum += band_p00
        block_count += 1
        if block_count >= int(args.block_records):
            block_sums.append(block_sum / block_count)
            block_counts.append(block_count)
            block_sum = np.zeros_like(block_sum)
            block_count = 0

        key = record_key(path0)
        path1 = map1.get(key)
        if path1 is None:
            continue
        try:
            values1 = read_record(path1, samples)
        except ValueError:
            continue
        if not auxiliary_ch1_ok(values1, args.ch1_max_range):
            continue
        pp00, pp11, pp01 = one_sided_spectra(
            values0,
            rate_hz,
            window,
            values1,
        )
        pair_p00 += pp00[band_mask]
        pair_p11 += pp11[band_mask]
        pair_p01 += pp01[band_mask]
        paired_keys.append(key)

    if block_count:
        block_sums.append(block_sum / block_count)
        block_counts.append(block_count)

    ch0_mean_psd = ch0_sum / len(accepted_paths0)
    ch0_asd = np.sqrt(np.maximum(ch0_mean_psd, 0.0))

    pair_count = len(paired_keys)
    if pair_count:
        pair_p00 /= pair_count
        pair_p11 /= pair_count
        pair_p01 /= pair_count
        ch1_asd = np.sqrt(np.maximum(pair_p11, 0.0))
        coherence = (
            np.abs(pair_p01) ** 2
            / np.maximum(pair_p00 * pair_p11, np.finfo(float).tiny)
        )
        phase = np.angle(pair_p01)
    else:
        ch1_asd = np.full_like(ch0_asd, np.nan)
        coherence = np.full_like(ch0_asd, np.nan)
        phase = np.full_like(ch0_asd, np.nan)

    lines, ch0_baseline_db, ch0_excess_db = detect_lines(
        band_frequency,
        ch0_asd,
        baseline_width_hz=args.baseline_width_hz,
        min_excess_db=args.min_excess_db,
        min_prominence_db=args.min_prominence_db,
        min_spacing_hz=args.min_spacing_hz,
        max_peaks=args.max_peaks,
    )

    if pair_count:
        kernel = odd_kernel_bins(args.baseline_width_hz, bin_hz)
        ch1_baseline_db = local_baseline_db(ch1_asd, kernel)
        ch1_db = 20.0 * np.log10(
            np.maximum(ch1_asd, np.finfo(float).tiny)
        )
        ch1_excess_db = ch1_db - ch1_baseline_db
    else:
        ch1_baseline_db = np.full_like(ch0_asd, np.nan)
        ch1_excess_db = np.full_like(ch0_asd, np.nan)

    block_psd = np.asarray(block_sums, dtype=float)
    for row in lines:
        index = int(row["band_index"])
        row["ch1_auxiliary_asd"] = (
            float(ch1_asd[index]) if pair_count else None
        )
        row["ch1_excess_over_local_baseline_dB"] = (
            float(ch1_excess_db[index]) if pair_count else None
        )
        row["ch0_ch1_coherence"] = (
            float(coherence[index]) if pair_count else None
        )
        row["ch0_ch1_cross_phase_rad"] = (
            float(phase[index]) if pair_count else None
        )
        if pair_count:
            if (
                coherence[index] >= 0.5
                and ch1_excess_db[index] >= args.min_excess_db
            ):
                classification = "common_mode_line_candidate"
            elif (
                coherence[index] < 0.1
                and ch1_excess_db[index] < args.min_excess_db
            ):
                classification = "ch0_specific_line_candidate"
            else:
                classification = "inconclusive"
        else:
            classification = "inconclusive_no_auxiliary_pairs"
        row["classification"] = classification

        block_rows = []
        for block_index, (psd, count) in enumerate(
            zip(block_psd, block_counts)
        ):
            baseline_psd = line_baseline_psd(
                band_frequency,
                psd,
                row["frequency_Hz"],
            )
            line_psd = float(psd[index])
            excess_db = (
                float(
                    10.0
                    * np.log10(
                        max(line_psd, np.finfo(float).tiny)
                        / max(baseline_psd, np.finfo(float).tiny)
                    )
                )
                if np.isfinite(baseline_psd)
                else None
            )
            block_rows.append(
                {
                    "block_index": block_index,
                    "records": int(count),
                    "line_asd": float(np.sqrt(max(line_psd, 0.0))),
                    "line_over_local_baseline_dB": excess_db,
                }
            )
        row["block_stability"] = block_rows
        values = [
            item["line_over_local_baseline_dB"]
            for item in block_rows
            if item["line_over_local_baseline_dB"] is not None
        ]
        row["block_stability_summary"] = (
            {
                "min_dB": float(np.min(values)),
                "median_dB": float(np.median(values)),
                "max_dB": float(np.max(values)),
                "std_dB": float(np.std(values)),
            }
            if values
            else None
        )

    return {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "no_notch_applied": True,
        "no_model_parameter_fit": True,
        "case": {
            "label": repeat_case["label"],
            "comparison_source": comparison_source,
            "experiment_path": str(experiment_path),
        },
        "acquisition": {
            "rate_Hz": rate_hz,
            "samples": samples,
            "native_fft_bin_Hz": bin_hz,
            "ch0_production_accepted_records": len(accepted_paths0),
            "ch0_ch1_auxiliary_exact_key_pairs": pair_count,
        },
        "frequency_band_Hz": {
            "min": float(args.min_hz),
            "max": float(args.max_hz),
            "native_binning_preserved": True,
        },
        "line_detection": {
            "baseline_width_Hz": float(args.baseline_width_hz),
            "minimum_excess_dB": float(args.min_excess_db),
            "minimum_prominence_dB": float(args.min_prominence_db),
            "minimum_spacing_Hz": float(args.min_spacing_hz),
            "max_reported_peaks": int(args.max_peaks),
            "lines": lines,
            "harmonic_relations": harmonic_relations(lines),
        },
        "ch1_auxiliary_semantics": {
            "role": "auxiliary_non_production",
            "pairing": "exact event key with CH0 production-accepted records",
            "acceptance": (
                "finite, exact sample length, raw peak-to-peak range below "
                f"{args.ch1_max_range:g}; this is not independent CH1 "
                "production acceptance"
            ),
        },
        "interpretation": {
            "common_mode_rule": (
                "coherence >= 0.5 and CH1 line excess >= configured CH0 "
                "detection threshold"
            ),
            "ch0_specific_rule": (
                "coherence < 0.1 and CH1 line excess below configured threshold"
            ),
            "otherwise": "inconclusive",
            "guardrail": (
                "Detected narrow lines must not be absorbed into TES/readout "
                "continuum parameters merely to improve a broadband fit. "
                "This diagnostic does not authorize removing or notching a "
                "line from production data; hardware/source provenance should "
                "be established first."
            ),
        },
        "_plot": {
            "frequency_Hz": band_frequency.tolist(),
            "ch0_asd": ch0_asd.tolist(),
            "ch1_auxiliary_asd": (
                ch1_asd.tolist() if pair_count else None
            ),
            "ch0_baseline_dB": ch0_baseline_db.tolist(),
            "ch0_excess_dB": ch0_excess_db.tolist(),
            "ch1_excess_dB": (
                ch1_excess_db.tolist() if pair_count else None
            ),
            "coherence": coherence.tolist() if pair_count else None,
        },
    }


def write_plot(result: dict, output: Path, show: bool):
    import matplotlib.pyplot as plt

    plot = result["_plot"]
    frequency_khz = np.asarray(plot["frequency_Hz"], dtype=float) / 1000.0
    ch0_asd = np.asarray(plot["ch0_asd"], dtype=float)
    ch0_excess = np.asarray(plot["ch0_excess_dB"], dtype=float)
    lines = result["line_detection"]["lines"]

    fig, axes = plt.subplots(
        4,
        1,
        figsize=(11, 10),
        sharex=False,
        gridspec_kw={"height_ratios": [1.4, 1.0, 1.0, 1.0]},
    )
    ax_asd, ax_excess, ax_coh, ax_block = axes

    ax_asd.plot(frequency_khz, ch0_asd, label="CH0 production mask")
    if plot["ch1_auxiliary_asd"] is not None:
        ax_asd.plot(
            frequency_khz,
            np.asarray(plot["ch1_auxiliary_asd"], dtype=float),
            label="CH1 auxiliary paired",
        )
    ax_asd.set_yscale("log")
    ax_asd.set_ylabel("ASD [raw / sqrt(Hz)]")
    ax_asd.set_title("High-frequency narrow-line diagnostic")
    ax_asd.grid(True, alpha=0.25)
    ax_asd.legend(frameon=False)

    ax_excess.plot(
        frequency_khz,
        ch0_excess,
        label="CH0 excess above local median",
    )
    if plot["ch1_excess_dB"] is not None:
        ax_excess.plot(
            frequency_khz,
            np.asarray(plot["ch1_excess_dB"], dtype=float),
            label="CH1 auxiliary excess",
        )
    ax_excess.axhline(
        result["line_detection"]["minimum_excess_dB"],
        linestyle="--",
        linewidth=1.0,
        label="detection threshold",
    )
    for row in lines:
        ax_excess.axvline(
            row["frequency_Hz"] / 1000.0,
            linestyle=":",
            linewidth=0.8,
        )
    ax_excess.set_ylabel("Line excess [dB]")
    ax_excess.grid(True, alpha=0.25)
    ax_excess.legend(frameon=False, fontsize=8)

    if plot["coherence"] is not None:
        ax_coh.plot(
            frequency_khz,
            np.asarray(plot["coherence"], dtype=float),
        )
        ax_coh.axhline(0.5, linestyle="--", linewidth=1.0)
        ax_coh.set_ylim(-0.02, 1.02)
    else:
        ax_coh.text(
            0.5,
            0.5,
            "No usable CH0/CH1 auxiliary pairs",
            transform=ax_coh.transAxes,
            ha="center",
            va="center",
        )
    ax_coh.set_ylabel("CH0/CH1 coherence")
    ax_coh.grid(True, alpha=0.25)

    for row in lines[:5]:
        blocks = row["block_stability"]
        x = [item["block_index"] for item in blocks]
        y = [item["line_over_local_baseline_dB"] for item in blocks]
        ax_block.plot(
            x,
            y,
            marker="o",
            markersize=3,
            linewidth=1.0,
            label=f"{row['frequency_Hz'] / 1000.0:.3f} kHz",
        )
    ax_block.set_xlabel(
        f"CH0 accepted-record block ({result['acquisition']['ch0_production_accepted_records']} records total)"
    )
    ax_block.set_ylabel("Line / local baseline [dB]")
    ax_block.grid(True, alpha=0.25)
    if lines:
        ax_block.legend(frameon=False, fontsize=8, ncol=2)

    for ax in (ax_asd, ax_excess, ax_coh):
        ax.set_xlim(
            result["frequency_band_Hz"]["min"] / 1000.0,
            result["frequency_band_Hz"]["max"] / 1000.0,
        )
        ax.set_xlabel("Frequency [kHz]")

    fig.suptitle(
        f"{result['case']['label']} | native bin "
        f"{result['acquisition']['native_fft_bin_Hz']:.3f} Hz | "
        f"{result['acquisition']['ch0_production_accepted_records']} CH0 records",
        fontsize=10,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result: dict):
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
    parser.add_argument("--min-hz", type=float, default=50_000.0)
    parser.add_argument("--max-hz", type=float, default=130_000.0)
    parser.add_argument("--baseline-width-hz", type=float, default=2_000.0)
    parser.add_argument("--min-excess-db", type=float, default=3.0)
    parser.add_argument("--min-prominence-db", type=float, default=2.0)
    parser.add_argument("--min-spacing-hz", type=float, default=500.0)
    parser.add_argument("--max-peaks", type=int, default=12)
    parser.add_argument("--block-records", type=int, default=100)
    parser.add_argument("--ch1-max-range", type=float, default=0.2)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    result = run(args)
    write_plot(result, args.figure, args.show)
    payload = cleaned_result(result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "native_fft_bin_Hz": result["acquisition"][
                    "native_fft_bin_Hz"
                ],
                "ch0_records": result["acquisition"][
                    "ch0_production_accepted_records"
                ],
                "auxiliary_pairs": result["acquisition"][
                    "ch0_ch1_auxiliary_exact_key_pairs"
                ],
                "detected_lines": [
                    {
                        "frequency_Hz": row["frequency_Hz"],
                        "excess_dB": row[
                            "excess_over_local_baseline_dB"
                        ],
                        "width_Hz": row[
                            "half_prominence_width_Hz"
                        ],
                        "coherence": row["ch0_ch1_coherence"],
                        "classification": row["classification"],
                    }
                    for row in result["line_detection"]["lines"]
                ],
                "harmonic_relations": result["line_detection"][
                    "harmonic_relations"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
