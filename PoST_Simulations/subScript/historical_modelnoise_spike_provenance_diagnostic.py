"""Trace why narrow high-frequency spikes disappear from historical modelnoise.

This diagnostic reconstructs several analysis paths from the same CH0 raw
records and compares them on the native FFT grid:

* all valid records, raw pre-analysis ASD
* current production accepted mask, raw ASD
* historical PeakThreshold accepted mask, raw ASD
* the same three masks after the historical 10 kHz second-order Bessel
  filtfilt stage
* records rejected by each mask
* stored CH0_noise/modelnoise.txt, when present

It is deliberately model-free.  No TES/readout simulation parameter is fitted
and no notch is applied.  The purpose is to determine whether a narrow line was
hidden by the digital filter, selectively removed by the historical pulse veto,
changed by record-mask provenance, or already differs in the stored target.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal, stats
from scipy.ndimage import median_filter

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_residual_dof_competition_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "historical_modelnoise_spike_provenance_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "historical_modelnoise_spike_provenance_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    one_sided_asd_from_power,
)
from subScript import high_frequency_line_noise_diagnostic as line_diag  # noqa: E402
from subScript import modelnoise_provenance_diagnostic as provenance  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


def natural_record_number(path: Path):
    try:
        return int(path.stem.split("_", 1)[1])
    except (IndexError, ValueError):
        return path.name


def natural_noise_paths(experiment_path: Path):
    root = experiment_path / "CH0_noise" / "rawdata"
    return sorted(root.glob("CH0_*.dat"), key=natural_record_number)


def comparison_order_paths(experiment_path: Path):
    """Ordering used by the tracked repeat comparison-spec implementation."""
    root = experiment_path / "CH0_noise" / "rawdata"
    return sorted(root.glob("CH0_*.dat"))


def read_record(path: Path, samples: int):
    values = np.fromfile(path, dtype=np.float64, offset=4)
    if values.shape != (int(samples),) or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid record: {path}")
    return values


def _float_or_none(value):
    if value is None:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if np.isfinite(result) else None


def read_setting_txt(path: Path):
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(lines) < 15:
        return {}
    labels = (
        "max_input",
        "minimum_input",
        "rate_Hz",
        "samples",
        "samples_to_log",
        "pretrigger_samples",
        "threshold",
    )
    values = {}
    for key, value in zip(labels, lines[8:15]):
        parsed = _float_or_none(value)
        values[key] = parsed if parsed is not None else value
    return values


def discover_historical_settings(
    experiment_path: Path,
    *,
    threshold_override=None,
    presamples_override=None,
    cutoff_override=None,
    acquisition_cutoff=None,
):
    """Resolve historical noise_main settings without silently inventing them."""
    result = {
        "threshold": None,
        "presamples": None,
        "cutoff_Hz": None,
        "eta_uA_per_V": None,
        "sources": {},
    }

    setting_json = experiment_path / "setting.json"
    if setting_json.exists():
        payload = json.loads(setting_json.read_text(encoding="utf-8"))
        config = payload.get("Config", {})
        main = payload.get("main", {})
        threshold = _float_or_none(config.get("threshold"))
        presamples = _float_or_none(config.get("presamples"))
        cutoff = _float_or_none(main.get("cutoff"))
        eta = _float_or_none(config.get("eta_uA_per_V"))
        if eta is None:
            eta = _float_or_none(config.get("eta"))
        if threshold is not None:
            result["threshold"] = threshold
            result["sources"]["threshold"] = str(setting_json)
        if presamples is not None:
            result["presamples"] = int(round(presamples))
            result["sources"]["presamples"] = str(setting_json)
        if cutoff is not None:
            result["cutoff_Hz"] = cutoff
            result["sources"]["cutoff_Hz"] = str(setting_json)
        if eta is not None:
            result["eta_uA_per_V"] = eta
            result["sources"]["eta_uA_per_V"] = str(setting_json)

    setting_txt = experiment_path / "Setting.txt"
    if setting_txt.exists():
        parsed = read_setting_txt(setting_txt)
        if result["threshold"] is None:
            value = _float_or_none(parsed.get("threshold"))
            if value is not None:
                result["threshold"] = value
                result["sources"]["threshold"] = str(setting_txt)
        if result["presamples"] is None:
            value = _float_or_none(parsed.get("pretrigger_samples"))
            if value is not None:
                result["presamples"] = int(round(value))
                result["sources"]["presamples"] = str(setting_txt)

    if threshold_override is not None:
        result["threshold"] = float(threshold_override)
        result["sources"]["threshold"] = "CLI --peak-threshold"
    if presamples_override is not None:
        result["presamples"] = int(presamples_override)
        result["sources"]["presamples"] = "CLI --presamples"
    if cutoff_override is not None:
        result["cutoff_Hz"] = float(cutoff_override)
        result["sources"]["cutoff_Hz"] = "CLI --cutoff-hz"

    if result["cutoff_Hz"] is None and acquisition_cutoff is not None:
        result["cutoff_Hz"] = float(acquisition_cutoff)
        result["sources"]["cutoff_Hz"] = (
            "tracked comparison acquisition fallback; not independently "
            "recovered from historical setting.json"
        )

    result["historical_peak_mask_available"] = bool(
        result["threshold"] is not None
        and result["presamples"] is not None
    )
    result["historical_bessel_available"] = bool(
        result["cutoff_Hz"] is not None
        and result["cutoff_Hz"] > 0.0
    )
    return result


def historical_peak_value(values: np.ndarray, presamples: int):
    start = int(presamples) - 1000
    stop = start + 500
    if start < 0 or stop > len(values):
        raise ValueError(
            "historical baseline window is outside the record; "
            f"presamples={presamples}, record_length={len(values)}"
        )
    baseline = float(np.mean(values[start:stop]))
    peak = float(np.max(values - baseline))
    return baseline, peak


def historical_peak_accept(values, presamples, threshold):
    baseline, peak = historical_peak_value(values, presamples)
    # noise_main.py contains (base <= -3 and base >= 3), which can never be
    # true.  Therefore the effective historical veto is peak >= threshold.
    return bool(peak < float(threshold)), baseline, peak


def bessel_coefficients(rate_hz: float, cutoff_hz: float):
    nyquist = float(rate_hz) / 2.0
    if not 0.0 < float(cutoff_hz) < nyquist:
        raise ValueError("Bessel cutoff must lie between 0 and Nyquist")
    return signal.bessel(2, float(cutoff_hz) / nyquist, "low")


def filtered_record(values, rate_hz, cutoff_hz):
    x = np.asarray(values, dtype=float) - float(np.mean(values))
    b, a = bessel_coefficients(rate_hz, cutoff_hz)
    return signal.filtfilt(b, a, x)


def periodogram_power(values, rate_hz, window):
    x = np.asarray(values, dtype=float)
    fft = np.fft.rfft(x * window)
    return np.abs(fft) ** 2


def pathway_accumulator(samples):
    return {
        "power_sum": np.zeros(int(samples) // 2 + 1, dtype=float),
        "count": 0,
    }


def add_power(accumulator, power):
    accumulator["power_sum"] += np.asarray(power, dtype=float)
    accumulator["count"] += 1


def asd_from_accumulator(accumulator, samples, rate_hz, window):
    if accumulator["count"] <= 0:
        return None
    window_power_gain = float(np.sqrt(np.mean(np.asarray(window) ** 2)))
    return one_sided_asd_from_power(
        accumulator["power_sum"] / accumulator["count"],
        int(samples),
        float(rate_hz),
        window_power_gain,
    )


def normalize_at(frequency, values, reference_hz=1000.0):
    frequency = np.asarray(frequency, dtype=float)
    values = np.asarray(values, dtype=float)
    index = int(np.argmin(np.abs(frequency - float(reference_hz))))
    reference = float(values[index])
    if not np.isfinite(reference) or reference <= 0.0:
        raise ValueError("invalid normalization reference")
    return values / reference


def local_line_metric(
    frequency,
    asd,
    center_hz,
    *,
    baseline_half_width_hz=1000.0,
    exclusion_half_width_hz=100.0,
):
    frequency = np.asarray(frequency, dtype=float)
    asd = np.asarray(asd, dtype=float)
    index = int(np.argmin(np.abs(frequency - float(center_hz))))
    distance = np.abs(frequency - float(frequency[index]))
    mask = (
        (distance <= float(baseline_half_width_hz))
        & (distance >= float(exclusion_half_width_hz))
    )
    if not np.any(mask):
        baseline_asd = float("nan")
        excess_db = float("nan")
    else:
        baseline_psd = float(np.median(asd[mask] ** 2))
        baseline_asd = float(np.sqrt(max(baseline_psd, 0.0)))
        excess_db = float(
            20.0
            * np.log10(
                max(float(asd[index]), np.finfo(float).tiny)
                / max(baseline_asd, np.finfo(float).tiny)
            )
        )
    return {
        "frequency_Hz": float(frequency[index]),
        "asd": float(asd[index]),
        "local_baseline_asd": baseline_asd,
        "excess_over_local_baseline_dB": excess_db,
    }


def theoretical_filtfilt_transfer_db(rate_hz, cutoff_hz, frequency_hz):
    b, a = bessel_coefficients(rate_hz, cutoff_hz)
    w = 2.0 * np.pi * np.asarray(frequency_hz, dtype=float) / float(rate_hz)
    _, response = signal.freqz(b, a, worN=w)
    # filtfilt applies the single-pass response forward and backward, so the
    # amplitude transfer is |H|^2.
    amplitude = np.abs(response) ** 2
    return 20.0 * np.log10(np.maximum(amplitude, np.finfo(float).tiny))


def mask_overlap(a, b):
    a = set(a)
    b = set(b)
    union = a | b
    return {
        "a_count": len(a),
        "b_count": len(b),
        "intersection_count": len(a & b),
        "a_only_count": len(a - b),
        "b_only_count": len(b - a),
        "union_count": len(union),
        "jaccard": float(len(a & b) / len(union)) if union else 1.0,
    }


def shape_audit(stored, fresh, frequency, smooth_width_hz):
    return provenance.shape_difference_analysis(
        np.asarray(stored, dtype=float),
        np.asarray(fresh, dtype=float),
        np.asarray(frequency, dtype=float),
        smooth_width_hz=float(smooth_width_hz),
    )


def load_stored_modelnoise(path, frequency):
    if not path.exists():
        return None
    stored = np.asarray(np.loadtxt(path), dtype=float)
    if stored.ndim != 1 or stored.size < 2:
        raise ValueError(f"invalid stored modelnoise: {path}")
    frequency = np.asarray(frequency, dtype=float)
    if stored.shape == frequency.shape:
        return stored
    source_frequency = np.linspace(
        float(frequency[0]),
        float(frequency[-1]),
        stored.size,
    )
    return np.interp(frequency, source_frequency, stored)


def _safe_corr(x, y, method):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    valid = np.isfinite(x) & np.isfinite(y)
    x = x[valid]
    y = y[valid]
    if x.size < 3 or np.std(x) == 0.0 or np.std(y) == 0.0:
        return None
    if method == "pearson":
        return float(np.corrcoef(x, y)[0, 1])
    if method == "spearman":
        return float(stats.spearmanr(x, y).statistic)
    raise ValueError(method)


def detect_candidate_lines(
    frequency,
    asd,
    min_hz,
    max_hz,
    baseline_width_hz,
    min_excess_db,
    min_prominence_db,
    min_spacing_hz,
    max_peaks,
):
    mask = (
        (frequency >= float(min_hz))
        & (frequency <= float(max_hz))
    )
    rows, _baseline, _excess = line_diag.detect_lines(
        np.asarray(frequency)[mask],
        np.asarray(asd)[mask],
        baseline_width_hz=baseline_width_hz,
        min_excess_db=min_excess_db,
        min_prominence_db=min_prominence_db,
        min_spacing_hz=min_spacing_hz,
        max_peaks=max_peaks,
    )
    return [float(row["frequency_Hz"]) for row in rows]


def current_accepted_names(
    comparison,
    experiment_path,
    comparison_source,
):
    indices = [
        int(value)
        for value in comparison["acquisition"]["accepted_record_indices"]
    ]
    if comparison_source == "tracked_comparison_spec_recomputed_mask":
        ordered = comparison_order_paths(experiment_path)
        ordering = "lexicographic, matching build_comparison_from_spec"
    else:
        # Existing comparison summaries in this project were generated using
        # sorted glob ordering.  Preserve that contract for index lookup.
        ordered = comparison_order_paths(experiment_path)
        ordering = "lexicographic, matching tracked comparison-summary indices"
    if indices and max(indices) >= len(ordered):
        raise IndexError("current accepted index exceeds raw record list")
    return {ordered[index].name for index in indices}, ordering, ordered


def convert_derived_units(asd, eta_uA_per_V):
    if asd is None:
        return None, "unavailable"
    if eta_uA_per_V is None:
        return np.asarray(asd, dtype=float), "raw_input_unit/rtHz"
    return (
        np.asarray(asd, dtype=float) * float(eta_uA_per_V) * 1.0e6,
        "pA/rtHz",
    )


def run(args):
    config = json.loads(args.config.read_text(encoding="utf-8"))
    manifest_path = ident.resolve_config_path(config["manifest"], args.config)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = shared.normalize_manifest(manifest, manifest_path)
    case = ident.case_by_label(cases, config["repeat_case_label"])
    comparison, experiment_path, comparison_source = ident.comparison_for_case(
        case
    )
    if args.experiment_path is not None:
        experiment_path = args.experiment_path

    acquisition = comparison["acquisition"]
    rate_hz = float(acquisition["rate_Hz"])
    samples = int(acquisition["samples"])
    acquisition_cutoff = float(acquisition.get("cutoff_Hz", 0.0))
    frequency = np.fft.rfftfreq(samples, d=1.0 / rate_hz)
    bin_hz = float(rate_hz / samples)
    window = np.hanning(samples)

    settings = discover_historical_settings(
        experiment_path,
        threshold_override=args.peak_threshold,
        presamples_override=args.presamples,
        cutoff_override=args.cutoff_hz,
        acquisition_cutoff=acquisition_cutoff,
    )

    current_names, current_ordering, current_ordered_paths = (
        current_accepted_names(
            comparison,
            experiment_path,
            comparison_source,
        )
    )
    natural_paths = natural_noise_paths(experiment_path)
    if not natural_paths:
        raise FileNotFoundError(
            f"no CH0 raw records under {experiment_path}"
        )
    lexical_names = [path.name for path in current_ordered_paths]
    natural_names = [path.name for path in natural_paths]

    pathway_names = (
        "all_valid_raw",
        "current_mask_raw",
        "current_rejected_raw",
        "all_valid_bessel",
        "current_mask_bessel",
        "current_rejected_bessel",
        "historical_peak_raw",
        "historical_peak_rejected_raw",
        "historical_peak_bessel",
        "historical_peak_rejected_bessel",
    )
    accumulators = {
        name: pathway_accumulator(samples)
        for name in pathway_names
    }

    valid_names = set()
    historical_names = set()
    historical_rejected_names = set()
    record_peak_rows = []
    invalid_records = []

    bessel_available = settings["historical_bessel_available"]
    historical_mask_available = settings["historical_peak_mask_available"]

    for path in natural_paths:
        try:
            values = read_record(path, samples)
        except ValueError:
            invalid_records.append(path.name)
            continue

        valid_names.add(path.name)
        raw_centered = values - float(np.mean(values))
        raw_power = periodogram_power(raw_centered, rate_hz, window)
        add_power(accumulators["all_valid_raw"], raw_power)

        is_current = path.name in current_names
        if is_current:
            add_power(accumulators["current_mask_raw"], raw_power)
        else:
            add_power(accumulators["current_rejected_raw"], raw_power)

        hist_accept = None
        hist_baseline = None
        hist_peak = None
        if historical_mask_available:
            hist_accept, hist_baseline, hist_peak = historical_peak_accept(
                values,
                settings["presamples"],
                settings["threshold"],
            )
            if hist_accept:
                historical_names.add(path.name)
                add_power(accumulators["historical_peak_raw"], raw_power)
            else:
                historical_rejected_names.add(path.name)
                add_power(
                    accumulators["historical_peak_rejected_raw"],
                    raw_power,
                )

        record_peak_rows.append(
            {
                "path": str(path),
                "name": path.name,
                "current_accept": bool(is_current),
                "historical_accept": hist_accept,
                "historical_baseline": hist_baseline,
                "historical_peak": hist_peak,
            }
        )

        if bessel_available:
            filtered = filtered_record(
                values,
                rate_hz,
                settings["cutoff_Hz"],
            )
            filtered_power = periodogram_power(
                filtered,
                rate_hz,
                window,
            )
            add_power(accumulators["all_valid_bessel"], filtered_power)
            if is_current:
                add_power(
                    accumulators["current_mask_bessel"],
                    filtered_power,
                )
            else:
                add_power(
                    accumulators["current_rejected_bessel"],
                    filtered_power,
                )
            if historical_mask_available:
                if hist_accept:
                    add_power(
                        accumulators["historical_peak_bessel"],
                        filtered_power,
                    )
                else:
                    add_power(
                        accumulators[
                            "historical_peak_rejected_bessel"
                        ],
                        filtered_power,
                    )

    derived = {}
    for name, accumulator in accumulators.items():
        raw_asd = asd_from_accumulator(
            accumulator,
            samples,
            rate_hz,
            window,
        )
        converted, units = convert_derived_units(
            raw_asd,
            settings["eta_uA_per_V"],
        )
        derived[name] = {
            "asd": converted,
            "count": int(accumulator["count"]),
            "units": units,
        }

    primary_for_detection = derived["current_mask_raw"]["asd"]
    if primary_for_detection is None:
        primary_for_detection = derived["all_valid_raw"]["asd"]
    candidate_frequencies = detect_candidate_lines(
        frequency,
        primary_for_detection,
        args.min_hz,
        args.max_hz,
        args.baseline_width_hz,
        args.min_excess_db,
        args.min_prominence_db,
        args.min_spacing_hz,
        args.max_peaks,
    )

    # Second pass: direct record-wise relationship between the historical
    # time-domain peak veto and the power in each detected spectral line.
    line_indices = {
        f: int(np.argmin(np.abs(frequency - f)))
        for f in candidate_frequencies
    }
    recordwise = {f: [] for f in candidate_frequencies}
    if historical_mask_available and candidate_frequencies:
        row_by_name = {
            row["name"]: row
            for row in record_peak_rows
        }
        for path in natural_paths:
            row = row_by_name.get(path.name)
            if row is None:
                continue
            values = read_record(path, samples)
            centered = values - float(np.mean(values))
            power = periodogram_power(centered, rate_hz, window)
            for f, index in line_indices.items():
                recordwise[f].append(
                    {
                        "peak": float(row["historical_peak"]),
                        "line_power": float(power[index]),
                        "accepted": bool(row["historical_accept"]),
                    }
                )

    stored_path = (
        args.stored_modelnoise
        if args.stored_modelnoise is not None
        else experiment_path / "CH0_noise" / "modelnoise.txt"
    )
    stored = load_stored_modelnoise(stored_path, frequency)
    stored_units = "pA/rtHz" if stored is not None else "unavailable"

    pathway_output = {}
    for name, row in derived.items():
        asd = row["asd"]
        pathway_output[name] = {
            "record_count": row["count"],
            "units": row["units"],
            "available": bool(asd is not None),
            "line_metrics": (
                {
                    f"{f:g}": local_line_metric(
                        frequency,
                        asd,
                        f,
                        baseline_half_width_hz=args.line_baseline_half_width_hz,
                        exclusion_half_width_hz=args.line_exclusion_half_width_hz,
                    )
                    for f in candidate_frequencies
                }
                if asd is not None
                else {}
            ),
        }

    if stored is not None:
        pathway_output["stored_modelnoise"] = {
            "record_count": None,
            "units": stored_units,
            "available": True,
            "path": str(stored_path),
            "line_metrics": {
                f"{f:g}": local_line_metric(
                    frequency,
                    stored,
                    f,
                    baseline_half_width_hz=args.line_baseline_half_width_hz,
                    exclusion_half_width_hz=args.line_exclusion_half_width_hz,
                )
                for f in candidate_frequencies
            },
        }
    else:
        pathway_output["stored_modelnoise"] = {
            "record_count": None,
            "units": stored_units,
            "available": False,
            "path": str(stored_path),
            "line_metrics": {},
        }

    bessel_evidence = {}
    if bessel_available:
        theory_db = theoretical_filtfilt_transfer_db(
            rate_hz,
            settings["cutoff_Hz"],
            candidate_frequencies,
        )
        for i, f in enumerate(candidate_frequencies):
            key = f"{f:g}"
            entry = {
                "frequency_Hz": f,
                "theoretical_filtfilt_amplitude_suppression_dB": float(
                    theory_db[i]
                ),
            }
            for prefix in ("all_valid", "current_mask", "historical_peak"):
                raw_name = f"{prefix}_raw"
                filt_name = f"{prefix}_bessel"
                raw_asd = derived[raw_name]["asd"]
                filt_asd = derived[filt_name]["asd"]
                if raw_asd is None or filt_asd is None:
                    continue
                index = line_indices[f]
                empirical = float(
                    20.0
                    * np.log10(
                        max(float(filt_asd[index]), np.finfo(float).tiny)
                        / max(float(raw_asd[index]), np.finfo(float).tiny)
                    )
                )
                raw_excess = pathway_output[raw_name][
                    "line_metrics"
                ][key]["excess_over_local_baseline_dB"]
                filt_excess = pathway_output[filt_name][
                    "line_metrics"
                ][key]["excess_over_local_baseline_dB"]
                entry[prefix] = {
                    "empirical_absolute_suppression_dB": empirical,
                    "empirical_minus_theory_dB": float(
                        empirical - theory_db[i]
                    ),
                    "line_excess_before_dB": raw_excess,
                    "line_excess_after_dB": filt_excess,
                    "line_excess_change_dB": float(
                        filt_excess - raw_excess
                    ),
                }
            bessel_evidence[key] = entry

    selection_evidence = {}
    for f in candidate_frequencies:
        key = f"{f:g}"
        entry = {"frequency_Hz": f}
        all_metric = pathway_output["all_valid_raw"]["line_metrics"][key]
        current_metric = pathway_output["current_mask_raw"]["line_metrics"][key]
        entry["current_mask"] = {
            "absolute_asd_change_dB": float(
                20.0
                * np.log10(
                    max(current_metric["asd"], np.finfo(float).tiny)
                    / max(all_metric["asd"], np.finfo(float).tiny)
                )
            ),
            "line_excess_change_dB": float(
                current_metric["excess_over_local_baseline_dB"]
                - all_metric["excess_over_local_baseline_dB"]
            ),
        }
        if historical_mask_available:
            hist_metric = pathway_output["historical_peak_raw"][
                "line_metrics"
            ][key]
            reject_metric = pathway_output[
                "historical_peak_rejected_raw"
            ]["line_metrics"].get(key)
            hist_entry = {
                "absolute_asd_change_dB": float(
                    20.0
                    * np.log10(
                        max(hist_metric["asd"], np.finfo(float).tiny)
                        / max(all_metric["asd"], np.finfo(float).tiny)
                    )
                ),
                "line_excess_change_dB": float(
                    hist_metric["excess_over_local_baseline_dB"]
                    - all_metric["excess_over_local_baseline_dB"]
                ),
            }
            if reject_metric is not None:
                hist_entry["rejected_over_accepted_line_asd_dB"] = float(
                    20.0
                    * np.log10(
                        max(reject_metric["asd"], np.finfo(float).tiny)
                        / max(hist_metric["asd"], np.finfo(float).tiny)
                    )
                )
            rows = recordwise.get(f, [])
            peaks = [row["peak"] for row in rows]
            powers = [row["line_power"] for row in rows]
            log_powers = np.log10(
                np.maximum(powers, np.finfo(float).tiny)
            )
            hist_entry["recordwise_peak_line_power_correlation"] = {
                "pearson_peak_vs_log10_line_power": _safe_corr(
                    peaks,
                    log_powers,
                    "pearson",
                ),
                "spearman_peak_vs_line_power": _safe_corr(
                    peaks,
                    powers,
                    "spearman",
                ),
                "record_count": len(rows),
            }
            selection_evidence_key = hist_entry
            entry["historical_peak_threshold"] = selection_evidence_key
        selection_evidence[key] = entry

    stored_shape_audits = {}
    if stored is not None:
        for name in (
            "all_valid_bessel",
            "current_mask_bessel",
            "historical_peak_bessel",
        ):
            asd = derived[name]["asd"]
            if asd is None:
                continue
            stored_shape_audits[name] = shape_audit(
                stored,
                asd,
                frequency,
                args.smooth_width_hz,
            )

    plot_curves = {}
    for name in (
        "all_valid_raw",
        "current_mask_raw",
        "historical_peak_raw",
        "all_valid_bessel",
        "current_mask_bessel",
        "historical_peak_bessel",
    ):
        asd = derived[name]["asd"]
        if asd is not None:
            plot_curves[name] = normalize_at(frequency, asd).tolist()
    if stored is not None:
        plot_curves["stored_modelnoise"] = normalize_at(
            frequency,
            stored,
        ).tolist()

    correlation_plot = None
    if historical_mask_available and candidate_frequencies:
        strongest = candidate_frequencies[0]
        correlation_plot = {
            "frequency_Hz": strongest,
            "peak": [
                row["peak"]
                for row in recordwise[strongest]
            ],
            "line_power": [
                row["line_power"]
                for row in recordwise[strongest]
            ],
            "accepted": [
                row["accepted"]
                for row in recordwise[strongest]
            ],
        }

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "no_notch_applied": True,
        "no_model_parameter_fit": True,
        "purpose": (
            "Separate historical PeakThreshold selection, the 10 kHz "
            "digital Bessel filtfilt, record-mask provenance, and stored "
            "modelnoise provenance as possible reasons narrow high-frequency "
            "spikes were not visible in modelnoise.txt."
        ),
        "case": {
            "label": case["label"],
            "experiment_path": str(experiment_path),
            "comparison_source": comparison_source,
        },
        "acquisition": {
            "rate_Hz": rate_hz,
            "samples": samples,
            "native_fft_bin_Hz": bin_hz,
            "valid_raw_records": len(valid_names),
            "invalid_raw_records": len(invalid_records),
            "invalid_record_names": invalid_records,
        },
        "historical_noise_main_reconstruction": {
            "settings": settings,
            "effective_record_veto": (
                "peak=max(raw-baseline); reject peak>=threshold. "
                "The historical expression "
                "(base <= -3 and base >= 3) is unsatisfiable and therefore "
                "does not veto records."
            ),
            "baseline_window": (
                "[presamples-1000 : presamples-500), matching "
                "gp.baseline(data, presamples, 1000, 500)"
            ),
            "fft_pipeline": (
                "record mean removal -> second-order digital Bessel low-pass "
                "-> scipy.signal.filtfilt -> Hann -> rFFT power average -> "
                "one-sided ASD"
            ),
        },
        "path_ordering_audit": {
            "current_mask_index_ordering": current_ordering,
            "historical_noise_main_ordering": "natural/natsorted",
            "lexicographic_and_natural_order_identical": bool(
                lexical_names == natural_names
            ),
            "note": (
                "Masks are compared by filename, not by raw list index, so "
                "different sorting conventions cannot silently select "
                "different records."
            ),
        },
        "masks": {
            "current_mask_count": len(current_names & valid_names),
            "historical_peak_mask_count": (
                len(historical_names)
                if historical_mask_available
                else None
            ),
            "historical_peak_rejected_count": (
                len(historical_rejected_names)
                if historical_mask_available
                else None
            ),
            "current_vs_historical": (
                mask_overlap(
                    current_names & valid_names,
                    historical_names,
                )
                if historical_mask_available
                else None
            ),
        },
        "candidate_lines_Hz": candidate_frequencies,
        "pathways": pathway_output,
        "selection_evidence": selection_evidence,
        "bessel_evidence": bessel_evidence,
        "stored_modelnoise": {
            "path": str(stored_path),
            "available": bool(stored is not None),
            "stored_vs_reconstructed_shape_audits": stored_shape_audits,
        },
        "interpretation_guardrails": [
            (
                "A large absolute Bessel suppression with little change in "
                "line-over-local-baseline dB means the line was attenuated "
                "with the surrounding continuum rather than selectively "
                "removed."
            ),
            (
                "A large reduction caused by the historical PeakThreshold "
                "mask, together with positive recordwise peak/line-power "
                "correlation and stronger lines in rejected records, is "
                "evidence that pulse veto selection also selected on the "
                "spectral line."
            ),
            (
                "A stored-vs-historical-reconstruction narrow discrepancy "
                "that remains after reproducing both mask and Bessel stage "
                "indicates additional target-file provenance differences."
            ),
            (
                "None of these diagnostics identifies the hardware source of "
                "the line and none authorizes notching it from production "
                "data."
            ),
        ],
        "_plot": {
            "frequency_Hz": frequency.tolist(),
            "curves_normalized_1kHz": plot_curves,
            "correlation": correlation_plot,
        },
    }
    return result


def write_plot(result, output, args):
    import matplotlib.pyplot as plt

    frequency = np.asarray(result["_plot"]["frequency_Hz"], dtype=float)
    band = (
        (frequency >= float(args.min_hz))
        & (frequency <= float(args.max_hz))
    )
    freq_khz = frequency[band] / 1000.0
    curves = result["_plot"]["curves_normalized_1kHz"]

    fig, axes = plt.subplots(
        5,
        1,
        figsize=(12, 13),
        gridspec_kw={"height_ratios": [1.2, 1.2, 1.0, 1.0, 1.0]},
    )
    ax_raw, ax_post, ax_excess, ax_corr, ax_bessel = axes

    raw_names = (
        ("all_valid_raw", "all valid raw"),
        ("current_mask_raw", "current mask raw"),
        ("historical_peak_raw", "historical PeakThreshold raw"),
    )
    for name, label in raw_names:
        if name in curves:
            ax_raw.plot(
                freq_khz,
                np.asarray(curves[name], dtype=float)[band],
                label=label,
            )
    ax_raw.set_yscale("log")
    ax_raw.set_ylabel("Normalized ASD")
    ax_raw.set_title("Pre-analysis paths")
    ax_raw.legend(frameon=False, fontsize=8)
    ax_raw.grid(True, alpha=0.25)

    post_names = (
        ("all_valid_bessel", "all valid + 10 kHz Bessel"),
        ("current_mask_bessel", "current mask + Bessel"),
        (
            "historical_peak_bessel",
            "historical PeakThreshold + Bessel",
        ),
        ("stored_modelnoise", "stored modelnoise.txt"),
    )
    for name, label in post_names:
        if name in curves:
            ax_post.plot(
                freq_khz,
                np.asarray(curves[name], dtype=float)[band],
                label=label,
            )
    ax_post.set_yscale("log")
    ax_post.set_ylabel("Normalized ASD")
    ax_post.set_title("Post-analysis / stored target paths")
    ax_post.legend(frameon=False, fontsize=8)
    ax_post.grid(True, alpha=0.25)

    candidate = result["candidate_lines_Hz"]
    pathway_names = (
        "all_valid_raw",
        "current_mask_raw",
        "historical_peak_raw",
        "current_mask_bessel",
        "historical_peak_bessel",
        "stored_modelnoise",
    )
    x = np.arange(len(candidate), dtype=float)
    offsets = np.linspace(-0.3, 0.3, len(pathway_names))
    for offset, name in zip(offsets, pathway_names):
        row = result["pathways"].get(name)
        if not row or not row.get("available"):
            continue
        values = [
            row["line_metrics"][f"{f:g}"][
                "excess_over_local_baseline_dB"
            ]
            for f in candidate
        ]
        ax_excess.plot(
            x + offset,
            values,
            marker="o",
            linestyle="none",
            label=name,
        )
    ax_excess.set_xticks(x)
    ax_excess.set_xticklabels(
        [f"{f / 1000.0:.3f}" for f in candidate],
        rotation=45,
        ha="right",
    )
    ax_excess.set_ylabel("Line excess [dB]")
    ax_excess.set_title("What actually removes the local spike?")
    ax_excess.legend(frameon=False, fontsize=7, ncol=2)
    ax_excess.grid(True, alpha=0.25)

    correlation = result["_plot"]["correlation"]
    if correlation is not None:
        peak = np.asarray(correlation["peak"], dtype=float)
        power = np.asarray(correlation["line_power"], dtype=float)
        accepted = np.asarray(correlation["accepted"], dtype=bool)
        if np.any(accepted):
            ax_corr.scatter(
                peak[accepted],
                power[accepted],
                s=10,
                alpha=0.5,
                label="historical accepted",
            )
        if np.any(~accepted):
            ax_corr.scatter(
                peak[~accepted],
                power[~accepted],
                s=10,
                alpha=0.5,
                label="historical rejected",
            )
        threshold = result["historical_noise_main_reconstruction"][
            "settings"
        ]["threshold"]
        if threshold is not None:
            ax_corr.axvline(
                threshold,
                linestyle="--",
                linewidth=1.0,
                label="PeakThreshold",
            )
        ax_corr.set_yscale("log")
        ax_corr.set_xlabel("Historical time-domain peak")
        ax_corr.set_ylabel("Native-bin FFT power")
        ax_corr.set_title(
            f"PeakThreshold coupling to "
            f"{correlation['frequency_Hz'] / 1000.0:.3f} kHz line"
        )
        ax_corr.legend(frameon=False, fontsize=8)
    else:
        ax_corr.text(
            0.5,
            0.5,
            "Historical PeakThreshold settings unavailable",
            transform=ax_corr.transAxes,
            ha="center",
            va="center",
        )
    ax_corr.grid(True, alpha=0.25)

    bessel = result["bessel_evidence"]
    if bessel:
        line_freq = np.asarray(
            [float(row["frequency_Hz"]) for row in bessel.values()],
            dtype=float,
        )
        theory = np.asarray(
            [
                float(row["theoretical_filtfilt_amplitude_suppression_dB"])
                for row in bessel.values()
            ],
            dtype=float,
        )
        order = np.argsort(line_freq)
        ax_bessel.plot(
            line_freq[order] / 1000.0,
            theory[order],
            marker="o",
            label="theoretical filtfilt",
        )
        for prefix in ("current_mask", "historical_peak"):
            points_x = []
            points_y = []
            for row in bessel.values():
                if prefix not in row:
                    continue
                points_x.append(row["frequency_Hz"] / 1000.0)
                points_y.append(
                    row[prefix]["empirical_absolute_suppression_dB"]
                )
            if points_x:
                ax_bessel.plot(
                    points_x,
                    points_y,
                    marker="x",
                    linestyle="none",
                    label=f"empirical {prefix}",
                )
        ax_bessel.legend(frameon=False, fontsize=8)
    ax_bessel.set_xlabel("Frequency [kHz]")
    ax_bessel.set_ylabel("Filtered/raw ASD [dB]")
    ax_bessel.set_title("10 kHz Bessel: theory vs measured suppression")
    ax_bessel.grid(True, alpha=0.25)

    for ax in (ax_raw, ax_post):
        ax.set_xlim(args.min_hz / 1000.0, args.max_hz / 1000.0)
        ax.set_xlabel("Frequency [kHz]")

    fig.suptitle(
        (
            f"{result['case']['label']} | "
            f"{result['acquisition']['native_fft_bin_Hz']:.3f} Hz/bin"
        ),
        fontsize=10,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if args.show:
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
    parser.add_argument("--stored-modelnoise", type=Path)
    parser.add_argument("--peak-threshold", type=float)
    parser.add_argument("--presamples", type=int)
    parser.add_argument("--cutoff-hz", type=float)
    parser.add_argument("--min-hz", type=float, default=50_000.0)
    parser.add_argument("--max-hz", type=float, default=130_000.0)
    parser.add_argument("--baseline-width-hz", type=float, default=2_000.0)
    parser.add_argument("--min-excess-db", type=float, default=3.0)
    parser.add_argument("--min-prominence-db", type=float, default=2.0)
    parser.add_argument("--min-spacing-hz", type=float, default=500.0)
    parser.add_argument("--max-peaks", type=int, default=12)
    parser.add_argument(
        "--line-baseline-half-width-hz",
        type=float,
        default=1000.0,
    )
    parser.add_argument(
        "--line-exclusion-half-width-hz",
        type=float,
        default=100.0,
    )
    parser.add_argument("--smooth-width-hz", type=float, default=1000.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    result = run(args)
    write_plot(result, args.figure, args)
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
                "historical_settings": result[
                    "historical_noise_main_reconstruction"
                ]["settings"],
                "mask_counts": result["masks"],
                "candidate_lines_Hz": result["candidate_lines_Hz"],
                "stored_modelnoise_available": result[
                    "stored_modelnoise"
                ]["available"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
