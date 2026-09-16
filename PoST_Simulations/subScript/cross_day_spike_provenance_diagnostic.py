"""Cross-day test of narrow high-frequency line provenance.

Compare, on the same native FFT grid:
  1. 2024-12-05 repeat raw records using its tracked current acceptance mask
  2. 2024-12-06 reference raw records using its tracked accepted mask
  3. 2024-12-06 stored CH0_noise/modelnoise.txt

For each day the raw records are reconstructed both before and after the same
10 kHz second-order Bessel filtfilt stage. Candidate narrow lines are detected
from both days and the stored target, clustered across nearby frequencies, and
then re-localized independently inside each dataset. This separates:
  * acquisition/day variation,
  * small line-frequency drift,
  * digital-filter attenuation,
  * stored-modelnoise generation/provenance differences,
  * and a line that is actually present in modelnoise but visually hidden on a
    broad log-log plot.

No TES/readout model is fitted and no notch is applied.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import median_filter

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = CONFIG_DIR / "readout_residual_dof_competition_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "cross_day_spike_provenance_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "cross_day_spike_provenance_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import high_frequency_line_noise_diagnostic as line_diag  # noqa: E402
from subScript import historical_modelnoise_spike_provenance_diagnostic as hist  # noqa: E402
from subScript import modelnoise_provenance_diagnostic as provenance  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


def _manifest_cases(config_path: Path):
    config = json.loads(config_path.read_text(encoding="utf-8"))
    manifest_path = ident.resolve_config_path(config["manifest"], config_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cases = shared.normalize_manifest(manifest, manifest_path)
    return config, cases


def _case_by_role(cases, role):
    rows = [case for case in cases if case.get("role") == role]
    if len(rows) != 1:
        raise ValueError(f"expected exactly one {role!r} case, got {len(rows)}")
    return rows[0]


def _accumulate_case(case, experiment_override=None):
    comparison, experiment_path, comparison_source = ident.comparison_for_case(case)
    if experiment_override is not None:
        experiment_path = experiment_override

    acquisition = comparison["acquisition"]
    rate_hz = float(acquisition["rate_Hz"])
    samples = int(acquisition["samples"])
    cutoff_hz = float(acquisition["cutoff_Hz"])
    frequency = np.fft.rfftfreq(samples, d=1.0 / rate_hz)
    window = np.hanning(samples)

    accepted_names, ordering, _ordered = hist.current_accepted_names(
        comparison,
        experiment_path,
        comparison_source,
    )
    paths = hist.natural_noise_paths(experiment_path)
    if not paths:
        raise FileNotFoundError(f"no CH0 raw records found under {experiment_path}")

    raw_acc = hist.pathway_accumulator(samples)
    post_acc = hist.pathway_accumulator(samples)
    all_raw_acc = hist.pathway_accumulator(samples)
    all_post_acc = hist.pathway_accumulator(samples)
    invalid = []

    for path in paths:
        try:
            values = hist.read_record(path, samples)
        except ValueError:
            invalid.append(path.name)
            continue

        centered = values - float(np.mean(values))
        raw_power = hist.periodogram_power(centered, rate_hz, window)
        post = hist.filtered_record(values, rate_hz, cutoff_hz)
        post_power = hist.periodogram_power(post, rate_hz, window)
        hist.add_power(all_raw_acc, raw_power)
        hist.add_power(all_post_acc, post_power)

        if path.name in accepted_names:
            hist.add_power(raw_acc, raw_power)
            hist.add_power(post_acc, post_power)

    raw_asd = hist.asd_from_accumulator(raw_acc, samples, rate_hz, window)
    post_asd = hist.asd_from_accumulator(post_acc, samples, rate_hz, window)
    all_raw_asd = hist.asd_from_accumulator(
        all_raw_acc, samples, rate_hz, window
    )
    all_post_asd = hist.asd_from_accumulator(
        all_post_acc, samples, rate_hz, window
    )
    if raw_asd is None or post_asd is None:
        raise RuntimeError(f"case {case['label']} retained no accepted records")

    settings = hist.discover_historical_settings(
        experiment_path,
        acquisition_cutoff=cutoff_hz,
    )
    return {
        "label": case["label"],
        "role": case["role"],
        "comparison_source": comparison_source,
        "experiment_path": experiment_path,
        "rate_Hz": rate_hz,
        "samples": samples,
        "cutoff_Hz": cutoff_hz,
        "frequency_Hz": frequency,
        "accepted_record_count": int(raw_acc["count"]),
        "all_valid_record_count": int(all_raw_acc["count"]),
        "invalid_records": invalid,
        "mask_ordering": ordering,
        "raw_asd": raw_asd,
        "post_bessel_asd": post_asd,
        "all_raw_asd": all_raw_asd,
        "all_post_bessel_asd": all_post_asd,
        "historical_settings": settings,
    }


def _detected_frequencies(
    frequency,
    asd,
    *,
    min_hz,
    max_hz,
    baseline_width_hz,
    min_excess_db,
    min_prominence_db,
    min_spacing_hz,
    max_peaks,
):
    mask = (frequency >= min_hz) & (frequency <= max_hz)
    rows, _base, _excess = line_diag.detect_lines(
        frequency[mask],
        asd[mask],
        baseline_width_hz=baseline_width_hz,
        min_excess_db=min_excess_db,
        min_prominence_db=min_prominence_db,
        min_spacing_hz=min_spacing_hz,
        max_peaks=max_peaks,
    )
    return [float(row["frequency_Hz"]) for row in rows]


def cluster_frequencies(frequencies, tolerance_hz):
    values = sorted(float(value) for value in frequencies)
    if not values:
        return []
    clusters = [[values[0]]]
    for value in values[1:]:
        if value - clusters[-1][-1] <= float(tolerance_hz):
            clusters[-1].append(value)
        else:
            clusters.append([value])
    return [
        {
            "anchor_Hz": float(np.mean(cluster)),
            "members_Hz": cluster,
            "min_Hz": float(min(cluster)),
            "max_Hz": float(max(cluster)),
        }
        for cluster in clusters
    ]


def _local_excess_curve(frequency, asd, baseline_width_hz):
    bin_hz = float(np.median(np.diff(frequency)))
    kernel = line_diag.odd_kernel_bins(baseline_width_hz, bin_hz)
    db = 20.0 * np.log10(np.maximum(asd, np.finfo(float).tiny))
    baseline = median_filter(db, size=kernel, mode="nearest")
    return db - baseline


def localize_line(
    frequency,
    asd,
    anchor_hz,
    search_half_width_hz,
    baseline_width_hz,
):
    frequency = np.asarray(frequency, dtype=float)
    asd = np.asarray(asd, dtype=float)
    excess = _local_excess_curve(frequency, asd, baseline_width_hz)
    mask = np.abs(frequency - float(anchor_hz)) <= float(search_half_width_hz)
    if not np.any(mask):
        return None
    candidates = np.flatnonzero(mask)
    index = int(candidates[np.argmax(excess[mask])])
    return {
        "frequency_Hz": float(frequency[index]),
        "asd": float(asd[index]),
        "excess_over_local_baseline_dB": float(excess[index]),
        "offset_from_anchor_Hz": float(frequency[index] - anchor_hz),
    }


def classify_line(
    repeat_post,
    reference_post,
    stored,
    *,
    material_db=3.0,
    stored_match_db=1.5,
    line_present_db=3.0,
):
    if repeat_post is None or reference_post is None:
        return "insufficient_raw_data"

    day_delta = (
        reference_post["excess_over_local_baseline_dB"]
        - repeat_post["excess_over_local_baseline_dB"]
    )
    stored_delta = None
    if stored is not None:
        stored_delta = (
            stored["excess_over_local_baseline_dB"]
            - reference_post["excess_over_local_baseline_dB"]
        )

    day_material = abs(day_delta) >= material_db
    stored_material = stored_delta is not None and abs(stored_delta) >= material_db
    stored_matches_reference = (
        stored_delta is not None and abs(stored_delta) <= stored_match_db
    )
    stored_present = (
        stored is not None
        and stored["excess_over_local_baseline_dB"] >= line_present_db
    )

    if day_material and stored_matches_reference:
        return "day_or_acquisition_variation_supported"
    if stored_material and not day_material:
        return "stored_generation_difference_supported"
    if day_material and stored_material:
        return "mixed_day_and_stored_difference"
    if stored_present and stored_matches_reference:
        return "line_persists_in_stored_target_visual_concealment_candidate"
    if stored is None:
        return "stored_target_unavailable"
    return "no_material_difference_at_this_line"


def _normalize(values, frequency):
    return hist.normalize_at(frequency, values, reference_hz=1000.0)


def run(args):
    config, cases = _manifest_cases(args.config)
    reference_case = _case_by_role(cases, "reference")
    repeat_case = ident.case_by_label(cases, config["repeat_case_label"])

    repeat = _accumulate_case(repeat_case, args.repeat_experiment_path)
    reference = _accumulate_case(reference_case, args.reference_experiment_path)

    for key in ("rate_Hz", "samples"):
        if repeat[key] != reference[key]:
            raise ValueError(
                f"cross-day comparison requires equal {key}: "
                f"{repeat[key]} vs {reference[key]}"
            )
    if repeat["frequency_Hz"].shape != reference["frequency_Hz"].shape:
        raise ValueError("cross-day FFT grids differ")
    frequency = reference["frequency_Hz"]

    stored_path, stored_discovery = hist.discover_stored_modelnoise_path(
        reference["experiment_path"],
        reference["historical_settings"],
        args.stored_modelnoise,
    )
    stored = (
        hist.load_stored_modelnoise(stored_path, frequency)
        if stored_path is not None
        else None
    )

    detected = []
    detection_sources = {
        "repeat_raw": repeat["raw_asd"],
        "reference_raw": reference["raw_asd"],
        "repeat_post_bessel": repeat["post_bessel_asd"],
        "reference_post_bessel": reference["post_bessel_asd"],
    }
    if stored is not None:
        detection_sources["stored_modelnoise"] = stored

    per_source_detected = {}
    for name, asd in detection_sources.items():
        rows = _detected_frequencies(
            frequency,
            asd,
            min_hz=args.min_hz,
            max_hz=args.max_hz,
            baseline_width_hz=args.baseline_width_hz,
            min_excess_db=args.min_excess_db,
            min_prominence_db=args.min_prominence_db,
            min_spacing_hz=args.min_spacing_hz,
            max_peaks=args.max_peaks_per_source,
        )
        per_source_detected[name] = rows
        detected.extend(rows)

    clusters = cluster_frequencies(detected, args.cluster_tolerance_hz)

    line_rows = []
    for cluster in clusters:
        anchor = cluster["anchor_Hz"]
        localized = {
            "repeat_raw": localize_line(
                frequency,
                repeat["raw_asd"],
                anchor,
                args.search_half_width_hz,
                args.baseline_width_hz,
            ),
            "repeat_post_bessel": localize_line(
                frequency,
                repeat["post_bessel_asd"],
                anchor,
                args.search_half_width_hz,
                args.baseline_width_hz,
            ),
            "reference_raw": localize_line(
                frequency,
                reference["raw_asd"],
                anchor,
                args.search_half_width_hz,
                args.baseline_width_hz,
            ),
            "reference_post_bessel": localize_line(
                frequency,
                reference["post_bessel_asd"],
                anchor,
                args.search_half_width_hz,
                args.baseline_width_hz,
            ),
            "stored_modelnoise": (
                localize_line(
                    frequency,
                    stored,
                    anchor,
                    args.search_half_width_hz,
                    args.baseline_width_hz,
                )
                if stored is not None
                else None
            ),
        }

        repeat_post = localized["repeat_post_bessel"]
        reference_post = localized["reference_post_bessel"]
        stored_line = localized["stored_modelnoise"]
        day_delta = (
            float(
                reference_post["excess_over_local_baseline_dB"]
                - repeat_post["excess_over_local_baseline_dB"]
            )
            if repeat_post is not None and reference_post is not None
            else None
        )
        stored_delta = (
            float(
                stored_line["excess_over_local_baseline_dB"]
                - reference_post["excess_over_local_baseline_dB"]
            )
            if stored_line is not None and reference_post is not None
            else None
        )
        frequency_shift = (
            float(
                reference_post["frequency_Hz"]
                - repeat_post["frequency_Hz"]
            )
            if repeat_post is not None and reference_post is not None
            else None
        )
        stored_frequency_shift = (
            float(
                stored_line["frequency_Hz"]
                - reference_post["frequency_Hz"]
            )
            if stored_line is not None and reference_post is not None
            else None
        )

        line_rows.append(
            {
                "cluster": cluster,
                "localized": localized,
                "reference_minus_repeat_line_excess_dB": day_delta,
                "reference_minus_repeat_peak_frequency_Hz": frequency_shift,
                "stored_minus_reference_reconstruction_line_excess_dB": stored_delta,
                "stored_minus_reference_peak_frequency_Hz": stored_frequency_shift,
                "classification": classify_line(
                    repeat_post,
                    reference_post,
                    stored_line,
                    material_db=args.material_difference_db,
                    stored_match_db=args.stored_match_db,
                    line_present_db=args.line_present_db,
                ),
            }
        )

    shape_audit = None
    if stored is not None:
        shape_audit = provenance.shape_difference_analysis(
            stored,
            reference["post_bessel_asd"],
            frequency,
            smooth_width_hz=args.smooth_width_hz,
        )

    # Explicit anchors guarantee that the previously discussed features are
    # reported even if automatic peak ranking changes.
    explicit_anchors = [
        59_640.0,
        80_000.0,
        88_000.0,
        96_000.0,
        97_410.0,
        100_000.0,
        100_580.0,
        117_330.0,
        119_290.0,
    ]
    explicit_rows = {}
    for anchor in explicit_anchors:
        if not args.min_hz <= anchor <= args.max_hz:
            continue
        explicit_rows[f"{anchor:g}"] = {
            name: (
                localize_line(
                    frequency,
                    asd,
                    anchor,
                    args.explicit_search_half_width_hz,
                    args.baseline_width_hz,
                )
                if asd is not None
                else None
            )
            for name, asd in {
                "repeat_raw": repeat["raw_asd"],
                "repeat_post_bessel": repeat["post_bessel_asd"],
                "reference_raw": reference["raw_asd"],
                "reference_post_bessel": reference["post_bessel_asd"],
                "stored_modelnoise": stored,
            }.items()
        }

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "no_notch_applied": True,
        "no_model_parameter_fit": True,
        "question": (
            "Was the apparent disappearance of narrow high-frequency spikes "
            "caused by acquisition/day variation, by the 10 kHz analysis "
            "filter, or by stored modelnoise generation?"
        ),
        "comparison": {
            "repeat": {
                key: value
                for key, value in repeat.items()
                if key
                not in {
                    "frequency_Hz",
                    "raw_asd",
                    "post_bessel_asd",
                    "all_raw_asd",
                    "all_post_bessel_asd",
                }
            },
            "reference": {
                key: value
                for key, value in reference.items()
                if key
                not in {
                    "frequency_Hz",
                    "raw_asd",
                    "post_bessel_asd",
                    "all_raw_asd",
                    "all_post_bessel_asd",
                }
            },
            "native_fft_bin_Hz": float(np.median(np.diff(frequency))),
        },
        "stored_modelnoise": {
            "available": bool(stored is not None),
            "path": str(stored_path) if stored_path is not None else None,
            "discovery": stored_discovery,
            "reference_reconstruction_shape_audit": shape_audit,
        },
        "detected_lines_by_source_Hz": per_source_detected,
        "line_clusters": line_rows,
        "explicit_anchor_checks": explicit_rows,
        "classification_thresholds_dB": {
            "material_difference": float(args.material_difference_db),
            "stored_match": float(args.stored_match_db),
            "line_present": float(args.line_present_db),
        },
        "interpretation_guardrails": [
            (
                "Day/acquisition variation is supported only when the two raw-"
                "derived post-Bessel spectra differ materially while the "
                "reference reconstruction and its stored modelnoise agree "
                "locally."
            ),
            (
                "Stored-generation difference is supported only when the "
                "reference raw reconstruction and the stored target differ "
                "materially at the same local line while the cross-day change "
                "is not material."
            ),
            (
                "A line with similar local excess in reference reconstruction "
                "and stored modelnoise did not disappear from the stored target; "
                "its absence from a broad plot is then a visualization-scale "
                "issue rather than filtering or record rejection."
            ),
            (
                "Frequency re-localization inside each cluster prevents a "
                "small day-to-day line drift from being mistaken for complete "
                "disappearance."
            ),
        ],
        "_plot": {
            "frequency_Hz": frequency.tolist(),
            "repeat_raw_norm": _normalize(
                repeat["raw_asd"], frequency
            ).tolist(),
            "reference_raw_norm": _normalize(
                reference["raw_asd"], frequency
            ).tolist(),
            "repeat_post_norm": _normalize(
                repeat["post_bessel_asd"], frequency
            ).tolist(),
            "reference_post_norm": _normalize(
                reference["post_bessel_asd"], frequency
            ).tolist(),
            "stored_norm": (
                _normalize(stored, frequency).tolist()
                if stored is not None
                else None
            ),
        },
    }
    return result


def write_plot(result, output, args):
    import matplotlib.pyplot as plt

    plot = result["_plot"]
    frequency = np.asarray(plot["frequency_Hz"], dtype=float)
    band = (
        (frequency >= float(args.min_hz))
        & (frequency <= float(args.max_hz))
    )
    x = frequency[band] / 1000.0

    fig, axes = plt.subplots(
        4,
        1,
        figsize=(12, 11),
        gridspec_kw={"height_ratios": [1.1, 1.1, 1.0, 1.0]},
    )
    raw_ax, post_ax, excess_ax, freq_ax = axes

    raw_ax.plot(
        x,
        np.asarray(plot["repeat_raw_norm"])[band],
        label="2024-12-05 repeat raw",
    )
    raw_ax.plot(
        x,
        np.asarray(plot["reference_raw_norm"])[band],
        label="2024-12-06 reference raw",
    )
    raw_ax.set_yscale("log")
    raw_ax.set_ylabel("ASD / ASD(1 kHz)")
    raw_ax.set_title("Same estimator before 10 kHz analysis filter")
    raw_ax.legend(frameon=False)
    raw_ax.grid(True, alpha=0.25)

    post_ax.plot(
        x,
        np.asarray(plot["repeat_post_norm"])[band],
        label="2024-12-05 + 10 kHz Bessel",
    )
    post_ax.plot(
        x,
        np.asarray(plot["reference_post_norm"])[band],
        label="2024-12-06 + 10 kHz Bessel",
    )
    if plot["stored_norm"] is not None:
        post_ax.plot(
            x,
            np.asarray(plot["stored_norm"])[band],
            label="2024-12-06 stored modelnoise.txt",
        )
    post_ax.set_yscale("log")
    post_ax.set_ylabel("ASD / ASD(1 kHz)")
    post_ax.set_title("Post-analysis reconstruction vs stored target")
    post_ax.legend(frameon=False)
    post_ax.grid(True, alpha=0.25)

    rows = result["line_clusters"]
    if rows:
        labels = [
            f"{row['cluster']['anchor_Hz']/1000.0:.3f}"
            for row in rows
        ]
        pos = np.arange(len(rows), dtype=float)
        repeat_excess = [
            row["localized"]["repeat_post_bessel"][
                "excess_over_local_baseline_dB"
            ]
            for row in rows
        ]
        ref_excess = [
            row["localized"]["reference_post_bessel"][
                "excess_over_local_baseline_dB"
            ]
            for row in rows
        ]
        excess_ax.plot(
            pos,
            repeat_excess,
            marker="o",
            label="12/05 reconstructed",
        )
        excess_ax.plot(
            pos,
            ref_excess,
            marker="o",
            label="12/06 reconstructed",
        )
        if result["stored_modelnoise"]["available"]:
            stored_excess = [
                (
                    row["localized"]["stored_modelnoise"][
                        "excess_over_local_baseline_dB"
                    ]
                    if row["localized"]["stored_modelnoise"] is not None
                    else np.nan
                )
                for row in rows
            ]
            excess_ax.plot(
                pos,
                stored_excess,
                marker="o",
                label="12/06 stored",
            )
        excess_ax.set_xticks(pos)
        excess_ax.set_xticklabels(labels, rotation=45, ha="right")
        excess_ax.legend(frameon=False, fontsize=8)
        excess_ax.set_ylabel("Local line excess [dB]")
        excess_ax.set_title("Line strength: day change vs stored-generation change")
        excess_ax.grid(True, alpha=0.25)

        repeat_freq = [
            row["localized"]["repeat_post_bessel"]["frequency_Hz"] / 1000.0
            for row in rows
        ]
        ref_freq = [
            row["localized"]["reference_post_bessel"]["frequency_Hz"] / 1000.0
            for row in rows
        ]
        freq_ax.plot(pos, repeat_freq, marker="o", label="12/05")
        freq_ax.plot(pos, ref_freq, marker="o", label="12/06 reconstructed")
        if result["stored_modelnoise"]["available"]:
            stored_freq = [
                (
                    row["localized"]["stored_modelnoise"]["frequency_Hz"]
                    / 1000.0
                    if row["localized"]["stored_modelnoise"] is not None
                    else np.nan
                )
                for row in rows
            ]
            freq_ax.plot(
                pos,
                stored_freq,
                marker="o",
                label="12/06 stored",
            )
        freq_ax.set_xticks(pos)
        freq_ax.set_xticklabels(labels, rotation=45, ha="right")
        freq_ax.set_ylabel("Localized peak [kHz]")
        freq_ax.set_title("Did the line disappear, or move?")
        freq_ax.legend(frameon=False, fontsize=8)
        freq_ax.grid(True, alpha=0.25)

    for ax in (raw_ax, post_ax):
        ax.set_xlim(args.min_hz / 1000.0, args.max_hz / 1000.0)
        ax.set_xlabel("Frequency [kHz]")
    freq_ax.set_xlabel("Cluster anchor [kHz]")

    fig.suptitle(
        "Cross-day narrow-line provenance: 2024-12-05 vs 2024-12-06",
        fontsize=11,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if args.show:
        plt.show()
    plt.close(fig)


def cleaned_result(result):
    return {key: value for key, value in result.items() if key != "_plot"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--repeat-experiment-path", type=Path)
    parser.add_argument("--reference-experiment-path", type=Path)
    parser.add_argument("--stored-modelnoise", type=Path)
    parser.add_argument("--min-hz", type=float, default=50_000.0)
    parser.add_argument("--max-hz", type=float, default=130_000.0)
    parser.add_argument("--baseline-width-hz", type=float, default=2_000.0)
    parser.add_argument("--min-excess-db", type=float, default=3.0)
    parser.add_argument("--min-prominence-db", type=float, default=2.0)
    parser.add_argument("--min-spacing-hz", type=float, default=500.0)
    parser.add_argument("--max-peaks-per-source", type=int, default=12)
    parser.add_argument("--cluster-tolerance-hz", type=float, default=300.0)
    parser.add_argument("--search-half-width-hz", type=float, default=500.0)
    parser.add_argument(
        "--explicit-search-half-width-hz",
        type=float,
        default=150.0,
    )
    parser.add_argument("--material-difference-db", type=float, default=3.0)
    parser.add_argument("--stored-match-db", type=float, default=1.5)
    parser.add_argument("--line-present-db", type=float, default=3.0)
    parser.add_argument("--smooth-width-hz", type=float, default=1000.0)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    result = run(args)
    write_plot(result, args.figure, args)
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
                "repeat_records": result["comparison"]["repeat"][
                    "accepted_record_count"
                ],
                "reference_records": result["comparison"]["reference"][
                    "accepted_record_count"
                ],
                "stored_modelnoise_available": result[
                    "stored_modelnoise"
                ]["available"],
                "line_summary": [
                    {
                        "anchor_Hz": row["cluster"]["anchor_Hz"],
                        "day_delta_dB": row[
                            "reference_minus_repeat_line_excess_dB"
                        ],
                        "stored_delta_dB": row[
                            "stored_minus_reference_reconstruction_line_excess_dB"
                        ],
                        "day_frequency_shift_Hz": row[
                            "reference_minus_repeat_peak_frequency_Hz"
                        ],
                        "classification": row["classification"],
                    }
                    for row in result["line_clusters"]
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
