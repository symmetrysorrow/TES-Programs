"""Model-free adjacent-day pre-analysis ASD drift and block-repeatability diagnostic.

The empirical comparison uses only accepted raw CH0 records from the reference
and repeat runs.  Detector/noise-model parameters are not used to construct the
measured day-to-day ASD ratio.

A separately tracked diagnostic transfer ratio (repeat-local / reference-shared)
is overlaid only after the empirical ratio is formed.  The repeat accepted
records are also split into deterministic non-overlapping blocks with the same
record count as the reference run, so day-to-day shape drift can be compared
against within-run sampling/stationarity variation.
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
DEFAULT_CONFIG = CONFIG_DIR / "readout_day_to_day_drift_config.json"
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_day_to_day_drift_diagnostic.json"
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    estimate_one_sided_asd,
)
from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import preanalysis_hybrid_pole_profile_diagnostic as profile  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def case_by_label(cases, label):
    matches = [case for case in cases if case["label"] == label]
    if len(matches) != 1:
        raise ValueError(
            f"expected exactly one manifest case named {label!r}; "
            f"found {len(matches)}"
        )
    return matches[0]


def comparison_for_case(case):
    if case["comparison_summary"] is not None:
        comparison = json.loads(
            case["comparison_summary"].read_text(encoding="utf-8")
        )
        source = "tracked_comparison_summary"
    else:
        comparison = shared.build_comparison_from_spec(
            case["comparison_spec"]
        )
        source = "tracked_comparison_spec_recomputed_mask"
    experiment_path = (
        case["experiment_path"]
        if case["experiment_path"] is not None
        else Path(comparison["experiment_path"])
    )
    return comparison, experiment_path, source


def pre_analysis_asd(
    comparison,
    experiment_path: Path,
    accepted_indices=None,
):
    acquisition = comparison["acquisition"]
    rate = float(acquisition["rate_Hz"])
    samples = int(acquisition["samples"])
    channel = str(acquisition.get("channel", "CH0"))
    raw_dir = experiment_path / f"{channel}_noise" / "rawdata"
    paths = sorted(raw_dir.glob(f"{channel}_*.dat"))
    if not paths:
        raise FileNotFoundError(
            f"no {channel} raw records found in {raw_dir}"
        )

    if accepted_indices is None:
        accepted_indices = [
            int(value)
            for value in acquisition["accepted_record_indices"]
        ]
    else:
        accepted_indices = [int(value) for value in accepted_indices]

    if not accepted_indices:
        raise ValueError("accepted record index list is empty")
    if max(accepted_indices) >= len(paths):
        raise IndexError("accepted record index exceeds raw record list")

    records = [
        base.read_record(paths[index])
        for index in accepted_indices
    ]
    asd, count = estimate_one_sided_asd(
        records,
        samples,
        rate,
        cutoff=0.0,
        remove_mean=True,
    )
    if count != len(accepted_indices):
        raise RuntimeError(
            "pre-analysis estimator changed the supplied record mask"
        )
    frequency = np.fft.rfftfreq(samples, d=1.0 / rate)
    return {
        "frequency_Hz": frequency,
        "asd": np.asarray(asd, dtype=float),
        "normalized": opt.normalize_at(
            frequency,
            asd,
            reference_hz=1000.0,
        ),
        "accepted_indices": accepted_indices,
        "accepted_records": int(count),
        "rate_Hz": rate,
        "samples": samples,
    }


def require_matching_frequency(reference, repeat):
    if reference["samples"] != repeat["samples"]:
        raise ValueError("reference/repeat sample counts differ")
    if not np.isclose(
        reference["rate_Hz"],
        repeat["rate_Hz"],
        rtol=0.0,
        atol=max(reference["rate_Hz"], 1.0) * 1e-12,
    ):
        raise ValueError("reference/repeat sample rates differ")
    np.testing.assert_allclose(
        reference["frequency_Hz"],
        repeat["frequency_Hz"],
        rtol=0.0,
        atol=max(reference["rate_Hz"], 1.0) * 1e-12,
    )


def normalized_transfer(
    frequency,
    parameters,
    fixed_pole_hz,
    reference_hz,
):
    magnitude = profile.hybrid_transfer_magnitude(
        frequency,
        parameters,
        fixed_pole_hz,
    )
    return opt.normalize_at(
        frequency,
        magnitude,
        reference_hz=reference_hz,
    )


def load_repeat_transfer(path: Path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    transfer = payload["transfer"]
    params = transfer["parameters"]
    for name in shared.TRANSFER_PARAMETER_NAMES:
        if float(params[name]) <= 0.0:
            raise ValueError(f"repeat transfer {name} must be positive")
    fixed = transfer.get("fixed_leadlag_pole_Hz")
    if transfer.get("leadlag_pole_is_infinite", fixed is None):
        fixed = None
    elif fixed is None:
        raise ValueError("finite repeat lead/lag pole requires a value")
    else:
        fixed = float(fixed)
    return {
        "parameters": {
            name: float(params[name])
            for name in shared.TRANSFER_PARAMETER_NAMES
        },
        "fixed_leadlag_pole_Hz": fixed,
        "source": payload.get("provenance", {}),
    }


def make_blocks(indices, block_size):
    block_size = int(block_size)
    if block_size <= 0:
        raise ValueError("block size must be positive")
    indices = list(indices)
    count = len(indices) // block_size
    return [
        indices[i * block_size : (i + 1) * block_size]
        for i in range(count)
    ], len(indices) - count * block_size


def interp_at(frequency, values, anchors):
    return [
        {
            "frequency_Hz": float(anchor),
            "value": float(np.interp(anchor, frequency, values)),
        }
        for anchor in anchors
    ]


def band_stats(frequency, empirical_db, transfer_db, bands):
    rows = {}
    residual = empirical_db - transfer_db
    for band in bands:
        name = str(band["name"])
        low = float(band["min"])
        high = float(band["max"])
        mask = (frequency >= low) & (frequency <= high)
        if not np.any(mask):
            continue
        rows[name] = {
            "empirical_day_ratio_mean_dB": float(
                np.mean(empirical_db[mask])
            ),
            "empirical_day_ratio_rms_dB": float(
                np.sqrt(np.mean(empirical_db[mask] ** 2))
            ),
            "transfer_ratio_mean_dB": float(
                np.mean(transfer_db[mask])
            ),
            "transfer_ratio_rms_dB": float(
                np.sqrt(np.mean(transfer_db[mask] ** 2))
            ),
            "empirical_minus_transfer_mean_dB": float(
                np.mean(residual[mask])
            ),
            "empirical_minus_transfer_rms_dB": float(
                np.sqrt(np.mean(residual[mask] ** 2))
            ),
        }
    return rows


def block_envelope(
    frequency,
    repeat_full_normalized,
    block_normalized,
    anchors,
    fit_mask,
):
    block_db = np.asarray(
        [
            20.0
            * np.log10(
                np.asarray(values, dtype=float)
                / repeat_full_normalized
            )
            for values in block_normalized
        ],
        dtype=float,
    )
    anchor_rows = []
    for anchor in anchors:
        index = int(np.argmin(np.abs(frequency - float(anchor))))
        values = block_db[:, index]
        anchor_rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "min_dB": float(np.min(values)),
                "p16_dB": float(np.percentile(values, 16.0)),
                "median_dB": float(np.median(values)),
                "p84_dB": float(np.percentile(values, 84.0)),
                "max_dB": float(np.max(values)),
                "max_abs_dB": float(np.max(np.abs(values))),
            }
        )

    rms_rows = [
        float(np.sqrt(np.mean(values[fit_mask] ** 2)))
        for values in block_db
    ]
    return {
        "block_ratio_dB": block_db,
        "anchor_envelope": anchor_rows,
        "fit_band_rms_dB": rms_rows,
        "fit_band_rms_max_dB": float(np.max(rms_rows)),
        "fit_band_rms_median_dB": float(np.median(rms_rows)),
    }


def compare_day_to_block(
    frequency,
    empirical_day_db,
    anchor_envelope,
    anchors,
    fit_mask,
    block_rms_max,
):
    rows = []
    exceeds = 0
    nonzero_anchors = 0
    for anchor, envelope in zip(anchors, anchor_envelope):
        index = int(np.argmin(np.abs(frequency - float(anchor))))
        day_value = float(empirical_day_db[index])
        max_abs = float(envelope["max_abs_dB"])
        ratio = (
            abs(day_value) / max_abs
            if max_abs > 0.0
            else (float("inf") if abs(day_value) > 0.0 else 0.0)
        )
        is_reference = np.isclose(
            float(anchor),
            1000.0,
            rtol=0.0,
            atol=1e-9,
        )
        exceeds_here = bool(abs(day_value) > max_abs)
        if not is_reference:
            nonzero_anchors += 1
            exceeds += int(exceeds_here)
        rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "empirical_day_ratio_dB": day_value,
                "within_repeat_max_abs_block_ratio_dB": max_abs,
                "day_over_block_max_abs_ratio": float(ratio),
                "day_exceeds_all_block_variation": exceeds_here,
            }
        )

    day_rms = float(
        np.sqrt(np.mean(empirical_day_db[fit_mask] ** 2))
    )
    return {
        "anchors": rows,
        "n_nonreference_anchors": int(nonzero_anchors),
        "n_anchors_day_exceeds_all_blocks": int(exceeds),
        "fraction_anchors_day_exceeds_all_blocks": float(
            exceeds / nonzero_anchors
            if nonzero_anchors
            else 0.0
        ),
        "empirical_day_fit_band_rms_dB": day_rms,
        "max_within_repeat_block_fit_band_rms_dB": float(
            block_rms_max
        ),
        "day_rms_over_block_rms_max_ratio": float(
            day_rms / block_rms_max
            if block_rms_max > 0.0
            else float("inf")
        ),
        "day_fit_band_rms_exceeds_all_blocks": bool(
            day_rms > block_rms_max
        ),
    }


def run(config, config_path: Path):
    manifest_path = resolve_config_path(
        config["manifest"],
        config_path,
    )
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    cases = shared.normalize_manifest(manifest, manifest_path)
    reference_case = case_by_label(
        cases,
        config["reference_case_label"],
    )
    repeat_case = case_by_label(
        cases,
        config["repeat_case_label"],
    )

    reference_comparison, reference_path, reference_source = (
        comparison_for_case(reference_case)
    )
    repeat_comparison, repeat_path, repeat_source = (
        comparison_for_case(repeat_case)
    )
    reference = pre_analysis_asd(
        reference_comparison,
        reference_path,
    )
    repeat = pre_analysis_asd(
        repeat_comparison,
        repeat_path,
    )
    require_matching_frequency(reference, repeat)

    comparison_cfg = config["comparison"]
    reference_hz = float(
        comparison_cfg["normalization_reference_Hz"]
    )
    frequency = reference["frequency_Hz"]
    reference_norm = opt.normalize_at(
        frequency,
        reference["asd"],
        reference_hz=reference_hz,
    )
    repeat_norm = opt.normalize_at(
        frequency,
        repeat["asd"],
        reference_hz=reference_hz,
    )
    empirical_ratio = repeat_norm / reference_norm
    empirical_db = 20.0 * np.log10(empirical_ratio)

    reference_transfer_path = resolve_config_path(
        config["reference_transfer"],
        config_path,
    )
    reference_transfer_payload = json.loads(
        reference_transfer_path.read_text(encoding="utf-8")
    )
    reference_transfer = shared.load_shared_transfer(
        reference_transfer_payload
    )
    repeat_transfer_path = resolve_config_path(
        config["repeat_local_transfer"],
        config_path,
    )
    repeat_transfer = load_repeat_transfer(
        repeat_transfer_path
    )
    reference_h = normalized_transfer(
        frequency,
        reference_transfer["parameters"],
        reference_transfer["fixed_leadlag_pole_Hz"],
        reference_hz,
    )
    repeat_h = normalized_transfer(
        frequency,
        repeat_transfer["parameters"],
        repeat_transfer["fixed_leadlag_pole_Hz"],
        reference_hz,
    )
    transfer_ratio = repeat_h / reference_h
    transfer_db = 20.0 * np.log10(transfer_ratio)
    unexplained_db = empirical_db - transfer_db

    fit_min = float(comparison_cfg["fit_min_Hz"])
    fit_max = float(comparison_cfg["fit_max_Hz"])
    fit_mask = (frequency >= fit_min) & (frequency <= fit_max)
    anchors = [
        float(value) for value in comparison_cfg["anchors_Hz"]
    ]

    block_cfg = config["block_size"]
    if block_cfg["mode"] == "reference_accepted_record_count":
        block_size = int(reference["accepted_records"])
    else:
        block_size = int(block_cfg["fallback_records"])
    blocks, remainder = make_blocks(
        repeat["accepted_indices"],
        block_size,
    )
    if len(blocks) < 2:
        raise ValueError(
            "repeat run does not contain at least two full blocks"
        )
    block_normalized = []
    for indices in blocks:
        block = pre_analysis_asd(
            repeat_comparison,
            repeat_path,
            accepted_indices=indices,
        )
        block_normalized.append(
            opt.normalize_at(
                frequency,
                block["asd"],
                reference_hz=reference_hz,
            )
        )
    envelope = block_envelope(
        frequency,
        repeat_norm,
        block_normalized,
        anchors,
        fit_mask,
    )
    day_vs_block = compare_day_to_block(
        frequency,
        empirical_db,
        envelope["anchor_envelope"],
        anchors,
        fit_mask,
        envelope["fit_band_rms_max_dB"],
    )

    corr = float(
        np.corrcoef(
            empirical_db[fit_mask],
            transfer_db[fit_mask],
        )[0, 1]
    )
    model_mismatch_rms = float(
        np.sqrt(np.mean(unexplained_db[fit_mask] ** 2))
    )

    anchor_rows = []
    for anchor in anchors:
        index = int(np.argmin(np.abs(frequency - anchor)))
        anchor_rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "empirical_repeat_over_reference_dB": float(
                    empirical_db[index]
                ),
                "repeat_local_over_reference_shared_transfer_dB": float(
                    transfer_db[index]
                ),
                "empirical_minus_transfer_dB": float(
                    unexplained_db[index]
                ),
            }
        )

    return {
        "diagnostic_only": True,
        "production_optimizer_unchanged": True,
        "production_noise_model_unchanged": True,
        "empirical_day_ratio_is_model_free": True,
        "digital_analysis_filter_applied_to_empirical_asd": False,
        "reference_case": {
            "label": reference_case["label"],
            "comparison_source": reference_source,
            "experiment_path": str(reference_path),
            "accepted_records": reference["accepted_records"],
        },
        "repeat_case": {
            "label": repeat_case["label"],
            "comparison_source": repeat_source,
            "experiment_path": str(repeat_path),
            "accepted_records": repeat["accepted_records"],
        },
        "normalization_reference_Hz": reference_hz,
        "empirical_semantics": (
            "each run: exact accepted raw CH0 mask -> mean removal -> "
            "Hann power average -> one-sided pre-analysis ASD -> independent "
            "normalization at 1 kHz; empirical ratio = repeat/reference"
        ),
        "transfer_comparison_semantics": (
            "diagnostic-only overlay: tracked repeat-local infinity-pole "
            "hybrid magnitude divided by tracked reference-shared hybrid "
            "magnitude, each independently normalized at 1 kHz"
        ),
        "block_repeatability": {
            "block_size_records": block_size,
            "n_full_nonoverlapping_blocks": int(len(blocks)),
            "unused_remainder_records": int(remainder),
            "anchor_envelope": envelope["anchor_envelope"],
            "fit_band_rms_dB_per_block": envelope[
                "fit_band_rms_dB"
            ],
            "fit_band_rms_max_dB": envelope[
                "fit_band_rms_max_dB"
            ],
            "fit_band_rms_median_dB": envelope[
                "fit_band_rms_median_dB"
            ],
        },
        "day_vs_within_repeat": day_vs_block,
        "empirical_vs_transfer_shape": {
            "fit_min_Hz": fit_min,
            "fit_max_Hz": fit_max,
            "correlation_coefficient_dB_shape": corr,
            "empirical_minus_transfer_rms_dB": model_mismatch_rms,
            "bands": band_stats(
                frequency,
                empirical_db,
                transfer_db,
                comparison_cfg["bands_Hz"],
            ),
            "anchors": anchor_rows,
        },
        "interpretation_flags": {
            "day_shape_exceeds_within_repeat_block_variation": bool(
                day_vs_block[
                    "day_fit_band_rms_exceeds_all_blocks"
                ]
            ),
            "transfer_ratio_tracks_day_shape_broadly": bool(
                corr >= 0.8 and model_mismatch_rms <= 1.0
            ),
            "strong_electronics_drift_claim_allowed": False,
        },
        "guardrail": config["guardrail"],
        "inputs": {
            "config": str(config_path),
            "manifest": str(manifest_path),
            "reference_transfer": str(reference_transfer_path),
            "repeat_local_transfer": str(repeat_transfer_path),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    config = json.loads(
        args.config.read_text(encoding="utf-8")
    )
    result = run(config, args.config)
    output = args.output or DEFAULT_OUTPUT
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(output),
                "reference_records": result["reference_case"][
                    "accepted_records"
                ],
                "repeat_records": result["repeat_case"][
                    "accepted_records"
                ],
                "blocks": result["block_repeatability"][
                    "n_full_nonoverlapping_blocks"
                ],
                "day_rms_over_block_rms_max": result[
                    "day_vs_within_repeat"
                ]["day_rms_over_block_rms_max_ratio"],
                "empirical_transfer_shape_correlation": result[
                    "empirical_vs_transfer_shape"
                ]["correlation_coefficient_dB_shape"],
                "empirical_minus_transfer_rms_dB": result[
                    "empirical_vs_transfer_shape"
                ]["empirical_minus_transfer_rms_dB"],
                "day_exceeds_within_repeat": result[
                    "interpretation_flags"
                ][
                    "day_shape_exceeds_within_repeat_block_variation"
                ],
                "transfer_tracks_day_shape": result[
                    "interpretation_flags"
                ]["transfer_ratio_tracks_day_shape_broadly"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
