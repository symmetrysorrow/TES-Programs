"""Direct cross-day stability test for the fitted Magnicon+c2 readout transfer.

This diagnostic deliberately removes the escape route present in the preceding
cross-application test.  It does not refit detector nuisance parameters or the
white floor while comparing the two readout shapes.  Instead it takes each
day's locally fitted Magnicon Bessel + c2 readout coordinates and compares the
normalized transfer magnitudes directly.

For c4=0 the readout magnitude is

    |H(f)| = sqrt(1 + c2 * (f/f_ref)^2)
             / sqrt((1-(f/fp)^2)^2 + (f/(Q*fp))^2).

The total cross-day difference is decomposed into:
  * Bessel denominator contribution (c2=0), and
  * residual numerator contribution sqrt(1+c2*x^2).

This is a repeatability diagnostic only.  It does not identify c2 with a
physical zero, nor does failure prove electronics drift.
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
DEFAULT_CONFIG = (
    CONFIG_DIR / "magnicon_xxf1_lpf_c2_transfer_shape_stability_config.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c2_transfer_shape_stability_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_lpf_c2_transfer_shape_stability_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import readout_effective_numerator_diagnostic as effective  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _normalized(values, frequency, reference_hz):
    return opt.normalize_at(
        np.asarray(frequency, dtype=float),
        np.asarray(values, dtype=float),
        reference_hz=float(reference_hz),
    )


def readout_components(frequency, readout, scale_hz, reference_hz):
    frequency = np.asarray(frequency, dtype=float)
    c2 = float(readout["c2"])
    if c2 < 0.0:
        raise ValueError("c2 must be non-negative")
    latent = {"v": float(np.sqrt(c2))}
    total = effective.effective_transfer_magnitude(
        frequency,
        float(readout["pole_Hz"]),
        float(readout["pole_Q"]),
        1,
        latent,
        float(scale_hz),
    )
    denominator_only = effective.effective_transfer_magnitude(
        frequency,
        float(readout["pole_Hz"]),
        float(readout["pole_Q"]),
        0,
        {},
        float(scale_hz),
    )
    x = frequency / float(scale_hz)
    numerator_only = np.sqrt(1.0 + c2 * x**2)

    return {
        "total": _normalized(total, frequency, reference_hz),
        "denominator_only": _normalized(
            denominator_only, frequency, reference_hz
        ),
        "numerator_only": _normalized(
            numerator_only, frequency, reference_hz
        ),
    }


def difference_db(repeat, reference):
    repeat = np.asarray(repeat, dtype=float)
    reference = np.asarray(reference, dtype=float)
    return 20.0 * np.log10(repeat / reference)


def band_metrics(frequency, values_db, bands):
    frequency = np.asarray(frequency, dtype=float)
    values_db = np.asarray(values_db, dtype=float)
    result = {}
    for band in bands:
        low = float(band["min"])
        high = float(band["max"])
        mask = (frequency >= low) & (frequency <= high)
        if not np.any(mask):
            continue
        values = values_db[mask]
        result[str(band["name"])] = {
            "mean_difference_dB": float(np.mean(values)),
            "rms_difference_dB": float(
                np.sqrt(np.mean(values**2))
            ),
            "max_abs_difference_dB": float(
                np.max(np.abs(values))
            ),
        }
    return result


def anchor_metrics(frequency, total_db, denominator_db, numerator_db, anchors):
    frequency = np.asarray(frequency, dtype=float)
    rows = []
    for anchor in anchors:
        index = int(np.argmin(np.abs(frequency - float(anchor))))
        rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "total_repeat_over_reference_dB": float(total_db[index]),
                "bessel_denominator_repeat_over_reference_dB": float(
                    denominator_db[index]
                ),
                "numerator_repeat_over_reference_dB": float(
                    numerator_db[index]
                ),
            }
        )
    return rows


def global_metrics(values_db):
    values_db = np.asarray(values_db, dtype=float)
    return {
        "mean_difference_dB": float(np.mean(values_db)),
        "rms_difference_dB": float(
            np.sqrt(np.mean(values_db**2))
        ),
        "max_abs_difference_dB": float(
            np.max(np.abs(values_db))
        ),
    }


def classify(total_metrics, bands, screen, cross_application_passes):
    max_band_rms = max(
        float(row["rms_difference_dB"])
        for row in bands.values()
    )
    passes = bool(
        float(total_metrics["rms_difference_dB"])
        <= float(screen["max_full_band_rms_difference_dB"])
        and max_band_rms
        <= float(screen["max_any_band_rms_difference_dB"])
        and float(total_metrics["max_abs_difference_dB"])
        <= float(screen["max_abs_difference_dB"])
    )
    if passes:
        classification = "direct_readout_transfer_shape_repeatable"
    elif cross_application_passes:
        classification = (
            "direct_readout_transfer_differs_but_detector_refit_absorbs_difference"
        )
    else:
        classification = "direct_readout_transfer_not_repeatable"
    return {
        "classification": classification,
        "direct_transfer_shape_passes_screen": passes,
        "bidirectional_cross_application_passes": bool(
            cross_application_passes
        ),
        "max_band_rms_difference_dB": float(max_band_rms),
        "stability_screen": screen,
    }


def _local_readouts(base_result):
    # In the cross-application payload, each source_readout is the source day's
    # local Magnicon+c2 solution frozen before applying it to the other day.
    repeat_readout = base_result["cross_application"][
        "repeat_readout_on_reference"
    ]["source_readout"]
    reference_readout = base_result["cross_application"][
        "reference_readout_on_repeat"
    ]["source_readout"]
    return reference_readout, repeat_readout


def run(config: dict, config_path: Path):
    base_config_path = resolve_config_path(
        config["base_cross_day_config"],
        config_path,
    )
    base_config = json.loads(
        base_config_path.read_text(encoding="utf-8")
    )
    base_result = crossday.run(base_config, base_config_path)

    reference_readout, repeat_readout = _local_readouts(base_result)
    scale_hz = float(
        base_result["reference_day"]["readout_reference_scale_Hz"]
    )
    repeat_scale_hz = float(
        base_result["repeat_day"]["readout_reference_scale_Hz"]
    )
    if not np.isclose(scale_hz, repeat_scale_hz, rtol=0.0, atol=0.0):
        raise ValueError("readout reference scales differ across days")

    p = base_result["_plot"]
    frequency = np.asarray(p["frequency_Hz"], dtype=float)
    reference_hz = float(config["normalization_reference_Hz"])

    reference_components = readout_components(
        frequency,
        reference_readout,
        scale_hz,
        reference_hz,
    )
    repeat_components = readout_components(
        frequency,
        repeat_readout,
        scale_hz,
        reference_hz,
    )

    total_db = difference_db(
        repeat_components["total"],
        reference_components["total"],
    )
    denominator_db = difference_db(
        repeat_components["denominator_only"],
        reference_components["denominator_only"],
    )
    numerator_db = difference_db(
        repeat_components["numerator_only"],
        reference_components["numerator_only"],
    )

    reconstruction_error = total_db - (
        denominator_db + numerator_db
    )
    max_reconstruction_error = float(
        np.max(np.abs(reconstruction_error))
    )
    if max_reconstruction_error > 1e-10:
        raise RuntimeError(
            "transfer decomposition failed: total dB difference does not "
            "equal denominator plus numerator contributions"
        )

    full = global_metrics(total_db)
    denominator_full = global_metrics(denominator_db)
    numerator_full = global_metrics(numerator_db)
    bands = band_metrics(
        frequency,
        total_db,
        config["bands_Hz"],
    )
    denominator_bands = band_metrics(
        frequency,
        denominator_db,
        config["bands_Hz"],
    )
    numerator_bands = band_metrics(
        frequency,
        numerator_db,
        config["bands_Hz"],
    )

    cross_application_passes = bool(
        base_result["interpretation"][
            "bidirectional_cross_application_passes"
        ]
    )
    interpretation = classify(
        full,
        bands,
        config["stability_screen"],
        cross_application_passes,
    )

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "tested_question": (
            "Are the locally fitted Magnicon+c2 readout transfer magnitudes "
            "themselves stable across 2024-12-05 and 2024-12-06, before "
            "detector nuisance or white-floor refits can absorb differences?"
        ),
        "normalization_reference_Hz": reference_hz,
        "readout_reference_scale_Hz": scale_hz,
        "transfer_definition": (
            "|H|=sqrt(1+c2*(f/f_ref)^2) / "
            "sqrt((1-(f/fp)^2)^2 + (f/(Q*fp))^2), with c4=0"
        ),
        "reference_day": {
            "label": base_result["reference_day"]["label"],
            "readout": {
                key: float(value)
                for key, value in reference_readout.items()
            },
            "profiled_cutoff_Hz": float(
                base_result["reference_day"]["profiled_cutoff_Hz"]
            ),
            "profiled_c2": float(
                base_result["reference_day"]["profiled_c2"]
            ),
        },
        "repeat_day": {
            "label": base_result["repeat_day"]["label"],
            "readout": {
                key: float(value)
                for key, value in repeat_readout.items()
            },
            "profiled_cutoff_Hz": float(
                base_result["repeat_day"]["profiled_cutoff_Hz"]
            ),
            "profiled_c2": float(
                base_result["repeat_day"]["profiled_c2"]
            ),
        },
        "direct_transfer_difference": {
            "semantics": "repeat_day / reference_day, each normalized at 1 kHz",
            "full_1_200k": full,
            "bands": bands,
            "anchors": anchor_metrics(
                frequency,
                total_db,
                denominator_db,
                numerator_db,
                config["anchors_Hz"],
            ),
        },
        "decomposition": {
            "bessel_denominator_contribution": {
                "full_1_200k": denominator_full,
                "bands": denominator_bands,
            },
            "residual_numerator_contribution": {
                "full_1_200k": numerator_full,
                "bands": numerator_bands,
            },
            "max_abs_reconstruction_error_dB": max_reconstruction_error,
        },
        "preceding_cross_application": {
            "classification": base_result["interpretation"][
                "classification"
            ],
            "bidirectional_cross_application_passes": (
                cross_application_passes
            ),
            "repeat_readout_on_reference": base_result[
                "cross_application"
            ]["repeat_readout_on_reference"][
                "comparison_to_reference_local_best"
            ],
            "reference_readout_on_repeat": base_result[
                "cross_application"
            ]["reference_readout_on_repeat"][
                "comparison_to_repeat_local_best"
            ],
        },
        "interpretation": {
            **interpretation,
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_cross_day_config": str(base_config_path),
        },
        "_plot": {
            "frequency_Hz": frequency.tolist(),
            "reference_total": reference_components["total"].tolist(),
            "repeat_total": repeat_components["total"].tolist(),
            "total_difference_dB": total_db.tolist(),
            "denominator_difference_dB": denominator_db.tolist(),
            "numerator_difference_dB": numerator_db.tolist(),
        },
    }
    return result


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    frequency = np.asarray(p["frequency_Hz"], dtype=float)

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.0), sharex=True)
    top, middle, bottom = axes

    top.loglog(
        frequency,
        np.asarray(p["reference_total"], dtype=float),
        label="12/06 readout transfer",
    )
    top.loglog(
        frequency,
        np.asarray(p["repeat_total"], dtype=float),
        label="12/05 readout transfer",
    )
    top.set_ylabel("Normalized |H|")
    top.legend(frameon=False)
    top.grid(True, which="both", alpha=0.2)

    middle.semilogx(
        frequency,
        np.asarray(p["total_difference_dB"], dtype=float),
        label="total",
    )
    middle.semilogx(
        frequency,
        np.asarray(p["denominator_difference_dB"], dtype=float),
        label="Bessel denominator",
    )
    middle.semilogx(
        frequency,
        np.asarray(p["numerator_difference_dB"], dtype=float),
        label="c2 numerator",
    )
    middle.axhline(0.0, linewidth=1.0)
    middle.set_ylabel("12/05 / 12/06 [dB]")
    middle.legend(frameon=False)
    middle.grid(True, which="both", alpha=0.2)

    names = list(result["direct_transfer_difference"]["bands"])
    rms = [
        result["direct_transfer_difference"]["bands"][name][
            "rms_difference_dB"
        ]
        for name in names
    ]
    x = np.arange(len(names), dtype=float)
    bottom.plot(x, rms, marker="o")
    bottom.set_xticks(x, names, rotation=20)
    bottom.set_ylabel("Band RMS difference [dB]")
    bottom.set_xlabel("Frequency band")
    bottom.grid(True, alpha=0.2)

    fig.suptitle(result["interpretation"]["classification"], fontsize=10)
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
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(config, args.config)
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
                "classification": result["interpretation"][
                    "classification"
                ],
                "full_rms_difference_dB": result[
                    "direct_transfer_difference"
                ]["full_1_200k"]["rms_difference_dB"],
                "max_abs_difference_dB": result[
                    "direct_transfer_difference"
                ]["full_1_200k"]["max_abs_difference_dB"],
                "max_band_rms_difference_dB": result[
                    "interpretation"
                ]["max_band_rms_difference_dB"],
                "bidirectional_cross_application_passes": result[
                    "preceding_cross_application"
                ]["bidirectional_cross_application_passes"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
