"""Quantify C_tes-anchor sensitivity of the fitted Magnicon+c2 transfer.

The physical-L+C_tes anchor diagnostic showed that c2 remains material and
cross-day repeatable for every fixed C_tes row, while the fitted c2 coordinate
moves strongly with the assumed C_tes.  This diagnostic asks whether that large
coordinate motion corresponds to an equally large motion of the actual
normalized readout transfer.

It reuses the preceding anchored fits, reconstructs the Magnicon+c2 transfer
for every C_tes multiplier, normalizes each transfer at 1 kHz, and compares:
  * anchor-to-anchor transfer shapes within 2024-12-06;
  * anchor-to-anchor transfer shapes within 2024-12-05; and
  * the geometric-mean transfer of the two days for each anchor.

Rows where c2 hits a configured bound are reported as censored and are excluded
from the primary finite-shape classification.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = (
    CONFIG_DIR
    / "magnicon_xxf1_Ctes_anchor_transfer_shape_sensitivity_config.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_Ctes_anchor_transfer_shape_sensitivity_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_Ctes_anchor_transfer_shape_sensitivity_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import magnicon_xxf1_lpf_c2_transfer_shape_stability_diagnostic as transferdiag  # noqa: E402
from subScript import magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic as anchored  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def _anchor_values(frequency, values_db, anchors):
    frequency = np.asarray(frequency, dtype=float)
    values_db = np.asarray(values_db, dtype=float)
    rows = []
    for anchor in anchors:
        index = int(np.argmin(np.abs(frequency - float(anchor))))
        rows.append(
            {
                "frequency_Hz": float(frequency[index]),
                "difference_dB": float(values_db[index]),
            }
        )
    return rows


def difference_summary(
    frequency,
    candidate,
    reference,
    bands,
    anchors,
):
    difference = transferdiag.difference_db(candidate, reference)
    return {
        "full_1_200k": transferdiag.global_metrics(difference),
        "bands": transferdiag.band_metrics(
            frequency,
            difference,
            bands,
        ),
        "anchors": _anchor_values(
            frequency,
            difference,
            anchors,
        ),
        "_difference_dB": np.asarray(difference, dtype=float),
    }


def geometric_mean_transfer(reference_day, repeat_day):
    reference_day = np.asarray(reference_day, dtype=float)
    repeat_day = np.asarray(repeat_day, dtype=float)
    if reference_day.shape != repeat_day.shape:
        raise ValueError("cross-day transfer arrays differ in shape")
    if np.any(reference_day <= 0.0) or np.any(repeat_day <= 0.0):
        raise ValueError("transfer values must be positive")
    return np.sqrt(reference_day * repeat_day)


def c2_is_censored(day_row: dict) -> bool:
    hits = day_row["c2_free"].get("readout_boundary_hits", {}).get(
        "c2", {}
    )
    return bool(hits.get("at_lower", False) or hits.get("at_upper", False))


def spread_summary(transfers: dict):
    if not transfers:
        raise ValueError("at least one transfer is required")
    stack = np.vstack(
        [np.asarray(value, dtype=float) for value in transfers.values()]
    )
    if np.any(stack <= 0.0):
        raise ValueError("transfer values must be positive")
    db = 20.0 * np.log10(stack)
    spread = np.max(db, axis=0) - np.min(db, axis=0)
    return {
        "rms_peak_to_peak_spread_dB": float(
            np.sqrt(np.mean(spread**2))
        ),
        "max_peak_to_peak_spread_dB": float(np.max(spread)),
        "_spread_dB": spread,
    }


def classify_pairwise(pairwise: dict, screen: dict) -> dict:
    if not pairwise:
        return {
            "classification": "insufficient_uncensored_anchor_pairs",
            "passes_anchor_sensitivity_screen": None,
            "max_pairwise_full_band_rms_difference_dB": None,
            "max_pairwise_abs_difference_dB": None,
        }

    max_rms = max(
        float(row["full_1_200k"]["rms_difference_dB"])
        for row in pairwise.values()
    )
    max_abs = max(
        float(row["full_1_200k"]["max_abs_difference_dB"])
        for row in pairwise.values()
    )
    passes = bool(
        max_rms
        <= float(screen["max_pairwise_full_band_rms_difference_dB"])
        and max_abs
        <= float(screen["max_pairwise_abs_difference_dB"])
    )
    return {
        "classification": (
            "total_readout_transfer_stable_across_Ctes_anchors"
            if passes
            else "required_total_readout_transfer_is_Ctes_anchor_sensitive"
        ),
        "passes_anchor_sensitivity_screen": passes,
        "max_pairwise_full_band_rms_difference_dB": float(max_rms),
        "max_pairwise_abs_difference_dB": float(max_abs),
    }


def _readout_for_day(base_result: dict, row: dict, day_key: str):
    fixed = base_result["magnicon_fixed_model"]
    return {
        "pole_Hz": float(fixed["pole_Hz"]),
        "pole_Q": float(fixed["pole_Q"]),
        "c2": float(row[day_key]["c2_free"]["c2"]),
        "c4": float(fixed["c4"]),
    }


def _pairwise_summaries(
    *,
    frequency,
    transfers,
    included_keys,
    bands,
    anchors,
):
    rows = {}
    for first, second in itertools.combinations(included_keys, 2):
        key = f"{first}_vs_{second}"
        rows[key] = difference_summary(
            frequency,
            transfers[second],
            transfers[first],
            bands,
            anchors,
        )
        rows[key]["semantics"] = f"{second} / {first}"
    return rows


def _clean_difference(row: dict):
    return {
        key: value
        for key, value in row.items()
        if key != "_difference_dB"
    }


def run(config: dict, config_path: Path):
    base_config_path = resolve_config_path(
        config["base_anchor_config"],
        config_path,
    )
    base_config = json.loads(
        base_config_path.read_text(encoding="utf-8")
    )
    base_result = anchored.run(base_config, base_config_path)

    frequency = np.asarray(
        base_result["_plot"]["frequency_Hz"],
        dtype=float,
    )
    scale_hz = 40000.0
    normalization_hz = float(
        config["normalization_reference_Hz"]
    )
    row_items = list(base_result["sensitivity_rows"].items())

    reference_transfers = {}
    repeat_transfers = {}
    mean_transfers = {}
    anchor_rows = {}
    censored = []

    for key, row in row_items:
        reference_readout = _readout_for_day(
            base_result, row, "reference_day"
        )
        repeat_readout = _readout_for_day(
            base_result, row, "repeat_day"
        )

        reference_components = transferdiag.readout_components(
            frequency,
            reference_readout,
            scale_hz,
            normalization_hz,
        )
        repeat_components = transferdiag.readout_components(
            frequency,
            repeat_readout,
            scale_hz,
            normalization_hz,
        )
        reference_transfer = reference_components["total"]
        repeat_transfer = repeat_components["total"]

        reference_transfers[key] = reference_transfer
        repeat_transfers[key] = repeat_transfer
        mean_transfers[key] = geometric_mean_transfer(
            reference_transfer,
            repeat_transfer,
        )

        reference_censored = c2_is_censored(row["reference_day"])
        repeat_censored = c2_is_censored(row["repeat_day"])
        row_censored = bool(reference_censored or repeat_censored)
        if row_censored:
            censored.append(key)

        anchor_rows[key] = {
            "C_tes_multiplier": float(row["C_tes_multiplier"]),
            "fixed_C_tes_J_per_K": float(
                row["fixed_C_tes_J_per_K"]
            ),
            "reference_c2": float(
                row["reference_day"]["c2_free"]["c2"]
            ),
            "repeat_c2": float(
                row["repeat_day"]["c2_free"]["c2"]
            ),
            "reference_continuum_rms_dB": float(
                row["reference_day"]["c2_free"][
                    "continuum_rms_dB"
                ]
            ),
            "repeat_continuum_rms_dB": float(
                row["repeat_day"]["c2_free"][
                    "continuum_rms_dB"
                ]
            ),
            "reference_c2_censored": reference_censored,
            "repeat_c2_censored": repeat_censored,
            "row_censored": row_censored,
            "cross_day_transfer_repeatability": row[
                "cross_day_c2_free_transfer"
            ],
        }

    all_keys = [key for key, _ in row_items]
    included_keys = [
        key
        for key in all_keys
        if not (
            bool(
                config[
                    "exclude_c2_boundary_rows_from_primary_classification"
                ]
            )
            and key in censored
        )
    ]

    reference_pairwise = _pairwise_summaries(
        frequency=frequency,
        transfers=reference_transfers,
        included_keys=included_keys,
        bands=config["bands_Hz"],
        anchors=config["anchors_Hz"],
    )
    repeat_pairwise = _pairwise_summaries(
        frequency=frequency,
        transfers=repeat_transfers,
        included_keys=included_keys,
        bands=config["bands_Hz"],
        anchors=config["anchors_Hz"],
    )
    mean_pairwise = _pairwise_summaries(
        frequency=frequency,
        transfers=mean_transfers,
        included_keys=included_keys,
        bands=config["bands_Hz"],
        anchors=config["anchors_Hz"],
    )

    reference_spread = spread_summary(
        {key: reference_transfers[key] for key in included_keys}
    )
    repeat_spread = spread_summary(
        {key: repeat_transfers[key] for key in included_keys}
    )
    mean_spread = spread_summary(
        {key: mean_transfers[key] for key in included_keys}
    )

    classification = classify_pairwise(
        mean_pairwise,
        config["anchor_sensitivity_screen"],
    )

    comparison_reference_key = anchored._multiplier_key(
        float(config["comparison_reference_multiplier"])
    )
    primary_material_key = anchored._multiplier_key(
        float(config["primary_material_multiplier"])
    )
    if comparison_reference_key not in mean_transfers:
        raise ValueError("comparison reference multiplier missing")
    if primary_material_key not in mean_transfers:
        raise ValueError("primary material multiplier missing")

    primary_vs_reference = difference_summary(
        frequency,
        mean_transfers[primary_material_key],
        mean_transfers[comparison_reference_key],
        config["bands_Hz"],
        config["anchors_Hz"],
    )
    primary_vs_reference["semantics"] = (
        f"{primary_material_key} / {comparison_reference_key} "
        "using geometric-mean cross-day transfer"
    )

    best_fit_key = min(
        all_keys,
        key=lambda key: 0.5
        * (
            anchor_rows[key]["reference_continuum_rms_dB"]
            + anchor_rows[key]["repeat_continuum_rms_dB"]
        ),
    )

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "tested_question": (
            "Do the large c2 coordinate changes induced by fixed C_tes "
            "anchors correspond to materially different normalized "
            "Magnicon+c2 readout transfer shapes?"
        ),
        "normalization_reference_Hz": normalization_hz,
        "readout_reference_scale_Hz": scale_hz,
        "anchor_rows": anchor_rows,
        "primary_classification_anchor_keys": included_keys,
        "censored_anchor_keys": censored,
        "reference_day_anchor_sensitivity": {
            "pairwise": {
                key: _clean_difference(row)
                for key, row in reference_pairwise.items()
            },
            "spread": {
                key: value
                for key, value in reference_spread.items()
                if key != "_spread_dB"
            },
        },
        "repeat_day_anchor_sensitivity": {
            "pairwise": {
                key: _clean_difference(row)
                for key, row in repeat_pairwise.items()
            },
            "spread": {
                key: value
                for key, value in repeat_spread.items()
                if key != "_spread_dB"
            },
        },
        "cross_day_geometric_mean_anchor_sensitivity": {
            "semantics": (
                "For each C_tes anchor, geometric mean of the normalized "
                "12/06 and 12/05 readout transfers, then anchor-to-anchor "
                "comparison."
            ),
            "pairwise": {
                key: _clean_difference(row)
                for key, row in mean_pairwise.items()
            },
            "spread": {
                key: value
                for key, value in mean_spread.items()
                if key != "_spread_dB"
            },
        },
        "primary_material_vs_reference_anchor": _clean_difference(
            primary_vs_reference
        ),
        "fit_quality_context": {
            "lowest_mean_c2_free_rms_anchor": best_fit_key,
            "mean_c2_free_rms_dB_by_anchor": {
                key: float(
                    0.5
                    * (
                        row["reference_continuum_rms_dB"]
                        + row["repeat_continuum_rms_dB"]
                    )
                )
                for key, row in anchor_rows.items()
            },
            "semantics": (
                "Descriptive fit-quality context only; the lowest-noise-fit "
                "RMS anchor is not treated as a physical C_tes measurement."
            ),
        },
        "interpretation": {
            **classification,
            "c2_bound_censoring_present": bool(censored),
            "excluded_censored_anchor_keys": censored,
            "comparison_reference_anchor": comparison_reference_key,
            "primary_material_anchor": primary_material_key,
            "primary_material_vs_reference_full_band_rms_difference_dB": (
                float(
                    primary_vs_reference["full_1_200k"][
                        "rms_difference_dB"
                    ]
                )
            ),
            "primary_material_vs_reference_max_abs_difference_dB": (
                float(
                    primary_vs_reference["full_1_200k"][
                        "max_abs_difference_dB"
                    ]
                )
            ),
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_anchor_config": str(base_config_path),
        },
        "_plot": {
            "frequency_Hz": frequency.tolist(),
            "reference_transfers": {
                key: value.tolist()
                for key, value in reference_transfers.items()
            },
            "repeat_transfers": {
                key: value.tolist()
                for key, value in repeat_transfers.items()
            },
            "mean_transfers": {
                key: value.tolist()
                for key, value in mean_transfers.items()
            },
            "comparison_reference_key": comparison_reference_key,
            "censored_anchor_keys": censored,
        },
    }
    return result


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    frequency = np.asarray(p["frequency_Hz"], dtype=float)
    ref_key = p["comparison_reference_key"]
    censored = set(p["censored_anchor_keys"])

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.0), sharex=True)
    top, middle, bottom = axes

    for key, values in p["reference_transfers"].items():
        suffix = " (censored)" if key in censored else ""
        top.loglog(
            frequency,
            np.asarray(values, dtype=float),
            label=f"Ctes {key}{suffix}",
        )
    top.set_ylabel("12/06 normalized |H|")
    top.grid(True, which="both", alpha=0.2)
    top.legend(frameon=False, fontsize=8)

    for key, values in p["repeat_transfers"].items():
        suffix = " (censored)" if key in censored else ""
        middle.loglog(
            frequency,
            np.asarray(values, dtype=float),
            label=f"Ctes {key}{suffix}",
        )
    middle.set_ylabel("12/05 normalized |H|")
    middle.grid(True, which="both", alpha=0.2)
    middle.legend(frameon=False, fontsize=8)

    reference = np.asarray(
        p["mean_transfers"][ref_key],
        dtype=float,
    )
    for key, values in p["mean_transfers"].items():
        if key == ref_key:
            continue
        difference = transferdiag.difference_db(
            np.asarray(values, dtype=float),
            reference,
        )
        suffix = " (censored)" if key in censored else ""
        bottom.semilogx(
            frequency,
            difference,
            label=f"{key}/{ref_key}{suffix}",
        )
    bottom.axhline(0.0, linewidth=1.0)
    bottom.set_xlabel("Frequency [Hz]")
    bottom.set_ylabel("Cross-day mean transfer diff [dB]")
    bottom.grid(True, which="both", alpha=0.2)
    bottom.legend(frameon=False, fontsize=8)

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
        json.dumps(cleaned_result(result), indent=2, allow_nan=False)
        + "\n",
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
                "uncensored_anchors": result[
                    "primary_classification_anchor_keys"
                ],
                "censored_anchors": result["censored_anchor_keys"],
                "max_pairwise_rms_difference_dB": result[
                    "interpretation"
                ][
                    "max_pairwise_full_band_rms_difference_dB"
                ],
                "max_pairwise_abs_difference_dB": result[
                    "interpretation"
                ]["max_pairwise_abs_difference_dB"],
                "primary_material_vs_reference_rms_dB": result[
                    "interpretation"
                ][
                    "primary_material_vs_reference_full_band_rms_difference_dB"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
