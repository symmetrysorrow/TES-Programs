"""Compare observable fitted ASD sensitivity to C_tes anchor choice.

The preceding C_tes-anchor transfer diagnostic showed that the inferred
Magnicon+c2 readout transfer changes by up to about 10 dB across uncensored
C_tes anchors.  This diagnostic asks whether that large decomposition change
survives in the actual fitted model ASD.

For each C_tes anchor and each adjacent-day dataset, the same anchored c2-free
fit is reconstructed in three layers:

  1. full fitted ASD, including the profiled post-filter white floor;
  2. the same fitted detector/readout model with post-filter white set to zero;
  3. the normalized Magnicon+c2 readout transfer alone.

All three are normalized at the same reference frequency and compared across
C_tes anchors.  Rows where c2 is at a configured bound are censored from the
primary finite-shape classification.

A stable full fitted ASD together with a strongly anchor-sensitive readout
transfer would demonstrate decomposition non-identifiability/cancellation.  It
would not identify which C_tes or readout correction is physically correct.
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
    / "magnicon_xxf1_Ctes_anchor_composite_ASD_sensitivity_config.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_Ctes_anchor_composite_ASD_sensitivity_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_Ctes_anchor_composite_ASD_sensitivity_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_Ctes_anchor_transfer_shape_sensitivity_diagnostic as transfer_sensitivity  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_transfer_shape_stability_diagnostic as transferdiag  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_physical_L_Ctes_anchor_c2_diagnostic as anchored  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def normalized_shape(frequency, values, reference_hz):
    frequency = np.asarray(frequency, dtype=float)
    values = np.asarray(values, dtype=float)
    return opt.normalize_at(
        frequency,
        values,
        reference_hz=float(reference_hz),
    )


def raw_c2_is_censored(fit_row: dict) -> bool:
    hits = fit_row.get("readout_boundary_hits", {}).get("c2", {})
    return bool(hits.get("at_lower", False) or hits.get("at_upper", False))


def zero_white_model(problem: dict, fit_row: dict):
    model, point = residual.model_for_candidate(
        fit_row["_detector_full"],
        fit_row["_readout_full"],
        problem["frequency"],
        problem["scale_hz"],
        0.0,
    )
    if model is None:
        raise RuntimeError(
            "zero-white model reconstruction failed: "
            + json.dumps(point, default=float)
        )
    return np.asarray(model, dtype=float)


def pairwise_summaries(
    *,
    frequency,
    shapes,
    included_keys,
    bands,
    anchors,
):
    rows = {}
    for first, second in itertools.combinations(included_keys, 2):
        key = f"{first}_vs_{second}"
        row = transfer_sensitivity.difference_summary(
            frequency,
            shapes[second],
            shapes[first],
            bands,
            anchors,
        )
        row["semantics"] = f"{second} / {first}"
        rows[key] = row
    return rows


def clean_difference(row: dict):
    return {
        key: value
        for key, value in row.items()
        if key != "_difference_dB"
    }


def sensitivity_ratio(observable_rms_dB, readout_rms_dB):
    observable = float(observable_rms_dB)
    readout = float(readout_rms_dB)
    if readout <= 0.0:
        return None
    return float(observable / readout)


def classify_decomposition(
    *,
    observable_classification: dict,
    readout_classification: dict,
):
    observable_pass = observable_classification[
        "passes_anchor_sensitivity_screen"
    ]
    readout_pass = readout_classification[
        "passes_anchor_sensitivity_screen"
    ]

    if observable_pass is True and readout_pass is False:
        return (
            "Ctes_readout_decomposition_nonidentifiable_"
            "but_observable_ASD_stable"
        )
    if observable_pass is False:
        return "observable_ASD_is_Ctes_anchor_sensitive"
    if observable_pass is True and readout_pass is True:
        return "observable_and_readout_shapes_stable_across_Ctes_anchors"
    return "Ctes_anchor_decomposition_classification_unresolved"


def _shape_triplet(
    *,
    problem: dict,
    fit_row: dict,
    normalization_hz: float,
):
    frequency = problem["frequency"]
    full_model = normalized_shape(
        frequency,
        fit_row["_model_full"],
        normalization_hz,
    )
    no_white = normalized_shape(
        frequency,
        zero_white_model(problem, fit_row),
        normalization_hz,
    )
    readout = transferdiag.readout_components(
        frequency,
        fit_row["_readout_full"],
        problem["scale_hz"],
        normalization_hz,
    )["total"]
    return {
        "full_fitted_ASD": np.asarray(full_model, dtype=float),
        "zero_post_filter_white_ASD": np.asarray(no_white, dtype=float),
        "readout_transfer": np.asarray(readout, dtype=float),
    }


def _geometric_mean_shapes(reference_shapes, repeat_shapes):
    return {
        name: transfer_sensitivity.geometric_mean_transfer(
            reference_shapes[name],
            repeat_shapes[name],
        )
        for name in reference_shapes
    }


def _layer_summary(
    *,
    frequency,
    layer_shapes,
    included_keys,
    bands,
    anchors,
    screen,
):
    pairwise = pairwise_summaries(
        frequency=frequency,
        shapes=layer_shapes,
        included_keys=included_keys,
        bands=bands,
        anchors=anchors,
    )
    spread = transfer_sensitivity.spread_summary(
        {key: layer_shapes[key] for key in included_keys}
    )
    classification = transfer_sensitivity.classify_pairwise(
        pairwise,
        screen,
    )
    return {
        "pairwise": {
            key: clean_difference(row)
            for key, row in pairwise.items()
        },
        "spread": {
            key: value
            for key, value in spread.items()
            if key != "_spread_dB"
        },
        "classification": classification,
        "_pairwise_raw": pairwise,
    }


def _clean_layer_summary(row: dict):
    return {
        key: value
        for key, value in row.items()
        if key != "_pairwise_raw"
    }


def run(config: dict, config_path: Path):
    base_config_path = resolve_config_path(
        config["base_anchor_config"],
        config_path,
    )
    base_config = json.loads(
        base_config_path.read_text(encoding="utf-8")
    )

    stack = compensation._load_stack(
        base_config,
        base_config_path,
    )
    base_cross = stack["base_cross"]
    c2_config = stack["c2_config"]
    magnicon_config = stack["magnicon_config"]

    reference_problem = crossday._problem_for_case(
        magnicon_config,
        stack["magnicon_config_path"],
        str(base_cross["reference_case_label"]),
    )
    repeat_problem = crossday._problem_for_case(
        magnicon_config,
        stack["magnicon_config_path"],
        str(base_cross["repeat_case_label"]),
    )
    if not np.allclose(
        reference_problem["frequency"],
        repeat_problem["frequency"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference/repeat frequency grids differ")
    if not np.isclose(
        reference_problem["scale_hz"],
        repeat_problem["scale_hz"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference/repeat readout scales differ")

    filter_cfg = base_config["magnicon_filter"]
    canonical = xxf1.second_order_bessel_canonical(
        float(filter_cfg["cutoff_Hz"]),
        str(filter_cfg["normalization"]),
    )
    readout_reference = {
        "pole_Hz": float(canonical["pole_Hz"]),
        "pole_Q": float(canonical["pole_Q"]),
        "c2": 0.0,
        "c4": float(filter_cfg["c4_fixed"]),
    }
    if readout_reference["c4"] != 0.0:
        raise ValueError("this diagnostic requires c4=0")

    fixed_L = float(base_config["external_L_H"])
    material_ctes = float(opt.C_TES_MATERIAL_J_PER_K)
    ctes_cfg = base_config["C_tes_anchor"]
    multipliers = [float(v) for v in ctes_cfg["multipliers"]]
    normalization_hz = float(
        config["normalization_reference_Hz"]
    )

    reference_layers = {
        "full_fitted_ASD": {},
        "zero_post_filter_white_ASD": {},
        "readout_transfer": {},
    }
    repeat_layers = {
        "full_fitted_ASD": {},
        "zero_post_filter_white_ASD": {},
        "readout_transfer": {},
    }
    mean_layers = {
        "full_fitted_ASD": {},
        "zero_post_filter_white_ASD": {},
        "readout_transfer": {},
    }
    anchor_rows = {}
    censored_keys = []

    for index, multiplier in enumerate(multipliers):
        key = anchored._multiplier_key(multiplier)
        fixed_ctes = material_ctes * multiplier

        reference_fit = anchored._fit_case(
            problem=reference_problem,
            fixed_L_H=fixed_L,
            fixed_C_tes_J_per_K=fixed_ctes,
            readout_reference=readout_reference,
            c2_upper_bound=float(c2_config["c2_upper_bound"]),
            materiality_screen=c2_config["materiality_screen"],
            seed_offset=21100 + 600 * index,
        )
        repeat_fit = anchored._fit_case(
            problem=repeat_problem,
            fixed_L_H=fixed_L,
            fixed_C_tes_J_per_K=fixed_ctes,
            readout_reference=readout_reference,
            c2_upper_bound=float(c2_config["c2_upper_bound"]),
            materiality_screen=c2_config["materiality_screen"],
            seed_offset=21400 + 600 * index,
        )

        reference_free = reference_fit["c2_free"]
        repeat_free = repeat_fit["c2_free"]
        reference_shape = _shape_triplet(
            problem=reference_problem,
            fit_row=reference_free,
            normalization_hz=normalization_hz,
        )
        repeat_shape = _shape_triplet(
            problem=repeat_problem,
            fit_row=repeat_free,
            normalization_hz=normalization_hz,
        )
        mean_shape = _geometric_mean_shapes(
            reference_shape,
            repeat_shape,
        )

        for layer in reference_layers:
            reference_layers[layer][key] = reference_shape[layer]
            repeat_layers[layer][key] = repeat_shape[layer]
            mean_layers[layer][key] = mean_shape[layer]

        reference_censored = raw_c2_is_censored(reference_free)
        repeat_censored = raw_c2_is_censored(repeat_free)
        row_censored = bool(reference_censored or repeat_censored)
        if row_censored:
            censored_keys.append(key)

        anchor_rows[key] = {
            "C_tes_multiplier": multiplier,
            "fixed_C_tes_J_per_K": fixed_ctes,
            "fixed_L_H": fixed_L,
            "reference_day": {
                "c2": float(reference_free["readout"]["c2"]),
                "white_asd_A_rtHz": float(
                    reference_free["profiled_white_asd_A_rtHz"]
                ),
                "continuum_rms_dB": xxf1._rms(reference_free),
                "c2_censored": reference_censored,
            },
            "repeat_day": {
                "c2": float(repeat_free["readout"]["c2"]),
                "white_asd_A_rtHz": float(
                    repeat_free["profiled_white_asd_A_rtHz"]
                ),
                "continuum_rms_dB": xxf1._rms(repeat_free),
                "c2_censored": repeat_censored,
            },
            "row_censored": row_censored,
        }

    all_keys = [anchored._multiplier_key(v) for v in multipliers]
    exclude_censored = bool(
        config["exclude_c2_boundary_rows_from_primary_classification"]
    )
    included_keys = [
        key
        for key in all_keys
        if not (exclude_censored and key in censored_keys)
    ]

    screen = config["observable_shape_screen"]
    frequency = reference_problem["frequency"]

    reference_summary = {}
    repeat_summary = {}
    mean_summary = {}
    for layer in reference_layers:
        reference_summary[layer] = _layer_summary(
            frequency=frequency,
            layer_shapes=reference_layers[layer],
            included_keys=included_keys,
            bands=config["bands_Hz"],
            anchors=config["anchors_Hz"],
            screen=screen,
        )
        repeat_summary[layer] = _layer_summary(
            frequency=frequency,
            layer_shapes=repeat_layers[layer],
            included_keys=included_keys,
            bands=config["bands_Hz"],
            anchors=config["anchors_Hz"],
            screen=screen,
        )
        mean_summary[layer] = _layer_summary(
            frequency=frequency,
            layer_shapes=mean_layers[layer],
            included_keys=included_keys,
            bands=config["bands_Hz"],
            anchors=config["anchors_Hz"],
            screen=screen,
        )

    observable_classification = mean_summary[
        "full_fitted_ASD"
    ]["classification"]
    no_white_classification = mean_summary[
        "zero_post_filter_white_ASD"
    ]["classification"]
    readout_classification = mean_summary[
        "readout_transfer"
    ]["classification"]

    comparison_key = anchored._multiplier_key(
        float(config["comparison_reference_multiplier"])
    )
    material_key = anchored._multiplier_key(
        float(config["primary_material_multiplier"])
    )
    if comparison_key not in mean_layers["full_fitted_ASD"]:
        raise ValueError("comparison reference multiplier missing")
    if material_key not in mean_layers["full_fitted_ASD"]:
        raise ValueError("primary material multiplier missing")

    primary_comparisons = {}
    for layer in mean_layers:
        row = transfer_sensitivity.difference_summary(
            frequency,
            mean_layers[layer][material_key],
            mean_layers[layer][comparison_key],
            config["bands_Hz"],
            config["anchors_Hz"],
        )
        row["semantics"] = (
            f"{material_key} / {comparison_key} using geometric-mean "
            f"cross-day {layer}"
        )
        primary_comparisons[layer] = row

    readout_rms = float(
        primary_comparisons["readout_transfer"]["full_1_200k"][
            "rms_difference_dB"
        ]
    )
    observable_rms = float(
        primary_comparisons["full_fitted_ASD"]["full_1_200k"][
            "rms_difference_dB"
        ]
    )
    no_white_rms = float(
        primary_comparisons["zero_post_filter_white_ASD"][
            "full_1_200k"
        ]["rms_difference_dB"]
    )

    classification = classify_decomposition(
        observable_classification=observable_classification,
        readout_classification=readout_classification,
    )

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "tested_question": (
            "When fixed C_tes assumptions drive large changes in the fitted "
            "Magnicon+c2 transfer, is the final observable fitted ASD shape "
            "also anchor-sensitive, or do detector dynamics/white compensate?"
        ),
        "normalization_reference_Hz": normalization_hz,
        "physical_anchor_context": {
            "fixed_L_H": fixed_L,
            "C_tes_material_J_per_K": material_ctes,
            "C_tes_multipliers": multipliers,
            "C_tes_source_semantics": ctes_cfg["source_semantics"],
        },
        "layer_semantics": {
            "full_fitted_ASD": (
                "best-fit pre-analysis ASD including the independently "
                "profiled post-filter white floor"
            ),
            "zero_post_filter_white_ASD": (
                "same best-fit detector/readout state with only the profiled "
                "post-filter white term set to zero"
            ),
            "readout_transfer": (
                "Magnicon 10 kHz Bessel times phenomenological c2 numerator"
            ),
        },
        "anchor_rows": anchor_rows,
        "primary_classification_anchor_keys": included_keys,
        "censored_anchor_keys": censored_keys,
        "reference_day_anchor_sensitivity": {
            key: _clean_layer_summary(value)
            for key, value in reference_summary.items()
        },
        "repeat_day_anchor_sensitivity": {
            key: _clean_layer_summary(value)
            for key, value in repeat_summary.items()
        },
        "cross_day_geometric_mean_anchor_sensitivity": {
            "semantics": (
                "For each C_tes anchor, geometric mean of the separately "
                "normalized 12/06 and 12/05 shapes, followed by anchor-to-"
                "anchor comparison."
            ),
            **{
                key: _clean_layer_summary(value)
                for key, value in mean_summary.items()
            },
        },
        "primary_material_vs_reference_anchor": {
            key: clean_difference(value)
            for key, value in primary_comparisons.items()
        },
        "interpretation": {
            "classification": classification,
            "observable_full_ASD_classification": (
                observable_classification["classification"]
            ),
            "zero_post_filter_white_ASD_classification": (
                no_white_classification["classification"]
            ),
            "readout_transfer_classification": (
                readout_classification["classification"]
            ),
            "observable_full_ASD_passes_screen": (
                observable_classification[
                    "passes_anchor_sensitivity_screen"
                ]
            ),
            "zero_post_filter_white_ASD_passes_screen": (
                no_white_classification[
                    "passes_anchor_sensitivity_screen"
                ]
            ),
            "readout_transfer_passes_screen": (
                readout_classification[
                    "passes_anchor_sensitivity_screen"
                ]
            ),
            "primary_material_vs_reference": {
                "observable_full_ASD_rms_difference_dB": observable_rms,
                "zero_post_filter_white_ASD_rms_difference_dB": (
                    no_white_rms
                ),
                "readout_transfer_rms_difference_dB": readout_rms,
                "observable_to_readout_sensitivity_ratio": (
                    sensitivity_ratio(observable_rms, readout_rms)
                ),
                "zero_white_to_readout_sensitivity_ratio": (
                    sensitivity_ratio(no_white_rms, readout_rms)
                ),
            },
            "c2_bound_censoring_present": bool(censored_keys),
            "excluded_censored_anchor_keys": censored_keys,
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_anchor_config": str(base_config_path),
        },
        "_plot": {
            "frequency_Hz": frequency.tolist(),
            "comparison_reference_key": comparison_key,
            "censored_anchor_keys": censored_keys,
            "mean_layers": {
                layer: {
                    key: values.tolist()
                    for key, values in rows.items()
                }
                for layer, rows in mean_layers.items()
            },
        },
    }
    return result


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    p = result["_plot"]
    frequency = np.asarray(p["frequency_Hz"], dtype=float)
    comparison_key = p["comparison_reference_key"]
    censored = set(p["censored_anchor_keys"])

    fig, axes = plt.subplots(3, 1, figsize=(10.0, 9.0), sharex=True)
    top, middle, bottom = axes

    for key, values in p["mean_layers"]["full_fitted_ASD"].items():
        suffix = " (censored)" if key in censored else ""
        top.loglog(
            frequency,
            np.asarray(values, dtype=float),
            label=f"Ctes {key}{suffix}",
        )
    top.set_ylabel("Cross-day mean full ASD")
    top.grid(True, which="both", alpha=0.2)
    top.legend(frameon=False, fontsize=8)

    for key, values in p["mean_layers"][
        "zero_post_filter_white_ASD"
    ].items():
        suffix = " (censored)" if key in censored else ""
        middle.loglog(
            frequency,
            np.asarray(values, dtype=float),
            label=f"Ctes {key}{suffix}",
        )
    middle.set_ylabel("Cross-day mean zero-white ASD")
    middle.grid(True, which="both", alpha=0.2)
    middle.legend(frameon=False, fontsize=8)

    full_reference = np.asarray(
        p["mean_layers"]["full_fitted_ASD"][comparison_key],
        dtype=float,
    )
    readout_reference = np.asarray(
        p["mean_layers"]["readout_transfer"][comparison_key],
        dtype=float,
    )
    for key in p["mean_layers"]["full_fitted_ASD"]:
        if key == comparison_key:
            continue
        suffix = " (censored)" if key in censored else ""
        full_difference = transferdiag.difference_db(
            np.asarray(
                p["mean_layers"]["full_fitted_ASD"][key],
                dtype=float,
            ),
            full_reference,
        )
        readout_difference = transferdiag.difference_db(
            np.asarray(
                p["mean_layers"]["readout_transfer"][key],
                dtype=float,
            ),
            readout_reference,
        )
        bottom.semilogx(
            frequency,
            full_difference,
            label=f"full {key}/{comparison_key}{suffix}",
        )
        bottom.semilogx(
            frequency,
            readout_difference,
            linestyle="--",
            label=f"readout {key}/{comparison_key}{suffix}",
        )
    bottom.axhline(0.0, linewidth=1.0)
    bottom.set_xlabel("Frequency [Hz]")
    bottom.set_ylabel("Anchor difference [dB]")
    bottom.grid(True, which="both", alpha=0.2)
    bottom.legend(frameon=False, fontsize=7)

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
                "primary_material_vs_reference": result[
                    "interpretation"
                ]["primary_material_vs_reference"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
