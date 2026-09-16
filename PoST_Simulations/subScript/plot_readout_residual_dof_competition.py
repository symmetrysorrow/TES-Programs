"""Plot the tracked residual-readout fit against the experimental spectrum.

This utility consumes the JSON written by
readout_residual_dof_competition_diagnostic.py and reconstructs the exact
pre-analysis experimental target plus selected fitted models.  It does not
rerun any optimizer.

The default figure emphasizes the current generalization candidate
(pole_section_plus_c2), while also showing the fixed-reference baseline and
the full_order2 fit so its holdout degradation remains visible.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import preanalysis_readout_biquad_diagnostic as base  # noqa: E402
from subScript import readout_lowmid_identifiability_diagnostic as ident  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as diag  # noqa: E402
from subScript import shared_readout_cross_dataset_diagnostic as shared  # noqa: E402


DEFAULT_RESULT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_residual_dof_competition_diagnostic.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "readout_residual_dof_competition.png"
)
DEFAULT_FAMILIES = (
    "fixed_reference_readout",
    "pole_section_plus_c2",
    "full_order2",
)


def row_by_name(result: dict, name: str) -> dict:
    for row in result["residual_readout_fits"]:
        if row["name"] == name:
            return row
    raise KeyError(f"fit family not found in diagnostic JSON: {name}")


def reconstruct_curves(
    result: dict,
    config: dict,
    config_path: Path,
    experiment_path_override: Path | None = None,
):
    manifest_path = diag.resolve_config_path(
        config["manifest"],
        config_path,
    )
    manifest = json.loads(
        manifest_path.read_text(encoding="utf-8")
    )
    cases = shared.normalize_manifest(
        manifest,
        manifest_path,
    )
    repeat_case = ident.case_by_label(
        cases,
        config["repeat_case_label"],
    )
    summary = json.loads(
        repeat_case["summary"].read_text(encoding="utf-8")
    )
    comparison, experiment_path, _ = ident.comparison_for_case(
        repeat_case
    )
    if experiment_path_override is not None:
        experiment_path = experiment_path_override

    full_args = base.fit_args(summary)
    frequency = np.geomspace(
        full_args.fit_min_hz,
        full_args.fit_max_hz,
        full_args.fit_points,
    )
    target_context = base.reconstruct_pre_analysis_target(
        comparison,
        experiment_path,
        frequency,
    )
    target = np.asarray(target_context["target"], dtype=float)

    inherited = dict(summary["best_case_parameters"])
    fixed_r = float(
        result["fixed_detector_dc_baseline"]["R_TES_Ohm"]
    )
    scale_hz = float(
        result["readout_reference"]["reference_scale_Hz"]
    )
    white_asd = float(
        result["fixed_white_floor"]["white_asd_A_rtHz"]
    )

    curves = {}
    for family_name in DEFAULT_FAMILIES:
        row = row_by_name(result, family_name)
        candidate = dict(inherited)
        candidate["R"] = fixed_r
        candidate.update(
            {
                key: float(value)
                for key, value in row["detector_candidate"].items()
            }
        )
        readout = {
            key: float(value)
            for key, value in row["readout"].items()
        }
        model, point = diag.model_for_candidate(
            candidate,
            readout,
            frequency,
            scale_hz,
            white_asd,
        )
        if model is None:
            raise RuntimeError(
                f"could not reconstruct stable model for {family_name}: "
                f"{point}"
            )
        normalized = diag.opt.normalize_at(
            frequency,
            model,
            reference_hz=base.REFERENCE_HZ,
        )
        curves[family_name] = {
            "normalized": np.asarray(normalized, dtype=float),
            "fit_rms_dB": float(
                row["residual_metrics"]["rms_residual_dB"]
            ),
            "holdout_rms_dB": float(
                row["holdout_metrics"]["residual_metrics"][
                    "rms_residual_dB"
                ]
            ),
        }

    return {
        "frequency_Hz": frequency,
        "target": target,
        "curves": curves,
        "accepted_records": int(
            target_context["accepted_records"]
        ),
        "experiment_path": str(experiment_path),
    }


def make_plot(
    result: dict,
    reconstructed: dict,
    output: Path,
    *,
    show: bool = False,
):
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "matplotlib is required to write the comparison plot"
        ) from exc

    frequency = reconstructed["frequency_Hz"]
    target = reconstructed["target"]
    curves = reconstructed["curves"]

    fit_min = float(result["fit_region"]["min_Hz"])
    fit_max = float(result["fit_region"]["max_Hz_exclusive"])
    hold_min = float(result["holdout_region"]["min_Hz"])
    hold_max = float(result["holdout_region"]["max_Hz"])

    labels = {
        "fixed_reference_readout": "Fixed reference",
        "pole_section_plus_c2": "Current best: pole + Q + c2",
        "full_order2": "Full order-2 (+c4)",
    }

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10.5, 7.5),
        sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.0]},
    )
    top, bottom = axes

    top.plot(
        frequency,
        target,
        linewidth=2.1,
        label="Experiment",
    )
    for family_name in DEFAULT_FAMILIES:
        curve = curves[family_name]
        top.plot(
            frequency,
            curve["normalized"],
            linewidth=1.5,
            label=labels[family_name],
        )

    top.set_xscale("log")
    top.set_yscale("log")
    top.set_ylabel("Normalized ASD (1 kHz = 1)")
    top.grid(True, which="both", alpha=0.25)
    top.legend(frameon=False, fontsize=9)
    top.set_title(
        "Residual readout competition: experiment vs fitted simulation"
    )

    for family_name in DEFAULT_FAMILIES:
        curve = curves[family_name]
        residual_db = 20.0 * np.log10(
            curve["normalized"] / target
        )
        label = (
            f"{labels[family_name]}  "
            f"(fit {curve['fit_rms_dB']:.3f} dB, "
            f"holdout {curve['holdout_rms_dB']:.3f} dB)"
        )
        bottom.plot(
            frequency,
            residual_db,
            linewidth=1.4,
            label=label,
        )

    bottom.axhline(0.0, linewidth=1.0)
    bottom.axhline(1.0, linewidth=0.8, linestyle="--")
    bottom.axhline(-1.0, linewidth=0.8, linestyle="--")
    bottom.axvline(fit_max, linewidth=1.0, linestyle=":")
    bottom.set_xscale("log")
    bottom.set_xlabel("Frequency [Hz]")
    bottom.set_ylabel("Model / experiment [dB]")
    bottom.grid(True, which="both", alpha=0.25)
    bottom.legend(frameon=False, fontsize=8)

    bottom.text(
        np.sqrt(fit_min * fit_max),
        bottom.get_ylim()[1] * 0.82,
        "fit region",
        ha="center",
        va="top",
        fontsize=9,
    )
    bottom.text(
        np.sqrt(hold_min * hold_max),
        bottom.get_ylim()[1] * 0.82,
        "strict holdout",
        ha="center",
        va="top",
        fontsize=9,
    )

    fig.suptitle(
        f"{result['repeat_case']['label']}  "
        f"({reconstructed['accepted_records']} accepted records)",
        y=0.995,
        fontsize=10,
    )
    fig.tight_layout()

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
    )
    if show:
        plt.show()
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot the saved residual-readout fit without rerunning optimization."
        )
    )
    parser.add_argument(
        "--result",
        type=Path,
        default=DEFAULT_RESULT,
        help="Diagnostic JSON to plot.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=diag.DEFAULT_CONFIG,
        help="Tracked diagnostic config used to resolve the repeat dataset.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help="PNG output path.",
    )
    parser.add_argument(
        "--experiment-path",
        type=Path,
        help=(
            "Optional override for the raw experiment directory when the "
            "tracked acquisition path is mounted elsewhere."
        ),
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure after writing it.",
    )
    args = parser.parse_args()

    result = json.loads(
        args.result.read_text(encoding="utf-8")
    )
    config = json.loads(
        args.config.read_text(encoding="utf-8")
    )

    reconstructed = reconstruct_curves(
        result,
        config,
        args.config,
        args.experiment_path,
    )
    make_plot(
        result,
        reconstructed,
        args.output,
        show=args.show,
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "result": str(args.result),
                "experiment_path": reconstructed[
                    "experiment_path"
                ],
                "accepted_records": reconstructed[
                    "accepted_records"
                ],
                "families": {
                    name: {
                        "fit_rms_dB": values["fit_rms_dB"],
                        "holdout_rms_dB": values[
                            "holdout_rms_dB"
                        ],
                    }
                    for name, values in reconstructed[
                        "curves"
                    ].items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
