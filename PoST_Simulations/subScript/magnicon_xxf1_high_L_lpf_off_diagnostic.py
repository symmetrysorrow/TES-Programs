"""Refit the high-L grid with the Magnicon 10 kHz LPF OFF (and ON as comparator).

Each fixed-L row is fitted by the shared filter-order fitter with order=0
(unity, LPF OFF) and order=2 (documented Bessel, LPF ON), c2 fixed to zero.
Rows are independent, so one row can be run per process (--L-nH) and merged
afterwards (--merge).
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
WORK_DIR = ROOT / ".noise_optimization_work_rsh_sweep"
DEFAULT_CONFIG = (
    ROOT / "config" / "magnicon_xxf1_high_L_lpf_off_diagnostic_config.json"
)
ROW_DIR = WORK_DIR / "high_L_lpf_off_rows"
DEFAULT_OUTPUT = WORK_DIR / "magnicon_xxf1_high_L_lpf_off_diagnostic.json"
DEFAULT_FIGURE = WORK_DIR / "magnicon_xxf1_high_L_lpf_off_diagnostic.png"

for path in (ROOT, REPOSITORY_ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_filter_order_thermal_backbone_diagnostic as orderdiag  # noqa: E402
from subScript import magnicon_xxf1_high_L_thermal_backbone_c2_diagnostic as highL  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_thermal_backbone_pure_zero_diagnostic as thermal  # noqa: E402

ORDER = {"off": 0, "on": 2}


def _row_key(L_nH: float) -> str:
    return f"{L_nH:g}_nH"


def run_row(config: dict, config_path: Path, L_nH: float) -> dict:
    grid = [float(v) for v in config["L_grid_nH"]]
    index = next(i for i, v in enumerate(grid) if np.isclose(v, L_nH))
    thermal_path = thermal.resolve_config_path(
        config["base_thermal_backbone_config"], config_path
    )
    thermal_config = json.loads(thermal_path.read_text(encoding="utf-8"))
    separated_path = thermal.resolve_config_path(
        thermal_config["base_separated_substrate_config"], thermal_path
    )
    separated_config = json.loads(separated_path.read_text(encoding="utf-8"))
    stack = compensation._load_stack(separated_config, separated_path)
    base_cross = stack["base_cross"]
    problems = {
        prefix: crossday._problem_for_case(
            stack["magnicon_config"],
            stack["magnicon_config_path"],
            str(base_cross[f"{prefix}_case_label"]),
        )
        for prefix in ("reference", "repeat")
    }
    fixed_base = highL._base_config_at_fixed_L(
        separated_config,
        L_nH * 1e-9,
        float(config["numerical_fixed_L_relative_half_width"]),
    )
    row_config = copy.deepcopy(config)
    row_config["optimizer"]["seed_offset"] = (
        int(config["optimizer"]["seed_offset"]) + 3000 * index
    )
    row = {"nominal_fixed_L_nH": L_nH, "states": {}}
    plot = {
        "frequency_Hz": problems["reference"]["frequency"].tolist(),
        "reference_target": problems["reference"]["target"].tolist(),
        "repeat_target": problems["repeat"]["target"].tolist(),
        "reference_raw": problems["reference"]["raw_target"].tolist(),
        "repeat_raw": problems["repeat"]["raw_target"].tolist(),
    }
    for state in config["states"]:
        fit = orderdiag.fit_order(
            order=ORDER[state],
            config=row_config,
            thermal_config=thermal_config,
            separated_config=fixed_base,
            material_C_tes=float(opt.C_TES_MATERIAL_J_PER_K),
            reference_problem=problems["reference"],
            repeat_problem=problems["repeat"],
        )
        print("fit", L_nH, state, fit["joint_continuum_rms_dB"], flush=True)
        row["states"][state] = {
            k: v for k, v in fit.items() if not k.startswith("_")
        }
        for prefix in ("reference", "repeat"):
            plot[f"{prefix}_{state}"] = fit["_models"][prefix].tolist()
    row["_plot"] = plot
    return row


def merge(config: dict) -> dict:
    rows = {}
    for L_nH in config["L_grid_nH"]:
        path = ROW_DIR / f"{_row_key(float(L_nH))}.json"
        if path.exists():
            rows[_row_key(float(L_nH))] = json.loads(
                path.read_text(encoding="utf-8")
            )
    return {
        "diagnostic_only": True,
        "purpose": config["purpose"],
        "guardrail": config["guardrail"],
        "rows": {
            key: {
                "nominal_fixed_L_nH": r["nominal_fixed_L_nH"],
                "joint_rms_dB": {
                    s: r["states"][s]["joint_continuum_rms_dB"]
                    for s in r["states"]
                },
                "off_minus_on_dB": r["states"]["off"][
                    "joint_continuum_rms_dB"
                ]
                - r["states"]["on"]["joint_continuum_rms_dB"],
                "states": r["states"],
            }
            for key, r in rows.items()
        },
        "_rows_with_plot": rows,
    }


def make_plot(result: dict, output: Path) -> None:
    import matplotlib.pyplot as plt

    rows = result["_rows_with_plot"]
    fig, axes = plt.subplots(
        2, 2, figsize=(13.8, 9.0), gridspec_kw={"height_ratios": [1, 1.6]}
    )
    keys = sorted(rows, key=lambda k: rows[k]["nominal_fixed_L_nH"])
    Ls = [rows[k]["nominal_fixed_L_nH"] for k in keys]
    ax = axes[0, 0]
    for state, label in (
        ("on", "LPF ON (2nd-order Bessel)"),
        ("off", "LPF OFF"),
    ):
        ax.semilogx(
            Ls,
            [rows[k]["states"][state]["joint_continuum_rms_dB"] for k in keys],
            "o-",
            label=label,
        )
    ax.set_xlabel("Fixed L [nH]")
    ax.set_ylabel("Joint RMS [dB], c2=0")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(frameon=False)
    axes[0, 1].axis("off")

    pick = min(
        keys, key=lambda k: rows[k]["states"]["off"]["joint_continuum_rms_dB"]
    )
    p = rows[pick]["_plot"]
    f = np.asarray(p["frequency_Hz"])
    for col, (prefix, title) in enumerate(
        (("reference", "2024-12-06"), ("repeat", "2024-12-05"))
    ):
        a = axes[1, col]
        t = np.asarray(p[f"{prefix}_target"])
        for state, label in (("on", "LPF ON"), ("off", "LPF OFF")):
            a.semilogx(
                f,
                20 * np.log10(np.asarray(p[f"{prefix}_{state}"]) / t),
                lw=1.1,
                label=label,
            )
        a.axhline(0, color="k", lw=0.9)
        a.set_title(
            f"{title}, L={rows[pick]['nominal_fixed_L_nH']:g} nH "
            "(best LPF-OFF row)",
            fontsize=10,
        )
        a.set_xlabel("Frequency [Hz]")
        a.set_ylabel("Model / measured [dB]")
        a.grid(True, which="both", alpha=0.2)
        a.legend(frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--L-nH", type=float, help="fit only this grid row and save it"
    )
    parser.add_argument("--merge", action="store_true")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    if args.L_nH is not None:
        row = run_row(config, args.config, args.L_nH)
        ROW_DIR.mkdir(parents=True, exist_ok=True)
        (ROW_DIR / f"{_row_key(args.L_nH)}.json").write_text(
            json.dumps(row, allow_nan=False), encoding="utf-8"
        )
    if args.merge:
        result = merge(config)
        make_plot(result, args.figure)
        cleaned = {k: v for k, v in result.items() if k != "_rows_with_plot"}
        args.output.write_text(
            json.dumps(cleaned, indent=2, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        print(
            json.dumps(
                {k: v["joint_rms_dB"] for k, v in cleaned["rows"].items()},
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
