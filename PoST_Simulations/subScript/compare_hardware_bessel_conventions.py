"""Diagnose the confirmed 100 kHz hardware Bessel convention.

This script holds the selected correlated-search balanced TES parameters fixed
and changes only the analog hardware transfer convention before the ADC:

A. current repository convention: 4th-order 100 kHz Bessel, ``norm="phase"``;
B. diagnostic bypass: no analog hardware attenuation on main or alias terms;
C. magnitude convention: 4th-order 100 kHz Bessel, ``norm="mag"`` (-3 dB at
   100 kHz under SciPy's magnitude normalization).

Every case then uses the same sampled/aliased ASD, the same finite-record RNG
seed and record count, and the same measured 10 kHz software Bessel analysis.
No TES parameter is re-ranked or re-fitted here.  The bypass is diagnostic only
and is not a proposed replacement for the confirmed physical hardware filter.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "PoST_Simulations"))
sys.path.insert(0, str(ROOT / "PoST_Simulations" / "subScript"))

from high_frequency_noise_comparison import (  # noqa: E402
    experiment_asd,
    log_interp,
    normalized,
)
from noise_measurement_model import (  # noqa: E402
    DEFAULT_FINITE_RECORD_SEED,
    DEFAULT_HARDWARE_BESSEL_ORDER,
    HARDWARE_BESSEL_CUTOFF_HZ,
    finite_record_post_analysis_asd,
    hardware_filter_magnitude,
    hardware_sampled_asd,
)


CASES = (
    {
        "key": "phase",
        "label": '100 kHz Bessel, phase norm (current)',
        "norm": "phase",
        "bypass": False,
    },
    {
        "key": "bypass",
        "label": "100 kHz hardware Bessel bypass (diagnostic)",
        "norm": "phase",
        "bypass": True,
    },
    {
        "key": "mag",
        "label": '100 kHz Bessel, mag norm (-3 dB @ 100 kHz)',
        "norm": "mag",
        "bypass": False,
    },
)


def dump(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def load_balanced_parameters(path: Path) -> tuple[dict, str, str]:
    document = json.loads(path.read_text(encoding="utf-8"))
    for section in ("best_finite", "best_deterministic"):
        row = document.get(section, {}).get("all")
        if isinstance(row, dict) and isinstance(row.get("parameters"), dict):
            identifier = str(row.get("trial_id", row.get("scenario_id", "unknown")))
            return dict(row["parameters"]), section, identifier
    raise RuntimeError(
        "search JSON does not contain best_finite['all'] or "
        "best_deterministic['all'] parameters"
    )


def response_db(value: float) -> float:
    return float(20.0 * np.log10(max(float(value), np.finfo(float).tiny)))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-root", type=Path, required=True)
    ap.add_argument("--case-dir", type=Path, required=True)
    ap.add_argument("--search-json", type=Path)
    ap.add_argument("--finite-records", type=int, default=0)
    ap.add_argument("--finite-seed", type=int, default=DEFAULT_FINITE_RECORD_SEED)
    args = ap.parse_args()

    search_json = args.search_json or (
        args.case_dir / "correlated_high_frequency_parameter_search.json"
    )
    parameters, parameter_section, parameter_id = load_balanced_parameters(search_json)

    exp_freq, exp_asd, exp_paths, acquisition = experiment_asd(args.target_root)
    rate = float(acquisition["rate_Hz"])
    sample = int(acquisition["samples"])
    analysis_cutoff = float(acquisition["analysis_bessel_cutoff_Hz"])
    record_count = int(args.finite_records) if args.finite_records > 0 else len(exp_paths)
    if record_count <= 0:
        raise RuntimeError("finite-record comparison requires at least one record")

    frequency = np.fft.rfftfreq(sample, d=1.0 / rate)
    exp_norm = normalized(exp_asd, exp_freq)
    simulations: dict[str, dict] = {}

    for case in CASES:
        pre_analysis = hardware_sampled_asd(
            parameters,
            frequency,
            rate_hz=rate,
            cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
            order=DEFAULT_HARDWARE_BESSEL_ORDER,
            norm=case["norm"],
            bypass=bool(case["bypass"]),
        )
        post_analysis = finite_record_post_analysis_asd(
            pre_analysis,
            sample,
            rate,
            analysis_cutoff_hz=analysis_cutoff,
            records=record_count,
            seed=int(args.finite_seed),
        )
        simulations[case["key"]] = {
            "case": case,
            "pre_analysis_asd": pre_analysis,
            "pre_analysis_normalized": normalized(pre_analysis, frequency),
            "post_analysis_asd": post_analysis,
            "post_analysis_normalized": normalized(post_analysis, frequency),
        }

    compare_points = np.asarray(
        [f for f in (10_000.0, 50_000.0, 100_000.0, 150_000.0, 200_000.0)
         if f <= rate / 2.0],
        dtype=float,
    )
    exp_at_points = log_interp(exp_freq[1:], exp_norm[1:], compare_points)
    metrics = {}
    for case in CASES:
        values = simulations[case["key"]]["post_analysis_normalized"]
        sim_at_points = log_interp(frequency[1:], values[1:], compare_points)
        pre_values = simulations[case["key"]]["pre_analysis_normalized"]
        pre_at_points = log_interp(frequency[1:], pre_values[1:], compare_points)
        metrics[case["key"]] = {
            "label": case["label"],
            "post_analysis_simulation_over_experiment": {
                f"{float(f):g}_Hz": float(sim / exp)
                for f, sim, exp in zip(compare_points, sim_at_points, exp_at_points)
            },
            "pre_analysis_normalized_asd": {
                f"{float(f):g}_Hz": float(value)
                for f, value in zip(compare_points, pre_at_points)
            },
        }

    response_at_cutoff = {}
    for case in CASES:
        response = hardware_filter_magnitude(
            np.asarray([HARDWARE_BESSEL_CUTOFF_HZ]),
            cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
            order=DEFAULT_HARDWARE_BESSEL_ORDER,
            norm=case["norm"],
            bypass=bool(case["bypass"]),
        )[0]
        response_at_cutoff[case["key"]] = {
            "magnitude": float(response),
            "dB": response_db(response),
        }

    result = {
        "stage": "hardware_bessel_convention_diagnostic",
        "purpose": (
            "Hold the correlated balanced TES parameters fixed and isolate the "
            "effect of the 100 kHz analog hardware Bessel convention."
        ),
        "parameter_source": str(search_json),
        "parameter_section": parameter_section,
        "parameter_id": parameter_id,
        "parameters_refit_for_this_diagnostic": False,
        "hardware_confirmed_physical_configuration": {
            "order": DEFAULT_HARDWARE_BESSEL_ORDER,
            "nominal_cutoff_Hz": HARDWARE_BESSEL_CUTOFF_HZ,
        },
        "bypass_is_diagnostic_only": True,
        "analysis_chain_fixed": {
            "rate_Hz": rate,
            "samples": sample,
            "analysis_bessel_cutoff_Hz": analysis_cutoff,
            "finite_record_count": record_count,
            "finite_record_seed": int(args.finite_seed),
        },
        "same_finite_record_seed_for_all_cases": True,
        "cases": [dict(case) for case in CASES],
        "hardware_magnitude_at_100_kHz": response_at_cutoff,
        "comparison_metrics": metrics,
    }
    dump(args.case_dir / "hardware_bessel_convention_diagnostic.json", result)

    import matplotlib.pyplot as plt

    plot_frequency = np.logspace(
        np.log10(rate / sample),
        np.log10(rate / 2.0),
        1200,
    )
    experiment_plot = log_interp(exp_freq[1:], exp_norm[1:], plot_frequency)

    plt.figure(figsize=(9.2, 5.7))
    plt.plot(
        plot_frequency,
        experiment_plot,
        color="black",
        lw=1.5,
        label=f"Experiment ({len(exp_paths)} accepted records)",
    )
    for case in CASES:
        values = simulations[case["key"]]["post_analysis_normalized"]
        curve = log_interp(frequency[1:], values[1:], plot_frequency)
        plt.plot(plot_frequency, curve, lw=1.15, label=case["label"])
    plt.axvline(HARDWARE_BESSEL_CUTOFF_HZ, color="black", lw=0.8, ls=":")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(rate / sample, rate / 2.0)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Normalized ASD (ASD / ASD at 1 kHz)")
    plt.title("Hardware Bessel convention diagnostic — fixed balanced TES parameters")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.case_dir / "hardware_bessel_convention_comparison.png", dpi=180)
    plt.close()

    ratio_mask = plot_frequency >= 1_000.0
    plt.figure(figsize=(9.2, 5.2))
    for case in CASES:
        values = simulations[case["key"]]["post_analysis_normalized"]
        curve = log_interp(frequency[1:], values[1:], plot_frequency)
        plt.plot(
            plot_frequency[ratio_mask],
            curve[ratio_mask] / experiment_plot[ratio_mask],
            lw=1.15,
            label=case["label"],
        )
    plt.axhline(1.0, color="black", lw=1.0)
    plt.axvline(HARDWARE_BESSEL_CUTOFF_HZ, color="black", lw=0.8, ls=":")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1_000.0, rate / 2.0)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Simulation / experiment")
    plt.title("Hardware Bessel convention diagnostic — residual ratio")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.case_dir / "hardware_bessel_convention_ratio.png", dpi=180)
    plt.close()

    plt.figure(figsize=(9.2, 5.2))
    for case in CASES:
        values = simulations[case["key"]]["pre_analysis_normalized"]
        curve = log_interp(frequency[1:], values[1:], plot_frequency)
        plt.plot(plot_frequency, curve, lw=1.15, label=case["label"])
    plt.axvline(HARDWARE_BESSEL_CUTOFF_HZ, color="black", lw=0.8, ls=":")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(rate / sample, rate / 2.0)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Pre-analysis normalized ASD")
    plt.title("Simulation before the 10 kHz software filter")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.case_dir / "hardware_bessel_convention_pre_analysis.png", dpi=180)
    plt.close()

    response_frequency = np.logspace(3.0, np.log10(rate / 2.0), 800)
    plt.figure(figsize=(8.6, 4.9))
    for case in CASES:
        response = hardware_filter_magnitude(
            response_frequency,
            cutoff_hz=HARDWARE_BESSEL_CUTOFF_HZ,
            order=DEFAULT_HARDWARE_BESSEL_ORDER,
            norm=case["norm"],
            bypass=bool(case["bypass"]),
        )
        plt.plot(response_frequency, response, lw=1.2, label=case["label"])
    plt.axvline(HARDWARE_BESSEL_CUTOFF_HZ, color="black", lw=0.8, ls=":")
    plt.axhline(1.0 / np.sqrt(2.0), color="black", lw=0.8, ls="--")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1_000.0, rate / 2.0)
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Analog hardware magnitude")
    plt.title("4th-order 100 kHz Bessel normalization conventions")
    plt.grid(True, which="both", alpha=0.25)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(args.case_dir / "hardware_bessel_convention_response.png", dpi=180)
    plt.close()

    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
