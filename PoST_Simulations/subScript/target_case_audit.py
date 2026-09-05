"""Capability-level readiness audit for a versioned target case.

This gate never fills values from the generic simulation input and never uses
the noise spectrum as a parameter source.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


OPERATING_POINT = ("T_c", "T_bath", "R", "G_tes-bath", "n")
PYTHON_STABILITY = OPERATING_POINT + (
    "C_abs", "C_tes", "G_abs-abs", "G_abs-tes", "R_l", "alpha", "beta", "L",
)
REDUCED_NOISE = PYTHON_STABILITY + ("rate", "samples")
CPP_PARITY = PYTHON_STABILITY + ("n_abs",)


def _missing(parameters: dict, names: tuple[str, ...]) -> list[str]:
    return [name for name in names if parameters.get(name) is None]


def _capability(parameters: dict, names: tuple[str, ...], reason: str) -> dict:
    missing = _missing(parameters, names)
    return {
        "ready": not missing,
        "required_parameters": list(names),
        "missing_parameters": missing,
        "blocking_reason": None if not missing else reason,
    }


def audit(case_dir: Path) -> dict:
    input_path = case_dir / "input.json"
    provenance_path = case_dir / "provenance.json"
    parameters = json.loads(input_path.read_text(encoding="utf-8"))
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    observed = provenance.get("independently_observed", {})
    calibration = provenance.get("parameter_status", {}).get("readout_calibration_A_per_V", {})
    acquisition_ready = all(
        key in observed for key in (
            "readout_rate_Hz", "readout_samples", "analysis_cutoff_Hz", "fft_bin_width_Hz"
        )
    )
    calibration_ready = calibration.get("value") is not None and calibration.get("status") in {
        "independently_measured", "independently_derived"
    }
    operating = _capability(parameters, OPERATING_POINT, "TES bath-power equilibrium cannot be evaluated.")
    stability = _capability(parameters, PYTHON_STABILITY, "Jacobian/stability requires the full reduced electrothermal parameter set.")
    reduced = _capability(parameters, REDUCED_NOISE, "MakeNoise() requires the reduced electrothermal set and FFT rate/sample settings.")
    cpp = _capability(parameters, CPP_PARITY, "C++ linearization inspection additionally requires n_abs; E is not required.")
    normalized = {
        "ready": reduced["ready"] and acquisition_ready,
        "required_parameters": list(REDUCED_NOISE) + ["target CH0 raw records"],
        "missing_parameters": reduced["missing_parameters"] + ([] if acquisition_ready else ["target acquisition metadata"]),
        "blocking_reason": None if reduced["ready"] and acquisition_ready else "A physical reduced spectrum and target acquisition grid are required; calibration is not required.",
        "normalization_is_not_calibration_evidence": True,
    }
    absolute_missing = [] if calibration_ready else ["readout_calibration_A_per_V"]
    absolute = {
        "ready": normalized["ready"] and not absolute_missing,
        "required_parameters": ["CH0 voltage/current calibration", "readout gain/transfer", *list(REDUCED_NOISE)],
        "missing_parameters": normalized["missing_parameters"] + absolute_missing,
        "blocking_reason": None if normalized["ready"] and not absolute_missing else "Absolute ASD comparison requires independently established CH0 calibration and readout gain.",
    }
    return {
        "status": "ready" if normalized["ready"] else "blocked",
        "runnable": False if not reduced["ready"] else True,
        "case_dir": case_dir.as_posix(),
        "source_commit": provenance.get("source_commit"),
        "capabilities": {
            "operating_point_ready": operating,
            "python_stability_ready": stability,
            "reduced_noise_ready": reduced,
            "cpp_parity_ready": cpp,
            "normalized_comparison_ready": normalized,
            "absolute_comparison_ready": absolute,
        },
        "noise_semantics": {
            "tes_johnson": "enabled physical source",
            "load_johnson": "enabled physical source",
            "tes_bath_tfn": "enabled physical source",
            "tes_absorber_tfn": "enabled physical source",
            "post_filter_white_asd_A_rtHz": parameters.get("post_filter_white_asd_A_rtHz", 0.0),
            "readout_white_asd_A_rtHz": parameters.get("readout_white_asd_A_rtHz", 0.0),
            "tes_resistance_fluctuation_model": parameters.get("tes_resistance_fluctuation_model", "none"),
            "tes_internal_model": parameters.get("tes_internal_model", "none"),
        },
        "generic_input_used": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit(args.case_dir), indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
