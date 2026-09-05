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

INDEPENDENT_STATUSES = {"independently_measured", "independently_derived"}
CONDITIONAL_STATUSES = {"setpoint_only", "same_device_prior"}


def _provenance_status(entry: dict | None) -> str:
    status = (entry or {}).get("status", "unresolved")
    if status in INDEPENDENT_STATUSES:
        return "independent_ready"
    if status in CONDITIONAL_STATUSES:
        return "conditional"
    return "not_admissible"


def _entry_for(parameters: dict, provenance: dict, name: str) -> dict:
    if name == "rate":
        return {"status": "independently_measured"} if parameters.get(name) is not None else {}
    if name == "samples":
        return {"status": "independently_measured"} if parameters.get(name) is not None else {}
    return provenance.get("parameter_status", {}).get(name, {})


def _capability(parameters: dict, provenance: dict, names: tuple[str, ...], reason: str) -> dict:
    missing = [name for name in names if parameters.get(name) is None]
    provenance_status = {
        name: _provenance_status(_entry_for(parameters, provenance, name))
        for name in names
        if name not in missing
    }
    inadmissible = [
        name for name, status in provenance_status.items() if status == "not_admissible"
    ]
    has_conditional = any(status == "conditional" for status in provenance_status.values())
    independent_ready = not missing and not inadmissible and not has_conditional
    return {
        "ready": independent_ready,
        "conditionally_ready": not missing and not inadmissible and has_conditional,
        "required_parameters": list(names),
        "missing_parameters": missing,
        "inadmissible_parameters": inadmissible,
        "provenance_status": provenance_status,
        "blocking_reason": (
            None
            if independent_ready
            else reason
            if missing
            else "One or more parameters have conditional or inadmissible provenance."
        ),
    }


def audit(case_dir: Path) -> dict:
    input_path = case_dir / "input.json"
    provenance_path = case_dir / "provenance.json"
    parameters = json.loads(input_path.read_text(encoding="utf-8"))
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    observed = provenance.get("independently_observed", {})
    calibration = provenance.get("parameter_status", {}).get("readout_calibration_A_per_V", {})
    transfer = provenance.get("parameter_status", {}).get("readout_gain_transfer", {})
    acquisition_ready = all(
        key in observed for key in (
            "readout_rate_Hz", "readout_samples", "analysis_cutoff_Hz", "fft_bin_width_Hz"
        )
    )
    operating = _capability(parameters, provenance, OPERATING_POINT, "TES bath-power equilibrium cannot be evaluated.")
    stability = _capability(parameters, provenance, PYTHON_STABILITY, "Jacobian/stability requires the full reduced electrothermal parameter set.")
    reduced = _capability(parameters, provenance, REDUCED_NOISE, "MakeNoise() requires the reduced electrothermal set and FFT rate/sample settings.")
    cpp = _capability(parameters, provenance, CPP_PARITY, "C++ linearization inspection additionally requires n_abs; E is not required.")
    acquisition_status = "independent_ready" if acquisition_ready else "not_admissible"
    normalized = {
        "ready": reduced["ready"] and acquisition_ready,
        "conditionally_ready": reduced["conditionally_ready"] and acquisition_ready,
        "required_parameters": list(REDUCED_NOISE) + ["target CH0 raw records"],
        "missing_parameters": reduced["missing_parameters"] + ([] if acquisition_ready else ["target acquisition metadata"]),
        "inadmissible_parameters": reduced["inadmissible_parameters"],
        "provenance_status": {**reduced["provenance_status"], "target CH0 raw records": acquisition_status},
        "blocking_reason": None if reduced["ready"] and acquisition_ready else "A physical reduced spectrum and target acquisition grid are required; calibration is not required.",
        "normalization_is_not_calibration_evidence": True,
    }
    absolute_missing = []
    if calibration.get("value") is None:
        absolute_missing.append("readout_calibration_A_per_V")
    if transfer.get("value") is None:
        absolute_missing.append("readout_gain_transfer")
    absolute_inadmissible = []
    for name, entry in (("readout_calibration_A_per_V", calibration), ("readout_gain_transfer", transfer)):
        if entry.get("value") is not None and _provenance_status(entry) == "not_admissible":
            absolute_inadmissible.append(name)
    absolute = {
        "ready": normalized["ready"] and not absolute_missing and not absolute_inadmissible,
        "conditionally_ready": normalized["conditionally_ready"] and not absolute_missing and not absolute_inadmissible,
        "required_parameters": ["CH0 voltage/current calibration", "readout gain/transfer", *list(REDUCED_NOISE)],
        "missing_parameters": normalized["missing_parameters"] + absolute_missing,
        "inadmissible_parameters": normalized["inadmissible_parameters"] + absolute_inadmissible,
        "provenance_status": {**normalized["provenance_status"], "readout_calibration_A_per_V": _provenance_status(calibration), "readout_gain_transfer": _provenance_status(transfer)},
        "blocking_reason": None if normalized["ready"] and not absolute_missing else "Absolute ASD comparison requires independently established CH0 calibration and readout gain.",
    }
    search = provenance.get("existing_data_search", {})
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
        "provenance_audit": {
            "detector_channel_linkage": provenance.get("detector_channel_linkage", {}).get("result", "unknown"),
            "existing_data_search": search.get("status", "not_recorded"),
            "target_physics_classification": search.get("conclusion", "not_recorded"),
            "exploratory_proxy_conclusion": provenance.get("exploratory_proxy_phase", {}).get("conclusion", "not_recorded"),
            "generic_input_used": False,
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
