"""Audit whether a versioned target case has enough independent physics data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


REQUIRED_PHYSICS = (
    "T_c", "T_bath", "R", "R_l", "alpha", "beta", "L", "n",
    "C_abs", "C_tes", "G_abs-abs", "G_abs-tes", "G_tes-bath", "n_abs", "E",
)


def audit(case_dir: Path) -> dict:
    input_path = case_dir / "input.json"
    provenance_path = case_dir / "provenance.json"
    parameters = json.loads(input_path.read_text(encoding="utf-8"))
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    missing = [name for name in REQUIRED_PHYSICS if parameters.get(name) is None]
    return {
        "status": "blocked_unresolved_physics_parameters" if missing else "ready_for_physics_audit",
        "runnable": not missing,
        "case_dir": case_dir.as_posix(),
        "source_commit": provenance.get("source_commit"),
        "required_physics_parameters": list(REQUIRED_PHYSICS),
        "missing_or_unresolved": missing,
        "independent_acquisition_metadata_present": all(
            key in provenance.get("independently_observed", {})
            for key in ("readout_rate_Hz", "readout_samples", "analysis_cutoff_Hz", "fft_bin_width_Hz")
        ),
        "operating_point_jacobian_stability": (
            "not run: target physics values are unresolved" if missing else "run target_case_audit and noise_model_audit"
        ),
        "comparison_policy": {
            "fit_parameters_to_noise": False,
            "use_1kHz_normalization_as_calibration": False,
            "absolute_comparison_without_independent_calibration": False,
        },
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
