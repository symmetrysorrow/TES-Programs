"""Write a small, non-fitting target-case comparison summary.

The report is intentionally able to finish with null simulation columns.  It
must not turn a directory label or the generic simulation input into a target
physical operating point.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))
sys.path.insert(0, str(REPOSITORY_ROOT / "PoST_Simulations"))
from Analyze_Experimental_Data.tes_analysis.noise_utils import (  # noqa: E402
    estimate_one_sided_asd,
)


POINTS_HZ = (10, 100, 1000, 3000, 5000, 7000, 10000)


def build_report(experiment_path: Path, case_dir: Path) -> dict:
    config = json.loads((experiment_path / "PulseConfig.json").read_text())
    rate = float(config["Readout"]["Rate"])
    samples = int(config["Readout"]["Sample"])
    cutoff = float(config["Analysis"]["CutoffFrequency"])
    rawdata_path = experiment_path / "CH0_noise" / "rawdata"

    def records():
        for path in sorted(rawdata_path.glob("CH0_*.dat")):
            yield np.frombuffer(path.read_bytes()[4:], dtype=np.float64).copy()

    def range_ok(values):
        return values.max() - values.min() <= 0.04

    spectrum, accepted_records = estimate_one_sided_asd(
        records(),
        samples,
        rate,
        cutoff=cutoff,
        remove_mean=True,
        accept_raw=range_ok,
        accept_processed=range_ok,
    )
    frequencies = np.arange(spectrum.size, dtype=float) * rate / samples
    reference_index = int(np.argmin(np.abs(frequencies - 1000.0)))
    reference = float(spectrum[reference_index])
    rows = []
    for point in POINTS_HZ:
        index = int(np.argmin(np.abs(frequencies - point)))
        experimental = float(spectrum[index])
        rows.append(
            {
                "frequency_Hz": float(frequencies[index]),
                "experimental_asd_native_units": experimental,
                "experimental_normalized_to_1kHz": experimental / reference,
                "simulation_asd_native_units": None,
                "simulation_normalized_to_1kHz": None,
                "simulation_over_experiment": None,
            }
        )
    return {
        "status": "blocked_unresolved_physics_parameters",
        "comparison_kind": "target_case_no_added_noise",
        "absolute_comparison_allowed": False,
        "reason": "No independently sourced target operating point and readout calibration are available.",
        "experiment_path": experiment_path.as_posix(),
        "case_input": (case_dir / "input.json").as_posix(),
        "source_files": [
            (experiment_path / "PulseConfig.json").as_posix(),
            (experiment_path / "Setting.txt").as_posix(),
            (experiment_path / "CH0_noise" / "rawdata").as_posix(),
        ],
        "acquisition": {
            "rate_Hz": rate,
            "samples": samples,
            "cutoff_Hz": cutoff,
            "fft_bin_width_Hz": rate / samples,
            "channel": "CH0",
            "accepted_records": int(accepted_records),
            "experimental_units": "raw CH0 voltage ASD units per sqrt(Hz); no target voltage-to-current calibration applied",
        },
        "normalization": "experimental and simulation columns would be normalized independently to the 1 kHz bin; this is shape-only and is not calibration evidence",
        "points": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiment-path", type=Path, required=True)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report(args.experiment_path, args.case_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
