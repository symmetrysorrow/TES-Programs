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
    accepted_noise_indices,
    estimate_one_sided_asd,
)


POINTS_HZ = (10, 100, 1000, 3000, 5000, 7000, 10000)


def build_report(experiment_path: Path, case_dir: Path) -> dict:
    config = json.loads((experiment_path / "PulseConfig.json").read_text())
    rate = float(config["Readout"]["Rate"])
    samples = int(config["Readout"]["Sample"])
    cutoff = float(config["Analysis"]["CutoffFrequency"])
    rawdata_path = experiment_path / "CH0_noise" / "rawdata"

    record_paths = sorted(rawdata_path.glob("CH0_*.dat"))

    def records(paths=record_paths):
        for path in paths:
            yield np.frombuffer(path.read_bytes()[4:], dtype=np.float64).copy()

    def range_ok(values):
        return values.max() - values.min() <= 0.04

    accepted_indices = accepted_noise_indices(
        records(),
        samples,
        rate,
        cutoff=cutoff,
        remove_mean=True,
        accept_raw=range_ok,
        accept_processed=range_ok,
    )
    accepted_paths = [record_paths[index] for index in accepted_indices]
    spectrum, accepted_records = estimate_one_sided_asd(
        records(accepted_paths), samples, rate, cutoff=cutoff, remove_mean=True
    )
    pre_spectrum, pre_accepted_records = estimate_one_sided_asd(
        records(accepted_paths), samples, rate, cutoff=0.0, remove_mean=True
    )
    frequencies = np.arange(spectrum.size, dtype=float) * rate / samples
    reference_index = int(np.argmin(np.abs(frequencies - 1000.0)))
    reference = float(spectrum[reference_index])
    pre_reference = float(pre_spectrum[reference_index])
    rows = []
    for point in POINTS_HZ:
        index = int(np.argmin(np.abs(frequencies - point)))
        experimental = float(spectrum[index])
        rows.append(
            {
                "frequency_Hz": float(frequencies[index]),
                "experimental_post_analysis_asd": experimental,
                "experimental_post_analysis_normalized": experimental / reference,
                "experimental_pre_analysis_asd": float(pre_spectrum[index]),
                "experimental_pre_analysis_normalized": float(pre_spectrum[index] / pre_reference),
                "simulation_post_analysis_asd": None,
                "simulation_post_analysis_normalized": None,
                "simulation_over_experiment": None,
            }
        )
    return {
        "status": "blocked_unresolved_physics_parameters",
        "comparison_kind": "target_case_intrinsic_physical_noise_only",
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
            "accepted_record_indices": [int(index) for index in accepted_indices],
            "experimental_units": "raw CH0 voltage ASD units per sqrt(Hz); no target voltage-to-current calibration applied",
        },
        "normalization": "experimental and simulation post-analysis columns would be normalized independently to the 1 kHz bin; this is shape-only and is not calibration evidence",
        "spectrum_semantics": {
            "experimental_post_analysis_asd": "shared estimator after per-record mean removal, 10 kHz digital Bessel filtfilt, Hann power average, and one-sided ASD normalization",
            "experimental_pre_analysis_asd": "same accepted raw CH0 records after per-record mean removal, Hann power average, and one-sided ASD normalization; no inverse filter",
            "experimental_pre_analysis_normalized": "normalized to the 1 kHz pre-analysis bin independently from post-analysis",
            "pre_analysis_status": "directly reconstructed from the same accepted record mask",
            "simulation_post_analysis_asd": "would include intrinsic TES Johnson, load Johnson, TES-bath TFN and TES-absorber TFN, with empirical white/readout floors set to zero",
        },
        "points": rows,
        "pre_analysis": {
            "accepted_records": int(pre_accepted_records),
            "reference_frequency_Hz": float(frequencies[reference_index]),
            "normalization_reference_asd": pre_reference,
            "filter": "none; direct raw-record analysis",
        },
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
