"""Summarize and compare canonical TES pulse waveforms."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [{key: float(value) for key, value in row.items() if value} for row in csv.DictReader(handle)]


def metrics(rows: list[dict[str, float]]) -> dict[str, float | int | None]:
    currents = [row["tes_current_A"] for row in rows]
    times = [row["time_s"] for row in rows]
    baseline = currents[0]
    peak = max(currents)
    trough = min(currents)
    amplitude = peak - trough
    def crossing(fraction: float) -> float | None:
        target = trough + fraction * amplitude
        for index in range(1, len(rows)):
            before = currents[index - 1] - target
            after = currents[index] - target
            if before == 0:
                return times[index - 1]
            if before * after <= 0 and currents[index] != currents[index - 1]:
                ratio = (target - currents[index - 1]) / (currents[index] - currents[index - 1])
                return times[index - 1] + ratio * (times[index] - times[index - 1])
        return None
    return {
        "row_count": len(rows),
        "baseline_current_A": baseline,
        "min_current_A": trough,
        "max_current_A": peak,
        "current_amplitude_A": amplitude,
        "relative_amplitude": amplitude / max(abs(baseline), 1.0e-300),
        "t10_s": crossing(0.1),
        "t50_s": crossing(0.5),
        "t90_s": crossing(0.9),
        "final_current_A": currents[-1],
        "final_temperature_K": rows[-1]["tes_temperature_K"],
    }


def compare(left: list[dict[str, float]], right: list[dict[str, float]]) -> dict[str, object]:
    right_by_step = {int(row["time_step"]): row for row in right}
    pairs = [(row, right_by_step[int(row["time_step"])]) for row in left if int(row["time_step"]) in right_by_step]
    report: dict[str, object] = {"status": "PASS" if pairs else "NOT_AVAILABLE", "common_timesteps": len(pairs)}
    if not pairs:
        return report
    for field in ("tes_current_A", "tes_temperature_K", "tes_resistance_ohm", "tes_power_W"):
        errors = [abs(a[field] - b[field]) for a, b in pairs]
        report[field] = {"max_absolute_difference": max(errors), "rmse": math.sqrt(sum(value * value for value in errors) / len(errors))}
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--waveform", action="append", required=True, help="label=CSV")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reference = read(args.reference)
    report: dict[str, object] = {"reference": str(args.reference.resolve()), "reference_metrics": metrics(reference), "waveforms": {}}
    for item in args.waveform:
        label, raw_path = item.split("=", 1)
        rows = read(Path(raw_path))
        report["waveforms"][label] = {"path": str(Path(raw_path).resolve()), "metrics": metrics(rows), "vs_reference": compare(reference, rows)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
