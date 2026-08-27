"""Rank material parameters by their measured tail sensitivity.

The input is the aligned trace output of ``post_multivariate_search.py``.
This is deliberately a proxy analysis: the current completed traces are the
dual-TES PoST geometry, so the result is a screening order for the shared
material parameters, not a claim about the final single-pixel geometry.

Example::

    python scripts/analysis/analyze_tail_parameter_sensitivity.py
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SEARCH_ROOT = ROOT / "artifacts" / "search" / "post_multivariate"
DEFAULT_OUTPUT_ROOT = ROOT / "artifacts" / "analysis" / "tail_parameter_sensitivity"
WINDOWS_MS = ((5.0, 20.0), (20.0, 40.0), (40.0, 70.0))


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_candidates(path: Path) -> dict[str, dict[str, float]]:
    candidates: dict[str, dict[str, float]] = {}
    with path.open("r", newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file):
            candidates[str(row["candidate_id"])] = {
                str(name): float(value)
                for name, value in json.loads(row["factors"]).items()
            }
    return candidates


def read_trace(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    names = set(data.dtype.names or ())
    time_name = "time_after_pulse_ms"
    residual_name = "residual"
    if time_name not in names or residual_name not in names:
        raise ValueError(f"{path}: expected {time_name} and {residual_name}")
    time_ms = np.asarray(data[time_name], dtype=float)
    residual = np.asarray(data[residual_name], dtype=float)
    mask = np.isfinite(time_ms) & np.isfinite(residual)
    if np.count_nonzero(mask) < 3:
        raise ValueError(f"{path}: too few finite trace samples")
    return time_ms[mask], residual[mask]


def trace_window_rmse(path: Path, start_ms: float, end_ms: float) -> float:
    time_ms, residual = read_trace(path)
    mask = (time_ms >= start_ms) & (time_ms < end_ms)
    if np.count_nonzero(mask) == 0:
        raise ValueError(f"{path}: no samples in {start_ms:g}-{end_ms:g} ms")
    return float(np.sqrt(np.mean(residual[mask] ** 2)))


def candidate_window_scores(
    search_root: Path, candidate_id: str
) -> dict[str, float]:
    directory = search_root / "candidates" / candidate_id
    trace_paths = sorted(directory.glob("aligned_ch0_extrema_*.csv"))
    if not trace_paths:
        raise FileNotFoundError(f"{candidate_id}: no aligned ch0 extrema traces")
    scores: dict[str, float] = {}
    for start_ms, end_ms in WINDOWS_MS:
        values = [trace_window_rmse(path, start_ms, end_ms) for path in trace_paths]
        scores[f"rmse_{start_ms:g}_{end_ms:g}_ms"] = float(np.mean(values))
    return scores


def correlation(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) == 0.0 or np.std(y) == 0.0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def analyze(search_root: Path, output_root: Path) -> dict[str, Path]:
    factors_by_id = read_candidates(search_root / "candidates.csv")
    rows: list[dict[str, Any]] = []
    for candidate_id, factors in factors_by_id.items():
        try:
            scores = candidate_window_scores(search_root, candidate_id)
        except FileNotFoundError:
            continue
        rows.append({"candidate_id": candidate_id, **factors, **scores})
    if len(rows) < 3:
        raise RuntimeError(
            "tail sensitivity needs at least three completed candidates with traces"
        )

    variable_names = sorted(factors_by_id[next(iter(factors_by_id))])
    sensitivity_rows: list[dict[str, Any]] = []
    for start_ms, end_ms in WINDOWS_MS:
        metric = f"rmse_{start_ms:g}_{end_ms:g}_ms"
        y = np.asarray([float(row[metric]) for row in rows], dtype=float)
        for variable_name in variable_names:
            x = np.asarray(
                [math.log(float(row[variable_name])) for row in rows], dtype=float
            )
            design = np.column_stack([np.ones(len(x)), x])
            coefficient = float(np.linalg.lstsq(design, y, rcond=None)[0][1])
            corr = correlation(x, y)
            sensitivity_rows.append(
                {
                    "window_ms": f"{start_ms:g}-{end_ms:g}",
                    "metric": metric,
                    "parameter": variable_name,
                    "coefficient_per_log_factor": coefficient,
                    "pearson_r": corr,
                    "abs_pearson_r": abs(corr) if math.isfinite(corr) else float("nan"),
                }
            )

    for window in {row["window_ms"] for row in sensitivity_rows}:
        subset = [row for row in sensitivity_rows if row["window_ms"] == window]
        subset.sort(key=lambda row: row["abs_pearson_r"], reverse=True)
        for rank, row in enumerate(subset, start=1):
            row["rank_by_abs_pearson_r"] = rank

    output_root.mkdir(parents=True, exist_ok=True)
    scores_path = output_root / "tail_candidate_scores.csv"
    with scores_path.open("w", newline="", encoding="utf-8") as file:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    sensitivity_path = output_root / "tail_parameter_sensitivity.csv"
    with sensitivity_path.open("w", newline="", encoding="utf-8") as file:
        fieldnames = list(sensitivity_rows[0].keys())
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sensitivity_rows)

    report = {
        "search_root": str(search_root.relative_to(ROOT)),
        "completed_candidates": [row["candidate_id"] for row in rows],
        "baseline_factor_vector": {name: 1.0 for name in variable_names},
        "baseline_trace_available": False,
        "interpretation": (
            "Proxy ranking from completed dual-TES PoST traces. The factor=1 baseline "
            "is the physical reference point but is not present as a completed tail trace. "
            "Use the ranking to choose single-pixel OFAT runs, then re-fit against the "
            "actual single-pixel baseline."
        ),
        "windows_ms": [list(window) for window in WINDOWS_MS],
        "candidate_score_file": str(scores_path.relative_to(ROOT)),
        "sensitivity_file": str(sensitivity_path.relative_to(ROOT)),
        "top_by_window": {},
    }
    for start_ms, end_ms in WINDOWS_MS:
        window = f"{start_ms:g}-{end_ms:g}"
        subset = [row for row in sensitivity_rows if row["window_ms"] == window]
        subset.sort(key=lambda row: row["abs_pearson_r"], reverse=True)
        report["top_by_window"][window] = [
            {
                "parameter": row["parameter"],
                "pearson_r": row["pearson_r"],
                "coefficient_per_log_factor": row["coefficient_per_log_factor"],
            }
            for row in subset
        ]
    report_path = output_root / "tail_parameter_sensitivity.json"
    report_path.write_text(
        json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return {
        "candidate_scores": scores_path,
        "sensitivity": sensitivity_path,
        "report": report_path,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--search-root", type=Path, default=DEFAULT_SEARCH_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    paths = analyze(
        args.search_root if args.search_root.is_absolute() else ROOT / args.search_root,
        args.output_root if args.output_root.is_absolute() else ROOT / args.output_root,
    )
    for name, path in paths.items():
        print(f"{name}: {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
