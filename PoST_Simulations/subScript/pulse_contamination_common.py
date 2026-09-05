"""Pulse-only template and fixed-threshold utilities for contamination audits."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from scipy.signal import fftconvolve


RATE_HZ = 500000.0
SAMPLES = 100000
PRETRIGGER = 5000
FULL_START = -256
FULL_LENGTH = 4096
TAIL_LENGTH = 4096
SLOW_DECIMATION = 10
SLOW_LENGTH = 4096
MAX_PULSE_RECORDS = 300


def pulse_paths(target_root: Path, channel: str) -> list[Path]:
    paths = sorted((target_root / f"{channel}_pulse/rawdata").glob(f"{channel}_*.dat"), key=lambda p: int(p.stem.split("_")[1]))
    return paths[:MAX_PULSE_RECORDS]


def noise_paths(target_root: Path, channel: str) -> list[Path]:
    return sorted((target_root / f"{channel}_noise/rawdata").glob(f"{channel}_*.dat"), key=lambda p: int(p.stem.split("_")[1]))


def read_record(path: Path) -> np.ndarray:
    values = np.fromfile(path, dtype=np.float64, offset=4)
    if values.size != SAMPLES or not np.all(np.isfinite(values)):
        raise ValueError(f"invalid record: {path}")
    return values


def record_key(path: Path) -> str:
    return path.stem.split("_", 1)[1]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def orient_pulse(values: np.ndarray) -> tuple[np.ndarray, int, float, float, float]:
    baseline = float(np.mean(values[: PRETRIGGER - 500]))
    baseline_std = float(np.std(values[: PRETRIGGER - 500]))
    oriented = values - baseline
    post = oriented[PRETRIGGER:]
    if abs(float(np.min(post))) > abs(float(np.max(post))):
        oriented = -oriented
    peak = int(PRETRIGGER + np.argmax(oriented[PRETRIGGER:]))
    amplitude = float(oriented[peak])
    if peak <= PRETRIGGER or amplitude <= 10.0 * max(baseline_std, np.finfo(float).tiny):
        raise ValueError("pulse peak is not above the pulse-record baseline")
    return oriented, peak, amplitude, baseline, baseline_std


def normalized_pulse(values: np.ndarray) -> dict:
    oriented, peak, amplitude, baseline, baseline_std = orient_pulse(values)
    full = oriented[peak + FULL_START: peak + FULL_START + FULL_LENGTH]
    tail = oriented[peak: peak + TAIL_LENGTH]
    slow = oriented[peak: peak + SLOW_DECIMATION * SLOW_LENGTH: SLOW_DECIMATION]
    if len(full) != FULL_LENGTH or len(tail) != TAIL_LENGTH or len(slow) != SLOW_LENGTH:
        raise ValueError("pulse does not contain the required tail window")
    normalized = full / amplitude
    tail_normalized = tail / amplitude
    slow_normalized = slow / amplitude
    return {
        "peak_index": peak,
        "amplitude_raw": amplitude,
        "baseline_raw": baseline,
        "baseline_std_raw": baseline_std,
        "full_template": normalized,
        "tail_template": tail_normalized,
        "slow_tail_template": slow_normalized,
    }


def matched_score(record: np.ndarray, template: np.ndarray, step: int = 1) -> tuple[float, int, float]:
    """Return normalized correlation score, lag, and signed amplitude."""
    values = np.asarray(record, dtype=float)
    kernel = np.asarray(template, dtype=float)
    if step > 1:
        values = values[::step]
    kernel = kernel - np.mean(kernel)
    norm = float(np.linalg.norm(kernel))
    if norm == 0.0 or values.size < kernel.size:
        return 0.0, 0, 0.0
    correlation = fftconvolve(values, kernel[::-1], mode="valid")
    starts = np.arange(correlation.size)
    # Normalize by local energy to prevent a large broadband excursion from
    # looking like a template match.
    energy = np.sqrt(fftconvolve(values**2, np.ones(kernel.size), mode="valid"))
    scores = correlation / np.maximum(norm * energy, np.finfo(float).tiny)
    index = int(np.argmax(scores))
    segment = values[index:index + kernel.size]
    amplitude = float(np.dot(segment, kernel) / max(np.dot(kernel, kernel), np.finfo(float).tiny))
    return float(scores[index]), int(index * step), amplitude


def scan_templates(record: np.ndarray, templates: dict) -> dict:
    rows = []
    for name, row in templates.items():
        score, lag, amplitude = matched_score(record, np.asarray(row["values"], dtype=float), int(row.get("sample_step", 1)))
        rows.append({"template": name, "score": score, "lag_samples": lag, "lag_s": lag / RATE_HZ, "amplitude_normalized": amplitude})
    return max(rows, key=lambda row: row["score"]) if rows else {"template": None, "score": 0.0, "lag_samples": 0, "lag_s": 0.0, "amplitude_normalized": 0.0}
