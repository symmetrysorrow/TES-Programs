"""Map HDF5 v2 physical deposits to one-dimensional TES absorber blocks.

This is a *pulse-model approximation*, not a PHITS prediction of the spatial
stopping-power density.  A SEGMENT's reconstructed energy is distributed
uniformly in path length along its straight transport chord.  The separate
PHITS score coordinate is not substituted for the segment geometry.

Events and deposits are streamed using the format-validated public readers;
incomplete/unresolved events or energy outside the simulated absorber fail
closed rather than producing a deceptively complete pulse.
"""

from __future__ import annotations

import math

import h5py
import numpy as np

from .event_hdf5 import FORMAT_NAME, iter_deposits, iter_events, require_legacy_point_events


SEGMENT_ALLOCATION_MODEL = "uniform_along_straight_chord_v1"


def event_format_version(path) -> int:
    """Return 2 for physical v2, 1 for a validated legacy point-only file."""
    with h5py.File(path, "r") as source:
        name = source.attrs.get("format", "")
        if isinstance(name, bytes):
            name = name.decode("utf-8")
        if name != FORMAT_NAME:
            raise ValueError(f"{path} is not a {FORMAT_NAME} event file")
        version = source.attrs.get("format_version")
        if version is not None and int(version) == 2:
            return 2
    require_legacy_point_events(path)
    return 1


def _coordinate(x: float, half_length_cm: float) -> float:
    x = float(x)
    if not math.isfinite(x):
        raise ValueError("deposit x coordinate is not finite")
    # Round-off at a PHITS/geometry boundary is not physical energy escape.
    tol = max(1e-12, half_length_cm * 1e-10)
    if x < -half_length_cm - tol or x > half_length_cm + tol:
        raise ValueError(
            f"deposit x={x} cm is outside absorber [-{half_length_cm}, {half_length_cm}] cm; "
            "check input.json length, coordinate origin, and PHITS geometry"
        )
    return min(half_length_cm, max(-half_length_cm, x))


def _index(x: float, half_length_cm: float, n_abs: int) -> int:
    return min(n_abs - 1, int((x + half_length_cm) * n_abs / (2.0 * half_length_cm)))


def allocate_deposit(energy_by_block, deposit, *, length_mm: float, mirror_x: bool = False):
    """Allocate one POINT or SEGMENT without losing energy at block boundaries.

    The total energy of an entirely in-absorber chord is conserved; a segment
    parallel to the absorber's x axis (including a zero-length chord) uses a
    single block.  A partly out-of-absorber chord is rejected, not clipped.
    """
    energy = deposit.get("energy_MeV")
    kind = deposit.get("kind")
    if kind not in ("POINT", "SEGMENT") or energy is None:
        raise ValueError(f"unresolved or unknown v2 deposit kind: {kind!r}")
    energy = float(energy)
    if not math.isfinite(energy) or energy < 0:
        raise ValueError("deposit energy must be finite and nonnegative")
    if energy == 0:
        return
    n_abs = len(energy_by_block)
    half_length_cm = float(length_mm) / 20.0
    if n_abs < 1 or not math.isfinite(half_length_cm) or half_length_cm <= 0:
        raise ValueError("absorber length and block count must be positive and finite")

    if kind == "POINT":
        xyz = deposit.get("score_xyz_cm")
        if xyz is None or len(xyz) != 3 or not all(math.isfinite(v) for v in xyz):
            raise ValueError("POINT deposit lacks a valid score position")
        x = -xyz[0] if mirror_x else xyz[0]
        energy_by_block[_index(_coordinate(x, half_length_cm), half_length_cm, n_abs)] += energy
        return

    start = deposit.get("start_xyz_cm")
    end = deposit.get("end_xyz_cm")
    if (start is None or end is None or len(start) != 3 or len(end) != 3
            or not all(math.isfinite(v) for v in (*start, *end))):
        raise ValueError("SEGMENT deposit lacks valid start/end coordinates")
    a = _coordinate(-start[0] if mirror_x else start[0], half_length_cm)
    b = _coordinate(-end[0] if mirror_x else end[0], half_length_cm)
    low, high = sorted((a, b))
    if low == high:  # parallel to x planes: the projected x is fixed
        energy_by_block[_index(low, half_length_cm, n_abs)] += energy
        return

    first = _index(low, half_length_cm, n_abs)
    last = _index(high, half_length_cm, n_abs)
    block_width = 2.0 * half_length_cm / n_abs
    portions = []
    for i in range(first, last + 1):
        left = -half_length_cm + i * block_width
        right = left + block_width
        overlap = max(0.0, min(high, right) - max(low, left))
        if overlap > 0:
            portions.append((i, energy * overlap / (high - low)))
    if not portions:
        raise ValueError("SEGMENT has no overlap with absorber blocks")
    # Put tiny floating-point remainder in the final intersected block so that
    # the reconstructed physical energy is conserved across all blocks.
    distributed = math.fsum(portion for _, portion in portions)
    if not math.isclose(distributed, energy, rel_tol=1e-9, abs_tol=1e-12):
        raise ValueError("SEGMENT block intersections do not conserve energy")
    for i, portion in portions:
        energy_by_block[i] += portion
    energy_by_block[portions[-1][0]] += energy - distributed


def iter_v2_block_energies(path, n_abs: int, length_mm: float, *, mirror_x: bool = False):
    """Yield (event_id, energies_MeV[n_abs]) for complete physical v2 events.

    The HDF5 event offsets are checked against the streamed deposit event IDs;
    no full-run array is materialized.  'all_absorption=undetermined' is allowed:
    a measured deposited-energy pulse does not imply full-energy absorption.
    """
    if int(n_abs) != n_abs or n_abs < 1:
        raise ValueError("n_abs must be a positive integer")
    if not math.isfinite(float(length_mm)) or float(length_mm) <= 0:
        raise ValueError("length_mm must be positive and finite")
    if event_format_version(path) != 2:
        raise ValueError("iter_v2_block_energies requires a HDF5 v2 file")
    deposits = iter(iter_deposits(path))
    for event_id, summary in iter_events(path):
        if not summary["complete"] or summary["unresolved_count"]:
            raise ValueError(f"event {event_id} is incomplete or has unresolved deposits")
        d0, d1 = summary["deposit_range"]
        if d0 < 0 or d1 < d0:
            raise ValueError(f"event {event_id} has an invalid deposit range")
        energy_by_block = np.zeros(n_abs, dtype=np.float64)
        raw_energies = []
        for _ in range(d0, d1):
            try:
                deposit_event_id, deposit = next(deposits)
            except StopIteration as exc:
                raise ValueError(f"event {event_id} is missing deposit rows") from exc
            if deposit_event_id != event_id:
                raise ValueError(f"event/deposit ID mismatch: {event_id} != {deposit_event_id}")
            if deposit["energy_MeV"] is None:
                raise ValueError(f"event {event_id} has unknown deposited energy")
            raw_energies.append(float(deposit["energy_MeV"]))
            try:
                allocate_deposit(energy_by_block, deposit, length_mm=length_mm, mirror_x=mirror_x)
            except ValueError as exc:
                raise ValueError(f"event {event_id} addition {deposit['addition_id']}: {exc}") from exc
        original = math.fsum(raw_energies)
        reported = summary["deposit_total_MeV"]
        if not math.isfinite(reported) or not math.isclose(original, reported, rel_tol=1e-9, abs_tol=1e-10):
            raise ValueError(f"event {event_id} has inconsistent deposited-energy total")
        if not math.isclose(float(energy_by_block.sum()), original, rel_tol=1e-9, abs_tol=1e-10):
            raise ValueError(f"event {event_id} lost energy during spatial allocation")
        yield event_id, energy_by_block
    if next(deposits, None) is not None:
        raise ValueError("v2 file has deposits not accounted for by event offsets")
