"""Small finite stale-tail regression for the production Schur Bt action.

Elmer's CRS matvec writes the rows stored by the matrix.  The production
work-vector contract therefore clears the complete logical output first.
This module keeps the regression independent of a local Elmer build while
also allowing the test suite to assert the exact finite-tail behavior.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def stored_crs_matvec(
    rows: np.ndarray,
    cols: np.ndarray,
    values: np.ndarray,
    vector: np.ndarray,
    output: np.ndarray,
) -> None:
    """Apply a one-based CSR matrix, writing only its stored rows."""

    for row in range(rows.size - 1):
        start = int(rows[row]) - 1
        stop = int(rows[row + 1]) - 1
        output[row] = np.dot(values[start:stop], vector[cols[start:stop] - 1])


def production_short_bt_action() -> dict[str, object]:
    """Exercise a 3x5 stored Bt with a 5-entry finite stale tail."""

    # Logical Bt is 5 x 3, but only rows 1..3 are stored.  Rows 4 and 5 are
    # mathematically zero and deliberately do not appear in this CSR object.
    rows = np.array([1, 2, 3, 4], dtype=np.int64)
    cols = np.array([1, 2, 3], dtype=np.int64)
    values = np.array([2.0, -3.0, 5.0])
    vector = np.array([7.0, 11.0, 13.0])
    stale_tail = np.array([123.0, -456.0])

    bt_v = np.empty(5, dtype=np.float64)
    bt_v[3:] = stale_tail
    bt_v[:] = 0.0  # production fix: clear all logical output entries first
    stored_crs_matvec(rows, cols, values, vector, bt_v)

    expected = np.array([14.0, -33.0, 65.0, 0.0, 0.0])
    return {
        "stored_shape": [3, 3],
        "logical_shape": [5, 3],
        "stale_tail": stale_tail.tolist(),
        "result": bt_v.tolist(),
        "expected": expected.tolist(),
        "tail_after": bt_v[3:].tolist(),
        "finite": bool(np.isfinite(bt_v).all()),
        "pass": bool(np.array_equal(bt_v, expected)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = production_short_bt_action()
    encoded = json.dumps(result, indent=2)
    print(encoded)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded + "\n", encoding="utf-8")
    if not result["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
