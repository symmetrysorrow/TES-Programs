"""Compare Elmer ``Linear System Save`` matrix/RHS dumps.

The dump is intentionally produced at ``Linear System Save Slot = linear
solve``: this is the matrix actually passed to MUMPS/HYPRE after mortar
restriction, multiplier augmentation, or penalty/elimination.
"""
from __future__ import annotations

import argparse
import hashlib
import math
from pathlib import Path


def records(path: Path, columns: int) -> list[tuple[int, ...] | tuple[int, int, float]]:
    result: list[tuple[int, ...] | tuple[int, int, float]] = []
    for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        fields = line.split()
        if len(fields) < columns:
            raise ValueError(f"{path}:{line_number}: expected {columns} fields")
        if columns == 3:
            result.append((int(fields[0]), int(fields[1]), float(fields[2])))
        else:
            result.append((int(fields[0]), float(fields[1])))
    return result


def signature(items: list[tuple[object, ...]]) -> str:
    digest = hashlib.sha256()
    for item in items:
        digest.update(repr(item).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def compare(left: Path, right: Path, columns: int) -> dict[str, object]:
    a = records(left, columns)
    b = records(right, columns)
    if columns == 3:
        am = {(int(i), int(j)): float(v) for i, j, v in a}
        bm = {(int(i), int(j)): float(v) for i, j, v in b}
    else:
        am = {int(i): float(v) for i, v in a}
        bm = {int(i): float(v) for i, v in b}
    keys = set(am) | set(bm)
    max_abs = max((abs(am.get(k, 0.0) - bm.get(k, 0.0)) for k in keys), default=0.0)
    mismatch = sum(not math.isclose(am.get(k, 0.0), bm.get(k, 0.0), rel_tol=0.0, abs_tol=0.0) for k in keys)
    return {
        "left_records": len(a),
        "right_records": len(b),
        "union_records": len(keys),
        "exact_mismatches": mismatch,
        "max_abs_difference": max_abs,
        "left_sha256": signature(a),
        "right_sha256": signature(b),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left_prefix", type=Path)
    parser.add_argument("right_prefix", type=Path)
    args = parser.parse_args()
    for suffix, columns in (("_a.dat", 3), ("_b.dat", 2)):
        report = compare(args.left_prefix.with_name(args.left_prefix.name + suffix), args.right_prefix.with_name(args.right_prefix.name + suffix), columns)
        print(f"{suffix}: {report}")


if __name__ == "__main__":
    main()
