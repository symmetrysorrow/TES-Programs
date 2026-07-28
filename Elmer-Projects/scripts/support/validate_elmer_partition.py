"""Validate rank counts and hybrid element conservation after ElmerGrid partitioning."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def header(path: Path) -> tuple[list[int], dict[int, int]]:
    lines = [line.split() for line in path.read_text(encoding="ascii").splitlines() if line.split()]
    totals = [int(value) for value in lines[0][:3]]
    counts: dict[int, int] = {}
    for fields in lines[2:]:
        if len(fields) >= 2:
            counts[int(fields[0])] = int(fields[1])
    return totals, counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("serial_mesh", type=Path)
    parser.add_argument("partition_mesh", type=Path)
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    serial_totals, serial_types = header(args.serial_mesh / "mesh.header")
    ranks = []
    type_sums: dict[int, int] = {}
    for rank in range(1, args.ranks + 1):
        totals, types = header(args.partition_mesh / f"partitioning.{args.ranks}" / f"part.{rank}.header")
        ranks.append({"rank": rank - 1, "nodes": totals[0], "elements": totals[1], "boundaries": totals[2], "types": types})
        for key, value in types.items():
            if key in (504, 706):
                type_sums[key] = type_sums.get(key, 0) + value
    payload = {
        "serial": {"nodes": serial_totals[0], "elements": serial_totals[1], "boundaries": serial_totals[2], "types": serial_types},
        "ranks": ranks,
        "hybrid_type_sums": type_sums,
        "hybrid_types_conserved": all(type_sums.get(kind) == serial_types.get(kind) for kind in (504, 706)),
        "volume_elements_conserved": sum(type_sums.get(kind, 0) for kind in (504, 706)) == serial_totals[1],
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    if not payload["hybrid_types_conserved"] or not payload["volume_elements_conserved"]:
        raise SystemExit("hybrid element counts were not conserved")


if __name__ == "__main__":
    main()
