"""Compare isolated custom steady output with the retained standard reference."""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path


NUMBER = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?")


def state_values(path: Path) -> list[float]:
    return [float(item.replace("D", "E")) for item in NUMBER.findall(path.read_text())]


def temperatures(path: Path) -> list[float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    perm_line = next(i for i, line in enumerate(lines) if line.startswith("Perm:"))
    count = int(lines[perm_line].split()[1])
    start = perm_line + 1 + count
    values = [float(lines[start + i].strip().replace("D", "E")) for i in range(count)]
    return values


def comparison(reference: list[float], candidate: list[float]) -> dict:
    if len(reference) != len(candidate):
        raise ValueError(f"length mismatch: {len(reference)} != {len(candidate)}")
    delta = [b - a for a, b in zip(reference, candidate)]
    return {
        "count": len(reference),
        "reference_min": min(reference), "reference_max": max(reference),
        "candidate_min": min(candidate), "candidate_max": max(candidate),
        "max_abs_difference": max(abs(value) for value in delta),
        "l2_difference": math.sqrt(sum(value * value for value in delta)),
        "reference_l2": math.sqrt(sum(value * value for value in reference)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-result", type=Path, required=True)
    parser.add_argument("--candidate-result", type=Path, required=True)
    parser.add_argument("--reference-state", type=Path, required=True)
    parser.add_argument("--candidate-state", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    ref_state, candidate_state = state_values(args.reference_state), state_values(args.candidate_state)
    payload = {
        "reference_result": str(args.reference_result.resolve()),
        "candidate_result": str(args.candidate_result.resolve()),
        "state_reference": ref_state,
        "state_candidate": candidate_state,
        "state_difference": [b - a for a, b in zip(ref_state, candidate_state)],
        "temperature": comparison(temperatures(args.reference_result), temperatures(args.candidate_result)),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
