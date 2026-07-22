"""Extract the converged TES current from an inner-circuit Elmer solver log."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


TIME_RE = re.compile(r"MAIN: Time:\s+(\d+)/(\d+):\s+([0-9.Ee+-]+)")
CURRENT_RE = re.compile(r"TESInnerCircuit:\s+T=\s*([0-9.Ee+-]+)\s+I=\s*([0-9.Ee+-]+)")


def extract(log_path: Path) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    step: tuple[int, int, float] | None = None
    currents: list[float] = []
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        tm = TIME_RE.search(line)
        if tm:
            if step is not None and currents:
                rows.append({"step": step[0], "time_s": step[2], "tes_current_A": currents[-1]})
            step = (int(tm.group(1)), int(tm.group(2)), float(tm.group(3)))
            currents = []
            continue
        cm = CURRENT_RE.search(line)
        if cm and step is not None:
            currents.append(float(cm.group(2)))
    if step is not None and currents:
        rows.append({"step": step[0], "time_s": step[2], "tes_current_A": currents[-1]})
    if not rows:
        raise RuntimeError(f"No MAIN/TESInnerCircuit records found in {log_path}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    rows = extract(args.log)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["step", "time_s", "tes_current_A"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Extracted {len(rows)} timesteps -> {args.output}")


if __name__ == "__main__":
    main()
