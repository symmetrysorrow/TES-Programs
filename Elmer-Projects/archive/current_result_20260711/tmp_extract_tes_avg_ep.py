from pathlib import Path
import sys

import numpy as np


DEFAULT_EP_PATH = Path("mesh_fragment_centered_test/case_constant_power_fragment_test.ep")


def main() -> None:
    ep_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_EP_PATH
    lines = ep_path.read_text(encoding="utf-8", errors="replace").splitlines()
    header = lines[0].split()
    node_count = int(header[0])
    element_count = int(header[1])

    group_line_index = 1 + 1 + node_count
    if lines[group_line_index] != "#group all":
        raise RuntimeError(f"Unexpected group marker: {lines[group_line_index]!r}")

    element_start = group_line_index + 1
    element_end = element_start + element_count
    value_lines = lines[element_end:]
    time_markers = [i for i, line in enumerate(value_lines) if line.startswith("#time")]
    if not time_markers:
        raise RuntimeError("No time marker found in ep file")

    completed_steps: list[list[str]] = []
    for marker_index in time_markers:
        block = value_lines[marker_index + 1 : marker_index + 1 + node_count]
        if len(block) == node_count and all(block):
            completed_steps.append(block)

    if not completed_steps:
        raise RuntimeError("No completed output step found in ep file")

    step_count = len(completed_steps)
    last_step = completed_steps[-1]

    tes_nodes: set[int] = set()
    for line in lines[element_start:element_end]:
        parts = line.split()
        if parts and parts[0].lower() == "tes":
            tes_nodes.update(int(node_id) - 1 for node_id in parts[2:])

    temperature = np.fromiter(
        (float(line.split()[0]) for line in last_step),
        dtype=float,
        count=node_count,
    )
    tes_idx = np.array(sorted(tes_nodes), dtype=int)
    tes_temperature = temperature[tes_idx]

    print(f"node_count={node_count}")
    print(f"element_count={element_count}")
    print(f"step_count={step_count}")
    next_expected_lines = time_markers[-1] + 1 + node_count
    print(f"truncated_tail_lines={max(len(value_lines) - next_expected_lines, 0)}")
    print(f"tes_node_count={len(tes_idx)}")
    print(f"temperature_min={temperature.min():.15f}")
    print(f"temperature_max={temperature.max():.15f}")
    print(f"tes_mean={tes_temperature.mean():.15f}")
    print(f"tes_min={tes_temperature.min():.15f}")
    print(f"tes_max={tes_temperature.max():.15f}")


if __name__ == "__main__":
    main()
