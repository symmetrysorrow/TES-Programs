"""Map a serial single-temperature Elmer restart onto an existing partition."""

from __future__ import annotations

import argparse
from pathlib import Path


def result_layout(lines: list[str]) -> tuple[int, int, list[int], list[float]]:
    node_line = next(line for line in lines if line.startswith(" Number Of Nodes:"))
    node_count = int(node_line.split()[-1])
    perm_header = next(i for i, line in enumerate(lines) if line.startswith("Perm:"))
    perm_start = perm_header + 1
    value_start = perm_start + node_count
    permutation = [0] * node_count
    for line in lines[perm_start:value_start]:
        node, position = (int(value) for value in line.split())
        permutation[node - 1] = position - 1
    values = [float(line) for line in lines[value_start : value_start + node_count]]
    return node_count, value_start, permutation, values


def global_temperature(serial_result: Path) -> list[float]:
    lines = serial_result.read_text(encoding="ascii").splitlines()
    _, _, _, values = result_layout(lines)
    # Elmer's ASCII result values are written in mesh-node order.  ``Perm``
    # describes the solver's internal ordering and is metadata for restoring
    # the variable; applying it a second time scrambles the physical field.
    return values


def partition_node_ids(path: Path) -> list[int]:
    return [int(line.split()[0]) for line in path.read_text(encoding="ascii").splitlines()]


def write_partition(
    template: Path,
    node_file: Path,
    output: Path,
    temperature: list[float],
) -> None:
    lines = template.read_text(encoding="ascii").splitlines()
    node_count, value_start, _, _ = result_layout(lines)
    global_ids = partition_node_ids(node_file)
    if len(global_ids) != node_count:
        raise ValueError(
            f"{node_file}: {len(global_ids)} nodes, template expects {node_count}"
        )
    partition_values = [temperature[global_id - 1] for global_id in global_ids]
    replacement = [f"  {value:.17g}" for value in partition_values]
    output.write_text(
        "\n".join(lines[:value_start] + replacement + lines[value_start + node_count :])
        + "\n",
        encoding="ascii",
    )


def map_serial_restart(
    serial_result: Path,
    mesh_dir: Path,
    template_base: str,
    output_base: str,
    partitions: int = 4,
    rank_suffix_base: int = 1,
) -> list[Path]:
    """Write rank-local restart files and return their paths.

    ElmerGrid numbers ``part.*.nodes`` from one, while MPI result files are
    normally numbered from zero.  Keeping these two conventions explicit
    prevents an otherwise silent rank rotation.
    """
    if rank_suffix_base not in (0, 1):
        raise ValueError("rank_suffix_base must be 0 or 1")
    temperature = global_temperature(serial_result)
    partition_dir = mesh_dir / f"partitioning.{partitions}"
    outputs: list[Path] = []
    for rank in range(partitions):
        output = mesh_dir / f"{output_base}.{rank + rank_suffix_base}"
        write_partition(
            mesh_dir / f"{template_base}.{rank}",
            partition_dir / f"part.{rank + 1}.nodes",
            output,
            temperature,
        )
        outputs.append(output)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("serial_result", type=Path)
    parser.add_argument("mesh_dir", type=Path)
    parser.add_argument("template_base", help="partitioned template result basename")
    parser.add_argument("output_base", help="output result basename")
    parser.add_argument("--partitions", type=int, default=4)
    parser.add_argument("--rank-suffix-base", type=int, default=1, choices=(0, 1),
                        help="output result suffix base: 1 preserves legacy .1..N; 0 writes MPI .0..N-1")
    args = parser.parse_args()

    map_serial_restart(
        args.serial_result,
        args.mesh_dir,
        args.template_base,
        args.output_base,
        args.partitions,
        args.rank_suffix_base,
    )
    print(f"Wrote {args.partitions} partitioned restart files in {args.mesh_dir}")


if __name__ == "__main__":
    main()
