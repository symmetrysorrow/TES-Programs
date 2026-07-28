from pathlib import Path

from scripts.prep.partition_elmer_restart import map_serial_restart, result_layout


def _result(path: Path, values: list[float]) -> None:
    count = len(values)
    lines = [
        f" Number Of Nodes: {count}",
        "Perm:",
        *[f"{node} {node}" for node in range(1, count + 1)],
        *[str(value) for value in values],
    ]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")


def _mapped_values(path: Path) -> list[float]:
    lines = path.read_text(encoding="ascii").splitlines()
    count, value_start, _, _ = result_layout(lines)
    return [float(value) for value in lines[value_start : value_start + count]]


def _mesh(tmp_path: Path) -> tuple[Path, Path]:
    serial = tmp_path / "serial.result"
    _result(serial, [10.0, 20.0])
    mesh = tmp_path / "mesh"
    parts = mesh / "partitioning.2"
    parts.mkdir(parents=True)
    _result(mesh / "template.0", [0.0])
    _result(mesh / "template.1", [0.0])
    # Rank zero owns global node 2; rank one owns global node 1.  This catches
    # both suffix rotation and accidental use of part.0.nodes.
    (parts / "part.1.nodes").write_text("2 0 0 0\n", encoding="ascii")
    (parts / "part.2.nodes").write_text("1 0 0 0\n", encoding="ascii")
    return serial, mesh


def test_mpi_zero_based_outputs_use_one_based_partition_nodes(tmp_path: Path) -> None:
    serial, mesh = _mesh(tmp_path)
    outputs = map_serial_restart(
        serial, mesh, "template", "mapped", partitions=2, rank_suffix_base=0
    )
    assert [path.name for path in outputs] == ["mapped.0", "mapped.1"]
    assert _mapped_values(mesh / "mapped.0") == [20.0]
    assert _mapped_values(mesh / "mapped.1") == [10.0]


def test_legacy_one_based_output_suffix_is_preserved(tmp_path: Path) -> None:
    serial, mesh = _mesh(tmp_path)
    outputs = map_serial_restart(
        serial, mesh, "template", "mapped", partitions=2, rank_suffix_base=1
    )
    assert [path.name for path in outputs] == ["mapped.1", "mapped.2"]
    assert all(path.is_file() for path in outputs)
