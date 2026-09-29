import struct
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
VTU_PATH = ROOT / "legacy" / "meshes" / "mesh_shifted" / "case_constant_power_shifted_t0001.vtu"
NODES_PATH = ROOT / "legacy" / "meshes" / "mesh_shifted" / "mesh.nodes"
ELEMENTS_PATH = ROOT / "legacy" / "meshes" / "mesh_shifted" / "mesh.elements"


def load_temperature_array() -> np.ndarray:
    raw = VTU_PATH.read_bytes()
    marker = raw.index(b"<AppendedData")
    header = raw[:marker].decode("utf-8", errors="ignore")
    line = next(line for line in header.splitlines() if 'Name="temperature"' in line)
    offset = int(line.split('offset="', 1)[1].split('"', 1)[0])
    start = raw.index(b"_", marker) + 1
    block_start = start + offset
    nbytes = struct.unpack("<I", raw[block_start:block_start + 4])[0]
    return np.frombuffer(raw, dtype="<f8", count=nbytes // 8, offset=block_start + 4).copy()


def load_nodes() -> dict[int, np.ndarray]:
    nodes: dict[int, np.ndarray] = {}
    with NODES_PATH.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            nodes[int(parts[0])] = np.array(
                [float(parts[2]), float(parts[3]), float(parts[4])],
                dtype=float,
            )
    return nodes


def main() -> None:
    temperature = load_temperature_array()
    nodes = load_nodes()

    tes_nodes: set[int] = set()
    tes_element_count = 0
    volume_sum = 0.0
    temp_volume_sum = 0.0

    with ELEMENTS_PATH.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            body = int(parts[1])
            etype = int(parts[2])
            if body != 101:
                continue
            conn = [int(x) for x in parts[3:]]
            tes_nodes.update(conn)
            tes_element_count += 1
            if etype != 504 or len(conn) != 4:
                raise RuntimeError(f"Unexpected TES element type: {etype} ({len(conn)} nodes)")
            p0, p1, p2, p3 = (nodes[nid] for nid in conn)
            volume = abs(np.linalg.det(np.column_stack((p1 - p0, p2 - p0, p3 - p0)))) / 6.0
            tavg = float(np.mean([temperature[nid - 1] for nid in conn]))
            volume_sum += volume
            temp_volume_sum += volume * tavg

    tes_nodal_average = float(np.mean([temperature[nid - 1] for nid in sorted(tes_nodes)]))
    tes_volume_average = temp_volume_sum / volume_sum

    print(f"tes_element_count={tes_element_count}")
    print(f"tes_node_count={len(tes_nodes)}")
    print(f"temperature_min={temperature.min():.15f}")
    print(f"temperature_max={temperature.max():.15f}")
    print(f"tes_nodal_average={tes_nodal_average:.15f}")
    print(f"tes_volume_average={tes_volume_average:.15f}")
    print(f"tes_volume_sum={volume_sum:.18e}")


if __name__ == "__main__":
    main()
