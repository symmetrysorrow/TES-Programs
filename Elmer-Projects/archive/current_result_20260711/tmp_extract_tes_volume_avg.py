from pathlib import Path
import sys

import meshio
import numpy as np


def main() -> None:
    mesh_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("mesh_shifted_merged")
    vtu_path = Path(sys.argv[2]) if len(sys.argv) > 2 else mesh_dir / "case_constant_power_t0001.vtu"

    mesh = meshio.read(vtu_path)
    temperature = mesh.point_data["temperature"]

    nodes: dict[int, np.ndarray] = {}
    with (mesh_dir / "mesh.nodes").open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            nodes[int(parts[0])] = np.array(
                [float(parts[2]), float(parts[3]), float(parts[4])],
                dtype=float,
            )

    tes_nodes: set[int] = set()
    tes_element_count = 0
    volume_sum = 0.0
    temp_volume_sum = 0.0

    with (mesh_dir / "mesh.elements").open("r", encoding="utf-8") as f:
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

    print(f"mesh_dir={mesh_dir}")
    print(f"vtu_path={vtu_path}")
    print(f"tes_element_count={tes_element_count}")
    print(f"tes_node_count={len(tes_nodes)}")
    print(f"tes_nodal_average={tes_nodal_average:.15f}")
    print(f"tes_volume_average={tes_volume_average:.15f}")
    print(f"tes_volume_sum={volume_sum:.18e}")


if __name__ == "__main__":
    main()
