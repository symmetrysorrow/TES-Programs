from pathlib import Path

import meshio
import numpy as np


I_BIAS = 7.15e-4
R_SH = 3.9e-3
R0 = 1.5527e-2
R_MIN = 1.0e-6
ALPHA = 256.46
BETA = 5.03
I0 = 1.4353734493231072e-4
TC = 0.16857
T0 = 0.16857


def tes_node_indices(mesh_elements: Path, tes_body_id: int = 101) -> list[int]:
    nodes: set[int] = set()
    with mesh_elements.open() as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 7 and int(parts[1]) == tes_body_id:
                nodes.update(int(node_id) - 1 for node_id in parts[3:])
    return sorted(nodes)


def load_nodes(mesh_nodes: Path) -> dict[int, np.ndarray]:
    nodes: dict[int, np.ndarray] = {}
    with mesh_nodes.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            nodes[int(parts[0])] = np.array(
                [float(parts[2]), float(parts[3]), float(parts[4])],
                dtype=float,
            )
    return nodes


def tes_volume_average(mesh_dir: Path, temperature: np.ndarray, tes_body_id: int = 101) -> float:
    nodes = load_nodes(mesh_dir / "mesh.nodes")
    volume_sum = 0.0
    temp_volume_sum = 0.0
    with (mesh_dir / "mesh.elements").open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts or int(parts[1]) != tes_body_id:
                continue
            etype = int(parts[2])
            conn = [int(x) for x in parts[3:]]
            if etype != 504 or len(conn) != 4:
                raise RuntimeError(f"Unexpected TES element type: {etype} ({len(conn)} nodes)")
            p0, p1, p2, p3 = (nodes[nid] for nid in conn)
            volume = abs(np.linalg.det(np.column_stack((p1 - p0, p2 - p0, p3 - p0)))) / 6.0
            tavg = float(np.mean([temperature[nid - 1] for nid in conn]))
            volume_sum += volume
            temp_volume_sum += volume * tavg
    return temp_volume_sum / volume_sum


def circuit_state(temperature: float) -> tuple[float, float, float]:
    a = R0 * (1.0 + ALPHA * (temperature - TC) / TC - BETA)
    b = R0 * BETA / I0
    discriminant = max((R_SH + a) ** 2 + 4.0 * b * I_BIAS * R_SH, 0.0)
    current = (discriminant**0.5 - (R_SH + a)) / (2.0 * b)
    current = min(max(current, 0.0), I_BIAS)
    resistance = a + b * current
    if resistance < R_MIN:
        resistance = R_MIN
        current = I_BIAS * R_SH / (R_SH + resistance)
    power = current * current * resistance
    return resistance, current, power


def summarize(mesh_dir: Path, vtu_path: Path, node_indices: list[int]) -> None:
    mesh = meshio.read(vtu_path)
    temperature = np.asarray(mesh.point_data["temperature"]).reshape(-1)
    tes_temperature = temperature[node_indices]
    nodal_average = float(tes_temperature.mean())
    volume_average = tes_volume_average(mesh_dir, temperature)
    print(vtu_path)
    print(
        "  TES temperature min/max = "
        f"{tes_temperature.min():.12g} "
        f"{tes_temperature.max():.12g} "
    )
    print(f"  TES nodal average = {nodal_average:.12g} K")
    print(f"  TES volume average = {volume_average:.12g} K")
    print(
        "  TES volume average - T0 = "
        f"{volume_average - T0:.12g} K "
        f"({(volume_average - T0) / T0 * 100.0:.12g} %)"
    )
    resistance, current, power = circuit_state(volume_average)
    print(
        "  R_TES / I_TES / P_TES = "
        f"{resistance:.12g} ohm "
        f"{current:.12g} A "
        f"{power:.12g} W"
    )


if __name__ == "__main__":
    mesh_dir = Path("mesh_shifted_merged")
    indices = tes_node_indices(mesh_dir / "mesh.elements")
    print(f"TES nodes: {len(indices)}")
    for path in sorted(mesh_dir.glob("case_*_t0001.vtu")):
        summarize(mesh_dir, path, indices)
