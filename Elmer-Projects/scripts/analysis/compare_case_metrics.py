from pathlib import Path
import json

import meshio
import numpy as np


ROOT = Path(__file__).resolve().parents[2]


def load_case_summary(json_path: Path) -> dict:
    data = json.loads(json_path.read_text(encoding="utf-8"))
    components = {}
    if "geometry_tree" in data:
        tree = data["geometry_tree"]["children"]
        for group in tree:
            if "children" not in group or not group["children"]:
                continue
            child = group["children"][0]
            components[group["name"]] = child
    elif "geometry" in data:
        tree = data["geometry"]["children"]
        for group in tree:
            if "children" not in group or not group["children"]:
                continue
            child = group["children"][0]
            components[group["name"]] = child
    elif "components" in data:
        for component in data["components"]:
            components[component["name"]] = component["position"]
    else:
        raise KeyError(f"Unsupported JSON structure in {json_path}")
    elmer = data.get("elmer_overrides", {})
    return {
        "path": str(json_path),
        "stycast_shape": elmer.get("stycast_shape", data.get("stycast_shape")),
        "fragment_mortar_interfaces": elmer.get("fragment_mortar_interfaces"),
        "membrane_k_constant": elmer.get("membrane_k_constant"),
        "y_positions": {
            name: components[name]["y"]
            for name in ["SiO2_1", "Si_1", "SiNx", "Si_2", "SiO2_2"]
        },
    }


def bbox_from_nodes(nodes_path: Path) -> dict[str, tuple[float, float]]:
    mins = [float("inf")] * 3
    maxs = [float("-inf")] * 3
    with nodes_path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            xyz = [float(parts[2]), float(parts[3]), float(parts[4])]
            for i, value in enumerate(xyz):
                mins[i] = min(mins[i], value)
                maxs[i] = max(maxs[i], value)
    return {
        "x": (mins[0], maxs[0]),
        "y": (mins[1], maxs[1]),
        "z": (mins[2], maxs[2]),
    }


def tes_node_indices(mesh_elements: Path, tes_body_id: int = 101) -> list[int]:
    nodes: set[int] = set()
    with mesh_elements.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 7 and int(parts[1]) == tes_body_id:
                nodes.update(int(node_id) - 1 for node_id in parts[3:])
    return sorted(nodes)


def load_nodes(nodes_path: Path) -> dict[int, np.ndarray]:
    nodes: dict[int, np.ndarray] = {}
    with nodes_path.open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            nodes[int(parts[0])] = np.array(
                [float(parts[2]), float(parts[3]), float(parts[4])],
                dtype=float,
            )
    return nodes


def tes_temperature_metrics(mesh_dir: Path, vtu_path: Path) -> dict[str, float]:
    mesh = meshio.read(vtu_path)
    temperature = np.asarray(mesh.point_data["temperature"]).reshape(-1)
    tes_idx = tes_node_indices(mesh_dir / "mesh.elements")
    tes = temperature[tes_idx]
    nodes = load_nodes(mesh_dir / "mesh.nodes")
    volume_sum = 0.0
    temp_volume_sum = 0.0
    with (mesh_dir / "mesh.elements").open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts or int(parts[1]) != 101:
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
    return {
        "tes_nodal_average": float(tes.mean()),
        "tes_volume_average": float(temp_volume_sum / volume_sum),
        "tes_min": float(tes.min()),
        "tes_max": float(tes.max()),
        "global_min": float(temperature.min()),
        "global_max": float(temperature.max()),
        "tes_node_count": float(len(tes_idx)),
    }


def element_and_boundary_counts(mesh_dir: Path) -> dict[str, int]:
    counts: dict[str, int] = {
        "tes_elements": 0,
        "bc_1004": 0,
        "bc_1104": 0,
        "bc_1105": 0,
        "bc_1204": 0,
        "bc_1205": 0,
        "bc_1305": 0,
        "bc_1306": 0,
    }
    with (mesh_dir / "mesh.elements").open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if parts and int(parts[1]) == 101:
                counts["tes_elements"] += 1
    with (mesh_dir / "mesh.boundary").open("r", encoding="utf-8") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            bc_id = int(parts[1])
            key = f"bc_{bc_id}"
            if key in counts:
                counts[key] += 1
    return counts


def print_case(title: str, json_path: Path | None, mesh_dir: Path, vtu_path: Path) -> None:
    print(f"[{title}]")
    if json_path is not None:
        summary = load_case_summary(json_path)
        print(f"json={summary['path']}")
        print(f"stycast_shape={summary['stycast_shape']}")
        print(f"fragment_mortar_interfaces={summary['fragment_mortar_interfaces']}")
        print(f"membrane_k_constant={summary['membrane_k_constant']}")
        for name, y in summary["y_positions"].items():
            print(f"{name}.y={y}")
    bbox = bbox_from_nodes(mesh_dir / "mesh.nodes")
    print(f"mesh_dir={mesh_dir}")
    print(f"bbox.x={bbox['x']}")
    print(f"bbox.y={bbox['y']}")
    print(f"bbox.z={bbox['z']}")
    metrics = tes_temperature_metrics(mesh_dir, vtu_path)
    for key, value in metrics.items():
        print(f"{key}={value}")
    counts = element_and_boundary_counts(mesh_dir)
    for key, value in counts.items():
        print(f"{key}={value}")
    print()


def main() -> None:
    print_case(
        "bad_current_project",
        ROOT / "elmer_project.json",
        ROOT / "mesh_shifted_merged",
        ROOT / "mesh_shifted_merged" / "case_constant_power_t0001.vtu",
    )
    print_case(
        "good_fragment_only",
        ROOT / "generated" / "tes_constant_power_compare_fragment.json",
        ROOT / "legacy" / "meshes" / "mesh_fragment_centered_test",
        ROOT / "legacy" / "meshes" / "mesh_fragment_centered_test" / "case_constant_power_fragment_test_t0001.vtu",
    )
    print_case(
        "mesh_shifted",
        None,
        ROOT / "legacy" / "meshes" / "mesh_shifted",
        ROOT / "legacy" / "meshes" / "mesh_shifted" / "case_constant_power_shifted_t0001.vtu",
    )
    print_case(
        "mesh_base",
        None,
        ROOT / "legacy" / "meshes" / "mesh",
        ROOT / "legacy" / "meshes" / "mesh" / "case_constant_power_t0001.vtu",
    )


if __name__ == "__main__":
    main()
