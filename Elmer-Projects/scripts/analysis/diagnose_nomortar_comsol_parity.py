"""Diagnose the no-mortar single-pixel model against the mortar reference.

This is an evidence collector, not a solver.  It deliberately separates:

* mesh and interface topology,
* mesh-integrated body volumes and pulse normalization,
* SIF material/circuit/pulse constants,
* current-waveform metrics against COMSOL, and
* optional VTU field jumps at shared interfaces.

The current repository does not retain VTU output for the long no-mortar
run, so the field-continuity section reports ``not_assessed`` unless VTU files
are supplied explicitly.

Usage from the repository root::

    python scripts/analysis/diagnose_nomortar_comsol_parity.py

The defaults use the existing no-mortar fine all-tetra result and the closest
existing mortar all-tetra reference (mesh_refined_3x).  These are not claimed
to be an identical-mesh A/B; the report makes that limitation explicit.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
DEFAULTS = {
    "nomortar_mesh": ROOT / "work/meshes/mesh_singlepixel_conformal_gpu_fine",
    "mortar_mesh": ROOT / "work/meshes/mesh_refined_3x",
    "nomortar_sif": ROOT / "generated/cases/case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step.sif",
    "mortar_sif": ROOT / "generated/cases/case_tes_mpi_comsol_grid_full_uniform_continuous.sif",
    "nomortar_series": ROOT / "results/case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step/case_tes_pulse_singlepixel_conformal_gpu_fine_hybrid_177step_series.csv",
    "mortar_series": ROOT / "artifacts/series/tes_pulse_20ms_3x_series.csv",
    "comsol": ROOT / "docs/Single-Pixel.txt",
    "same_mesh_probe_iterations": ROOT / "artifacts/comparison/same_mesh_mortar_steady_probe/case_tes_steady_singlepixel_conformal_same_mesh_mortar_probe_iterations.csv",
    "coarse_nomortar_iterations": ROOT / "results/case_tes_steady_singlepixel_conformal_gpu/case_tes_steady_singlepixel_conformal_gpu_iterations.csv",
    "outdir": ROOT / "artifacts/comparison/nomortar_parity_diagnostic",
}


BODY_NAMES = {
    100: "abs",
    101: "TES",
    102: "Stycast",
    103: "Membrane_SiNx",
    104: "SiO2_1",
    105: "Si_1",
    106: "SiNx",
    107: "Si_2",
    108: "SiO2_2",
    109: "Membrane_Si1",
}

FACE_DEFINITIONS = {
    504: ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)),
    706: ((0, 1, 2), (3, 4, 5), (0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5)),
}
ELEMENT_NODES = {504: 4, 706: 6}


def _tet_volume(points: np.ndarray) -> float:
    a, b, c, d = points
    return abs(float(np.linalg.det(np.column_stack((b - a, c - a, d - a))))) / 6.0


def _wedge_volume(points: np.ndarray) -> float:
    # Standard 6-node wedge split into three tetrahedra.
    return (
        _tet_volume(points[[0, 1, 2, 3]])
        + _tet_volume(points[[1, 2, 4, 3]])
        + _tet_volume(points[[2, 4, 5, 3]])
    )


def parse_names(path: Path) -> tuple[dict[str, int], dict[int, str]]:
    body_to_id: dict[str, int] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"\s*!?\s*\$\s*([^=]+?)\s*=\s*(-?\d+)", line)
        if match:
            name, value = match.group(1).strip(), int(match.group(2))
            if name in BODY_NAMES.values():
                body_to_id[name] = value
    return body_to_id, {value: key for key, value in body_to_id.items()}


def mesh_report(mesh_dir: Path) -> dict:
    nodes: dict[int, np.ndarray] = {}
    with (mesh_dir / "mesh.nodes").open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) >= 5:
                nodes[int(parts[0])] = np.asarray([float(x) for x in parts[2:5]])

    _, body_names = parse_names(mesh_dir / "mesh.names")
    body_volume: defaultdict[int, float] = defaultdict(float)
    body_element_counts: Counter[tuple[int, int]] = Counter()
    face_owners: defaultdict[tuple[int, ...], list[int]] = defaultdict(list)
    body_nodes: defaultdict[int, set[int]] = defaultdict(set)
    element_count = 0
    unsupported_types: Counter[int] = Counter()

    with (mesh_dir / "mesh.elements").open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 4:
                continue
            body, element_type = int(parts[1]), int(parts[2])
            connectivity = tuple(int(value) for value in parts[3:])
            expected = ELEMENT_NODES.get(element_type)
            if expected is None or len(connectivity) != expected:
                unsupported_types[element_type] += 1
                continue
            element_count += 1
            body_element_counts[(body, element_type)] += 1
            body_nodes[body].update(connectivity)
            points = np.asarray([nodes[node] for node in connectivity])
            if element_type == 504:
                body_volume[body] += _tet_volume(points)
            else:
                body_volume[body] += _wedge_volume(points)
            for face_definition in FACE_DEFINITIONS[element_type]:
                face = tuple(sorted(connectivity[index] for index in face_definition))
                face_owners[face].append(body)

    interfaces: Counter[tuple[int, int]] = Counter()
    interface_nodes: defaultdict[tuple[int, int], set[int]] = defaultdict(set)
    exterior_faces: Counter[int] = Counter()
    for face, owners in face_owners.items():
        unique_owners = sorted(set(owners))
        if len(unique_owners) >= 2:
            for left, right in zip(unique_owners, unique_owners[1:]):
                pair = (left, right)
                interfaces[pair] += 1
                interface_nodes[pair].update(face)
        elif len(unique_owners) == 1:
            exterior_faces[unique_owners[0]] += 1

    body_rows = []
    for body in sorted(set(body_volume) | set(body_names)):
        body_rows.append({
            "body_id": body,
            "body": body_names.get(body, BODY_NAMES.get(body, f"body_{body}")),
            "elements_tet": body_element_counts[(body, 504)],
            "elements_wedge": body_element_counts[(body, 706)],
            "nodes": len(body_nodes[body]),
            "volume_m3": body_volume.get(body, 0.0),
        })

    interface_rows = []
    for pair in sorted(interfaces):
        interface_rows.append({
            "left_body": body_names.get(pair[0], str(pair[0])),
            "right_body": body_names.get(pair[1], str(pair[1])),
            "left_body_id": pair[0],
            "right_body_id": pair[1],
            "faces": interfaces[pair],
            "unique_nodes": len(interface_nodes[pair]),
        })

    return {
        "path": str(mesh_dir),
        "nodes": len(nodes),
        "elements_supported": element_count,
        "element_types": {str(k): int(v) for k, v in body_element_counts.items()},
        "unsupported_element_types": dict(unsupported_types),
        "body_rows": body_rows,
        "interfaces": interface_rows,
        "exterior_faces_by_body": {
            body_names.get(body, str(body)): count for body, count in sorted(exterior_faces.items())
        },
    }


def parse_sif(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    wanted = {
        "Apply Mortar BCs", "Linear System Solver", "Linear System Direct Method",
        "Linear System Iterative Method", "TES Bias Current", "TES Shunt Resistance",
        "TES R0", "TES Rmin", "TES Alpha", "TES Beta", "TES I0", "TES Tc", "TES T0",
        "TES Volume", "Pulse Energy", "Pulse Start Time", "Pulse Duration", "Pulse Sigma",
        "Pulse Radius", "Pulse Center X", "Pulse Center Y", "Pulse Center Z",
        "Pulse Discrete Norm", "Lumped Mass Matrix", "BDF Order",
    }
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "=" not in line or line.lstrip().startswith("!"):
            continue
        key, value = line.split("=", 1)
        key = key.strip().strip('"')
        if key in wanted:
            values[key] = value.strip()
    return values


def _numeric(value: str) -> float | None:
    match = re.search(r"[-+]?\d+(?:\.\d*)?(?:[eEdD][-+]?\d+)?", value)
    if not match:
        return None
    return float(match.group(0).replace("D", "E").replace("d", "e"))


def sif_parity(nomortar: dict[str, str], mortar: dict[str, str]) -> dict:
    rows = []
    for key in sorted(set(nomortar) | set(mortar)):
        left, right = nomortar.get(key), mortar.get(key)
        left_num, right_num = _numeric(left or ""), _numeric(right or "")
        if left_num is not None and right_num is not None:
            delta = left_num - right_num
            equal = bool(np.isclose(left_num, right_num, rtol=1e-12, atol=1e-15))
        else:
            delta = None
            equal = left == right
        rows.append({"key": key, "nomortar": left, "mortar": right, "equal": equal, "delta": delta})
    return {"rows": rows, "all_equal": all(row["equal"] for row in rows)}


def load_series(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.genfromtxt(path, delimiter=",", names=True)
    return np.asarray(data["time_s"], dtype=float) * 1e3, np.asarray(data["tes_current_A"], dtype=float) * 1e6


def load_comsol(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, comments="%", encoding="utf-8")
    # Columns are time, absorber T, Stycast T, TES T, ammeter current,
    # and TES resistance.
    return data[:, 0], data[:, 4]


def crossing_time(time_ms: np.ndarray, response: np.ndarray, level: float) -> float | None:
    indices = np.flatnonzero((response[:-1] < level) & (response[1:] >= level))
    if len(indices) == 0:
        return None
    index = int(indices[0])
    fraction = (level - response[index]) / (response[index + 1] - response[index])
    return float(time_ms[index] + fraction * (time_ms[index + 1] - time_ms[index]))


def waveform_metrics(time_ms: np.ndarray, current_uA: np.ndarray, label: str) -> dict:
    pulse_start = 20.02
    baseline_mask = (time_ms >= 19.5) & (time_ms <= pulse_start)
    baseline = float(np.mean(current_uA[baseline_mask]))
    response = baseline - current_uA
    post = np.flatnonzero(time_ms >= pulse_start)
    peak_index = int(post[np.argmax(response[post])])
    peak = float(response[peak_index])
    t10 = crossing_time(time_ms[:peak_index + 1], response[:peak_index + 1], 0.1 * peak)
    t90 = crossing_time(time_ms[:peak_index + 1], response[:peak_index + 1], 0.9 * peak)
    return {
        "label": label,
        "samples": int(len(time_ms)),
        "baseline_uA": baseline,
        "minimum_uA": float(current_uA[peak_index]),
        "peak_drop_uA": peak,
        "peak_time_ms": float(time_ms[peak_index]),
        "rise_time_10_90_ms": None if t10 is None or t90 is None else t90 - t10,
        "t10_ms": t10,
        "t90_ms": t90,
    }


def compare_waveforms(reference: dict, candidate: dict) -> dict:
    fields = ("baseline_uA", "minimum_uA", "peak_drop_uA", "peak_time_ms", "rise_time_10_90_ms")
    result = {}
    for field in fields:
        left, right = reference.get(field), candidate.get(field)
        if left is None or right is None:
            result[field] = None
        else:
            result[field] = {"candidate_minus_reference": right - left, "relative_percent": 100.0 * (right - left) / left if left else None}
    return result


def final_iteration(path: Path, label: str) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty iteration file: {path}")
    row = rows[-1]
    return {
        "label": label,
        "nonlinear_iter": int(row["nonlinear_iter"]),
        "tes_temperature_K": float(row["tes_temperature_K"]),
        "previous_current_uA": float(row["previous_current_A"]) * 1e6,
        "raw_current_uA": float(row["raw_current_A"]) * 1e6,
        "tes_resistance_ohm": float(row["tes_resistance_ohm"]),
        "raw_power_W": float(row["raw_power_W"]),
        "relaxed_power_W": float(row["relaxed_power_W"]),
    }


def optional_vtu_report(_nomortar: Path | None, _mortar: Path | None) -> dict:
    if _nomortar is None or _mortar is None:
        return {"status": "not_assessed", "reason": "VTU paths were not supplied"}
    return {"status": "not_implemented", "reason": "VTU paths supplied; field extraction requires a mesh/field mapping pass"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("nomortar_mesh", "mortar_mesh", "nomortar_sif", "mortar_sif", "nomortar_series", "mortar_series", "comsol", "outdir"):
        parser.add_argument(f"--{key.replace('_', '-')}", type=Path, default=DEFAULTS[key])
    parser.add_argument("--same-mesh-probe-iterations", type=Path, default=DEFAULTS["same_mesh_probe_iterations"])
    parser.add_argument("--coarse-nomortar-iterations", type=Path, default=DEFAULTS["coarse_nomortar_iterations"])
    parser.add_argument("--nomortar-vtu", type=Path)
    parser.add_argument("--mortar-vtu", type=Path)
    args = parser.parse_args()
    outdir = args.outdir if args.outdir.is_absolute() else ROOT / args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    nomortar_mesh = args.nomortar_mesh if args.nomortar_mesh.is_absolute() else ROOT / args.nomortar_mesh
    mortar_mesh = args.mortar_mesh if args.mortar_mesh.is_absolute() else ROOT / args.mortar_mesh
    nomortar_sif = args.nomortar_sif if args.nomortar_sif.is_absolute() else ROOT / args.nomortar_sif
    mortar_sif = args.mortar_sif if args.mortar_sif.is_absolute() else ROOT / args.mortar_sif
    nomortar_series = args.nomortar_series if args.nomortar_series.is_absolute() else ROOT / args.nomortar_series
    mortar_series = args.mortar_series if args.mortar_series.is_absolute() else ROOT / args.mortar_series
    comsol_path = args.comsol if args.comsol.is_absolute() else ROOT / args.comsol
    probe_iterations = args.same_mesh_probe_iterations if args.same_mesh_probe_iterations.is_absolute() else ROOT / args.same_mesh_probe_iterations
    coarse_iterations = args.coarse_nomortar_iterations if args.coarse_nomortar_iterations.is_absolute() else ROOT / args.coarse_nomortar_iterations

    comsol_t, comsol_i = load_comsol(comsol_path)
    nomortar_t, nomortar_i = load_series(nomortar_series)
    mortar_t, mortar_i = load_series(mortar_series)
    waveforms = {
        "comsol": waveform_metrics(comsol_t, comsol_i, "COMSOL"),
        "nomortar": waveform_metrics(nomortar_t, nomortar_i, "Elmer no-mortar"),
        "mortar": waveform_metrics(mortar_t, mortar_i, "Elmer mortar reference"),
    }

    report = {
        "scope": {
            "same_mesh_ab": False,
            "limitation": "Existing no-mortar fine and mortar refined-3x meshes differ in node count, element type mix, and refinement."
        },
        "mesh": {
            "nomortar": mesh_report(nomortar_mesh),
            "mortar": mesh_report(mortar_mesh),
        },
        "sif": {
            "nomortar": parse_sif(nomortar_sif),
            "mortar": parse_sif(mortar_sif),
            "parity": sif_parity(parse_sif(nomortar_sif), parse_sif(mortar_sif)),
        },
        "waveforms": waveforms,
        "waveform_comparisons": {
            "nomortar_vs_comsol": compare_waveforms(waveforms["comsol"], waveforms["nomortar"]),
            "mortar_vs_comsol": compare_waveforms(waveforms["comsol"], waveforms["mortar"]),
            "nomortar_vs_mortar": compare_waveforms(waveforms["mortar"], waveforms["nomortar"]),
        },
        "same_mesh_mortar_probe": {
            "coarse_no_mortar": final_iteration(coarse_iterations, "coarse no-mortar steady"),
            "no_mortar": final_iteration(
                ROOT / "results/case_tes_steady_singlepixel_conformal_gpu_fine/case_tes_steady_singlepixel_conformal_gpu_fine_iterations.csv",
                "no-mortar fine steady",
            ),
            "mortar_enabled_same_mesh": final_iteration(probe_iterations, "mortar-enabled same fine mesh steady"),
        },
        "field_continuity": optional_vtu_report(args.nomortar_vtu, args.mortar_vtu),
    }
    (outdir / "diagnostic.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    nom = waveforms["nomortar"]
    mor = waveforms["mortar"]
    com = waveforms["comsol"]
    nom_mesh = report["mesh"]["nomortar"]
    mor_mesh = report["mesh"]["mortar"]
    same_nom = report["same_mesh_mortar_probe"]["no_mortar"]
    same_mor = report["same_mesh_mortar_probe"]["mortar_enabled_same_mesh"]
    coarse_nom = report["same_mesh_mortar_probe"]["coarse_no_mortar"]
    lines = [
        "# No-mortar / mortar COMSOL parity diagnostic",
        "",
        "## Scope",
        "",
        "The existing cases are not an identical-mesh A/B. The no-mortar case uses the fine all-tetra conformal mesh; the mortar reference uses the refined-3x all-tetra mesh. This report therefore identifies established differences but does not attribute every difference to mortar alone.",
        "",
        "## Mesh and interface topology",
        "",
        f"| quantity | no-mortar | mortar reference |",
        "|---|---:|---:|",
        f"| nodes | {nom_mesh['nodes']:,} | {mor_mesh['nodes']:,} |",
        f"| supported elements | {nom_mesh['elements_supported']:,} | {mor_mesh['elements_supported']:,} |",
        "| element mix | " + f"tet:{sum(row['elements_tet'] for row in nom_mesh['body_rows']):,}" + " | " + f"tet:{sum(row['elements_tet'] for row in mor_mesh['body_rows']):,}" + " |",
        "",
        "Interface face counts (shared finite-element faces; mortar contact faces are explicit boundary pairs and therefore do not appear as internal faces here):",
        "",
        "| interface | no-mortar faces | mortar reference faces |",
        "|---|---:|---:|",
    ]
    nom_interfaces = {(row["left_body"], row["right_body"]): row["faces"] for row in nom_mesh["interfaces"]}
    mor_interfaces = {(row["left_body"], row["right_body"]): row["faces"] for row in mor_mesh["interfaces"]}
    for key in sorted(set(nom_interfaces) | set(mor_interfaces)):
        lines.append(f"| {key[0]} / {key[1]} | {nom_interfaces.get(key, 0):,} | {mor_interfaces.get(key, 0):,} |")
    lines += [
        "",
        "## SIF parity",
        "",
        "The common circuit/material/pulse constants are compared in `diagnostic.json`. `Apply Mortar BCs` is intentionally different; the report also records the pulse discrete norm, which is mesh-dependent even though the source is normalized to the requested pulse energy.",
        "",
        "## Waveform evidence",
        "",
        "| metric | COMSOL | no-mortar | mortar reference |",
        "|---|---:|---:|---:|",
        f"| baseline current [µA] | {com['baseline_uA']:.6f} | {nom['baseline_uA']:.6f} | {mor['baseline_uA']:.6f} |",
        f"| minimum current [µA] | {com['minimum_uA']:.6f} | {nom['minimum_uA']:.6f} | {mor['minimum_uA']:.6f} |",
        f"| peak current drop [µA] | {com['peak_drop_uA']:.6f} | {nom['peak_drop_uA']:.6f} | {mor['peak_drop_uA']:.6f} |",
        f"| peak delay [ms] | {com['peak_time_ms'] - 20.02:.6f} | {nom['peak_time_ms'] - 20.02:.6f} | {mor['peak_time_ms'] - 20.02:.6f} |",
        f"| 10–90 rise [ms] | {com['rise_time_10_90_ms']:.6f} | {nom['rise_time_10_90_ms']:.6f} | {mor['rise_time_10_90_ms']:.6f} |",
        "",
        "## Same-mesh mortar probe",
        "",
        "| quantity | no-mortar coarse | no-mortar fine | mortar enabled on the same fine mesh | fine→same-mesh mortar |",
        "|---|---:|---:|---:|---:|",
        f"| final TES temperature [K] | {coarse_nom['tes_temperature_K']:.12g} | {same_nom['tes_temperature_K']:.12g} | {same_mor['tes_temperature_K']:.12g} | {same_mor['tes_temperature_K'] - same_nom['tes_temperature_K']:+.3g} |",
        f"| final raw current [µA] | {coarse_nom['raw_current_uA']:.9f} | {same_nom['raw_current_uA']:.9f} | {same_mor['raw_current_uA']:.9f} | {same_mor['raw_current_uA'] - same_nom['raw_current_uA']:+.9f} |",
        f"| final relaxed power [W] | {coarse_nom['relaxed_power_W']:.12g} | {same_nom['relaxed_power_W']:.12g} | {same_mor['relaxed_power_W']:.12g} | {same_mor['relaxed_power_W'] - same_nom['relaxed_power_W']:+.3g} |",
        "",
        "The same-mesh probe converged with MUMPS and changed the raw current by only about -0.0153 µA. The existing coarse→fine no-mortar steady change is much larger (about -36.8 µA), so the dominant observed sensitivity is spatial discretization. The same-mesh mortar flag alone does not reproduce the large improvement seen in the separate mortar reference.",
        "",
        "## Current conclusion",
        "",
        "The no-mortar waveform remains substantially different from COMSOL in the existing data. The same-mesh probe shows that simply adding mortar constraints does not fix it. The strongest current hypothesis is insufficient spatial resolution/topology in the thin TES/Stycast/interface stack: the no-mortar fine mesh has only 1,615 TES tetrahedra and 901 Stycast tetrahedra, while the mortar refined-3x reference has 7,929 and 2,865 respectively. The next test should refine the no-mortar interface stack itself, then compare the converged MUMPS waveform.",
        "",
        "This model solves a thermal PDE plus a lumped TES circuit; it does not solve an electric-potential/current-density PDE. Temperature/heat-flux continuity requires VTU field output and was not assessed here because the stored long runs do not include VTU fields. Electric potential/current-density continuity is not an available field in this model.",
        "",
    ]
    (outdir / "summary.md").write_text("\n".join(lines), encoding="utf-8")
    print(outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
