"""Extract absorber volume-average temperatures from Phase30 VTU output."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import meshio
import numpy as np

from scripts.support.mesh_names import parse_mesh_names

ROOT = Path(__file__).resolve().parents[2]
MESH = ROOT / "work/meshes/mesh_hybrid_abs_tet_layers_prism_stack17_abs35r50_noextend"
VTU_DIR = ROOT / "work/meshes"
CASE = "case_p19_pulse_phase30_abs_temperature"
PULSE_S = 20.020e-3
OUT = ROOT / "artifacts/hybrid_prism_diagnostics/phase30_abs_temperature"


def absorber_tets() -> tuple[np.ndarray, np.ndarray]:
    body = parse_mesh_names(MESH / "mesh.names").bodies["abs"]
    nodes = np.loadtxt(MESH / "mesh.nodes", usecols=(2, 3, 4))
    tets = []
    with (MESH / "mesh.elements").open() as f:
        for line in f:
            p = line.split()
            if len(p) >= 7 and p[1] == str(body) and p[2] == "504":
                tets.append([int(x) - 1 for x in p[3:7]])
    tet = np.asarray(tets, dtype=int)
    a, b, c, d = (nodes[tet[:, i]] for i in range(4))
    vol = np.abs(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a))) / 6.0
    return tet, vol


def output_times() -> np.ndarray:
    stages = [(18e-6, 1, 1), (1999.9995e-9, 1, 1), (0.1e-12, 10, 1), (999e-12, 1, 1),
              (0.1e-12, 10, 1), (10e-9, 10, 1), (100e-9, 9, 1), (1e-6, 9, 1),
              (10e-6, 9, 1), (5e-6, 80, 10)]
    times, t = [], PULSE_S - 20e-6
    for dt, count, interval in stages:
        for j in range(1, count + 1):
            t += dt
            if j % interval == 0:
                times.append(t)
    return np.asarray(times)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    tet, vol = absorber_tets()
    files = sorted(VTU_DIR.glob(f"{CASE}_t*.vtu"))
    times = output_times()
    if len(files) != len(times):
        raise RuntimeError(f"VTU/time mismatch: {len(files)} files vs {len(times)} output times")
    values = []
    for path in files:
        mesh = meshio.read(path)
        temp = np.asarray(mesh.point_data["temperature"]).reshape(-1)
        values.append(float(np.sum(vol * temp[tet].mean(axis=1)) / np.sum(vol)))
    elmer_t = (times - PULSE_S) * 1e6
    comsol = np.loadtxt(ROOT / "docs/Single-Pixel.txt", comments="%", usecols=(0, 1), encoding="utf-8")
    comsol_t = (comsol[:, 0] * 1e-3 - PULSE_S) * 1e6
    comsol_T = comsol[:, 1]
    # Compare absolute temperature and baseline-relative rise over the common window.
    base_e = float(np.mean(np.asarray(values)[elmer_t < 0]))
    base_c = float(np.mean(comsol_T[comsol_t < 0]))
    values = np.asarray(values)
    rows = [[t, et, et - base_e, ct, ct - base_c] for t, et, ct in zip(elmer_t, values, np.interp(elmer_t, comsol_t, comsol_T))]
    with (OUT / "absorber_temperature_comparison.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["time_from_pulse_us", "elmer_abs_T_K", "elmer_delta_T_K", "comsol_abs_T_K", "comsol_delta_T_K"]); w.writerows(rows)
    metrics = {"elmer_baseline_K": base_e, "comsol_baseline_K": base_c, "elmer_peak_K": float(values.max()), "comsol_peak_K": float(comsol_T.max()), "elmer_peak_time_us": float(elmer_t[values.argmax()]), "comsol_peak_time_us": float(comsol_t[comsol_T.argmax()])}
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True, constrained_layout=True)
    ax1.plot(comsol_t, comsol_T, label="COMSOL", color="#2166ac", linestyle=":", linewidth=2)
    ax1.plot(elmer_t, values, "o-", label="Elmer Phase30 absorber average", color="#6a3d9a", linewidth=1.8, markersize=3)
    ax2.plot(comsol_t, comsol_T - base_c, label="COMSOL", color="#2166ac", linestyle=":", linewidth=2)
    ax2.plot(elmer_t, values - base_e, "o-", label="Elmer Phase30", color="#6a3d9a", linewidth=1.8, markersize=3)
    for ax in (ax1, ax2):
        ax.set_xscale("symlog", linthresh=0.1, linscale=1); ax.set_xlim(-2, 500); ax.axvline(0, color="#555", linestyle="--", linewidth=.9); ax.grid(True, alpha=.25); ax.legend(frameon=False, loc="best")
    ax1.set_ylabel("Absorber temperature [K]"); ax1.set_title(f"Absorber temperature: COMSOL vs Elmer {CASE}")
    ax2.set_ylabel("Temperature rise [K]"); ax2.set_xlabel("Time from pulse [µs] (symlog)")
    fig.savefig(OUT / "absorber_temperature_comsol_vs_elmer.png", dpi=240); fig.savefig(OUT / "absorber_temperature_comsol_vs_elmer.svg"); plt.close(fig)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
