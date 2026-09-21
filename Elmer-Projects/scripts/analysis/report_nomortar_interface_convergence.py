"""Write the evidence report for the no-mortar interface/refinement study.

This report combines the existing parity diagnostic with the newly generated
20/40 um conformal mesh, its MUMPS steady solve, and a short post-pulse probe.
"""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/comparison/nomortar_interface_convergence"
DIAG = ROOT / "scripts/analysis/diagnose_nomortar_comsol_parity.py"


def load_diag_module():
    spec = importlib.util.spec_from_file_location("nomortar_diag", DIAG)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DIAG}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def final_row(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    row = rows[-1]
    return {
        "iterations": int(row["nonlinear_iter"]),
        "temperature_K": float(row["tes_temperature_K"]),
        "current_uA": float(row["raw_current_A"]) * 1e6,
        "resistance_ohm": float(row["tes_resistance_ohm"]),
        "raw_power_W": float(row["raw_power_W"]),
        "relaxed_power_W": float(row["relaxed_power_W"]),
    }


def transient_final_rows(path: Path) -> list[dict]:
    """Keep the converged nonlinear row for each saved timestep."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected: dict[int, dict] = {}
    for row in rows:
        selected[int(row["time_step"])] = row
    return [
        {
            "time_step": int(row["time_step"]),
            "time_s": float(row["time_s"]),
            "temperature_K": float(row["tes_temperature_K"]),
            "current_uA": float(row["raw_current_A"]) * 1e6,
            "resistance_ohm": float(row["tes_resistance_ohm"]),
            "raw_power_W": float(row["raw_power_W"]),
            "nonlinear_iter": int(row["nonlinear_iter"]),
        }
        for row in selected.values()
    ]


def main() -> int:
    diag = load_diag_module()
    meshes = {
        "coarse_no_mortar": ROOT / "work/meshes/mesh_singlepixel_conformal_gpu",
        "fine_no_mortar": ROOT / "work/meshes/mesh_singlepixel_conformal_gpu_fine",
        "refine20_no_mortar": ROOT / "work/meshes/mesh_singlepixel_conformal_gpu_refine20",
        "mortar_reference": ROOT / "work/meshes/mesh_refined_3x",
    }
    iterations = {
        "coarse_no_mortar": ROOT / "results/case_tes_steady_singlepixel_conformal_gpu/case_tes_steady_singlepixel_conformal_gpu_iterations.csv",
        "fine_no_mortar": ROOT / "results/case_tes_steady_singlepixel_conformal_gpu_fine/case_tes_steady_singlepixel_conformal_gpu_fine_iterations.csv",
        "refine20_no_mortar": ROOT / "results/case_tes_steady_singlepixel_conformal_gpu_refine20/case_tes_steady_singlepixel_conformal_gpu_refine20_iterations.csv",
    }
    mesh_report = {name: diag.mesh_report(path) for name, path in meshes.items()}
    steady = {name: final_row(path) for name, path in iterations.items()}
    transient_path = ROOT / "results/case_tes_pulse_singlepixel_conformal_gpu_refine20_2step/case_tes_pulse_singlepixel_conformal_gpu_refine20_2step_iterations.csv"
    comsol = diag.waveform_metrics(
        *diag.load_comsol(ROOT / "docs/Single-Pixel.txt"), "COMSOL"
    )
    existing = json.loads(
        (ROOT / "artifacts/comparison/nomortar_parity_diagnostic/diagnostic.json").read_text(
            encoding="utf-8"
        )
    )

    result = {
        "mesh": mesh_report,
        "steady": steady,
        "transient_probe": transient_final_rows(transient_path),
        "comsol_waveform": comsol,
        "existing_waveforms": existing["waveforms"],
        "sif_parity": existing["sif"]["parity"],
        "field_checks": {
            "temperature": {
                "status": "topologically_continuous_on_no_mortar_mesh",
                "basis": "no-mortar interfaces are shared finite-element faces/nodes; one temperature DOF is used at shared nodes",
                "numerical_jump": "not_measured_from_VTU",
            },
            "heat_flux": {
                "status": "not_assessed",
                "reason": "stored production runs have no VTU output; requires two-sided face-flux extraction",
            },
            "electric_potential": {
                "status": "not_applicable",
                "reason": "the model uses a lumped TES circuit and does not solve an electric-potential PDE",
            },
            "current_density": {
                "status": "not_applicable",
                "reason": "the model has no distributed electric-current-density field",
            },
        },
        "probe_status": {
            "refine20_steady": "completed_exit_0",
            "manifest": "results/case_tes_steady_singlepixel_conformal_gpu_refine20/manifest.json",
            "solver_log": "results/case_tes_steady_singlepixel_conformal_gpu_refine20/solver.log",
        },
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    c = comsol["baseline_uA"]
    lines = [
        "# No-mortar interface and COMSOL convergence study",
        "",
        "## Main result",
        "",
        "The same conformal mesh with mortar enabled changed the steady current by only -0.0153 uA.  The large mismatch was therefore not caused by the mortar flag alone.  Refining the no-mortar mesh moved the steady current from 191.673 uA (coarse) to 154.859 uA (fine), and to 144.764 uA at iteration 25 of the 20/40 um probe; COMSOL is 143.055 uA.  Spatial resolution in the TES/Stycast stack is the dominant demonstrated effect.",
        "",
        "## Mesh and steady convergence",
        "",
        "| case | nodes | tetrahedra | TES tets | Stycast tets | TES volume [m3] | final current [uA] | vs COMSOL |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name, label in [
        ("coarse_no_mortar", "coarse no-mortar"),
        ("fine_no_mortar", "fine no-mortar"),
        ("refine20_no_mortar", "20/40 um no-mortar probe"),
        ("mortar_reference", "mortar reference (different mesh)"),
    ]:
        mesh = mesh_report[name]
        tes = next(row for row in mesh["body_rows"] if row["body"] == "TES")
        current = steady.get(name, {}).get("current_uA")
        current_text = "n/a" if current is None else f"{current:.6f} ({100*(current-c)/c:+.3f}%)"
        lines.append(
            f"| {label} | {mesh['nodes']:,} | {mesh['elements_supported']:,} | "
            f"{tes['elements_tet']:,} | "
            f"{next(row for row in mesh['body_rows'] if row['body'] == 'Stycast')['elements_tet']:,} | "
            f"{tes['volume_m3']:.6e} | {current_text} |"
        )
    lines += [
        "",
        "The refine20 steady solve completed with exit 0 at nonlinear iteration 21 under the 1e-7 convergence setting and produced a restart `.result`. The earlier 25-iteration interrupted probe is retained separately as additional convergence evidence.",
        "",
        "## Interface checks",
        "",
        "- No-mortar temperature is topologically continuous because the interface uses shared nodes and therefore a single temperature degree of freedom. A numerical two-sided heat-flux jump was not measured because the stored runs have no VTU field output.",
        "- Electric potential and current-density continuity are not applicable to this Elmer model: electrical behavior is represented by the lumped TES circuit, not a distributed electric PDE.",
        "- The common COMSOL/Elmer circuit, material, pulse, and initial-state constants are already parity-checked in `nomortar_parity_diagnostic/diagnostic.json`; the mesh-dependent pulse discrete norm remains a separate discretization quantity.",
        "",
        "## Transient convergence",
        "",
        "The existing full waveform comparison remains confounded by different meshes: no-mortar fine has baseline 154.851 uA, peak drop 6.266 uA, and delay 0.100 ms; the mortar reference has 148.130 uA, 7.725 uA, and 0.500 ms; COMSOL has 143.055 uA, 7.775 uA, and 0.428 ms.",
        "",
        "The completed refine20 short probe gives 144.711034 uA at 1 ns and 144.725127 uA at 101 ns after the event (about +0.000253 and +0.014346 uA relative to its steady 144.710781 uA). This confirms a stable, very small immediate response, but it does not determine the full 100 us-scale amplitude or time constant. A full transient refine20 waveform remains computationally expensive with MUMPS.",
        "",
        "## Files",
        "",
        "- `report.json` contains machine-readable mesh, steady, parity, and field-check status.",
        "- `results/case_tes_steady_singlepixel_conformal_gpu_refine20/solver.log` and `.../case_tes_steady_singlepixel_conformal_gpu_refine20_iterations.csv` contain the probe evidence.",
    ]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(OUT / "summary.md")
    print(OUT / "report.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
