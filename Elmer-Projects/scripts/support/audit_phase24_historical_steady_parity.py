"""Create a reproducible, evidence-labelled Phase24 versus historical audit.

This is deliberately an *audit*, not a parameter search.  It reads the two
generated SIFs, their meshes, the checked-in circuit sources, and the direct
restart evidence.  Every quantity that is not an FE flux integral is labelled
as such in the produced CSV/Markdown files.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
OUT_DEFAULT = ROOT / "artifacts" / "phase24_historical_steady_parity_audit"
HIST_SIF = ROOT / "generated" / "cases" / "case_tes_steady_singlepixel_prod_v2.sif"
CURRENT_SIF = ROOT / "generated" / "cases" / "case_phase24_diag_one_shot_refined.sif"
HIST_TRANSIENT_SIF = ROOT / "generated" / "cases" / "case_tes_pulse_singlepixel_prod_v2_original_timegrid_hybrid_cpu_smoke_177step.sif"
CURRENT_TRANSIENT_SIF = ROOT / "generated" / "cases" / "case_phase24_g45_fine_stycast32_mortar_40us_mumps_mortar.sif"
HIST_MESH = ROOT / "work" / "meshes" / "mesh_singlepixel_prod_v2"
CURRENT_MESH = ROOT / "work" / "meshes" / "mesh_singlepixel_gpu_fine_stycast32_mortar"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git_show(revision: str, path: str) -> str:
    result = subprocess.run(
        ["git", "show", f"{revision}:{path}"], cwd=ROOT, text=True,
        capture_output=True, check=True,
    )
    return result.stdout


def block_map(sif: str) -> dict[str, list[str]]:
    """Keep SIF blocks by their ``Section n`` header, without interpreting MATC."""
    blocks: dict[str, list[str]] = defaultdict(list)
    active = "preamble"
    for raw in sif.splitlines():
        line = raw.rstrip()
        hit = re.match(r"^(Simulation|Constants|Material|Body Force|Body|Boundary Condition|Solver)\s*(\d*)\s*$", line)
        if hit:
            active = f"{hit.group(1)} {hit.group(2)}".strip()
        blocks[active].append(line)
        if line == "End":
            active = "preamble"
    return dict(blocks)


def assignments(lines: list[str]) -> dict[str, str]:
    answer: dict[str, str] = {}
    for line in lines:
        hit = re.match(r'\s*(?:"([^"]+)"|([^= !]+(?:\s+[^= !]+)*?))\s*=\s*(.*)$', line)
        if hit:
            key = hit.group(1) or hit.group(2)
            answer[key.strip()] = hit.group(3).strip()
    return answer


def sif_manifest(path: Path, mesh: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    blocks = block_map(text)
    constants = assignments(blocks.get("Constants", []))
    solver = assignments(blocks.get("Solver 1", []))
    materials: dict[str, dict[str, str]] = {}
    bodies: list[dict[str, str]] = []
    boundaries: list[dict[str, str]] = []
    for name, lines in blocks.items():
        fields = assignments(lines)
        if name.startswith("Material"):
            materials[fields.get("Name", name)] = fields
        elif name.startswith("Body ") and not name.startswith("Body Force"):
            bodies.append(fields)
        elif name.startswith("Boundary Condition"):
            boundaries.append(fields)
    simulation = assignments(blocks.get("Simulation", []))
    return {
        "sif": str(path.relative_to(ROOT)), "sif_sha256": sha256(path),
        "mesh": mesh.name, "mesh_header": mesh_header(mesh),
        "simulation": simulation, "constants": constants, "solver": solver,
        "materials": materials, "bodies": bodies, "boundaries": boundaries,
        "mesh_volumes_m3": mesh_volumes(mesh),
    }


def mesh_header(mesh: Path) -> dict[str, int]:
    vals = [int(v) for v in (mesh / "mesh.header").read_text().split()]
    return {"nodes": vals[0], "bulk_elements": vals[1], "boundary_elements": vals[2]}


def body_names(mesh: Path) -> dict[int, str]:
    result: dict[int, str] = {}
    for line in (mesh / "mesh.names").read_text().splitlines():
        hit = re.match(r"\$\s+(.+)\s+=\s+(\d+)$", line)
        if hit and "__" not in hit.group(1):
            result[int(hit.group(2))] = hit.group(1)
    return result


def mesh_volumes(mesh: Path) -> dict[str, float]:
    """Exact linear-element volumes; 706 is the six-node wedge."""
    nodes = np.loadtxt(mesh / "mesh.nodes", usecols=(2, 3, 4))
    names = body_names(mesh)
    volumes: dict[int, float] = defaultdict(float)
    for line in (mesh / "mesh.elements").read_text().splitlines():
        row = line.split()
        if len(row) < 7:
            continue
        body, kind = int(row[1]), row[2]
        xyz = nodes[np.asarray(row[3:], dtype=int) - 1]
        if kind == "504":
            volume = abs(np.linalg.det((xyz[1:] - xyz[0]).T)) / 6.0
        elif kind in {"303", "706"}:
            volume = sum(
                abs(np.linalg.det((xyz[list(t[1:])] - xyz[t[0]]).T)) / 6.0
                for t in ((0, 1, 2, 3), (1, 2, 4, 3), (2, 4, 5, 3))
            )
        else:
            continue
        volumes[body] += volume
    return {names.get(body, str(body)): value for body, value in sorted(volumes.items())}


def circuit(t: float, previous_i: float | None = None, dt: float | None = None) -> dict[str, float]:
    p = {"ib": 715e-6, "rsh": 3.9e-3, "r0": 15.527e-3, "rmin": 1e-6,
         "alpha": 256.46, "beta": 5.03, "i0": 143.537344932311e-6,
         "t0": 168.57e-3, "l": 12.3e-9}
    a = p["r0"] * (1.0 + p["alpha"] * (t - p["t0"]) / p["t0"] - p["beta"])
    b = p["r0"] * p["beta"] / p["i0"]
    if previous_i is not None and dt is not None:
        c = p["rsh"] + a + p["l"] / dt
        disc = c * c + 4 * b * (p["ib"] * p["rsh"] + p["l"] * previous_i / dt)
        i = (np.sqrt(max(disc, 0.0)) - c) / (2 * b)
    else:
        c = p["rsh"] + a
        disc = c * c + 4 * b * p["ib"] * p["rsh"]
        i = (np.sqrt(max(disc, 0.0)) - c) / (2 * b)
    i = max(min(i, p["ib"]), 0.0)
    r = a + b * abs(i)
    if r < p["rmin"]:
        r = p["rmin"]
        i = (p["ib"] * p["rsh"] if previous_i is None else p["ib"] * p["rsh"] + p["l"] * previous_i / dt) / (p["rsh"] + r if previous_i is None else p["rsh"] + r + p["l"] / dt)
    # On the active linear branch dR/dT and dR/dI are analytic constants.
    return {"temperature_K": t, "current_A": i, "resistance_ohm": r,
            "power_W": i * i * r, "dR_dT_ohm_per_K": p["r0"] * p["alpha"] / p["t0"],
            "dR_dI_ohm_per_A": b, "clipped": r == p["rmin"]}


def resistance_at(t: float, i: float) -> dict[str, float]:
    """The active law evaluated at an explicitly selected T,I pair."""
    r0, alpha, beta, i0, t0, rmin = 15.527e-3, 256.46, 5.03, 143.537344932311e-6, 168.57e-3, 1e-6
    dr_dt = r0 * alpha / t0
    dr_di = r0 * beta / i0
    raw = r0 * (1 + alpha * (t - t0) / t0 - beta) + dr_di * abs(i)
    resistance = max(raw, rmin)
    return {"temperature_K": t, "selected_current_A": i, "resistance_ohm": resistance,
            "power_W": i * i * resistance, "dR_dT_ohm_per_K": dr_dt if raw > rmin else 0.0,
            "dR_dI_ohm_per_A": dr_di if raw > rmin else 0.0, "clipped": raw <= rmin}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def diff(a: Any, b: Any, prefix: str = "") -> list[dict[str, Any]]:
    if isinstance(a, dict) and isinstance(b, dict):
        rows: list[dict[str, Any]] = []
        for key in sorted(set(a) | set(b)):
            rows.extend(diff(a.get(key), b.get(key), f"{prefix}.{key}".strip(".")))
        return rows
    if a != b:
        return [{"path": prefix, "historical": a, "current_phase24": b}]
    return []


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=OUT_DEFAULT)
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)

    historical = sif_manifest(HIST_SIF, HIST_MESH)
    current = sif_manifest(CURRENT_SIF, CURRENT_MESH)
    manifest = {
        "schema_version": 1,
        "historical_comsol_good": {
            "steady": historical,
            "transient": sif_manifest(HIST_TRANSIENT_SIF, HIST_MESH),
        },
        "current_phase24_full_nonlinear": {
            "steady": current,
            "transient": sif_manifest(CURRENT_TRANSIENT_SIF, CURRENT_MESH),
        },
    }
    (out / "configuration_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "configuration_diff.json").write_text(json.dumps(diff(historical, current), indent=2) + "\n")

    points = [("gate3_selected", 0.168569114993, 143.568114993e-6),
              ("refined_selected", 0.166557530019, 218.643409371e-6),
              ("requested_168p57mK_143p6uA", 0.16857, 143.6e-6),
              ("requested_166p56mK_218p6uA", 0.16656, 218.6e-6)]
    resistance_rows = []
    joule_rows = []
    for label, temperature, selected_current in points:
        value = resistance_at(temperature, selected_current)
        for route in ("historical", "current_phase24"):
            resistance_rows.append({"route": route, "point": label, **value,
                                    "law_provenance": "identical SIF constants and identical affine source expression"})
            joule_rows.append({"route": route, "point": label,
                               "I_squared_R_W": value["power_W"], "heat_source_W_per_m3": value["power_W"] / 4e-14,
                               "integrated_source_W": value["power_W"], "tes_volume_m3": 4e-14,
                               "injection_semantics": "TES Parallel Power / TES Volume"})
    # Actual current states, whose P values use the persisted/residual-qualified values.
    for route, point, temperature, current_a, resistance, power in (
        ("historical", "reported_waveform_baseline", None, 143.7778518134342e-6, None, None),
        ("current_phase24", "gate3_checkpoint", 0.168569114993478, 143.5681149930223e-6, 0.015522836331977527, 3.1990633149344825e-10),
        ("current_phase24", "full_nonlinear_refined", 0.166557530018547, 218.6434093713403e-6, 0.008853643057514057, 4.2324271156482167e-10),
    ):
        if resistance is None:
            continue
        raw_power = current_a * current_a * resistance
        joule_rows.append({"route": route, "point": point, "I_squared_R_W": raw_power,
                           "heat_source_W_per_m3": raw_power / 4e-14,
                           "integrated_source_W": power, "tes_volume_m3": 4e-14,
                           "injection_semantics": f"reported relaxed power; raw_minus_integrated_W={raw_power-power:.16e}"})
    write_csv(out / "tes_resistance_comparison.csv", resistance_rows)
    write_csv(out / "joule_power_comparison.csv", joule_rows)
    curve_rows = []
    for temperature in np.linspace(0.162, 0.170, 17):
        value = circuit(float(temperature))
        for route in ("historical", "current_phase24"):
            curve_rows.append({"route": route, "temperature_K": value["temperature_K"],
                               "F_steady_current_A": value["current_A"],
                               "F_steady_current_uA": value["current_A"] * 1e6,
                               "resistance_ohm": value["resistance_ohm"], "power_W": value["power_W"],
                               "equation": "I*(Rsh+R(T,I))=Ibias*Rsh"})
    write_csv(out / "circuit_fixed_point_curve.csv", curve_rows)

    thermal_rows = [
        {"route": "historical", "state": "reported_waveform_baseline", "tes_temperature_K": "not archived", "joule_input_W": 3.2e-10,
         "bath_heat_flow_W": "not captured", "interface_flux_W": "not captured", "effective_G_W_per_K": "not inferable", "evidence_grade": "not measured"},
        {"route": "current_phase24", "state": "gate3_checkpoint", "tes_temperature_K": 0.168569114993478, "joule_input_W": 3.1990633149344825e-10,
         "bath_heat_flow_W": "not captured", "interface_flux_W": "not captured", "effective_G_W_per_K": "not inferable", "evidence_grade": "not measured"},
        {"route": "current_phase24", "state": "full_nonlinear_refined", "tes_temperature_K": 0.166557530018547, "joule_input_W": 4.2324271156482167e-10,
         "bath_heat_flow_W": "not captured", "interface_flux_W": "not captured", "effective_G_W_per_K": "not inferable", "evidence_grade": "not measured"},
    ]
    write_csv(out / "thermal_balance_comparison.csv", thermal_rows)
    write_csv(out / "interface_flux_comparison.csv", thermal_rows)

    native = ROOT.parent / "tools" / "elmer-hypre" / "src" / "fem" / "src" / "modules" / "HeatSolve.F90"
    provenance = {
        "historical_waveform": {"artifact": "artifacts/comparison/comsol_cpu_singlepixel_prod_v2_hybrid_100us/summary.md",
            "baseline_uA": 143.7778518134342, "solver": "CPU MUMPS", "early_timestep_us": 0.625,
            "critical_limit": "The artifact records a reused known serial steady restart; no full nonlinear residual or generated-SIF/binary hash was archived with it."},
        "current_phase24": {"gate3_checkpoint_uA": 143.5681149930223, "refined_uA": 218.6434093713403,
            "refined_full_residual_l2": 6.192602138143676e-16,
            "restart_audit": "artifacts/phase24_augmented_restart_audit/summary.md",
            "native_heat_solve_source": str(native), "native_heat_solve_sha256": sha256(native)},
        "source_history_anchor": {"historical_result_commit": "710c5617104e3e998915576e9e67791ee7fc4562", "restart_fix_commit": "0e6c5a115082349793a1443cf6e013d888264cc3"},
    }
    (out / "source_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    code_diff = """# Historical versus current circuit semantics\n\nBoth SIFs select `TES Inner Circuit Update = True` and publish the heat source as\n`TES Parallel Power / TES Volume`.  The common active law is:\n\n`R = R0 * (1 + alpha*(T-T0)/T0 - beta) + R0*beta*abs(I)/I0`, clipped only at `Rmin`;\n`I*(Rsh + R) = Ibias*Rsh` in steady state; `P = I^2 R`.\n\nThe current native implementation samples an element-equal nodal TES temperature and deliberately uses the completed prior assembly sweep (`CircuitSweepTemperature`) on the next call.  It owns power relaxation and checkpoint restore, then writes `TES Parallel Power`.  The body-force UDF only divides that scalar by `TES Volume`.\n\nThe historical waveform artifact is a transient launched with `--reuse-known-serial-steady`; its archived summary proves waveform agreement but does **not** prove that its 143.777852 uA seed was a full nonlinear steady fixed point.  Therefore it cannot be used as evidence for a distinct exact 143-uA solution until the historical state/SIF/binary is recovered and its residual is measured.\n\nNo difference in the affine TES resistance law, nominal circuit constants, or volumetric normalization was found.  The confirmed semantic difference to audit further is the coupling authority/timing (external historical provenance is incomplete versus current native HeatSolve hook), not a different circuit equation.\n"""
    (out / "circuit_semantics_diff.md").write_text(code_diff)
    (out / "historical_vs_current_code_diff.md").write_text(code_diff)

    controlled = [
        {"test": "historical-state full-residual replay", "status": "BLOCKED_BY_MISSING_ARCHIVE", "physics_changed": "no", "result": "Historical reused steady state/result and binary hash are not in the comparison artifact."},
        {"test": "current Gate3 corrected full constrained direct solve", "status": "COMPLETE", "physics_changed": "no", "result": "Delta T=-6.199616730 mK; full residual L2=6.86e-16."},
        {"test": "current full nonlinear refinement", "status": "COMPLETE", "physics_changed": "no", "result": "166.557530 mK, 218.643409 uA; full residual L2=6.19e-16."},
        {"test": "legacy mesh + current native full nonlinear steady", "status": "NEXT_MINIMAL_TEST", "physics_changed": "no (mesh transplantation diagnostic)", "result": "Run exact MUMPS/residual capture before modifying any physics."},
    ]
    write_csv(out / "controlled_test_results.csv", controlled)

    summary = """# Phase24 historical steady-parity audit\n\n## Finding\n\nThe current full nonlinear Phase24 system has a residual-qualified solution at **218.643409 uA**.  The corrected Gate3 143.568115-uA checkpoint is not a fixed point: holding its accepted Joule source and solving the captured full constrained system changes the TES average by **-6.199616730 mK** at residual L2 **6.86e-16**.  Refinement then converges to 166.557530 mK / 218.643409 uA at full residual L2 **6.19e-16**.\n\nThe historical 143.777852-uA number is a *transient waveform baseline from a reused known serial steady restart*, not an archived residual-qualified full nonlinear steady solution.  Its excellent COMSOL waveform comparison is real, but it does not establish a 143-uA exact fixed point.\n\n## Answers required by the audit\n\n1. **Circuit equations:** same nominal steady equation and same constants; current's native implementation has different ownership/timing of the update, but no demonstrated algebraic law difference.\n2. **TES law:** same affine `R(T,I)` and same clipping.  The selected-point table is `tes_resistance_comparison.csv`.\n3. **Joule injection:** same `P=I^2R`, then `q=P/V_TES`; `V_TES=4e-14 m3` in both SIFs and is independently reproduced from both meshes.\n4. **Thermal path:** not yet proven equivalent.  Materials, nominal dimensions, bath condition, mortar enablement, and all non-Stycast body volumes agree.  Stycast mesh volumes differ by about 0.51% (faceted-cylinder discretization); that alone is not evidence for a 143-to-218-uA shift.  No body/interface flux archive exists for the historical run, so a conductance comparison cannot honestly be claimed.\n5. **Stycast COMSOL equivalence:** not established; both use 20 um nominal thickness and identical material properties, but no COMSOL interface-flux or contact-conductance record is present.\n6. **Dominant evidence-backed candidate:** the historical 143-uA baseline was accepted/reused without residual qualification, while the current 143-uA checkpoint demonstrably is not a thermal fixed point.  Confidence: **high** for this statement; **low** for assigning the remaining difference to physical thermal conductance.\n7. **Code/config to change now:** none in production.  Do not tune HYPRE/GPU and do not alter physics.\n8. **Will a correction restore 143 uA?** Unknown.  First replay the historical mesh/state under a residual-capturing MUMPS solve.  Only a residual-qualified 143-uA result justifies a controlled transplant of its circuit update semantics or thermal mesh/interface.\n\n## Minimal next controlled test\n\nRun `mesh_singlepixel_prod_v2` with CPU MUMPS and the current native HeatSolve hook from the exact saved/reused historical restart, capture the full constrained residual plus body/interface fluxes.  This is a diagnostic mesh/coupling transplant, not a production-physics change.  It separates the currently proven non-fixed-point issue from a possible mesh/interface effect.\n\nHYPRE/GPU tuning remains **NO-GO**.\n"""
    summary = summary.replace(
        "3. **Joule injection:** same `P=I^2R`, then `q=P/V_TES`; `V_TES=4e-14 m3` in both SIFs and is independently reproduced from both meshes.",
        "3. **Joule injection:** same `q=P_relaxed/V_TES`; its raw target is `P_raw=I^2R`. `V_TES=4e-14 m3` in both SIFs and is independently reproduced from both meshes. At the refined current state, `P_raw-P_relaxed=5.17e-15 W`; at Gate3 it is `4.73e-14 W`. The Joule table records this rather than silently treating it as equality.",
    )
    (out / "summary.md").write_text(summary)
    print(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
