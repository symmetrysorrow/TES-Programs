"""Run isolated native-MUMPS fixed-power coupling diagnostics.

The production mesh and SIFs are untouched.  Each case uses the same
material/geometry/TES constants as the Phase24 diagnostic, freezes only the
checkpoint TES power, and captures the native mortar saddle system.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts" / "phase24_native_coupling_diagnostics"
MESH = "mesh_phase24_native_coupling_diag"
MESH_ROOT = ROOT / "work" / "meshes" / MESH
TEMPLATE = ROOT / "generated" / "cases" / "case_phase24_diag_one_shot_refined.sif"
SOLVER = ROOT.parent / "tools" / "elmer-hypre" / "install-steady-full-capture" / "bin" / "ElmerSolver.exe"
RUNTIME = SOLVER.parent
TOOLCHAIN = Path(r"C:\msys64\ucrt64\bin")
BASE_POWER = 3.203004762115138e-10
FRACTIONAL_DELTA_POWER = 0.05


def write_state(path: Path, power: float) -> None:
    # Circuit current/resistance are retained from the historical checkpoint;
    # HeatSolve's FreezePower path uses only the checkpoint power as its TES
    # heat source and still reports the algebraic current for diagnostics.
    path.write_text(
        f"  1.6856317580942831E-01  1.4377467748651372E-04 "
        f"1.5494931351952193E-02  {power:.16E}  1.4353734493231094E-04\n",
        encoding="utf-8",
    )


def make_sif(name: str, state: Path, capture: Path, power: float) -> Path:
    text = TEMPLATE.read_text(encoding="utf-8")
    text = text.replace('"mesh_singlepixel_gpu_fine_stycast32_mortar"', f'"{MESH}"')
    text = re.sub(r"^  Restart File = .*\r?\n", "", text, flags=re.MULTILINE)
    text = text.replace(
        "../work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/case_phase24_diag_one_shot_refined.result",
        f"../work/meshes/{MESH}/{name}.result",
    )
    text = text.replace(
        "case_phase24_diag_one_shot_refined_series.csv",
        f"{name}_series.csv",
    ).replace(
        "case_phase24_diag_one_shot_refined_iterations.csv",
        f"{name}_iterations.csv",
    )
    text = text.replace(
        "work/meshes/mesh_singlepixel_gpu_fine_stycast32_mortar/case_phase24_diag_one_shot_refined.state",
        f"work/meshes/{MESH}/{name}.state",
    )
    text = text.replace(
        "  Nonlinear System Max Iterations = 84",
        "  Nonlinear System Max Iterations = 8",
    )
    text = re.sub(r"^  \"Phase24 Full Restriction Capture[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Restart State Audit[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^  \"Phase24 Restart Audit Fail Fast\"[^\n]*\r?\n", "", text, flags=re.MULTILINE)
    needle = '  "TES Inner Circuit Step Commit" = Logical True\n'
    capture_rel = capture.relative_to(ROOT).as_posix()
    text = text.replace(
        needle,
        needle
        + '  "TES Inner Circuit Freeze Power" = Logical True\n'
        + '  "Phase24 Full Restriction Capture" = Logical True\n'
        + f'  "Phase24 Full Restriction Capture Prefix" = String "{capture_rel}"\n'
        + '  "Phase24 Full Restriction Capture Timestep" = Integer 1\n'
        + '  "Phase24 Full Restriction Capture Max Iterations" = Integer 1\n',
        1,
    )
    text = text.replace(
        "case_phase24_diag_one_shot_refined.result",
        f"{name}.result",
    )
    path = OUT / f"{name}.sif"
    path.write_text(text, encoding="utf-8")
    return path


def run_case(name: str, power: float) -> dict[str, object]:
    capture = OUT / "capture" / name
    (capture / "ts0001_nl0001").mkdir(parents=True, exist_ok=True)
    state = MESH_ROOT / f"{name}.state"
    write_state(state, power)
    sif = make_sif(name, state, capture, power)
    log = OUT / f"{name}.solver.log"
    env = os.environ.copy()
    env["ELMER_HOME"] = str(SOLVER.parent.parent)
    env["PATH"] = os.pathsep.join(
        [str(TOOLCHAIN), str(SOLVER.parent), str(RUNTIME / "share" / "elmersolver" / "lib"), env.get("PATH", "")]
    )
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            [str(SOLVER), str(sif.relative_to(ROOT).as_posix())],
            cwd=ROOT,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    if completed.returncode:
        raise RuntimeError(f"{name} exited {completed.returncode}; see {log}")
    return {"name": name, "power_W": power, "sif": str(sif.relative_to(ROOT)), "log": str(log.relative_to(ROOT))}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repeat", action="store_true", help="rerun cases even when capture metadata exists")
    args = parser.parse_args()
    if not (MESH_ROOT / "mesh.header").is_file():
        raise SystemExit(f"diagnostic mesh missing: {MESH_ROOT}")
    OUT.mkdir(parents=True, exist_ok=True)
    cases = [
        ("fixed_power_minus5pct", BASE_POWER * (1.0 - FRACTIONAL_DELTA_POWER)),
        ("fixed_power_center", BASE_POWER),
        ("fixed_power_plus5pct", BASE_POWER * (1.0 + FRACTIONAL_DELTA_POWER)),
    ]
    records = []
    for name, power in cases:
        metadata = OUT / "capture" / name / "ts0001_nl0001" / "metadata.json"
        if metadata.is_file() and not args.repeat:
            records.append({"name": name, "power_W": power, "status": "existing"})
        else:
            records.append(run_case(name, power))
    (OUT / "run_manifest.json").write_text(json.dumps(records, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(records, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
