"""User-facing commands for the TOML scenario workflow."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

from scripts.support.scenario_config import CIRCUIT_KEYS, ROOT, CASE_DIR, compile_project, load_scenario, quantity, write_project
from scripts.support.reconcile_project import reconcile_project


def run_command(argv: list[str], *, cwd: Path = ROOT) -> None:
    subprocess.run(argv, cwd=cwd, check=True)


def ensure_circuit_udf() -> Path:
    source = ROOT / "tes_parallel_circuit.f90"
    target_dir = ROOT / "work" / "udf"
    target = target_dir / ("tes_parallel_circuit.dll" if os.name == "nt" else "tes_parallel_circuit.so")
    if not target.is_file() or target.stat().st_mtime_ns < source.stat().st_mtime_ns:
        compiler = shutil.which("elmerf90") or shutil.which("elmerf90.bat")
        if not compiler:
            raise RuntimeError("elmerf90 not found; add Elmer bin to PATH to build tes_parallel_circuit")
        target_dir.mkdir(parents=True, exist_ok=True)
        run_command([compiler, str(source), "-o", str(target)], cwd=target_dir)
    if not target.is_file():
        raise RuntimeError(f"elmerf90 did not create {target}")
    # Elmer's Windows loader can resolve a DLL in its working directory
    # before PATH, so keep the root copy in sync with the compiled UDF.
    runtime_copy = ROOT / target.name
    if not runtime_copy.is_file() or hashlib.sha256(runtime_copy.read_bytes()).digest() != hashlib.sha256(target.read_bytes()).digest():
        shutil.copy2(target, runtime_copy)
    print(f"circuit UDF ready: {target.relative_to(ROOT)}")
    return target_dir


def mesh_command(scenario: dict, project: Path, mesh_name: str) -> None:
    target = ROOT / "work" / "meshes" / mesh_name
    if (target / "mesh.names").exists():
        print(f"reusing {target.relative_to(ROOT)}")
        return
    if scenario["source_mesh"]:
        raise ValueError(f"registered mesh is missing: {target}; restore or build {mesh_name} before running")
    geometry, mesh = scenario["geometry"], scenario["mesh"]
    msh = ROOT / "generated" / "scenarios" / "meshes" / f"{mesh_name}.msh"
    msh.parent.mkdir(parents=True, exist_ok=True)
    positions = ",".join(f"{t['x']:.17g}:{t['y']:.17g}" for t in geometry["tes"])
    generator = [sys.executable, str(ROOT / "generate_hybrid_prism_geometry.py"), str(project),
                 "--output", str(msh), "--tes-positions", positions,
                 "--absorber-length", str(geometry["absorber_length"]),
                 "--absorber-width", str(geometry["absorber_width"]),
                 "--stycast-diameter", str(geometry["stycast_diameter"]),
                 "--global-mesh-size", str(mesh["global_size"]),
                 "--stack-local-size", str(mesh["tes_local_size"]),
                 "--absorber-local-size", str(mesh["absorber_local_size"]),
                 "--absorber-local-radius", str(mesh["absorber_local_radius"]),
                 "--stycast-layers", str(mesh["stycast_layers"]),
                 "--si-1-layers", str(mesh["si_1_layers"]),
                 "--sio2-1-layers", str(mesh["sio2_1_layers"]),
                 "--sinx-layers", str(mesh["sinx_layers"]),
                 "--disable-mesh-size-extend-from-boundary", "--conformal-tes-stack", "--conformal-abs"]
    run_command(generator)
    elmergrid = shutil.which("ElmerGrid") or shutil.which("ElmerGrid.exe")
    if not elmergrid:
        prefix = Path(os.environ.get("ELMER_HOME", r"C:\Program Files\Elmer 26.1-Release"))
        candidate = prefix / "bin" / ("ElmerGrid.exe" if os.name == "nt" else "ElmerGrid")
        if candidate.is_file():
            elmergrid = str(candidate)
    if not elmergrid:
        raise RuntimeError("ElmerGrid not found; set ELMER_HOME or add ElmerGrid to PATH")
    target.parent.mkdir(parents=True, exist_ok=True)
    run_command([elmergrid, "14", "2", str(msh), "-merge", "1e-10", "-out", str(target)])
    if not (target / "mesh.names").exists():
        raise RuntimeError("ElmerGrid did not create mesh.names")
    msh_digest = hashlib.sha256()
    with msh.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            msh_digest.update(block)
    provenance = {"mesh_hash": scenario["mesh_hash"], "source": str(scenario["case_file"]), "model": str(scenario["model_file"]),
                  "geometry": geometry, "mesh": mesh, "msh_sha256": msh_digest.hexdigest()}
    (target / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"mesh ready: {target.relative_to(ROOT)}")


def explain(scenario: dict, mesh_name: str, steady: str, pulse: str) -> None:
    print(f"case: {scenario['name']}")
    print(f"layout mesh: {mesh_name} (hash {scenario['mesh_hash']})")
    print(f"steady: {steady}; pulse: {pulse if scenario['physics']['pulse'] else '(none)'}")
    print(f"absorber: {scenario['geometry']['absorber_length']:.9g} m x {scenario['geometry']['absorber_width']:.9g} m")
    project, _, _, _ = compile_project(scenario)
    common = reconcile_project(project)["parameters"]
    circuits = project["cases"][steady]["tes_circuits"]
    for i, (tes, overrides) in enumerate(zip(scenario["tes"], circuits), 1):
        print(f"  TES {i}: {tes['id']} at ({tes['x']:.9g}, {tes['y']:.9g}) m")
        for key in sorted(CIRCUIT_KEYS):
            value = quantity(overrides[key], common, key) if key in overrides else common[key]
            source = tes["sources"].get(key, f"{scenario['model_file'].name} [parameters]")
            if key == "I_0" and key in overrides and key not in tes["circuit"]:
                source = "derived from circuit overrides"
            print(f"    {key} = {value:.9g}  [{source}]")


def summary(scenario: dict) -> None:
    project, _, steady, _ = compile_project(scenario)
    common = reconcile_project(project)["parameters"]
    print(f"Case: {scenario['name']} | Model: {scenario['model_file'].stem} | {len(scenario['tes'])} TES")
    print(f"Absorber: {scenario['geometry']['absorber_length'] * 1e3:.4g} x {scenario['geometry']['absorber_width'] * 1e3:.4g} mm")
    print("TES       x [mm]   y [mm]   I_bias [uA]   R_sh [mohm]   alpha")
    for tes, values in zip(scenario["tes"], project["cases"][steady]["tes_circuits"]):
        bias = quantity(values["I_bias"], common, "I_bias") if "I_bias" in values else common["I_bias"]
        shunt = quantity(values["R_sh"], common, "R_sh") if "R_sh" in values else common["R_sh"]
        alpha = quantity(values["alpha"], common, "alpha") if "alpha" in values else common["alpha"]
        print(f"{tes['id']:<9} {tes['x']*1e3:>6.2f}   {tes['y']*1e3:>6.2f}      {bias*1e6:>7.2f}        {shunt*1e3:>7.3f}      {alpha:>7.2f}")
    print(f"Edit cases/{scenario['name']}.toml for run values; {scenario['model_file'].relative_to(ROOT)} for geometry and mesh.")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="TOML component/layout/mesh/case workflow")
    parser.add_argument("action", choices=("validate", "summary", "explain", "compile", "mesh", "udf", "run"))
    parser.add_argument("case_file", help="case name or cases/<name>.toml")
    parser.add_argument("--dry-run", action="store_true", help="show run plan without building or running")
    parser.add_argument("--mpi-procs", type=int, default=1)
    args = parser.parse_args(argv)
    if args.case_file.lower().endswith(".toml"):
        case_file = Path(args.case_file)
        if not case_file.is_absolute():
            case_file = ROOT / case_file
    else:
        case_file = CASE_DIR / f"{args.case_file}.toml"
    try:
        scenario = load_scenario(case_file)
        _, mesh_name, steady, pulse = compile_project(scenario)
        if args.action == "validate":
            print(f"OK: {len(scenario['tes'])} TES; mesh {mesh_name}; case {scenario['name']}")
            return 0
        if args.action == "explain":
            explain(scenario, mesh_name, steady, pulse)
            return 0
        if args.action == "summary":
            summary(scenario)
            return 0
        if args.action == "udf":
            ensure_circuit_udf()
            return 0
        if args.action == "run" and args.dry_run:
            mesh_dir = ROOT / "work" / "meshes" / mesh_name
            mesh_status = "reuse existing" if (mesh_dir / "mesh.names").is_file() else ("required prebuilt mesh missing" if scenario["source_mesh"] else "generate with Gmsh and convert with ElmerGrid")
            print(f"Case: {scenario['name']} ({len(scenario['tes'])} TES)")
            print(f"Mesh: {mesh_status}; model: {scenario['model_file'].relative_to(ROOT)}")
            print("Solve: steady" + (" -> pulse" if scenario["physics"]["pulse"] else ""))
            print("Dry run: no generation or calculation. Remove --dry-run to execute.")
            return 0
        steady_project, _, _, _ = write_project(scenario, steady_only=True)
        if args.action == "compile":
            if scenario["physics"]["pulse"]:
                write_project(scenario)
            print(f"compiled: {steady_project.relative_to(ROOT)}")
            return 0
        mesh_command(scenario, steady_project, mesh_name)
        if args.action == "mesh":
            return 0
        udf_dir = ensure_circuit_udf()
        # Keep the steady project stable when only pulse parameters change.
        from run import restart_result_is_reusable
        result = ROOT / "work" / "meshes" / mesh_name / f"{steady}.result"
        if args.mpi_procs > 1:
            partition = result.parent / f"partitioning.{args.mpi_procs}"
            if not partition.is_dir():
                raise ValueError(f"MPI mesh partition missing: {partition}; partition the mesh before --mpi-procs {args.mpi_procs}")
            result = result.with_name(result.name + ".0")
        states = [result.parent / f"{steady}_{i}.state" for i in range(1, len(scenario["tes"]) + 1)]
        reusable = restart_result_is_reusable(result, steady, steady_project) and all(p.is_file() for p in states)
        if not reusable:
            run_command([sys.executable, str(ROOT / "run.py"), steady, "--project", str(steady_project), "--mpi-procs", str(args.mpi_procs), "--runtime-bin", str(udf_dir)])
        else:
            print(f"reusing steady: {steady}")
        if scenario["physics"]["pulse"]:
            pulse_project, _, _, _ = write_project(scenario)
            run_command([sys.executable, str(ROOT / "run.py"), pulse, "--project", str(pulse_project), "--mpi-procs", str(args.mpi_procs), "--runtime-bin", str(udf_dir)])
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
