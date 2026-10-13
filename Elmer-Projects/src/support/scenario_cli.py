"""User-facing commands for the project/model TOML workflow."""

from __future__ import annotations

import sys as _sys
from pathlib import Path as _Path
_src_path = str(_Path(__file__).resolve().parents[1])
if _src_path not in _sys.path:
    _sys.path.insert(0, _src_path)

import argparse
import hashlib
import json
import os
import shutil
import shlex
import subprocess
import sys
from pathlib import Path

from support.scenario_config import CIRCUIT_KEYS, ROOT, PROJECT_DIR, compile_project, load_scenario, quantity, write_project
from support.reconcile_project import reconcile_project
from support.elmer_wsl_backend import ELMER_PREFIX


def run_command(argv: list[str], *, cwd: Path = ROOT) -> None:
    subprocess.run(argv, cwd=cwd, check=True)


def solver_runtime_args(scenario: dict) -> list[str]:
    """Use the repository's HYPRE build for projects requesting HYPRE."""
    if os.name != "nt" or "hypre" not in scenario["project"].get("solver", {}).get("linear_system", "").lower():
        return []
    # This runtime includes the multi-TES HeatSolve module.
    solver = ROOT / "workspace/work/elmer-latest/bin/ElmerSolver.exe"
    if not solver.is_file():
        raise RuntimeError(f"HYPRE-enabled ElmerSolver not found: {solver}")
    candidates = []
    if os.environ.get("ELMER_HOME"):
        candidates.append(Path(os.environ["ELMER_HOME"]).parent / "msys64/ucrt64/bin")
    hypre_dll = shutil.which("libHYPRE.dll")
    if hypre_dll:
        candidates.append(Path(hypre_dll).parent)
    candidates.append(Path("C:/msys64/ucrt64/bin"))
    toolchain = next((p for p in candidates if (p / "libHYPRE.dll").is_file()), None)
    if toolchain is None:
        raise RuntimeError("HYPRE runtime DLLs not found; add MSYS2 ucrt64/bin to PATH")
    return ["--elmer-solver", str(solver), "--toolchain-bin", str(toolchain)]


def ensure_circuit_udf() -> Path:
    source = ROOT / "src" / "fortran" / "tes_parallel_circuit.f90"
    target_dir = ROOT / "workspace/work" / "udf"
    target = target_dir / ("tes_parallel_circuit.dll" if os.name == "nt" else "tes_parallel_circuit.so")
    if not target.is_file() or target.stat().st_mtime_ns < source.stat().st_mtime_ns:
        compiler = shutil.which("elmerf90") or shutil.which("elmerf90.bat")
        if not compiler:
            raise RuntimeError("elmerf90 not found; add Elmer bin to PATH to build tes_parallel_circuit")
        target_dir.mkdir(parents=True, exist_ok=True)
        run_command([compiler, str(source), "-o", str(target)], cwd=target_dir)
    if not target.is_file():
        raise RuntimeError(f"elmerf90 did not create {target}")
    print(f"circuit UDF ready: {target.relative_to(ROOT)}")
    return target_dir


def mesh_command(scenario: dict, project: Path, mesh_name: str) -> None:
    target = ROOT / "workspace/work" / "meshes" / mesh_name
    if (target / "mesh.names").exists():
        if scenario["source_mesh"] is None:
            provenance = target / "PROVENANCE.json"
            if not provenance.is_file():
                raise RuntimeError(
                    f"Generated mesh {target} is missing PROVENANCE.json; "
                    "move the old cache aside and regenerate before reuse"
                )
            try:
                recorded = json.loads(provenance.read_text(encoding="utf-8"))
            except (OSError, ValueError) as exc:
                raise RuntimeError(f"Generated mesh has invalid provenance: {provenance}") from exc
            if recorded.get("mesh_hash") != scenario["mesh_hash"]:
                raise RuntimeError(
                    f"Generated mesh {target} does not match its TOML configuration; "
                    "move this cache aside before rebuilding"
                )
            required = ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary")
            missing = [name for name in required if not (target / name).is_file()]
            if missing:
                raise RuntimeError(f"Incomplete cached generated mesh {target}: {', '.join(missing)}")
        print(f"reusing {target.relative_to(ROOT)}")
        return
    if scenario["source_mesh"]:
        raise ValueError(f"registered mesh is missing: {target}; restore or build {mesh_name} before running")
    geometry, mesh = scenario["geometry"], scenario["mesh"]
    msh = ROOT / "workspace/generated" / "scenarios" / "meshes" / f"{mesh_name}.msh"
    msh.parent.mkdir(parents=True, exist_ok=True)
    positions = ",".join(f"{t['x']:.17g}:{t['y']:.17g}" for t in geometry["tes"])
    generator = [sys.executable, str(ROOT / "src" / "generate_hybrid_prism_geometry.py"), str(project),
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
    target.parent.mkdir(parents=True, exist_ok=True)
    if elmergrid:
        run_command([elmergrid, "14", "2", str(msh), "-merge", "1e-10", "-out", str(target)])
    elif os.name == "nt":
        # A clean Windows checkout has no host ElmerGrid installation.
        # Convert the host-generated Gmsh mesh with our bootstrapped CPU Elmer in WSL.
        from support.elmer_wsl_backend import ELMER_PREFIX, WSL_DISTRO, _wsl_path
        shell = (
            shlex.quote(f"{ELMER_PREFIX}/elmer/bin/ElmerGrid") + " 14 2 "
            + shlex.quote(_wsl_path(msh)) + " -merge 1e-10 -out "
            + shlex.quote(_wsl_path(target))
        )
        run_command(["wsl.exe", "-d", WSL_DISTRO, "-u", "root", "--exec", "bash", "-lc", shell])
    else:
        raise RuntimeError("ElmerGrid not found; install ElmerGrid or add it to PATH")
    if not (target / "mesh.names").exists():
        raise RuntimeError("ElmerGrid did not create mesh.names")
    msh_digest = hashlib.sha256()
    with msh.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            msh_digest.update(block)
    provenance = {"mesh_hash": scenario["mesh_hash"], "source": str(scenario["project_file"]), "model": str(scenario["model_file"]),
                  "geometry": geometry, "mesh": mesh, "msh_sha256": msh_digest.hexdigest()}
    (target / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"mesh ready: {target.relative_to(ROOT)}")


def explain(scenario: dict, mesh_name: str, steady: str, pulse: str) -> None:
    print(f"project: {scenario['name']}")
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
                source = "derived from project circuit values"
            print(f"    {key} = {value:.9g}  [{source}]")


def summary(scenario: dict) -> None:
    project, _, steady, _ = compile_project(scenario)
    common = reconcile_project(project)["parameters"]
    print(f"Project: {scenario['name']} | Model: {scenario['model_file'].relative_to(ROOT)} | {len(scenario['tes'])} TES")
    print(f"Absorber: {scenario['geometry']['absorber_length'] * 1e3:.4g} x {scenario['geometry']['absorber_width'] * 1e3:.4g} mm")
    print("TES       x [mm]   y [mm]   I_bias [uA]   R_sh [mohm]   alpha")
    for tes, values in zip(scenario["tes"], project["cases"][steady]["tes_circuits"]):
        bias = quantity(values["I_bias"], common, "I_bias") if "I_bias" in values else common["I_bias"]
        shunt = quantity(values["R_sh"], common, "R_sh") if "R_sh" in values else common["R_sh"]
        alpha = quantity(values["alpha"], common, "alpha") if "alpha" in values else common["alpha"]
        print(f"{tes['id']:<9} {tes['x']*1e3:>6.2f}   {tes['y']*1e3:>6.2f}      {bias*1e6:>7.2f}        {shunt*1e3:>7.3f}      {alpha:>7.2f}")
    print(f"Edit {scenario['project_file'].relative_to(ROOT)} for run values; {scenario['model_file'].relative_to(ROOT)} for geometry and mesh.")


def check_mesh_files(scenario: dict, mesh_name: str) -> None:
    mesh_dir = ROOT / "workspace" / "work" / "meshes" / mesh_name
    required = ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary", "mesh.names")
    missing = [name for name in required if not (mesh_dir / name).is_file()]
    empty = [name for name in required if (mesh_dir / name).is_file() and (mesh_dir / name).stat().st_size == 0]
    if missing or empty:
        details = []
        if missing:
            details.append(f"missing: {', '.join(missing)}")
        if empty:
            details.append(f"empty: {', '.join(empty)}")
        raise ValueError(f"{mesh_dir}: {'; '.join(details)}")

    header = (mesh_dir / "mesh.header").read_text(encoding="utf-8").splitlines()
    try:
        counts = [int(value) for value in header[0].split()]
    except (IndexError, ValueError) as exc:
        raise ValueError(f"{mesh_dir / 'mesh.header'} has an invalid element-count header") from exc
    if len(counts) != 3 or any(count <= 0 for count in counts):
        raise ValueError(f"{mesh_dir / 'mesh.header'} must declare positive node, element, and boundary counts")

    from support.mesh_names import parse_mesh_names

    names = parse_mesh_names(mesh_dir / "mesh.names")
    tes_bodies = [name for name in names.bodies if name == "TES" or name.startswith("TES_T")]
    if not names.bodies or not names.boundaries:
        raise ValueError(f"{mesh_dir / 'mesh.names'} must define body and boundary names")
    if len(tes_bodies) != len(scenario["tes"]):
        raise ValueError(
            f"{mesh_dir / 'mesh.names'} has {len(tes_bodies)} TES bodies; expected {len(scenario['tes'])}"
        )


def check_report(scenario: dict, mesh_name: str) -> bool:
    for label in (
        "Project config",
        "Model config",
        "Units/dimensions",
        "TES/circuit mapping",
        "Mesh definition",
    ):
        print(f"{label:<21} OK")
    try:
        check_mesh_files(scenario, mesh_name)
    except (OSError, ValueError) as exc:
        print(f"{'Mesh files':<21} ERROR: {exc}")
        return False
    print(f"{'Mesh files':<21} OK")
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="project/model TOML workflow")
    parser.add_argument("action", choices=("validate", "summary", "explain", "compile", "mesh", "udf", "run"))
    parser.add_argument("project_file", help="project name or projects/<name>.toml")
    parser.add_argument("--dry-run", action="store_true", help="print the run plan without building or running")
    parser.add_argument("--mpi-procs", type=int, default=1)
    parser.add_argument(
        "--backend",
        choices=("elmer_wsl", "elmergpu", "native"),
        default="elmer_wsl",
        help="calculation runtime (default: elmer_wsl; native preserves the former local solver path)",
    )
    args = parser.parse_args(argv)
    # Only the CLI boundary accepts the historical backend name.
    from src.run import normalize_backend_name
    args.backend = normalize_backend_name(args.backend)
    if args.project_file.lower().endswith(".toml"):
        project_file = Path(args.project_file)
        if not project_file.is_absolute():
            project_file = ROOT / project_file
    else:
        project_file = PROJECT_DIR / f"{args.project_file}.toml"
    try:
        scenario = load_scenario(project_file)
        _, mesh_name, steady, pulse = compile_project(scenario)
        if args.action == "validate":
            return 0 if check_report(scenario, mesh_name) else 2
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
            runtime_args = solver_runtime_args(scenario) if args.backend == "native" else []
            mesh_dir = ROOT / "workspace/work" / "meshes" / mesh_name
            mesh_status = "reuse existing" if (mesh_dir / "mesh.names").is_file() else ("required prebuilt mesh missing" if scenario["source_mesh"] else "generate with Gmsh and convert with ElmerGrid")
            print(f"Execution plan: {scenario['name']}")
            print(f"Mesh: {mesh_status} ({mesh_name}); model: {scenario['model_file'].relative_to(ROOT)}")
            print("Solve: steady (reuse a compatible result or calculate)" + (" -> pulse" if scenario["physics"]["pulse"] else ""))
            print(f"MPI ranks: {args.mpi_procs}")
            print(f"Backend: {'Elmer WSL (ElmerGPU distro)' if args.backend == 'elmer_wsl' else 'native Elmer'}")
            linear_system = scenario["project"].get("solver", {}).get("linear_system", "")
            if linear_system:
                print(f"Linear solver: {linear_system} (HYPRE GPU: {'enabled' if linear_system.endswith('_gpu') else 'disabled'})")
            if runtime_args:
                print(f"Solver: {runtime_args[1]}")
                print(f"Toolchain DLLs: {runtime_args[3]}")
            print("Dry run: no mesh generation or calculation.")
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
        runtime_args = solver_runtime_args(scenario) if args.backend == "native" else []
        udf_dir = ensure_circuit_udf() if args.backend == "native" else None
        linear_system = scenario["project"].get("solver", {}).get("linear_system", "")
        print(f"Backend: {'Elmer WSL (ElmerGPU distro)' if args.backend == 'elmer_wsl' else 'native Elmer'}")
        print(f"Linear solver: {linear_system or '(solver default)'}; HYPRE GPU: {'enabled' if linear_system.endswith('_gpu') else 'disabled'}")
        if args.backend == "elmer_wsl":
            print(f"Solver runtime: /opt/elmer-gpu/elmer/bin/{'ElmerSolver_mpi' if args.mpi_procs > 1 else 'ElmerSolver'}")
        elif runtime_args:
            print(f"Solver runtime: {runtime_args[1]}")
        if args.backend == "elmer_wsl":
            from support.elmer_wsl_backend import configure_fast_case

            configure_fast_case(steady_project, steady)
        mesh_dir = ROOT / "workspace/work" / "meshes" / mesh_name
        if args.mpi_procs > 1 and args.backend == "elmer_wsl":
            from support.elmer_wsl_backend import ensure_partition

            ensure_partition(mesh_dir, args.mpi_procs)
        # Keep the steady project stable when only pulse parameters change.
        from src.run import restart_result_is_reusable
        result = mesh_dir / f"{steady}.result"
        if args.mpi_procs > 1:
            partition = result.parent / f"partitioning.{args.mpi_procs}"
            if args.backend == "native" and not partition.is_dir():
                raise ValueError(f"MPI mesh partition missing: {partition}; partition the mesh before --mpi-procs {args.mpi_procs}")
        states = [result.parent / f"{steady}_{i}.state" for i in range(1, len(scenario["tes"]) + 1)]
        solver = (
            f"{ELMER_PREFIX}/elmer/bin/{'ElmerSolver_mpi' if args.mpi_procs > 1 else 'ElmerSolver'}"
            if args.backend == "elmer_wsl" else (runtime_args[1] if runtime_args else None)
        )
        if solver is None:
            from src.run import ELMERSOLVER

            solver = ELMERSOLVER
        reusable = restart_result_is_reusable(
            result, steady, steady_project, mpi_procs=args.mpi_procs,
            backend=args.backend, solver=solver, required_paths=states,
        )
        if not reusable:
            if args.backend == "elmer_wsl":
                from support.elmer_wsl_backend import run_case as run_elmer_wsl_case

                run_elmer_wsl_case(steady_project, steady, mesh_name, args.mpi_procs)
            else:
                run_command([sys.executable, str(ROOT / "src" / "run.py"), steady, "--project", str(steady_project), "--mpi-procs", str(args.mpi_procs), "--runtime-bin", str(udf_dir), *runtime_args])
        else:
            print(f"reusing steady: {steady}")
        if scenario["physics"]["pulse"]:
            pulse_project, _, _, _ = write_project(scenario)
            if args.backend == "elmer_wsl":
                from support.elmer_wsl_backend import configure_fast_case, run_case as run_elmer_wsl_case

                configure_fast_case(pulse_project, pulse)
                run_elmer_wsl_case(pulse_project, pulse, mesh_name, args.mpi_procs)
            else:
                run_command([sys.executable, str(ROOT / "src" / "run.py"), pulse, "--project", str(pulse_project), "--mpi-procs", str(args.mpi_procs), "--runtime-bin", str(udf_dir), *runtime_args])
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
