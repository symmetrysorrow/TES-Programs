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
import shlex
import shutil
import subprocess
import sys
import uuid
from pathlib import Path

from support.scenario_config import CIRCUIT_KEYS, ROOT, PROJECT_DIR, compile_project, load_scenario, quantity, write_project
from support.reconcile_project import reconcile_project

WSL_DISTRO = os.environ.get("TES_WSL_DISTRO", "ElmerGPU")
WSL_ELMER_PREFIX = os.environ.get("TES_ELMER_PREFIX_WSL", "/opt/elmer-gpu/cpu/elmer")
WSL_HYPRE_LIB = os.environ.get("TES_HYPRE_LIB_WSL", "/opt/elmer-gpu/cpu/hypre/lib")
WSL_UDF_DIR = os.environ.get("TES_UDF_DIR_WSL", "/opt/elmer-gpu/udf")


def run_command(argv: list[str], *, cwd: Path = ROOT) -> None:
    subprocess.run(argv, cwd=cwd, check=True)


def _wsl_path(path: Path) -> str:
    resolved = path.resolve()
    if os.name != "nt" or not resolved.drive:
        return resolved.as_posix()
    return f"/mnt/{resolved.drive[0].lower()}{resolved.as_posix()[2:]}"


def _wsl_bash(script: str, *, check: bool = True) -> subprocess.CompletedProcess:
    proc = subprocess.run(
        ["wsl.exe", "-d", WSL_DISTRO, "-u", "root", "--exec", "bash", "-lc", script],
        cwd=ROOT, capture_output=True, text=True,
    )
    if check and proc.returncode:
        detail = (proc.stdout or "")[-3000:] + (proc.stderr or "")[-3000:]
        raise RuntimeError(f"WSL ({WSL_DISTRO}) command failed ({proc.returncode}):\n{detail}")
    return proc


def ensure_wsl_partition(mesh_name: str, mpi_procs: int) -> None:
    if mpi_procs <= 1:
        return
    mesh_dir = ROOT / "workspace/work/meshes" / mesh_name
    if (mesh_dir / f"partitioning.{mpi_procs}").is_dir():
        return
    prefix = shlex.quote(WSL_ELMER_PREFIX)
    hypre = shlex.quote(WSL_HYPRE_LIB)
    mesh_parent = shlex.quote(_wsl_path(mesh_dir.parent))
    mesh = shlex.quote(mesh_name)
    script = (
        "set -euo pipefail; "
        f"export ELMER_HOME={prefix}; "
        f"export PATH={prefix}/bin:$PATH; "
        f"export LD_LIBRARY_PATH={prefix}/lib/elmersolver:{hypre}:$LD_LIBRARY_PATH; "
        f"cd {mesh_parent}; ElmerGrid 2 2 {mesh} -metiskway {mpi_procs}"
    )
    _wsl_bash(script)
    if not (mesh_dir / f"partitioning.{mpi_procs}").is_dir():
        raise RuntimeError(f"WSL ElmerGrid did not create partitioning.{mpi_procs} for {mesh_name}")


def run_case_wsl(project: Path, case: str, mesh_name: str, mpi_procs: int) -> None:
    """Run a generated case with the WSL CPU Release solver on ext4 storage."""
    project = project.resolve()
    run_command([sys.executable, str(ROOT / "src" / "sync_elmer_parameters.py"), str(project)])
    mesh_dir = ROOT / "workspace/work/meshes" / mesh_name
    sif = ROOT / "workspace/generated/cases" / f"{case}.sif"
    if not sif.is_file():
        raise FileNotFoundError(f"generated SIF missing: {sif}")
    if not mesh_dir.is_dir():
        raise FileNotFoundError(f"mesh directory missing: {mesh_dir}")
    if mpi_procs > 1 and not (mesh_dir / f"partitioning.{mpi_procs}").is_dir():
        raise FileNotFoundError(f"MPI partition missing: {mesh_dir / f'partitioning.{mpi_procs}'}")

    runtime = f"/tmp/tes-main-{os.getpid()}-{uuid.uuid4().hex[:10]}"
    runtime_q = shlex.quote(runtime)
    root_wsl = shlex.quote(_wsl_path(ROOT))
    project_rel = project.relative_to(ROOT).as_posix()
    mesh_rel = f"workspace/work/meshes/{mesh_name}"
    project_parent = shlex.quote(f"{runtime}/{Path(project_rel).parent.as_posix()}")
    stage = "\n".join([
        "set -euo pipefail",
        f"RUNTIME={runtime_q}",
        f"mkdir -p \"$RUNTIME/workspace/generated/cases\" \"$RUNTIME/workspace/work/meshes\" {project_parent}",
        f"cp -a {shlex.quote(_wsl_path(ROOT / 'src'))} \"$RUNTIME/\"",
        f"cp -a {shlex.quote(_wsl_path(sif))} \"$RUNTIME/workspace/generated/cases/\"",
        f"cp -a {shlex.quote(_wsl_path(project))} \"$RUNTIME/{shlex.quote(project_rel)}\"",
        f"cp -a {shlex.quote(_wsl_path(mesh_dir))} \"$RUNTIME/workspace/work/meshes/\"",
        f"ln -s {shlex.quote(mesh_rel)} \"$RUNTIME/{shlex.quote(mesh_name)}\"",
    ])
    _wsl_bash(stage)

    prefix = shlex.quote(WSL_ELMER_PREFIX)
    hypre = shlex.quote(WSL_HYPRE_LIB)
    udf = shlex.quote(WSL_UDF_DIR)
    solver = shlex.quote(f"{WSL_ELMER_PREFIX}/bin/ElmerSolver_mpi")
    case_q = shlex.quote(case)
    project_q = shlex.quote(project_rel)
    mesh_q = shlex.quote(mesh_rel)
    execute = "\n".join([
        "set -euo pipefail",
        f"RUNTIME={runtime_q}",
        f"ROOT_WSL={root_wsl}",
        f"export ELMER_HOME={prefix}",
        f"export PATH={prefix}/bin:/usr/bin:/bin",
        f"export LD_LIBRARY_PATH={udf}:{hypre}:{prefix}/lib/elmersolver:/usr/lib/wsl/lib",
        f"export TES_UDF_DIR_WSL={udf}",
        f"export TES_HYPRE_LIB_WSL={hypre}",
        "export OMP_NUM_THREADS=1",
        "export HWLOC_COMPONENTS=-opencl,-rsmi,-levelzero",
        "export OCL_ICD_VENDORS=/opt/elmer-gpu/no-ocl-icd",
        "export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1",
        f"test -f {udf}/tes_parallel_circuit.so",
        f"test -f {udf}/tes_transient_heat_source_t0.so",
        f"test -x {solver}",
        "cd \"$RUNTIME\"",
        "set +e",
        f"python3 src/run.py {case_q} --project {project_q} --skip-sync --mpi-procs {mpi_procs} --elmer-solver {solver}",
        "RC=$?",
        "set -e",
        f"mkdir -p \"$ROOT_WSL/outputs/results/{case_q}\" \"$ROOT_WSL/{mesh_q}\"",
        f"if test -d \"$RUNTIME/outputs/results/{case_q}\"; then cp -a \"$RUNTIME/outputs/results/{case_q}/.\" \"$ROOT_WSL/outputs/results/{case_q}/\"; fi",
        f"cp -a \"$RUNTIME/{mesh_q}/.\" \"$ROOT_WSL/{mesh_q}/\"",
        "if test \"$RC\" -eq 0; then rm -rf -- \"$RUNTIME\"; else echo \"WSL runtime retained at $RUNTIME\"; fi",
        "exit \"$RC\"",
    ])
    proc = _wsl_bash(execute, check=False)
    if proc.stdout:
        print(proc.stdout, end="")
    if proc.returncode:
        detail = (proc.stderr or "")[-3000:]
        raise RuntimeError(f"WSL case {case} failed ({proc.returncode}). {detail}")


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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="project/model TOML workflow")
    parser.add_argument("action", choices=("validate", "summary", "explain", "compile", "mesh", "udf", "run"))
    parser.add_argument("project_file", help="project name or projects/<name>.toml")
    parser.add_argument("--dry-run", action="store_true", help="show run plan without building or running")
    parser.add_argument("--mpi-procs", type=int, default=6)
    parser.add_argument("--backend", choices=("wsl-fast", "windows"), default="wsl-fast")
    args = parser.parse_args(argv)
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
            print(f"OK: {len(scenario['tes'])} TES; mesh {mesh_name}; project {scenario['name']}")
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
            mesh_dir = ROOT / "workspace/work" / "meshes" / mesh_name
            mesh_status = "reuse existing" if (mesh_dir / "mesh.names").is_file() else ("required prebuilt mesh missing" if scenario["source_mesh"] else "generate with Gmsh and convert with ElmerGrid")
            print(f"Project: {scenario['name']} ({len(scenario['tes'])} TES)")
            print(f"Mesh: {mesh_status}; model: {scenario['model_file'].relative_to(ROOT)}")
            print("Solve: steady" + (" -> pulse" if scenario["physics"]["pulse"] else ""))
            print(f"Backend: {args.backend}; MPI ranks: {args.mpi_procs}; threads/rank: 1")
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
        udf_dir = ensure_circuit_udf() if args.backend == "windows" else None
        if args.backend == "wsl-fast":
            ensure_wsl_partition(mesh_name, args.mpi_procs)
        # Keep the steady project stable when only pulse parameters change.
        from src.run import restart_result_is_reusable
        result = ROOT / "workspace/work" / "meshes" / mesh_name / f"{steady}.result"
        if args.mpi_procs > 1:
            partition = result.parent / f"partitioning.{args.mpi_procs}"
            if not partition.is_dir():
                raise ValueError(f"MPI mesh partition missing: {partition}; partition the mesh before --mpi-procs {args.mpi_procs}")
            result = result.with_name(result.name + ".0")
        states = [result.parent / f"{steady}_{i}.state" for i in range(1, len(scenario["tes"]) + 1)]
        reusable = restart_result_is_reusable(result, steady, steady_project) and all(p.is_file() for p in states)
        if not reusable:
            if args.backend == "wsl-fast":
                run_case_wsl(steady_project, steady, mesh_name, args.mpi_procs)
            else:
                run_command([sys.executable, str(ROOT / "src" / "run.py"), steady, "--project", str(steady_project), "--mpi-procs", str(args.mpi_procs), "--runtime-bin", str(udf_dir)])
        else:
            print(f"reusing steady: {steady}")
        if scenario["physics"]["pulse"]:
            pulse_project, _, _, _ = write_project(scenario)
            if args.backend == "wsl-fast":
                run_case_wsl(pulse_project, pulse, mesh_name, args.mpi_procs)
            else:
                run_command([sys.executable, str(ROOT / "src" / "run.py"), pulse, "--project", str(pulse_project), "--mpi-procs", str(args.mpi_procs), "--runtime-bin", str(udf_dir)])
        return 0
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
