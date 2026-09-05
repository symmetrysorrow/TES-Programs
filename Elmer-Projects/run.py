"""Run an Elmer case defined in elmer_project.json (redesign plan, Phase 2).

    python run.py <case_name> [--dry-run] [--skip-sync] [--force-deps]

- Regenerates generated/cases/ first (unless --skip-sync).
- Resolves the restart chain: if the case has `restart_from` and the
  dependency's `.result` file is missing (or its recorded run failed), the
  dependency is run first (recursively). `--force-deps` reruns dependencies
  even when present.
- Runs ElmerSolver from the repository root, captures the log, and collects
  the case outputs (VTU/EP from the mesh directory, the series CSV from the
  root, the solver log) into results/<case>/. The `.result` file stays in the
  mesh directory because it is the restart interface between cases.
- Writes results/<case>/manifest.json with input hashes and the run outcome.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ELMERSOLVER = r"C:\Program Files\Elmer 26.1-Release\bin\ElmerSolver.exe"


def load_model(project_path: Path) -> dict:
    sys.path.insert(0, str(ROOT))
    from scripts.support.reconcile_project import reconcile_project

    return reconcile_project(json.loads(project_path.read_text(encoding="utf-8")))


def find_case_insensitive(path: Path) -> Path | None:
    """Return an existing path matching *path*, tolerating a different case
    on the final component. Windows/NTFS is case-preserving on create but
    case-insensitive on lookup; the dual-TES UDF's OPEN(..., STATUS='REPLACE')
    reuses whichever case an earlier run created for a series CSV, so the
    file that lands on disk does not always match the case written into the
    SIF (e.g. '..._l_series.csv' on disk vs '..._L_series.csv' in Constants).
    Returns None if no match is found."""
    if path.exists():
        return path
    if not path.parent.exists():
        return None
    target = path.name.lower()
    for candidate in path.parent.iterdir():
        if candidate.name.lower() == target:
            return candidate
    return None


ELECTRICAL_SERIES_FIELDS = (
    "time_s",
    "time_step",
    "nonlinear_iter",
    "tes_temperature_K",
    "tes_current_A",
    "tes_resistance_ohm",
    "tes_power_W",
    "bias_current_A",
    "shunt_resistance_ohm",
)


def _csv_float(row: dict[str, str], *names: str) -> float | None:
    for name in names:
        value = row.get(name)
        if value not in (None, ""):
            try:
                return float(value)
            except ValueError:
                return None
    return None


def _electrical_rows_from_iteration(
    iteration_file: Path, bias_current: float | None, shunt_resistance: float | None
) -> list[dict[str, object]]:
    with iteration_file.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    # One canonical row per timestep: retain the last nonlinear circuit update
    # so the artifact describes the converged electrical state rather than a
    # log-derived hand selection.
    latest: dict[str, dict[str, str]] = {}
    for row in rows:
        step = row.get("time_step", "0")
        latest[step] = row
    result: list[dict[str, object]] = []
    for row in sorted(latest.values(), key=lambda item: int(item.get("time_step", "0"))):
        result.append(
            {
                "time_s": _csv_float(row, "time_s") or 0.0,
                "time_step": int(float(row.get("time_step", "0"))),
                "nonlinear_iter": int(float(row.get("nonlinear_iter", "0"))),
                "tes_temperature_K": _csv_float(row, "tes_temperature_K") or 0.0,
                "tes_current_A": _csv_float(row, "raw_current_A", "tes_current_A") or 0.0,
                "tes_resistance_ohm": _csv_float(row, "tes_resistance_ohm") or 0.0,
                "tes_power_W": _csv_float(row, "raw_power_W", "tes_power_W") or 0.0,
                "bias_current_A": bias_current,
                "shunt_resistance_ohm": shunt_resistance,
            }
        )
    return result


def ensure_electrical_series(
    root: Path,
    out_dir: Path,
    spec: dict,
    model: dict,
    series_name: str,
    iteration_path: Path | None,
) -> str:
    """Create the canonical series artifact when a solver only emits iterations.

    Custom HeatSolve builds currently emit an iteration CSV for the first
    transient step but may not emit the legacy summary CSV until a timestep
    boundary.  This normalization is deterministic and keeps analysis out of
    solver-log scraping.  A steady state falls back to the final iteration row
    or the shared five-value state file.
    """
    destination = out_dir / series_name
    if destination.exists():
        with destination.open(encoding="utf-8", newline="") as handle:
            header = set(next(csv.reader(handle), []))
        if set(ELECTRICAL_SERIES_FIELDS).issubset(header):
            return "solver_series"
        # Older HeatSolve builds emit a five-column legacy series.  Prefer the
        # richer iteration stream below so every collected parity artifact has
        # the same timestep/nonlinear/bias schema.
        if iteration_path and iteration_path.exists():
            destination.unlink()
        else:
            return "legacy_solver_series"
    params = model.get("parameters", {})
    bias = params.get("I_bias")
    shunt = params.get("R_sh")
    rows: list[dict[str, object]] = []
    if iteration_path and iteration_path.exists():
        rows = _electrical_rows_from_iteration(iteration_path, bias, shunt)
    if not rows:
        state_file = spec.get("state_file")
        if state_file:
            state_path = root / state_file
            if state_path.exists():
                values = state_path.read_text(encoding="utf-8").split()
                if len(values) >= 4:
                    rows = [{
                        "time_s": 0.0,
                        "time_step": 0,
                        "nonlinear_iter": 0,
                        "tes_temperature_K": float(values[0]),
                        "tes_current_A": float(values[1]),
                        "tes_resistance_ohm": float(values[2]),
                        "tes_power_W": float(values[3]),
                        "bias_current_A": bias,
                        "shunt_resistance_ohm": shunt,
                    }]
    if not rows:
        return "not_available"
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ELECTRICAL_SERIES_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    return "normalized_from_iteration_or_state"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def runtime_environment(
    elmer_solver: str, runtime_bin: str | None
) -> tuple[dict[str, str], Path, Path]:
    """Return a loader environment pinned to one Elmer installation."""
    solver = Path(elmer_solver).resolve()
    if not solver.is_file():
        raise FileNotFoundError(f"ElmerSolver not found: {solver}")
    prefix = solver.parent.parent
    parts = [str(solver.parent), str(prefix / "share" / "elmersolver" / "lib")]
    if runtime_bin:
        runtime = Path(runtime_bin).resolve()
        if not runtime.is_dir():
            raise FileNotFoundError(f"runtime DLL directory not found: {runtime}")
        parts.append(str(runtime))
    env = os.environ.copy()
    env["ELMER_HOME"] = str(prefix)
    env["PATH"] = os.pathsep.join([*parts, env.get("PATH", "")])
    return env, solver, prefix


def configure_amgx_sif(
    text: str, config: Path, mesh_dir: Path | None = None,
    constraint_mode: str = "default",
    constraint_penalty: float = 1.0e4,
) -> str:
    """Switch every HeatSolve solver block in *text* to AMGX.

    The generated production SIFs intentionally keep their normal direct
    solver setting.  AMGX is injected only into the execution copy so a GPU
    experiment cannot silently change the reproducible source SIF or project
    JSON.  The config path is written in POSIX form because the GPU run is
    launched inside WSL.
    """
    valid_constraint_modes = {
        "default", "no-scaling", "slave", "master", "slave-transpose", "master-transpose",
        "dual-lagrange", "penalty", "schur", "stabilized",
    }
    if constraint_mode not in valid_constraint_modes:
        raise ValueError(f"invalid AMGX constraint mode: {constraint_mode}")
    if constraint_penalty <= 0:
        raise ValueError("AMGX constraint penalty must be positive")
    config_token = config.as_posix()
    lines = text.splitlines()
    out: list[str] = []
    in_solver = False
    amgx_block = False
    changed = False
    inserted_config = False
    iterative_present = False
    preconditioning_present = False
    mortar_present = False
    constraint_elimination_present = False
    scaling_method_present = False
    linear_system_scaling_present = False

    solver_header = re.compile(r"^\s*Solver\s+\d+\s*$", re.IGNORECASE)
    end_line = re.compile(r"^\s*End\s*$", re.IGNORECASE)
    linear_solver = re.compile(r"^(\s*)Linear System Solver\s*=", re.IGNORECASE)
    amgx_config = re.compile(r"^\s*AMGX Config\s*=", re.IGNORECASE)
    iterative_method = re.compile(
        r"^(\s*)Linear System Iterative Method\s*=", re.IGNORECASE
    )
    preconditioning = re.compile(
        r"^(\s*)Linear System Preconditioning\s*=", re.IGNORECASE
    )
    direct_method = re.compile(
        r"^\s*Linear System Direct Method\s*=", re.IGNORECASE
    )
    apply_mortar = re.compile(
        r"^\s*Apply Mortar BCs\s*=\s*(?:Logical\s+)?True\s*$", re.IGNORECASE
    )
    eliminate_constraints = re.compile(
        r"^(\s*)Eliminate Linear Constraints\s*=", re.IGNORECASE
    )
    scaling_method = re.compile(
        r"^(\s*)Linear System Scaling Method\s*=", re.IGNORECASE
    )
    mesh_file = re.compile(r"^(\s*)(Restart File|Output File)\s*=\s*(.+)$", re.IGNORECASE)

    for line in lines:
        if solver_header.match(line):
            in_solver = True
            amgx_block = False
            inserted_config = False
            iterative_present = False
            preconditioning_present = False
            mortar_present = False
            constraint_elimination_present = False
            scaling_method_present = False
            linear_system_scaling_present = False
        if in_solver and end_line.match(line):
            if amgx_block and not iterative_present:
                method = "GMRES" if constraint_mode == "schur" else (
                    "GCR" if constraint_mode == "stabilized" else "FGMRES"
                )
                out.append(f'  Linear System Iterative Method = "{method}"')
            if amgx_block and not preconditioning_present:
                out.append('  Linear System Preconditioning = "AMG"')
            if amgx_block and not inserted_config:
                out.append(f'  AMGX Config = String "{config_token}"')
                inserted_config = True
            # Mortar coupling otherwise appends zero-diagonal Lagrange
            # multiplier rows.  They make pointwise AMG smoothers singular.
            # Elmer's native constraint elimination gives AMGX the reduced
            # thermal system while preserving the same mortar continuity.
            if (amgx_block and mortar_present and not constraint_elimination_present
                    and constraint_mode not in {"penalty", "schur", "stabilized"}):
                out.append("  Eliminate Linear Constraints = Logical True")
            if amgx_block and mortar_present and constraint_mode == "penalty":
                out.append("  Penalty Linear Constraints = Logical True")
                out.append(f"  Linear Constraint Penalty = Real {constraint_penalty:.17g}")
            if amgx_block and mortar_present and constraint_mode == "schur":
                out.append(
                    f"  AMGX Schur Augmentation = Real {constraint_penalty:.17g}"
                )
                out.append("  Linear System Scaling = Logical True")
                out.append("  Linear System Refactorize = Logical False")
                out.append("  AMGX Allow Not Converged = Logical True")
                out.append("  Linear System Max Iterations = Integer 300")
                out.append("  Linear System Min Iterations = Integer 1")
                out.append("  Linear System Convergence Tolerance = Real 1.0e-6")
                out.append("  Linear System GMRES Restart = Integer 100")
                out.append("  Linear System Abort Not Converged = Logical True")
            if amgx_block and mortar_present and constraint_mode == "stabilized":
                out.append(
                    f"  AMGX Constraint Stabilization = Real {1.0 / constraint_penalty:.17g}"
                )
                out.append("  Linear System Scaling = Logical True")
                out.append("  Linear System Refactorize = Logical False")
                out.append("  AMGX Allow Not Converged = Logical True")
                out.append("  Linear System Max Iterations = Integer 300")
                out.append("  Linear System Min Iterations = Integer 1")
                out.append("  Linear System Convergence Tolerance = Real 1.0e-6")
                out.append("  Linear System Abort Not Converged = Logical True")
            if amgx_block and constraint_mode.startswith("slave"):
                out.append("  Eliminate Slave = Logical True")
            if amgx_block and constraint_mode.startswith("master"):
                out.append("  Eliminate From Master = Logical True")
            if amgx_block and constraint_mode.endswith("transpose"):
                out.append("  Use Transpose Values = Logical True")
            if amgx_block and not scaling_method_present and constraint_mode != "no-scaling":
                out.append('  Linear System Scaling Method = "row equilibration"')
            if amgx_block and constraint_mode == "no-scaling" and not linear_system_scaling_present:
                out.append("  Linear System Scaling = Logical False")
            out.append(line)
            in_solver = False
            continue
        if mesh_dir is not None:
            path_match = mesh_file.match(line)
            if path_match:
                indent, keyword, value = path_match.groups()
                filename = Path(value.strip().strip('"')).name
                line = f"{indent}{keyword} = {(mesh_dir / filename).as_posix()}"
        if in_solver:
            if apply_mortar.match(line):
                mortar_present = True
            match = linear_solver.match(line)
            if match:
                indent = match.group(1)
                if constraint_mode == "schur":
                    solver_name = "AMGX Schur"
                elif constraint_mode == "stabilized":
                    solver_name = "AMGX Stabilized"
                else:
                    solver_name = "AMGX"
                out.append(f'{indent}Linear System Solver = "{solver_name}"')
                amgx_block = True
                changed = True
                continue
            if amgx_block and direct_method.match(line):
                changed = True
                continue
            if amgx_block and iterative_method.match(line):
                indent = iterative_method.match(line).group(1)
                method = "GMRES" if constraint_mode == "schur" else (
                    "GCR" if constraint_mode == "stabilized" else "FGMRES"
                )
                out.append(f'{indent}Linear System Iterative Method = "{method}"')
                iterative_present = True
                changed = True
                continue
            if amgx_block and preconditioning.match(line):
                indent = preconditioning.match(line).group(1)
                out.append(f'{indent}Linear System Preconditioning = "AMG"')
                preconditioning_present = True
                changed = True
                continue
            if amgx_block and amgx_config.match(line):
                out.append(f'  AMGX Config = String "{config_token}"')
                inserted_config = True
                changed = True
                continue
            match = eliminate_constraints.match(line)
            if amgx_block and match:
                enabled = "False" if constraint_mode in {"penalty", "schur", "stabilized"} else "True"
                out.append(f"{match.group(1)}Eliminate Linear Constraints = Logical {enabled}")
                constraint_elimination_present = True
                changed = True
                continue
            match = scaling_method.match(line)
            if amgx_block and match:
                # The eliminated TES heat matrix is nonsymmetric.  Row
                # equilibration bounds its 1e7-scale coefficient spread and
                # makes AMGX's L1 smoother effective.
                out.append(
                    f'{match.group(1)}Linear System Scaling Method = "row equilibration"'
                )
                scaling_method_present = True
                changed = True
                continue
            if amgx_block and re.match(r"^(\s*)Linear System Scaling\s*=", line, re.IGNORECASE):
                indent = re.match(r"^(\s*)Linear System Scaling\s*=", line, re.IGNORECASE).group(1)
                value = "False" if constraint_mode == "no-scaling" else line.split("=", 1)[1].strip()
                out.append(f"{indent}Linear System Scaling = Logical {value}" if constraint_mode == "no-scaling" else line)
                linear_system_scaling_present = True
                changed = True
                continue
        out.append(line)

    if not amgx_block or not changed:
        raise ValueError("SIF does not contain a Linear System Solver block to configure for AMGX")
    result = "\n".join(out)
    if constraint_mode == "dual-lagrange":
        result = configure_mortar_dual_lagrange(result)
    return result + ("\n" if text.endswith(("\n", "\r")) else "")


def configure_mortar_dual_lagrange(text: str) -> str:
    """Use a dual Lagrange basis on Mortar BCs without changing field bases."""
    lines = text.splitlines()
    out: list[str] = []
    block: list[str] = []
    in_boundary = False
    boundary_header = re.compile(r"^\s*Boundary\s+Condition\s+\d+\s*$", re.IGNORECASE)
    end_line = re.compile(r"^\s*End\s*$", re.IGNORECASE)
    mortar = re.compile(r"^\s*Mortar\s+BC\s*=", re.IGNORECASE)
    dual_setting = re.compile(
        r"^\s*(?:Use Biorthogonal Basis|Biorthogonal Dual (?:Slave|Master|Lagrange Coefficients))\s*=",
        re.IGNORECASE,
    )

    def emit_boundary(items: list[str]) -> list[str]:
        if not any(mortar.match(item) for item in items):
            return items
        clean = [item for item in items[:-1] if not dual_setting.match(item)]
        clean.extend([
            "  Use Biorthogonal Basis = Logical True",
            "  Biorthogonal Dual Slave = Logical False",
            "  Biorthogonal Dual Master = Logical False",
            "  Biorthogonal Dual Lagrange Coefficients = Logical True",
            items[-1],
        ])
        return clean

    for line in lines:
        if not in_boundary and boundary_header.match(line):
            in_boundary = True
            block = [line]
            continue
        if in_boundary:
            block.append(line)
            if end_line.match(line):
                out.extend(emit_boundary(block))
                block = []
                in_boundary = False
            continue
        out.append(line)
    if block:
        out.extend(block)
    return "\n".join(out)


def write_runtime_sif(
    source: Path, out_dir: Path, udf_dll: str | None,
    amgx_config: str | None = None, mesh_dir: Path | None = None,
    amgx_constraint_mode: str = "default",
    amgx_constraint_penalty: float = 1.0e4,
) -> tuple[Path, bool]:
    """Create an execution-only SIF with optional UDF and AMGX overrides."""
    if udf_dll is None and amgx_config is None:
        return source, False

    text = source.read_text(encoding="utf-8")
    changed = False
    if udf_dll is not None:
        dll = Path(udf_dll).resolve()
        if not dll.is_file():
            raise FileNotFoundError(f"UDF DLL not found: {dll}")
        token = '"tes_transient_heat_source_t0"'
        target_functions = ("TESTransientHeatSource", "AbsorberWindowPulseHeatSource")
        uses_target = any(name in text for name in target_functions)
        if uses_target:
            if token not in text:
                raise ValueError(
                    f"{source}: expected Procedure library token {token} was not found"
                )
            text = text.replace(token, f'"{dll.as_posix()}"')
            changed = True

    config = None
    if amgx_config is not None:
        config = Path(amgx_config).resolve()
        if not config.is_file():
            raise FileNotFoundError(f"AMGX config not found: {config}")
        text = configure_amgx_sif(
            text, config, mesh_dir=mesh_dir,
            constraint_mode=amgx_constraint_mode,
            constraint_penalty=amgx_constraint_penalty,
        )
        changed = True

    if not changed:
        return source, False
    out_dir.mkdir(parents=True, exist_ok=True)
    runtime_sif = out_dir / "runtime.sif"
    runtime_sif.write_text(text, encoding="utf-8", newline="\n")
    return runtime_sif, True


def runtime_artifacts(
    solver: Path, prefix: Path, udf_dll: str | None
) -> dict[str, str]:
    """Best-effort file hashes for loader-provenance in a run manifest."""
    paths = {"solver": solver, "libelmersolver": solver.parent / "libelmersolver.dll"}
    heat = prefix / "share" / "elmersolver" / "lib" / "HeatSolve.dll"
    if heat.exists():
        paths["HeatSolve"] = heat
    if udf_dll:
        paths["udf"] = Path(udf_dll).resolve()
    return {name: sha256(path) for name, path in paths.items() if path.is_file()}


def mpi_launcher(env: dict[str, str]) -> str:
    """Resolve mpiexec against the pinned child environment, not the parent."""
    launcher = shutil.which("mpiexec", path=env["PATH"])
    if launcher is None:
        raise FileNotFoundError("mpiexec was not found in the selected runtime PATH")
    return str(Path(launcher).resolve())


def restart_chain(model: dict, case_name: str) -> list[str]:
    """Dependency-first list of cases to consider for *case_name*."""
    chain: list[str] = []
    seen: set[str] = set()

    def visit(name: str) -> None:
        if name in seen:
            raise ValueError(f"restart_from cycle involving '{name}'")
        seen.add(name)
        spec = model["cases"].get(name)
        if spec is None:
            raise ValueError(f"unknown case '{name}'")
        dep = spec.get("restart_from")
        if dep:
            visit(dep)
        chain.append(name)

    visit(case_name)
    return chain


def mesh_dir_of(model: dict, case_name: str) -> Path:
    spec = model["cases"][case_name]
    return ROOT / "work" / "meshes" / model["meshes"][spec["mesh"]]["dir"]


def result_file_of(model: dict, case_name: str) -> Path | None:
    spec = model["cases"][case_name]
    if not spec.get("output_result"):
        return None
    serial = mesh_dir_of(model, case_name) / f"{case_name}.result"
    # MPI ResultOutput writes one restart file per rank, conventionally using
    # ``.result.0``, ``.result.1``, ... .  Returning rank 0 here makes the
    # dependency resolver treat a completed parallel steady run as a valid
    # restart source while the SIF still references the common base name.
    parallel_rank0 = mesh_dir_of(model, case_name) / f"{case_name}.result.0"
    return serial if serial.exists() or not parallel_rank0.exists() else parallel_rank0


def restart_result_is_reusable(
    result: Path, case_name: str, project_path: Path
) -> bool:
    """Return whether a restart result may safely be reused.

    Elmer can create a ``.result`` before a later nonlinear iteration or a
    direct-solver allocation fails.  A restart also must not be reused after
    its project has been regenerated with different solver settings.  When a
    manifest is available, its status and project hash are authoritative;
    untracked legacy results remain reusable for backwards compatibility.
    """
    if not result.is_file():
        return False
    manifest_path = ROOT / "results" / case_name / "manifest.json"
    if not manifest_path.is_file():
        return True
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    inputs = manifest.get("inputs_sha256", {})
    return (
        manifest.get("exit_code") == 0
        and not manifest.get("errors")
        and inputs.get(project_path.name) == sha256(project_path)
    )


def preexisting_restart_paths(model: dict, case_name: str, mpi_procs: int) -> list[Path]:
    """Resolve the complete restart interface required by an explicit case."""
    spec = model["cases"][case_name]
    if not spec.get("preexisting_restart"):
        return []
    base = spec.get("restart_file_base")
    if not base:
        raise ValueError(f"{case_name}: preexisting_restart requires restart_file_base")
    state = spec.get("state_file")
    if not state:
        raise ValueError(f"{case_name}: preexisting_restart requires state_file")
    mesh = mesh_dir_of(model, case_name)
    results = (
        [mesh / f"{base}.result"]
        if mpi_procs == 1
        else [mesh / f"{base}.result.{rank}" for rank in range(mpi_procs)]
    )
    return [*results, ROOT / state]


def validate_preexisting_restart(
    model: dict, case_name: str, mpi_procs: int
) -> list[Path]:
    needed = preexisting_restart_paths(model, case_name, mpi_procs)
    missing = [str(path) for path in needed if not path.is_file()]
    if missing:
        raise FileNotFoundError("preexisting restart missing: " + ", ".join(missing))
    return needed


def run_case(
    model: dict,
    case_name: str,
    project_path: Path,
    elmer_solver: str,
    mpi_procs: int,
    udf_dll: str | None = None,
    runtime_bin: str | None = None,
    amgx_config: str | None = None,
    amgx_constraint_mode: str = "default",
    amgx_constraint_penalty: float = 1.0e4,
) -> int:
    spec = model["cases"][case_name]
    sif = ROOT / "generated" / "cases" / f"{case_name}.sif"
    mesh_dir = mesh_dir_of(model, case_name)
    out_dir = ROOT / "results" / case_name
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = out_dir / "solver.log"
    runtime_sif, runtime_override_applied = write_runtime_sif(
        sif, out_dir, udf_dll, amgx_config,
        amgx_constraint_mode=amgx_constraint_mode,
        amgx_constraint_penalty=amgx_constraint_penalty,
    )
    udf_applied = bool(
        udf_dll
        and any(
            name in sif.read_text(encoding="utf-8")
            for name in ("TESTransientHeatSource", "AbsorberWindowPulseHeatSource")
        )
    )
    env, solver_path, prefix = runtime_environment(elmer_solver, runtime_bin)
    if amgx_config or not Path(elmer_solver).suffix.lower() == ".exe":
        # Linux Elmer prefixes restart/output names with the mesh basename.
        # The historical ../work/meshes/... paths therefore need that
        # basename to exist in the launch CWD for ``mesh/../work`` to be
        # normalized by the POSIX filesystem.  The real mesh remains under
        # work/meshes; this is only an empty path anchor for WSL runs.
        mesh_anchor = ROOT / model["meshes"][spec["mesh"]]["dir"]
        mesh_anchor.mkdir(parents=True, exist_ok=True)

    inputs = {
        project_path.name: sha256(project_path),
        sif.name: sha256(sif),
        "mesh.header": sha256(mesh_dir / "mesh.header"),
    }
    preexisting_inputs = {
        str(path.relative_to(ROOT)): sha256(path)
        for path in validate_preexisting_restart(model, case_name, mpi_procs)
    }
    root_dlls = ("tes_heat_source_t0.dll",)
    if udf_dll is None:
        root_dlls = ("tes_transient_heat_source_t0.dll", *root_dlls)
    for dll in root_dlls:
        if (ROOT / dll).exists():
            inputs[dll] = sha256(ROOT / dll)
    amgx_path = Path(amgx_config).resolve() if amgx_config else None
    if amgx_path is not None:
        inputs["amgx_config"] = sha256(amgx_path)

    started = datetime.now().isoformat(timespec="seconds")
    display_sif = (
        runtime_sif.relative_to(ROOT)
        if runtime_sif.is_relative_to(ROOT)
        else runtime_sif
    )
    print(f"[{case_name}] ElmerSolver {display_sif} (log: {log_path.relative_to(ROOT)})")
    with log_path.open("w", encoding="utf-8") as log:
        command = [str(solver_path), str(runtime_sif)]
        if mpi_procs > 1:
            command = [mpi_launcher(env), "-n", str(mpi_procs), *command]
        proc = subprocess.run(
            command,
            cwd=ROOT,
            stdout=log,
            stderr=subprocess.STDOUT,
            env=env,
        )
    finished = datetime.now().isoformat(timespec="seconds")

    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    convergence = re.findall(r"ComputeChange:.*", log_text)[-4:]
    fatal = re.findall(r"(?:ERROR|Fatal).*", log_text)[:3]
    amgx_failures = re.findall(
        r".*(?:Caught amgx exception|AMGX did not converge).*",
        log_text,
        flags=re.IGNORECASE,
    )[:3]
    for line in amgx_failures:
        if line not in fatal:
            fatal.append(line)
    load_failures = re.findall(
        r".*(?:cannot open shared object file|failed to load|Unable to load).*",
        log_text,
        flags=re.IGNORECASE,
    )
    for line in load_failures:
        if line not in fatal:
            fatal.append(line)
    reported_done = "MAIN: *** Elmer Solver: ALL DONE ***" in log_text
    if not reported_done:
        fatal.append("solver did not report 'ALL DONE'")
    if proc.returncode != 0:
        fatal.append(f"solver process exited with code {proc.returncode}")
    completed = reported_done and proc.returncode == 0 and not fatal

    collected: list[str] = []
    electrical_series_source = "not_available"
    for pattern in (f"{case_name}_t*.vtu", f"{case_name}*.pvtu", f"{case_name}.ep"):
        for src in mesh_dir.glob(pattern):
            shutil.move(str(src), out_dir / src.name)
            collected.append(src.name)
    series = spec.get("series_file")
    if series:
        found = find_case_insensitive(ROOT / series)
        if found is not None:
            # Single-instance cases (single-pixel): one series CSV, exact name.
            shutil.move(str(found), out_dir / series)
            collected.append(series)
            electrical_series_source = "solver_series"
        else:
            # Dual-TES cases: the builder writes one series CSV per circuit
            # instance with an 'L'/'R' side tag inserted before the
            # '_series.csv' suffix (scripts/support/build_cases.py,
            # _side_series_file) instead of the single base name.
            sys.path.insert(0, str(ROOT))
            from scripts.support.build_cases import _side_series_file

            for side in ("L", "R"):
                side_name = _side_series_file(series, side)
                found_side = find_case_insensitive(ROOT / side_name)
                if found_side is not None:
                    shutil.move(str(found_side), out_dir / side_name)
                    collected.append(side_name)
    iteration_series = spec.get("iteration_series_file")
    iteration_output_path: Path | None = None
    if iteration_series:
        found_iteration = find_case_insensitive(ROOT / iteration_series)
        if found_iteration is not None:
            shutil.move(str(found_iteration), out_dir / iteration_series)
            collected.append(iteration_series)
            iteration_output_path = out_dir / iteration_series
    if series:
        electrical_series_source = ensure_electrical_series(
            ROOT, out_dir, spec, model, series, iteration_output_path
        )
        if (out_dir / series).exists() and series not in collected:
            collected.append(series)
    result_file = result_file_of(model, case_name)

    manifest = {
        "case": case_name,
        "started": started,
        "finished": finished,
        "exit_code": proc.returncode,
        "inputs_sha256": inputs,
        "preexisting_restart_inputs_sha256": preexisting_inputs,
        "runtime_sif": str(runtime_sif.relative_to(ROOT)) if runtime_sif != sif else None,
        "runtime_sif_sha256": sha256(runtime_sif),
        "runtime_override_applied": runtime_override_applied,
        "udf_applied": udf_applied,
        "udf_dll": str(Path(udf_dll).resolve()) if udf_applied else None,
        "udf_sha256": sha256(Path(udf_dll).resolve()) if udf_applied else None,
        "solver": str(solver_path),
        "elmer_prefix": str(prefix),
        "runtime_bin": str(Path(runtime_bin).resolve()) if runtime_bin else None,
        "amgx_config": str(amgx_path) if amgx_path else None,
        "amgx_config_sha256": sha256(amgx_path) if amgx_path else None,
        "runtime_artifacts_sha256": runtime_artifacts(
            solver_path, prefix, udf_dll if udf_applied else None
        ),
        "mesh": spec["mesh"],
        "restart_from": spec.get("restart_from"),
        "collected_outputs": sorted(collected),
        "electrical_series_source": electrical_series_source,
        "result_file": str(result_file.relative_to(ROOT)) if result_file else None,
        "convergence_tail": convergence,
        "solver_completed": completed,
        "errors": fatal,
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    status = "OK" if completed else "FAILED"
    print(f"[{case_name}] {status} (exit {proc.returncode}); outputs -> {out_dir.relative_to(ROOT)}")
    for line in convergence[-2:]:
        print(f"[{case_name}]   {line.strip()}")
    for line in fatal:
        print(f"[{case_name}]   {line.strip()}")
    return 0 if completed else (proc.returncode or 1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case")
    parser.add_argument("--dry-run", action="store_true", help="print the plan without running")
    parser.add_argument("--skip-sync", action="store_true", help="do not regenerate generated/cases first")
    parser.add_argument("--force-deps", action="store_true", help="rerun restart dependencies even if their .result exists")
    parser.add_argument(
        "--project",
        default="elmer_project.json",
        help="path to the project JSON (default: elmer_project.json). Use an "
        "alternate file to try different parameters without touching the "
        "main config; its case must use a distinct name so results/ and "
        "generated/cases/ outputs don't collide with the main project's.",
    )
    parser.add_argument(
        "--udf-dll",
        help="pin tes_transient_heat_source_t0 to this DLL in results/<case>/runtime.sif",
    )
    parser.add_argument(
        "--runtime-bin",
        help="runtime DLL directory appended after the selected Elmer install's bin/module paths",
    )
    parser.add_argument(
        "--amgx-config",
        help=(
            "AMGX JSON config to inject into the target case's execution-only "
            "SIF; restart dependencies retain their configured solver"
        ),
    )
    parser.add_argument(
        "--amgx-constraint-mode",
        choices=[
            "default", "no-scaling", "slave", "master", "slave-transpose", "master-transpose",
            "dual-lagrange", "penalty", "schur", "stabilized",
        ],
        default="default",
        help="Mortar constraint elimination variant used by the AMGX target",
    )
    parser.add_argument(
        "--amgx-constraint-penalty",
        type=float,
        default=1.0e4,
        help="normalized penalty strength for --amgx-constraint-mode penalty",
    )
    parser.add_argument(
        "--elmer-solver",
        default=ELMERSOLVER,
        help="path to ElmerSolver executable; use this to select an alternate build such as the HYPRE/MPI solver",
    )
    parser.add_argument(
        "--mpi-procs",
        type=int,
        default=1,
        help="number of MPI ranks; requires a mesh/partitioning.N directory when N > 1",
    )
    args = parser.parse_args()

    project_path = Path(args.project)
    if not project_path.is_absolute():
        project_path = ROOT / project_path

    if not args.skip_sync:
        subprocess.run(
            [sys.executable, str(ROOT / "sync_elmer_parameters.py"), str(project_path)],
            cwd=ROOT,
            check=True,
        )

    model = load_model(project_path)
    chain = restart_chain(model, args.case)
    validate_preexisting_restart(model, args.case, args.mpi_procs)
    target_preexisting = bool(model["cases"][args.case].get("preexisting_restart"))
    if target_preexisting and args.force_deps:
        raise ValueError("--force-deps cannot be used with preexisting_restart")

    plan: list[str] = []
    for name in chain[:-1]:
        result = result_file_of(model, name)
        if result is None:
            raise ValueError(f"'{name}' is a restart dependency but does not set output_result")
        if not target_preexisting and (
            args.force_deps or not restart_result_is_reusable(result, name, project_path)
        ):
            plan.append(name)
        else:
            print(f"[{name}] successful restart field {result.relative_to(ROOT)} exists - skipping (use --force-deps to rerun)")
    plan.append(args.case)

    print("run plan: " + " -> ".join(plan))
    if args.dry_run:
        return 0

    for name in plan:
        # A restart dependency is an input-generation step, not part of the
        # requested target override.  Keeping its configured direct solver is
        # also important for this TES model: eliminating mortar multipliers is
        # required by AMGX but changes the steady initial field measurably.
        target_amgx_config = args.amgx_config if name == args.case else None
        code = run_case(
            model, name, project_path, args.elmer_solver, args.mpi_procs,
            args.udf_dll, args.runtime_bin, target_amgx_config,
            args.amgx_constraint_mode, args.amgx_constraint_penalty,
        )
        if code != 0:
            print(f"aborting chain: {name} failed")
            return code
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
