"""TES simulation front end: one small case file instead of long command lines.

    python tes_sim.py run       cases/single_pixel_h8.toml   # partition -> steady -> transient -> summary
    python tes_sim.py steady    cases/post2_h8.toml
    python tes_sim.py transient cases/post2_h8.toml          # runs the steady first if needed
    python tes_sim.py show      cases/post2_h8.toml          # resolved settings, no run
    python tes_sim.py summary   cases/post2_h8.toml          # per-circuit table of the last run
    python tes_sim.py mesh-import gmsh/foo.msh --name mesh_foo [--roles roles.toml]

A case file (TOML) names a registered mesh, optionally the role table of its
bodies (which body is a TES circuit, which material, bath, pulse bodies),
per-TES circuit overrides, the pulse and the time window.  Solver settings
come from a preset (the validated fast Phase24 HYPRE route by default) and
can be overridden in [solver].  See docs/extending_tes_models.md.

The heavy lifting (SIF generation, WSL/MPI launch) is done by the validated
scripts/support/run_phase24_gpu_elmergpu.py, called as a subprocess.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
RUNNER = ROOT / "scripts" / "support" / "run_phase24_gpu_elmergpu.py"
MESH_REGISTRY = ROOT / "artifacts" / "phase24_conformal_hybrid" / "project.json"
RUNS = ROOT / "runs"
DISTRO = "ElmerGPU"
ELMER_PREFIX = "/opt/elmer-gpu"
PULSE_START = 0.02002  # s; the transient windows put the pulse here

# ---------------------------------------------------------------------------
# Presets: the validated solver routes.  "fast" is the best 75 ms setting
# (artifacts/phase24_conformal_hybrid/summary.md); linear_tolerance None means
# "auto": 3e-11 for one TES, 1e-11 for several (larger systems need the
# tighter tolerance to stay below the 1e-7 K coupled nonlinear criterion).
PRESETS: dict[str, dict] = {
    "fast": {
        "elmer": "cpu/elmer",
        "linear_system": "iterative_hypre_pcg_boomeramg",
        "linear_tolerance": None,
        "boomer_amg": {"strong_threshold": 0.25, "relax_type": 8, "coarsen_type": 8},
        "nonlinear_tolerance": 1.0e-6,
        "bdf_order": 2,
        "predictor": True,
        "predictor_max_delta": 1.0e-3,
        "window": "75ms_bdf2h15",
        "extra_solver_keywords": [
            "TES Inner Circuit Newton = Logical True",
            "Phase24 BDF1 Step Start Times(2) = Real 0.02002 0.020020001",
            "Phase24 Prism Lumped Preconditioner = Logical True",
            "Phase24 Precond Matrix Reuse = Logical True",
            "Phase24 Coupled Temperature Tolerance = Real 1e-7",
            "TES Inner Circuit Predict From Temperature = Logical True",
            "TES Inner Circuit Newton Slope Reuse = Logical True",
            "BoomerAMG Aggressive Levels = Integer 1",
        ],
        # Steady: the same circuit Newton and TES-temperature exit as the
        # transient (the legacy relaxed circuit update stops on small field
        # increments long before the circuit has converged).
        "steady_extra_solver_keywords": [
            "Phase24 Prism Lumped Preconditioner = Logical True",
            "Phase24 Precond Matrix Reuse = Logical True",
            "BoomerAMG Aggressive Levels = Integer 1",
            "TES Inner Circuit Newton = Logical True",
            "Phase24 Coupled Temperature Tolerance = Real 1e-7",
        ],
    },
}
PRESETS["accurate"] = {**PRESETS["fast"], "linear_tolerance": 1.0e-12, "window": "75ms_bdf2"}

DEFAULT_PULSE = {"energy": "1332[keV]", "start": "20.02[ms]", "duration": "1[ns]", "sigma": "50[um]", "center": "auto"}


# ---------------------------------------------------------------------------
def load_case(path: Path) -> dict:
    cfg = tomllib.loads(path.read_text(encoding="utf-8"))
    if "mesh" not in cfg:
        raise SystemExit(f"{path}: 'mesh' is required")
    cfg.setdefault("name", path.stem)
    cfg.setdefault("mpi", 6)
    cfg.setdefault("preset", "fast")
    if cfg["preset"] not in PRESETS:
        raise SystemExit(f"{path}: unknown preset '{cfg['preset']}' (known: {sorted(PRESETS)})")
    return cfg


def _registry() -> dict:
    return json.loads(MESH_REGISTRY.read_text(encoding="utf-8"))


def mesh_dir(cfg: dict) -> Path:
    meshes = _registry()["meshes"]
    if cfg["mesh"] not in meshes:
        raise SystemExit(f"mesh '{cfg['mesh']}' is not registered (use: tes_sim.py mesh-import)")
    return ROOT / "work" / "meshes" / meshes[cfg["mesh"]]["dir"]


def roles_of(cfg: dict) -> dict | None:
    return cfg.get("roles") or _registry()["meshes"][cfg["mesh"]].get("roles")


def tes_bodies(cfg: dict) -> list[str]:
    """TES circuit bodies (role table or naming convention), circuit order."""
    from generate_project_geometry import reconcile_project
    import scripts.support.build_cases as bc
    from scripts.support.mesh_names import parse_mesh_names

    model = reconcile_project(_registry())
    bc.set_active_roles(roles_of(cfg), model)
    return bc.resolve_tes_body_names(parse_mesh_names(mesh_dir(cfg) / "mesh.names"))


def resolved(cfg: dict) -> dict:
    """Preset + [solver] overrides, with the automatic linear tolerance."""
    s = dict(PRESETS[cfg["preset"]])
    over = dict(cfg.get("solver", {}))
    if "boomer_amg" in over:
        s["boomer_amg"] = {**s["boomer_amg"], **over.pop("boomer_amg")}
    for key in ("extra_solver_keywords", "steady_extra_solver_keywords"):
        s[key] = list(s[key]) + list(over.pop(key, []))
    s.update(over)
    n_tes = len(tes_bodies(cfg))
    if s["linear_tolerance"] is None:
        s["linear_tolerance"] = 3.0e-11 if n_tes <= 1 else 1.0e-11
    s["n_tes"] = n_tes
    s["window"] = cfg.get("transient", {}).get("window", s["window"])
    return s


def _tag(cfg: dict, extra: str = "") -> str:
    base = re.sub(r"[^A-Za-z0-9]", "", cfg["name"])[:10] or "case"
    # Includes the preset contents, so a changed preset never reuses old runs.
    blob = json.dumps({**{k: cfg.get(k) for k in ("mesh", "roles", "tes", "solver", "preset", "mpi")},
                       "preset_contents": PRESETS[cfg["preset"]]}, sort_keys=True) + extra
    return f"{base}{hashlib.sha1(blob.encode()).hexdigest()[:4]}"


def _common_args(cfg: dict, s: dict) -> list[str]:
    args = [
        "--mesh", cfg["mesh"], "--mpi", str(cfg["mpi"]), "--elmer", s["elmer"],
        "--env", f"LD_LIBRARY_PATH={ELMER_PREFIX}/udf:{ELMER_PREFIX}/cpu/hypre/lib:{ELMER_PREFIX}/cpu/elmer/lib/elmersolver",
        "--env", "PHASE24_HYPRE_REPORT=1",
        "--linear-system", s["linear_system"], "--linear-tolerance", repr(float(s["linear_tolerance"])),
        "--case-option", "membrane_k_udf=true", "--case-option", "phase24_vector_assembly=true",
    ]
    for key, value in s["boomer_amg"].items():
        args += ["--solver-option", f"boomer_amg_{key}={json.dumps(value)}"]
    roles = roles_of(cfg)
    if roles:
        args += ["--case-option", "roles=" + json.dumps(roles)]
    if cfg.get("tes"):
        args += ["--case-option", "tes_circuits=" + json.dumps(cfg["tes"])]
    if cfg.get("materials_temperature_tables"):
        args += ["--case-option", "temperature_tables=" + json.dumps(cfg["materials_temperature_tables"])]
    return args


def _run_runner(args: list[str]) -> dict:
    proc = subprocess.run([sys.executable, str(RUNNER), *args], cwd=ROOT, capture_output=True, text=True,
                          env={**__import__("os").environ, "MSYS_NO_PATHCONV": "1", "MSYS2_ARG_CONV_EXCL": "*"})
    records = [json.loads(line) for line in proc.stdout.splitlines() if line.startswith("{")]
    if proc.returncode or not records:
        sys.stderr.write(proc.stdout[-3000:] + proc.stderr[-3000:])
        raise SystemExit(f"runner failed (exit {proc.returncode})")
    merged = {}
    for r in records:
        merged.update(r)
    return merged


def _wsl(cmd: str) -> subprocess.CompletedProcess:
    return subprocess.run(["wsl.exe", "-d", DISTRO, "-u", "root", "--exec", "bash", "-lc", cmd],
                          capture_output=True, text=True)


def _wsl_path(p: Path) -> str:
    r = p.resolve()
    return f"/mnt/{r.drive[0].lower()}{r.as_posix()[2:]}"


def ensure_partition(cfg: dict) -> None:
    d = mesh_dir(cfg)
    n = int(cfg["mpi"])
    if n <= 1 or (d / f"partitioning.{n}").is_dir():
        return
    print(f"partitioning {d.name} into {n} parts ...", flush=True)
    env = (f"export ELMER_HOME={ELMER_PREFIX}/cpu/elmer LD_LIBRARY_PATH={ELMER_PREFIX}/cpu/elmer/lib/elmersolver:"
           f"{ELMER_PREFIX}/cpu/hypre/lib PATH={ELMER_PREFIX}/cpu/elmer/bin:$PATH")
    r = _wsl(f"{env}; cd '{_wsl_path(d.parent)}' && ElmerGrid 2 2 {d.name} -metiskway {n}")
    if not (d / f"partitioning.{n}").is_dir():
        raise SystemExit("ElmerGrid partitioning failed:\n" + r.stdout[-2000:] + r.stderr[-2000:])


def _state_path(cfg: dict) -> Path:
    return RUNS / cfg["name"] / "state.json"


def _load_state(cfg: dict) -> dict:
    p = _state_path(cfg)
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _save_state(cfg: dict, state: dict) -> None:
    p = _state_path(cfg)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(state, indent=2) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
def cmd_steady(cfg: dict, force: bool = False) -> str:
    s = resolved(cfg)
    tag = _tag(cfg, "steady")
    state = _load_state(cfg)
    if not force and state.get("steady_tag") == tag and state.get("steady_case"):
        print(f"steady: reusing {state['steady_case']}")
        return state["steady_case"]
    ensure_partition(cfg)
    args = ["steady", *_common_args(cfg, s), "--tag", tag,
            "--case-option", "extra_solver_keywords=" + json.dumps(s["steady_extra_solver_keywords"])]
    print(f"steady: {cfg['mesh']}, {s['n_tes']} TES, {cfg['mpi']} ranks ...", flush=True)
    rec = _run_runner(args)
    if rec.get("exit_code"):
        raise SystemExit(f"steady failed; see {rec.get('launcher_log')}")
    state.update(steady_tag=tag, steady_case=rec["case"])
    _save_state(cfg, state)
    print(f"steady: done in {rec['elapsed_s']:.0f} s ({rec['case']})")
    return rec["case"]


def cmd_transient(cfg: dict) -> str:
    s = resolved(cfg)
    steady_case = cmd_steady(cfg)
    tr = cfg.get("transient", {})
    pulse = {**DEFAULT_PULSE, **tr.get("pulse", {})}
    tag = _tag(cfg, "transient" + json.dumps(tr, sort_keys=True))
    args = ["transient", *_common_args(cfg, s), "--tag", tag, "--window", s["window"],
            "--reference-case", steady_case,
            "--solver-option", f"nonlinear_convergence_tolerance={s['nonlinear_tolerance']!r}",
            "--case-option", f"bdf_order={s['bdf_order']}",
            "--case-option", f"phase24_bdf2_predictor={json.dumps(bool(s['predictor']))}",
            "--case-option", "pulse=" + json.dumps(pulse),
            "--case-option", "extra_simulation_keywords=" + json.dumps(
                [f"Phase24 BDF2 Predictor Max Delta = Real {s['predictor_max_delta']!r}"]),
            "--case-option", "extra_solver_keywords=" + json.dumps(s["extra_solver_keywords"])]
    print(f"transient: window {s['window']}, pulse {pulse['center']}, {cfg['mpi']} ranks ...", flush=True)
    rec = _run_runner(args)
    if rec.get("exit_code"):
        raise SystemExit(f"transient failed; see {rec.get('launcher_log')}")
    state = _load_state(cfg)
    state.update(transient_case=rec["case"], transient_elapsed_s=rec["elapsed_s"])
    _save_state(cfg, state)
    print(f"transient: done in {rec['elapsed_s']:.0f} s ({rec['case']})")
    cmd_summary(cfg)
    return rec["case"]


def series_files(case: str) -> list[Path]:
    d = ROOT / "results" / case
    single = d / f"{case}_series.csv"
    if single.exists():
        return [single]
    return sorted(d.glob(f"{case}_[0-9]*_series.csv"), key=lambda p: int(p.stem.split("_")[-2]))


def cmd_summary(cfg: dict) -> None:
    import numpy as np

    state = _load_state(cfg)
    case = state.get("transient_case")
    if not case:
        raise SystemExit("no transient run recorded for this case file")
    rows = []
    names = tes_bodies(cfg)
    for k, path in enumerate(series_files(case), start=1):
        d = np.genfromtxt(path, delimiter=",", names=True)
        t, current = d["time_s"], d["tes_current_A"] * 1e6
        i0 = max(int(np.searchsorted(t, PULSE_START)) - 1, 0)
        ip = int(np.argmin(current))
        rows.append({
            "circuit": k, "body": names[k - 1] if k - 1 < len(names) else "?",
            "baseline_uA": float(current[i0]), "peak_drop_uA": float(current[i0] - current[ip]),
            "peak_time_us": float((t[ip] - PULSE_START) * 1e6), "series": str(path.relative_to(ROOT)),
        })
    print(f"{'k':>2} {'body':16s} {'baseline uA':>12} {'peak drop uA':>13} {'peak time us':>13}")
    for r in rows:
        print(f"{r['circuit']:>2} {r['body']:16s} {r['baseline_uA']:12.4f} {r['peak_drop_uA']:13.4f} {r['peak_time_us']:13.1f}")
    out = RUNS / cfg["name"] / "summary.json"
    out.write_text(json.dumps({"case": case, "elapsed_s": state.get("transient_elapsed_s"), "circuits": rows},
                              indent=2) + "\n", encoding="utf-8")
    print(f"summary: {out.relative_to(ROOT)}")


def cmd_show(cfg: dict) -> None:
    s = resolved(cfg)
    print(json.dumps({"case_file": cfg, "mesh_dir": str(mesh_dir(cfg).relative_to(ROOT)),
                      "tes_circuits": tes_bodies(cfg), "roles": roles_of(cfg), "resolved_solver": s},
                     indent=2, ensure_ascii=False))


def cmd_mesh_import(msh: Path, name: str, roles_file: Path | None, notes: str) -> None:
    """Convert a Gmsh mesh (ElmerGrid, coincident nodes merged) and register it."""
    target = ROOT / "work" / "meshes" / name
    env = (f"export ELMER_HOME={ELMER_PREFIX}/cpu/elmer LD_LIBRARY_PATH={ELMER_PREFIX}/cpu/elmer/lib/elmersolver:"
           f"{ELMER_PREFIX}/cpu/hypre/lib PATH={ELMER_PREFIX}/cpu/elmer/bin:$PATH")
    r = _wsl(f"{env}; cd '{_wsl_path(ROOT)}' && ElmerGrid 14 2 '{_wsl_path(msh)}' -merge 1e-10 -out '{_wsl_path(target)}'")
    if not (target / "mesh.header").exists():
        raise SystemExit("ElmerGrid conversion failed:\n" + r.stdout[-2000:] + r.stderr[-2000:])
    data = _registry()
    entry = {"geometry": "custom", "dir": name, "notes": notes or f"imported from {msh}",
             "recipe": {"commands": [f"ElmerGrid 14 2 {msh.as_posix()} -merge 1e-10 -out work/meshes/{name}"]}}
    if roles_file:
        entry["roles"] = tomllib.loads(roles_file.read_text(encoding="utf-8"))
    data["meshes"][name] = entry
    MESH_REGISTRY.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"registered {name} ({(target / 'mesh.header').read_text().split()[0]} nodes)")
    print((target / "mesh.names").read_text())


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "steady", "transient", "show", "summary"):
        p = sub.add_parser(name)
        p.add_argument("case_file", type=Path)
        if name == "steady":
            p.add_argument("--force", action="store_true", help="rerun even if a matching steady exists")
    p = sub.add_parser("mesh-import")
    p.add_argument("msh", type=Path)
    p.add_argument("--name", required=True)
    p.add_argument("--roles", type=Path, default=None, help="TOML role table stored with the mesh")
    p.add_argument("--notes", default="")
    a = ap.parse_args(argv)
    if a.cmd == "mesh-import":
        cmd_mesh_import(a.msh, a.name, a.roles, a.notes)
        return 0
    cfg = load_case(a.case_file)
    if a.cmd in ("run", "transient"):
        cmd_transient(cfg)
    elif a.cmd == "steady":
        cmd_steady(cfg, force=a.force)
    elif a.cmd == "show":
        cmd_show(cfg)
    else:
        cmd_summary(cfg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
