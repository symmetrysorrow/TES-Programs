"""Entry point for inspecting and running Elmer-Projects.

New work starts from a TOML file in projects/. Legacy JSON and old TOML
workflows remain available for reproducing older calculations.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def project_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def read_project(value: str) -> dict:
    path = project_path(value)
    if not path.is_file():
        raise ValueError(f"project JSON not found: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("project JSON must contain an object")
    return data


def check_project(data: dict) -> list[str]:
    from scripts.support.reconcile_project import reconcile_project

    errors = []
    model = reconcile_project(data)
    geometries = model.get("geometries", {})
    meshes = model.get("meshes", {})
    cases = model.get("cases", {})
    for name, mesh in meshes.items():
        if mesh.get("geometry") not in geometries:
            errors.append(f"mesh {name}: unknown geometry {mesh.get('geometry')!r}")
        if not mesh.get("dir"):
            errors.append(f"mesh {name}: missing output directory")
    for name, case in cases.items():
        if case.get("mesh") not in meshes:
            errors.append(f"case {name}: unknown mesh {case.get('mesh')!r}")
        dependency = case.get("restart_from")
        if dependency and dependency not in cases:
            errors.append(f"case {name}: unknown restart case {dependency!r}")
    return errors


def dispatch(script: str, args: list[str]) -> int:
    return subprocess.run([sys.executable, str(ROOT / script), *args], cwd=ROOT).returncode


def scenario_action(command: str, value: str, *, dry_run: bool = False, mpi_procs: int = 1) -> int:
    from scripts.support.scenario_cli import main as scenario_main

    action = {"show": "summary", "check": "validate", "mesh": "mesh", "run": "run"}[command]
    target = str(project_path(value)) if value.lower().endswith(".toml") else value
    args = [action, target, "--mpi-procs", str(mpi_procs)]
    if dry_run:
        args.append("--dry-run")
    return scenario_main(args)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="TES simulation: choose a project in projects/, inspect it, then run it.",
        epilog="Start with: python main.py show single_pixel_alpha240; python main.py run single_pixel_alpha240 --dry-run",
    )
    sub = parser.add_subparsers(dest="command")

    for name, help_text in (
        ("list", "list available projects"),
        ("check", "validate a named project without calculating"),
        ("show", "show model, TES positions and final run values"),
        ("mesh", "build or reuse a project mesh (legacy JSON mesh names also work)"),
        ("run", "run a project; build or reuse its mesh automatically"),
    ):
        p = sub.add_parser(name, help=help_text)
        p.add_argument("--project", default="elmer_project.json", help="legacy project JSON to use")
        if name == "list":
            p.add_argument("--all", action="store_true", help="include legacy TOML and JSON entries")
        if name in {"show", "mesh", "run"}:
            p.add_argument("name", help="project name from projects/, or a legacy JSON case/mesh name")
        if name == "check":
            p.add_argument("name", nargs="?", help="project name (omit for legacy project JSON)")
        if name == "mesh":
            p.add_argument("--record-only", action="store_true", help="legacy JSON only: record provenance for an existing mesh")
        if name == "run":
            p.add_argument("--dry-run", action="store_true", help="print the run plan only")
            p.add_argument("--force-deps", action="store_true", help="legacy JSON only: rerun restart dependencies")
            p.add_argument("--skip-sync", action="store_true", help="legacy JSON only: skip SIF regeneration")
            p.add_argument("--mpi-procs", type=int, default=1, help="number of MPI ranks")

    toml = sub.add_parser("toml", help="use the legacy cases/*.toml / tes_sim.py workflow")
    toml.add_argument("action", choices=("show", "run", "steady", "transient", "summary"))
    toml.add_argument("case_file", help="example: cases/single_pixel_h8.toml")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        print("TES simulation: edit a file in projects/, then run it.\n")
        print("  python main.py list")
        print("  python main.py show single_pixel_alpha240")
        print("  python main.py run single_pixel_alpha240 --dry-run")
        print("  python main.py run single_pixel_alpha240")
        return 0

    try:
        if args.command == "toml":
            path = project_path(args.case_file)
            if not path.is_file():
                raise ValueError(f"legacy TOML case not found: {path}")
            return dispatch("tes_sim.py", [args.action, str(path)])

        from scripts.support.scenario_config import PROJECT_DIR, read_toml

        if args.command == "list":
            print("Projects (use 'show <name>' or 'run <name>'):")
            for path in sorted(PROJECT_DIR.glob("*.toml")):
                if read_toml(path).get("schema_version") == 3:
                    print(f"  {path.stem}")
            if not args.all:
                print("Use --all to show legacy TOML and JSON entries.")
                return 0

            data = read_project(args.project)
            print("Legacy TOML cases:")
            for path in sorted((ROOT / "cases").glob("*.toml")):
                print(f"  {path.relative_to(ROOT)}")
            for label, key in (
                ("Legacy JSON geometries", "geometries"),
                ("Legacy JSON meshes", "meshes"),
                ("Legacy JSON cases", "cases"),
            ):
                print(f"{label} ({len(data.get(key, {}))}):")
                for name, entry in data.get(key, {}).items():
                    detail = entry.get("geometry", "") if key == "meshes" else entry.get("mesh", "") if key == "cases" else ""
                    print(f"  {name}" + (f"  [{detail}]" if detail else ""))
            return 0

        if args.command in {"show", "check", "mesh", "run"} and args.name and not args.name.lower().endswith(".toml"):
            candidate = PROJECT_DIR / f"{args.name}.toml"
            if candidate.is_file() and read_toml(candidate).get("schema_version") == 3:
                if args.command == "mesh" and getattr(args, "record_only", False):
                    raise ValueError("--record-only is only for legacy JSON meshes")
                return scenario_action(
                    args.command,
                    args.name,
                    dry_run=getattr(args, "dry_run", False),
                    mpi_procs=getattr(args, "mpi_procs", 1),
                )

        if args.command in {"show", "check", "mesh", "run"} and args.name and args.name.lower().endswith(".toml"):
            path = project_path(args.name)
            if not path.is_file():
                raise ValueError(f"TOML file not found: {path}")
            project = read_toml(path)
            if project.get("schema_version") == 3:
                if args.command == "mesh" and getattr(args, "record_only", False):
                    raise ValueError("--record-only is only for legacy JSON meshes")
                return scenario_action(
                    args.command,
                    args.name,
                    dry_run=getattr(args, "dry_run", False),
                    mpi_procs=getattr(args, "mpi_procs", 1),
                )
            if args.command == "check":
                raise ValueError(f"{args.name} is a legacy TOML; use: python main.py toml show {args.name}")
            if args.command == "mesh":
                raise ValueError("legacy TOML mesh operations use tes_sim.py, not main.py mesh")
            if getattr(args, "dry_run", False):
                raise ValueError("--dry-run is unavailable for legacy TOML; no calculation was started")
            return dispatch("tes_sim.py", [args.command, str(path)])

        data = read_project(args.project)
        if args.command == "check":
            errors = check_project(data)
            if errors:
                for error in errors:
                    print(error, file=sys.stderr)
                return 1
            print(
                f"OK: {len(data.get('geometries', {}))} geometries, "
                f"{len(data.get('meshes', {}))} meshes, {len(data.get('cases', {}))} cases"
            )
            return 0

        name = args.name
        cases, meshes = data.get("cases", {}), data.get("meshes", {})
        if args.command == "show":
            if name in cases:
                print(json.dumps({"case": name, **cases[name]}, ensure_ascii=False, indent=2))
            elif name in meshes:
                print(json.dumps({"mesh": name, **meshes[name]}, ensure_ascii=False, indent=2))
            else:
                raise ValueError(f"unknown project, legacy case or mesh: {name}; use 'list --all'")
            return 0

        if args.command == "mesh":
            if name not in meshes:
                raise ValueError(f"unknown project or legacy mesh: {name}; use 'list --all'")
            cmd = [name, "--project", str(project_path(args.project))]
            if args.record_only:
                cmd.append("--record-only")
            return dispatch("build_mesh.py", cmd)

        if name not in cases:
            raise ValueError(f"unknown project or legacy case: {name}; use 'list --all'")
        cmd = [name, "--project", str(project_path(args.project)), "--mpi-procs", str(args.mpi_procs)]
        for flag in ("dry_run", "force_deps", "skip_sync"):
            if getattr(args, flag):
                cmd.append("--" + flag.replace("_", "-"))
        if args.dry_run and not args.skip_sync:
            cmd.append("--skip-sync")
        return dispatch("run.py", cmd)
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
