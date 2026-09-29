"""Entry point for Elmer-Projects.

New simulations use projects/*.toml. Legacy JSON and old TOML workflows are
available only through the explicit legacy namespace.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))


def root_path(value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def read_legacy_project(value: str) -> dict:
    path = root_path(value)
    if not path.is_file():
        raise ValueError(f"legacy project JSON not found: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("legacy project JSON must contain an object")
    return data


def check_legacy_project(data: dict) -> list[str]:
    from support.reconcile_project import reconcile_project

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
    return subprocess.run([sys.executable, str(ROOT / "src" / script), *args], cwd=ROOT).returncode


def resolve_project_toml(value: str) -> Path:
    from support.scenario_config import PROJECT_DIR, read_toml

    if value.lower().endswith(".toml"):
        path = root_path(value)
    else:
        path = PROJECT_DIR / f"{value}.toml"

    if not path.is_file():
        raise ValueError(
            f"unknown project: {value!r}; use 'python main.py list'. "
            "Legacy inputs are under 'python main.py legacy ...'."
        )

    data = read_toml(path)
    if data.get("schema_version") != 3:
        display = path.relative_to(ROOT) if path.is_relative_to(ROOT) else path
        raise ValueError(
            f"{display} is not a schema_version = 3 project; "
            "use 'python main.py legacy ...' for old inputs."
        )
    return path


def project_action(command: str, value: str, *, dry_run: bool = False, mpi_procs: int = 1) -> int:
    from support.scenario_cli import main as scenario_main

    action = {"show": "summary", "check": "validate", "mesh": "mesh", "run": "run"}[command]
    project_file = resolve_project_toml(value)
    args = [action, str(project_file), "--mpi-procs", str(mpi_procs)]
    if dry_run:
        args.append("--dry-run")
    return scenario_main(args)


def add_legacy_json_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--project",
        default="projects/legacy/elmer_project.json",
        help="legacy project JSON to use",
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="TES simulation projects. New work starts in projects/*.toml.",
        epilog=(
            "Start with: python main.py list; python main.py show single_pixel_alpha240; "
            "python main.py run single_pixel_alpha240 --dry-run"
        ),
    )
    parser.add_argument(
        "-i",
        "--interactive",
        action="store_true",
        help="open the guided menu (also used when no command is given)",
    )
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("list", help="list schema-v3 projects")

    check = sub.add_parser("check", help="validate a project without calculating")
    check.add_argument("name", help="project name or projects/<name>.toml")

    show = sub.add_parser("show", help="show model, TES IDs/positions and final run values")
    show.add_argument("name", help="project name or projects/<name>.toml")

    mesh = sub.add_parser("mesh", help="build or reuse the mesh for a project")
    mesh.add_argument("name", help="project name or projects/<name>.toml")

    run = sub.add_parser("run", help="run a project; build or reuse its mesh automatically")
    run.add_argument("name", help="project name or projects/<name>.toml")
    run.add_argument("--dry-run", action="store_true", help="print the run plan only")
    run.add_argument("--mpi-procs", type=int, default=1, help="number of MPI ranks")

    legacy = sub.add_parser(
        "legacy",
        help="reproduce old JSON or projects/legacy/cases/*.toml workflows",
        description="Legacy workflows. New projects should not use these commands.",
    )
    legacy_sub = legacy.add_subparsers(dest="legacy_command")

    legacy_list = legacy_sub.add_parser("list", help="list legacy TOML and JSON entries")
    add_legacy_json_options(legacy_list)

    legacy_check = legacy_sub.add_parser("check", help="validate a legacy project JSON")
    add_legacy_json_options(legacy_check)

    legacy_show = legacy_sub.add_parser("show", help="show a legacy JSON case or mesh")
    legacy_show.add_argument("name")
    add_legacy_json_options(legacy_show)

    legacy_mesh = legacy_sub.add_parser("mesh", help="build a legacy JSON mesh")
    legacy_mesh.add_argument("name")
    legacy_mesh.add_argument("--record-only", action="store_true", help="record provenance for an existing mesh")
    add_legacy_json_options(legacy_mesh)

    legacy_run = legacy_sub.add_parser("run", help="run a legacy JSON case")
    legacy_run.add_argument("name")
    legacy_run.add_argument("--dry-run", action="store_true", help="print the legacy run plan only")
    legacy_run.add_argument("--force-deps", action="store_true", help="rerun restart dependencies")
    legacy_run.add_argument("--skip-sync", action="store_true", help="skip SIF regeneration")
    legacy_run.add_argument("--mpi-procs", type=int, default=1, help="number of MPI ranks")
    add_legacy_json_options(legacy_run)

    legacy_toml = legacy_sub.add_parser("toml", help="use the old projects/legacy/cases/*.toml / tes_sim.py workflow")
    legacy_toml.add_argument("action", choices=("show", "run", "steady", "transient", "summary"))
    legacy_toml.add_argument("case_file", help="example: projects/legacy/cases/single_pixel_h8.toml")
    return parser


def interactive_main() -> int:
    """Run a guided menu for users who do not want to enter CLI arguments."""
    try:
        import questionary
    except ImportError:
        install_command = f'"{sys.executable}" -m pip install questionary'
        if sys.platform == "win32":
            install_command = f"& {install_command}"
        print(
            "The interactive menu requires questionary for this Python interpreter.\n"
            f"Install it with: {install_command}",
            file=sys.stderr,
        )
        return 2

    from support.scenario_config import PROJECT_DIR, read_toml

    projects = [
        path.stem
        for path in sorted(PROJECT_DIR.glob("*.toml"))
        if read_toml(path).get("schema_version") == 3
    ]
    actions = {
        "プロジェクトの内容を見る": "show",
        "プロジェクトを検証する": "check",
        "計算を実行する": "run",
        "終了": "exit",
    }

    try:
        while True:
            selected_action = questionary.select(
                "何をしますか？",
                choices=list(actions),
            ).ask()
            command = actions.get(selected_action)
            if command in (None, "exit"):
                return 0

            if not projects:
                print("schema-v3 project が見つかりません。projects/ を確認してください。")
                return 1

            project = questionary.select(
                "プロジェクトを選んでください。",
                choices=projects,
            ).ask()
            if project is None:
                return 0

            dry_run = False
            if command == "run":
                run_mode = questionary.select(
                    "実行方法を選んでください。",
                    choices=[
                        questionary.Choice("実行計画だけ確認する（dry-run）", value="dry-run"),
                        questionary.Choice("計算を実行する", value="run"),
                        questionary.Choice("戻る", value="cancel"),
                    ],
                ).ask()
                if run_mode in (None, "cancel"):
                    continue
                dry_run = run_mode == "dry-run"

            try:
                result = project_action(command, project, dry_run=dry_run)
            except (ValueError, OSError, json.JSONDecodeError) as exc:
                print(f"error: {exc}", file=sys.stderr)
                result = 2

            if result != 0:
                return result
    except (KeyboardInterrupt, EOFError):
        print()
        return 0


def legacy_main(args: argparse.Namespace) -> int:
    if not args.legacy_command:
        print("Legacy workflows:")
        print("  python main.py legacy list")
        print("  python main.py legacy show <json-case-or-mesh>")
        print("  python main.py legacy run <json-case>")
        print("  python main.py legacy toml show projects/legacy/cases/<name>.toml")
        return 0

    if args.legacy_command == "toml":
        path = root_path(args.case_file)
        if not path.is_file():
            raise ValueError(f"legacy TOML case not found: {path}")
        return dispatch("tes_sim.py", [args.action, str(path)])

    data = read_legacy_project(args.project)

    if args.legacy_command == "list":
        print("Legacy TOML cases:")
        for path in sorted((ROOT / "projects" / "legacy" / "cases").glob("*.toml")):
            print(f"  {path.relative_to(ROOT)}")
        for label, key in (
            ("Legacy JSON geometries", "geometries"),
            ("Legacy JSON meshes", "meshes"),
            ("Legacy JSON cases", "cases"),
        ):
            print(f"{label} ({len(data.get(key, {}))}):")
            for name, entry in data.get(key, {}).items():
                detail = (
                    entry.get("geometry", "")
                    if key == "meshes"
                    else entry.get("mesh", "")
                    if key == "cases"
                    else ""
                )
                print(f"  {name}" + (f"  [{detail}]" if detail else ""))
        return 0

    if args.legacy_command == "check":
        errors = check_legacy_project(data)
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

    if args.legacy_command == "show":
        if name in cases:
            print(json.dumps({"case": name, **cases[name]}, ensure_ascii=False, indent=2))
        elif name in meshes:
            print(json.dumps({"mesh": name, **meshes[name]}, ensure_ascii=False, indent=2))
        else:
            raise ValueError(f"unknown legacy case or mesh: {name}; use 'python main.py legacy list'")
        return 0

    if args.legacy_command == "mesh":
        if name not in meshes:
            raise ValueError(f"unknown legacy mesh: {name}; use 'python main.py legacy list'")
        cmd = [name, "--project", str(root_path(args.project))]
        if args.record_only:
            cmd.append("--record-only")
        return dispatch("build_mesh.py", cmd)

    if name not in cases:
        raise ValueError(f"unknown legacy case: {name}; use 'python main.py legacy list'")
    cmd = [name, "--project", str(root_path(args.project)), "--mpi-procs", str(args.mpi_procs)]
    for flag in ("dry_run", "force_deps", "skip_sync"):
        if getattr(args, flag):
            cmd.append("--" + flag.replace("_", "-"))
    if args.dry_run and not args.skip_sync:
        cmd.append("--skip-sync")
    return dispatch("run.py", cmd)


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.interactive or (not args.command and sys.stdin.isatty()):
        return interactive_main()
    if not args.command:
        print("TES simulation: edit a file in projects/, then run it.\n")
        print("  python main.py list")
        print("  python main.py show single_pixel_alpha240")
        print("  python main.py run single_pixel_alpha240 --dry-run")
        print("  python main.py run single_pixel_alpha240")
        print("\nOld workflows: python main.py legacy ...")
        return 0

    try:
        if args.command == "legacy":
            return legacy_main(args)

        if args.command == "list":
            from support.scenario_config import PROJECT_DIR, read_toml

            print("Projects (use 'show <name>' or 'run <name>'):")
            for path in sorted(PROJECT_DIR.glob("*.toml")):
                if read_toml(path).get("schema_version") == 3:
                    print(f"  {path.stem}")
            return 0

        return project_action(
            args.command,
            args.name,
            dry_run=getattr(args, "dry_run", False),
            mpi_procs=getattr(args, "mpi_procs", 1),
        )
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
