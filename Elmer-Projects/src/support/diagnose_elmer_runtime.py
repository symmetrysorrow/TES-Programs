"""Record the Elmer runtime selected for a reproducible hybrid-mesh test.

This is deliberately read-only: it neither builds nor copies DLLs.  It records
the executable, core solver DLL and HeatSolve module actually selected from an
Elmer install prefix, plus PE import names and the DLLs resolvable from the
process search path.  The resulting JSON can be compared before every A/B run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import struct
import sys
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TOOLS_ELMER = ROOT.parent / "tools" / "elmer-hypre"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pe_imports(path: Path) -> list[str]:
    """Read PE import-library names without requiring Visual Studio tools."""
    data = path.read_bytes()
    if data[:2] != b"MZ":
        return []
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    if data[pe : pe + 4] != b"PE\0\0":
        return []
    sections = struct.unpack_from("<H", data, pe + 6)[0]
    optional = struct.unpack_from("<H", data, pe + 20)[0]
    opt = pe + 24
    magic = struct.unpack_from("<H", data, opt)[0]
    data_dir = opt + (112 if magic == 0x20B else 96)
    import_rva, _ = struct.unpack_from("<II", data, data_dir + 8)
    section_table = opt + optional
    table: list[tuple[int, int, int]] = []
    for index in range(sections):
        base = section_table + index * 40
        virtual_size, virtual_address, raw_size, raw_offset = struct.unpack_from(
            "<IIII", data, base + 8
        )
        table.append((virtual_address, max(virtual_size, raw_size), raw_offset))

    def offset(rva: int) -> int:
        for virtual_address, size, raw_offset in table:
            if virtual_address <= rva < virtual_address + size:
                return raw_offset + rva - virtual_address
        raise ValueError(f"RVA {rva:#x} is outside PE sections")

    if not import_rva:
        return []
    imports: list[str] = []
    cursor = offset(import_rva)
    while cursor + 20 <= len(data):
        original_first, _, _, name_rva, first_thunk = struct.unpack_from("<IIIII", data, cursor)
        if not any((original_first, name_rva, first_thunk)):
            break
        name_offset = offset(name_rva)
        end = data.find(b"\0", name_offset)
        if end < 0:
            break
        imports.append(data[name_offset:end].decode("ascii", errors="replace"))
        cursor += 20
    return imports


def runtime_paths(solver: Path, prefix: Path) -> list[Path]:
    paths = [solver.parent, prefix / "bin", prefix / "lib", prefix / "share" / "elmersolver" / "lib"]
    paths += [Path(part) for part in os.environ.get("PATH", "").split(os.pathsep) if part]
    unique: list[Path] = []
    seen: set[str] = set()
    for candidate in paths:
        try:
            key = str(candidate.resolve()).lower()
        except OSError:
            key = str(candidate).lower()
        if key not in seen:
            unique.append(candidate)
            seen.add(key)
    return unique


def resolve_dll(name: str, paths: list[Path]) -> str | None:
    for directory in paths:
        candidate = directory / name
        if candidate.is_file():
            return str(candidate.resolve())
    return None


def file_record(path: Path, search_paths: list[Path]) -> dict:
    stat = path.stat()
    imports = pe_imports(path)
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "size_bytes": stat.st_size,
        "modified_utc": datetime.fromtimestamp(stat.st_mtime, timezone.utc).isoformat(),
        "imports": imports,
        "resolved_imports": {name: resolve_dll(name, search_paths) for name in imports},
    }


def find_heatsolve(prefix: Path, solver: Path) -> Path | None:
    candidates = [
        prefix / "share" / "elmersolver" / "lib" / "HeatSolve.dll",
        prefix / "bin" / "HeatSolve.dll",
        solver.parent / "HeatSolve.dll",
    ]
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def known_build_matches(path: Path) -> list[str]:
    """Return build trees containing an identical named binary, if any."""
    if not TOOLS_ELMER.is_dir():
        return []
    expected = sha256(path)
    matches: list[str] = []
    for build in sorted(TOOLS_ELMER.glob("build*")):
        if not build.is_dir():
            continue
        for candidate in build.rglob(path.name):
            if candidate.is_file() and sha256(candidate) == expected:
                matches.append(str(build.resolve()))
                break
    return matches


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--prefix", type=Path, help="Elmer install prefix (defaults to solver/..)")
    parser.add_argument("--extra-dll", action="append", default=[], type=Path,
                        help="external UDF DLL to hash and inspect (excluded from install provenance gate)")
    args = parser.parse_args()
    solver = args.solver.resolve()
    if not solver.is_file():
        raise SystemExit(f"ElmerSolver not found: {solver}")
    prefix = (args.prefix or solver.parent.parent).resolve()
    paths = runtime_paths(solver, prefix)
    core = solver.parent / "libelmersolver.dll"
    if not core.is_file():
        resolved = resolve_dll("libelmersolver.dll", paths)
        core = Path(resolved) if resolved else core
    heat = find_heatsolve(prefix, solver)
    records = {"ElmerSolver": file_record(solver, paths)}
    if core.is_file():
        records["libelmersolver"] = file_record(core, paths)
    if heat:
        records["HeatSolve"] = file_record(heat, paths)
    extras: dict[str, dict] = {}
    for dll in args.extra_dll:
        if not dll.is_file():
            raise SystemExit(f"extra DLL not found: {dll}")
        extras[dll.name] = file_record(dll.resolve(), paths)
    for record in records.values():
        record["matching_build_trees"] = known_build_matches(Path(record["path"]))
    known = [set(record["matching_build_trees"]) for record in records.values()]
    known_nonempty = [item for item in known if item]
    shared_build = set.intersection(*known_nonempty) if len(known_nonempty) >= 2 else set()
    critical_external: dict[str, str] = {}
    for record in records.values():
        for name, resolved in record["resolved_imports"].items():
            if name.lower().startswith(("libgfortran", "libgcc", "libquadmath", "libgomp", "libmumps", "libhypre")) and resolved:
                if not str(resolved).lower().startswith(str(prefix).lower()):
                    critical_external[name] = resolved
    result = {
        "schema": 1,
        "label": args.label,
        "captured_utc": datetime.now(timezone.utc).isoformat(),
        "prefix": str(prefix),
        "runtime_search_paths": [str(path) for path in paths],
        "artifacts": records,
        "external_udfs": extras,
        "assessment": {
            "required_artifacts_present": sorted(records) == ["ElmerSolver", "HeatSolve", "libelmersolver"],
            "same_install_prefix": all(str(item["path"]).lower().startswith(str(prefix).lower()) for item in records.values()),
            "shared_known_build_trees": sorted(shared_build),
            "same_known_build_tree": bool(shared_build) if len(known_nonempty) >= 2 else None,
            "critical_runtime_dlls_outside_prefix": critical_external,
            "note": "same_install_prefix is a layout/provenance gate, not a binary ABI proof.",
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    print(json.dumps(result["assessment"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
