"""PHITS 3.341 / Pb EGS5 -> HDF5 v2 physical converter.

The old signed reaction-difference / normalization converter is removed.
Only the bounded physical score reconstruction is supported.  No automatic
TES position allocation or Pulse_Ms handoff is performed.
"""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from shutil import which


@dataclass(frozen=True)
class Dump2EventResult:
    output_path: Path


def _executable() -> str:
    configured = os.environ.get("TES_CPP_DUMP2EVENT_EXECUTABLE")
    if configured:
        return configured
    name = "dump2event.exe" if os.name == "nt" else "dump2event"
    bundled = files("tes_cpp").joinpath("bin", name)
    if bundled.is_file():
        return str(bundled)
    repo_root = Path(__file__).resolve().parents[3]
    for build_root in (repo_root / "tes_cpp" / "build", repo_root / "build"):
        if build_root.is_dir():
            matches = sorted(build_root.rglob(name))
            if matches:
                return str(matches[0])
    found = which("dump2event")
    if found:
        return found
    raise RuntimeError("dump2event executable was not found; set TES_CPP_DUMP2EVENT_EXECUTABLE")


def dump2event(
    dump_path: str | Path,
    output_path: str | Path,
    *,
    physics_config: str | Path | None = None,
    stopping_table: str | Path | None = None,
    physics: str = "pb-egs5",
    **removed_legacy_options: object,
) -> Dump2EventResult:
    """Convert a validated Pb/EGS5 Dumpall into an HDF5 v2 event ledger.

    Required input contract: PHITS 3.341, Pb region 1000, material 2,
    EGS5 medium 2, 1-keV charged cutoffs, unit weights and the matching
    private stopping-power table.  The native driver validates the PHITS
    input and fails closed for unsupported physics.

    Any former legacy option (including input_energy, save_all, and
    full_energy_only) raises instead of selecting an unphysical formula.
    """
    if removed_legacy_options:
        raise ValueError(
            "old dump2event options were removed: "
            + ", ".join(sorted(removed_legacy_options))
            + "; use physics_config and stopping_table for HDF5 v2"
        )
    if physics != "pb-egs5":
        raise ValueError("only PHITS 3.341 Pb/EGS5 physical reconstruction is supported")
    if physics_config is None:
        raise ValueError("Pb EGS5 requires physics_config")
    dump = Path(dump_path)
    config = Path(physics_config)
    table = Path(stopping_table) if stopping_table is not None else None
    output = Path(output_path)
    if output.suffix.lower() not in {".h5", ".hdf5"}:
        raise ValueError("physical dump2event output must be HDF5 v2 (.h5/.hdf5)")
    for name, path in (("Dumpall", dump), ("physics_config", config)):
        if not path.is_file():
            raise FileNotFoundError(f"{name} does not exist: {path}")
    if table is not None and not table.is_file():
        raise FileNotFoundError(f"explicit stopping_table does not exist: {table}")
    if output.exists():
        raise FileExistsError(f"output already exists; will not overwrite: {output}")
    command = [_executable(), str(dump), str(output), "--physics-config", str(config)]
    if table is not None:
        command += ["--stopping-table", str(table)]
    subprocess.run(command, check=True)
    return Dump2EventResult(output_path=output)
