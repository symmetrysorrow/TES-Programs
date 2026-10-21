"""Standalone wrapper around the Pb/EGS5 physical HDF5 v2 converter."""
from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from shutil import which


@dataclass(frozen=True)
class ConversionResult:
    output_path: Path


def _executable() -> str:
    configured = os.environ.get("TES_CPP_DUMP2EVENT_EXECUTABLE")
    if configured:
        return configured
    name = "dump2event.exe" if os.name == "nt" else "dump2event"
    bundled = files("tes_cpp_dump2event").joinpath("bin", name)
    if bundled.is_file():
        return str(bundled)
    found = which("dump2event")
    if found:
        return found
    raise RuntimeError("dump2event executable not found; set TES_CPP_DUMP2EVENT_EXECUTABLE")


def convert(
    dump_path: str | Path,
    output_path: str | Path,
    *,
    physics_config: str | Path | None = None,
    stopping_table: str | Path | None = None,
    physics: str = "pb-egs5",
    **removed_legacy_options: object,
) -> ConversionResult:
    """Convert PHITS 3.341 Pb/EGS5 Dumpall to HDF5 v2 (only supported mode)."""
    if removed_legacy_options:
        raise ValueError("obsolete converter options are unsupported: " +
                         ", ".join(sorted(removed_legacy_options)))
    if physics != "pb-egs5":
        raise ValueError("only PHITS 3.341 Pb EGS5 physics is supported")
    if physics_config is None or stopping_table is None:
        raise ValueError("Pb EGS5 requires physics_config and stopping_table")
    dump = Path(dump_path)
    config = Path(physics_config)
    table = Path(stopping_table)
    output = Path(output_path)
    if output.suffix.lower() not in {".h5", ".hdf5"}:
        raise ValueError("output must be HDF5 v2 (.h5/.hdf5)")
    for name, path in (("Dumpall", dump), ("physics_config", config), ("stopping_table", table)):
        if not path.is_file():
            raise FileNotFoundError(f"{name} does not exist: {path}")
    if output.exists():
        raise FileExistsError(f"will not overwrite existing HDF5: {output}")
    subprocess.run([_executable(), str(dump), str(output),
                    "--physics-config", str(config),
                    "--stopping-table", str(table)], check=True)
    return ConversionResult(output_path=output)
