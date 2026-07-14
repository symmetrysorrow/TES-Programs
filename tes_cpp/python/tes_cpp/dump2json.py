from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from shutil import which


@dataclass(frozen=True)
class Dump2JsonResult:
    output_path: Path
    full_energy_event_ids: tuple[int, ...]


def _executable() -> str:
    configured = os.environ.get("TES_CPP_DUMP2JSON_EXECUTABLE")
    if configured:
        return configured
    name = "dump2json.exe" if os.name == "nt" else "dump2json"
    bundled = files("tes_cpp").joinpath("bin", name)
    if bundled.is_file():
        return str(bundled)
    found = which("dump2json")
    if found:
        return found
    raise RuntimeError("dump2json executable was not found; set TES_CPP_DUMP2JSON_EXECUTABLE")


def dump2json(
    dump_path: str | Path,
    output_path: str | Path,
    *,
    input_energy: float,
    save_all: bool = False,
) -> Dump2JsonResult:
    """Convert ``dumpall.dat`` into the legacy-compatible ``batch.json`` schema."""
    output = Path(output_path)
    command = [_executable(), str(dump_path), str(output), "--input-energy", str(input_energy)]
    index_path = output.with_name("FullEnergyList.dat")
    if save_all:
        command += ["--save-all", "--full-energy-list", str(index_path)]
    subprocess.run(command, check=True)
    event_ids = ()
    if save_all and index_path.exists():
        event_ids = tuple(int(line) for line in index_path.read_text().splitlines() if line)
    return Dump2JsonResult(output_path=output, full_energy_event_ids=event_ids)
