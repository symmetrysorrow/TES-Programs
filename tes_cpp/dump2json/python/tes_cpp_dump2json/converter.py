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
    full_energy_event_ids: tuple[int, ...]


def _executable() -> str:
    configured = os.environ.get("TES_CPP_DUMP2JSON_EXECUTABLE")
    if configured:
        return configured
    bundled = files("tes_cpp_dump2json").joinpath("bin", "dump2json.exe" if os.name == "nt" else "dump2json")
    if bundled.is_file():
        return str(bundled)
    command = which("dump2json")
    if command:
        return command
    raise RuntimeError("dump2json executable was not found; set TES_CPP_DUMP2JSON_EXECUTABLE")


def convert(
    dump_path: str | Path,
    output_path: str | Path,
    *,
    input_energy: float,
    save_all: bool = False,
) -> ConversionResult:
    """Convert one PHITS ``dumpall.dat`` file to the legacy ``batch.json`` schema."""
    dump_path, output_path = Path(dump_path), Path(output_path)
    command = [_executable(), str(dump_path), str(output_path), "--input-energy", str(input_energy)]
    full_energy_path = output_path.with_name("FullEnergyList.dat")
    if save_all:
        command += ["--save-all", "--full-energy-list", str(full_energy_path)]
    subprocess.run(command, check=True)
    event_ids = ()
    if save_all and full_energy_path.exists():
        event_ids = tuple(int(line) for line in full_energy_path.read_text().splitlines() if line)
    return ConversionResult(output_path=output_path, full_energy_event_ids=event_ids)
