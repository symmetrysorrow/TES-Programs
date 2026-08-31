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
    full_energy_event_ids: tuple[int, ...]


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
    input_energy: float,
    save_all: bool = False,
    full_energy_only: bool = False,
) -> Dump2EventResult:
    """Convert ``dumpall.dat`` into event JSON or ``.h5`` data.

    HDF5 conversion is streaming, so Python never loads the full native JSON
    document into memory.
    """
    output = Path(output_path)
    command = [_executable(), str(dump_path), str(output), "--input-energy", str(input_energy)]
    index_path = output.with_name("FullEnergyList.dat")
    if save_all:
        command += ["--save-all"]
    if save_all or full_energy_only:
        command += ["--full-energy-list", str(index_path)]
    if full_energy_only:
        command += ["--full-energy-only"]
    subprocess.run(command, check=True)
    event_ids = ()
    if (save_all or full_energy_only) and index_path.exists():
        event_ids = tuple(int(line) for line in index_path.read_text().splitlines() if line)
    return Dump2EventResult(output_path=output, full_energy_event_ids=event_ids)
