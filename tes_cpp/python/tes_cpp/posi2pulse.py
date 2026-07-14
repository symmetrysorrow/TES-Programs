from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from shutil import which
from tempfile import NamedTemporaryFile


@dataclass(frozen=True)
class Pulse:
    position: int
    time: list[float]
    ch0: list[float]
    ch1: list[float]


def _executable() -> str:
    configured = os.environ.get("TES_CPP_POSI2PULSE_EXECUTABLE")
    if configured:
        return configured
    name = "posi2pulse.exe" if os.name == "nt" else "posi2pulse"
    bundled = files("tes_cpp").joinpath("bin", name)
    if bundled.is_file():
        return str(bundled)
    found = which("posi2pulse")
    if found:
        return found
    raise RuntimeError("posi2pulse executable was not found; set TES_CPP_POSI2PULSE_EXECUTABLE")


def _read(path: Path) -> list[Pulse]:
    document = json.loads(path.read_text(encoding="utf-8"))
    time = document["time"]
    return [
        Pulse(position=int(position), time=time, **waveform)
        for position, waveform in document["pulses"].items()
    ]


def posi2pulse(
    input_path: str | Path,
    positions: list[int],
    *,
    output_path: str | Path | None = None,
    threads: int | None = None,
) -> list[Pulse] | Path:
    """Generate CH0/CH1 pulses for one or more one-based absorber positions.

    Without ``output_path``, return pulse data. With it, save one JSON file and
    return its path, avoiding JSON re-reading in Python.
    """
    if not positions:
        raise ValueError("positions must not be empty")
    if threads is not None and threads < 1:
        raise ValueError("threads must be at least 1")
    target: Path
    temporary = output_path is None
    if temporary:
        with NamedTemporaryFile(suffix=".json", delete=False) as file:
            target = Path(file.name)
    else:
        target = Path(output_path)
    command = [_executable(), str(input_path), str(target), "--positions", ",".join(map(str, positions))]
    if threads is not None:
        command.extend(["--threads", str(threads)])
    try:
        subprocess.run(command, check=True)
        return _read(target) if temporary else target
    finally:
        if temporary:
            target.unlink(missing_ok=True)
