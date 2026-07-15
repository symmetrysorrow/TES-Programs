from __future__ import annotations

import json
import os
import subprocess
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path
from shutil import which
from tempfile import NamedTemporaryFile

import h5py
import numpy as np


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
    # When running from a source checkout, also look in common CMake build
    # directories (single- and multi-configuration generators).
    repo_root = Path(__file__).resolve().parents[3]
    for build_root in (repo_root / "tes_cpp" / "build", repo_root / "build"):
        if build_root.is_dir():
            matches = sorted(build_root.rglob(name))
            if matches:
                return str(matches[0])
    found = which("posi2pulse")
    if found:
        return found
    raise RuntimeError("posi2pulse executable was not found; set TES_CPP_POSI2PULSE_EXECUTABLE")


def _read(path: Path) -> list[Pulse]:
    if path.suffix.lower() in {".h5", ".hdf5"}:
        with h5py.File(path, "r") as file:
            time = file["time"][:].tolist()
            return [
                Pulse(int(position), time, file["ch0"][index].tolist(), file["ch1"][index].tolist())
                for index, position in enumerate(file["event_id"].asstr()[:])
            ]
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

    Without ``output_path``, return pulse data. With it, save one HDF5 file
    (``.h5``/``.hdf5``) or legacy JSON file and return its path.
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
    native_target = target
    hdf5_target = target.suffix.lower() in {".h5", ".hdf5"}
    if hdf5_target:
        with NamedTemporaryFile(suffix=".json", delete=False) as file:
            native_target = Path(file.name)
    command = [_executable(), str(input_path), str(native_target), "--positions", ",".join(map(str, positions))]
    if threads is not None:
        command.extend(["--threads", str(threads)])
    try:
        subprocess.run(command, check=True)
        if hdf5_target:
            document = json.loads(native_target.read_text(encoding="utf-8"))
            with h5py.File(target, "w") as file:
                file.attrs["format"] = "tes-pulses"
                file.attrs["format_version"] = 1
                file.attrs["input_json"] = json.dumps(document["input"], separators=(",", ":"))
                file.create_dataset("time", data=np.asarray(document["time"], dtype=np.float64))
                pulse_items = list(document["pulses"].items())
                file.create_dataset("event_id", data=np.asarray([key for key, _ in pulse_items], dtype=h5py.string_dtype("utf-8")))
                file.create_dataset("ch0", data=np.asarray([value["ch0"] for _, value in pulse_items]), compression="gzip")
                file.create_dataset("ch1", data=np.asarray([value["ch1"] for _, value in pulse_items]), compression="gzip")
        return _read(target) if temporary else target
    finally:
        if native_target != target:
            native_target.unlink(missing_ok=True)
        if temporary:
            target.unlink(missing_ok=True)
