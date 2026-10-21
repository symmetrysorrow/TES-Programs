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
    input_energy: float | None = None,
    save_all: bool = False,
    full_energy_only: bool = False,
    normalize_deposits: bool = False,
    reaction_audit: str | Path | None = None,
    reaction_audit_all: bool = False,
    physics: str | None = None,
    physics_config: str | Path | None = None,
    stopping_table: str | Path | None = None,
) -> Dump2EventResult:
    """Convert ``dumpall.dat`` into event JSON or ``.h5`` data.

    HDF5 conversion is streaming, so Python never loads the full native JSON
    document into memory. ``physics="pb-egs5"`` selects an experimental,
    fail-closed PHITS 3.341 / Pb / unit-weight HDF5 v2 path; it does not change
    the default conversion behavior.
    """
    output = Path(output_path)
    command = [_executable(), str(dump_path), str(output)]
    if physics is not None:
        if physics != "pb-egs5":
            raise ValueError("physics must be None or the experimental 'pb-egs5' mode")
        if physics_config is None or stopping_table is None:
            raise ValueError("pb-egs5 requires physics_config and stopping_table")
        config = Path(physics_config)
        table = Path(stopping_table)
        if not config.is_file():
            raise FileNotFoundError(f"physics_config does not exist: {config}")
        if not table.is_file():
            raise FileNotFoundError(f"stopping_table does not exist: {table}")
        if output.suffix.lower() not in {".h5", ".hdf5"}:
            raise ValueError("pb-egs5 mode writes HDF5 v2; output must end in .h5 or .hdf5")
        if save_all or full_energy_only or normalize_deposits or reaction_audit or reaction_audit_all:
            raise ValueError("legacy event-selection/audit options are not supported in pb-egs5 mode")
        command += ["--physics", physics, "--hdf5-version", "2",
                    "--physics-config", str(config), "--stopping-table", str(table)]
    else:
        if physics_config is not None or stopping_table is not None:
            raise ValueError("physics_config/stopping_table require explicit physics='pb-egs5'")
        if input_energy is None:
            raise ValueError("input_energy is required for the default converter")
        command += ["--input-energy", str(input_energy)]
    index_path = output.with_name("FullEnergyList.dat")
    if save_all:
        command += ["--save-all"]
    if save_all or full_energy_only:
        command += ["--full-energy-list", str(index_path)]
    if full_energy_only:
        command += ["--full-energy-only"]
    if normalize_deposits:
        if physics is not None:
            raise ValueError("--normalize-deposits is incompatible with pb-egs5 physical reconstruction")
        command += ["--normalize-deposits"]
    if reaction_audit_all and reaction_audit is None:
        raise ValueError("reaction_audit_all requires reaction_audit")
    if reaction_audit is not None:
        command += ["--reaction-audit", str(reaction_audit)]
    if reaction_audit_all:
        command += ["--reaction-audit-all"]
    subprocess.run(command, check=True)
    event_ids = ()
    if (save_all or full_energy_only) and index_path.exists():
        event_ids = tuple(int(line) for line in index_path.read_text().splitlines() if line)
    return Dump2EventResult(output_path=output, full_energy_event_ids=event_ids)
