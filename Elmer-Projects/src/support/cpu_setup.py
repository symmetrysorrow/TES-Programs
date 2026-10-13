"""Set up the CPU-only TES toolchain from a Windows host with WSL2 Ubuntu 24.04."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def wsl_mount_path(path: Path) -> str:
    path = path.resolve()
    if not path.drive:
        raise ValueError(f"Expected a Windows-drive checkout for WSL: {path}")
    return f"/mnt/{path.drive[0].lower()}/{path.relative_to(path.anchor).as_posix()}"


def setup_cpu() -> None:
    if sys.version_info < (3, 11):
        raise RuntimeError("Python 3.11 or later is required for TOML configuration.")
    if os.name != "nt":
        raise RuntimeError(
            "Automatic host setup currently targets Windows + WSL2 Ubuntu 24.04. "
            "On Ubuntu, start with 'sudo bash src/support/setup_cpu_ubuntu.sh'; "
            "native-Linux execution is not integrated with the CLI yet."
        )

    distro = os.environ.get("TES_WSL_DISTRO", "ElmerGPU")
    try:
        subprocess.run(
            ["wsl.exe", "-d", distro, "-u", "root", "--exec", "true"], check=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        raise RuntimeError(
            f"WSL distribution {distro!r} is unavailable. Set up an Ubuntu 24.04 "
            "WSL2 distribution first; select it using TES_WSL_DISTRO. "
            "The default distribution name is 'ElmerGPU'."
        ) from exc

    venv = ROOT / ".venv"
    if not venv.exists():
        subprocess.run([sys.executable, "-m", "venv", str(venv)], check=True)
    python = venv / "Scripts" / "python.exe"
    if not python.is_file():
        raise RuntimeError(f"Virtual environment is incomplete: {python}")
    subprocess.run(
        [str(python), "-m", "pip", "install", "-r", str(ROOT / "requirements-cpu.txt")],
        check=True,
    )
    script = wsl_mount_path(ROOT / "src" / "support" / "setup_cpu_ubuntu.sh")
    subprocess.run(
        ["wsl.exe", "-d", distro, "-u", "root", "--exec", "bash", script],
        check=True,
    )
    print("CPU setup completed. Keep this terminal open for the next command:")
    print(f'  "{python}" main.py run single_pixel_generated --dry-run')
