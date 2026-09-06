"""Write the reproducible HYPRE HIP build manifest after a successful install."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT.parent / "tools"
BUILD = TOOLS / "hypre-hip-v3-1-0-build"
INSTALL = TOOLS / "hypre-hip-v3-1-0-install"
SOURCE = TOOLS / "hypre-hip-v3-1-0-src"
OUT = ROOT / "artifacts" / "hypre_hip_build.json"


def command(args: list[str]) -> dict[str, object]:
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True, check=False)
    return {
        "command": args,
        "returncode": result.returncode,
        "stdout": result.stdout[-20000:],
        "stderr": result.stderr[-5000:],
    }


def cache_value(name: str) -> str | None:
    cache = BUILD / "CMakeCache.txt"
    if not cache.exists():
        return None
    prefix = f"{name}:"
    for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[-1]
    return None


def main() -> None:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    source_commit = command(["git", "-C", str(SOURCE), "rev-parse", "HEAD"])
    library = INSTALL / "lib" / "libHYPRE.so"
    manifest = {
        "schema": "tes.hypre_hip_build.v1",
        "status": "PASS" if library.exists() else "FAIL",
        "source": {
            "path": str(SOURCE),
            "commit": source_commit.get("stdout", "").strip() or None,
            "compatibility_patch": "docs/hypre_hip_rocm714_compat.patch",
        },
        "build": {
            "path": str(BUILD),
            "install_prefix": str(INSTALL),
            "cmake_build_type": cache_value("CMAKE_BUILD_TYPE"),
            "shared": cache_value("BUILD_SHARED_LIBS"),
            "mpi": cache_value("HYPRE_ENABLE_MPI"),
            "hip": cache_value("HYPRE_ENABLE_HIP"),
            "hip_architecture": cache_value("CMAKE_HIP_ARCHITECTURES"),
            "hip_flags": cache_value("CMAKE_HIP_FLAGS"),
            "rocsparse": cache_value("HYPRE_ENABLE_ROCSPARSE"),
            "rocrand": cache_value("HYPRE_ENABLE_ROCRAND"),
            "rocblas": cache_value("HYPRE_ENABLE_ROCBLAS"),
            "rocsolver": cache_value("HYPRE_ENABLE_ROCSOLVER"),
            "rocthrust": cache_value("HYPRE_ENABLE_ROCTHRUST"),
            "unified_memory": cache_value("HYPRE_ENABLE_UNIFIED_MEMORY"),
            "device_malloc_async": cache_value("HYPRE_ENABLE_DEVICE_MALLOC_ASYNC"),
            "gpu_aware_mpi": cache_value("HYPRE_ENABLE_GPU_AWARE_MPI"),
            "umpire": cache_value("HYPRE_ENABLE_UMPIRE"),
        },
        "installed_files": [
            str(path)
            for path in sorted(INSTALL.glob("lib/*"))
            if path.is_file() or path.is_symlink()
        ],
        "runtime_link_check": command(
            [
                "wsl.exe",
                "-d",
                "Ubuntu",
                "--",
                "bash",
                "-lc",
                "ldd /mnt/d/Github/TES-Programs/tools/hypre-hip-v3-1-0-install/lib/libHYPRE.so",
            ]
        )
        if library.exists()
        else None,
    }
    OUT.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(OUT)


if __name__ == "__main__":
    main()
