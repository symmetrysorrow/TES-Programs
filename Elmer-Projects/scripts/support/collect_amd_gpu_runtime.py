"""Collect the AMD/HIP runtime gate without changing the active Elmer build."""
from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
ARTIFACTS = ROOT / "artifacts"
HIP_SDK = Path(r"D:\AI\rocm10")


def run(
    command: list[str], *, timeout: int = 60, env: dict[str, str] | None = None
) -> dict[str, Any]:
    try:
        completed = subprocess.run(
            command,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
            env=env,
        )
        return {
            "command": command,
            "returncode": completed.returncode,
            "stdout": completed.stdout[-30000:],
            "stderr": completed.stderr[-10000:],
        }
    except FileNotFoundError as exc:
        return {"command": command, "returncode": None, "error": str(exc)}
    except subprocess.TimeoutExpired as exc:
        return {"command": command, "returncode": None, "timeout": timeout, "error": str(exc)}


def powershell(script: str) -> dict[str, Any]:
    return run(["powershell.exe", "-NoProfile", "-Command", script])


def wsl(script: str, *, timeout: int = 60) -> dict[str, Any]:
    return run(["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", script], timeout=timeout)


def parse_json_stdout(result: dict[str, Any]) -> Any:
    try:
        return json.loads(result.get("stdout", ""))
    except (TypeError, json.JSONDecodeError):
        return None


def main() -> None:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    gpu_query = powershell(
        "Get-CimInstance Win32_VideoController | "
        "Select-Object Name,AdapterRAM,DriverVersion,VideoProcessor,Status | "
        "ConvertTo-Json -Compress"
    )
    hip_info = run([str(HIP_SDK / "Scripts" / "hipInfo.exe")])
    arch = run([str(HIP_SDK / "Scripts" / "offload-arch.exe")])
    sdk_core = HIP_SDK / "Lib" / "site-packages" / "_rocm_sdk_core"
    sdk_libs_root = HIP_SDK / "Lib" / "site-packages" / "_rocm_sdk_libraries" / "bin"
    hip_env = os.environ.copy()
    hip_env.update(
        {
            "HIP_PATH": str(HIP_SDK),
            "ROCM_PATH": str(HIP_SDK),
            "HIP_PLATFORM": "amd",
            "PATH": os.pathsep.join(
                [
                    str(HIP_SDK / "Scripts"),
                    str(sdk_core / "bin"),
                    str(sdk_core / "lib" / "llvm" / "bin"),
                    str(sdk_libs_root),
                    hip_env.get("PATH", ""),
                ]
            ),
        }
    )
    hipconfig = run([str(sdk_core / "bin" / "hipconfig.exe"), "--full"], env=hip_env)
    sdk_libs = sorted(
        str(path)
        for path in HIP_SDK.rglob("*")
        if path.is_file()
        and path.name.lower() in {"hipsparse.dll", "rocsparse.dll", "rocrand.dll"}
    )
    wsl_rocm = wsl(
        "export PATH=/opt/rocm-wsl/bin:/opt/rocm/core-7.14/bin:/usr/bin; "
        "export HIP_PATH=/opt/rocm/core-7.14; "
        "export ROCM_PATH=/opt/rocm/core-7.14; "
        "export LD_LIBRARY_PATH=/opt/rocm-wsl/lib:/opt/rocm/core-7.14/lib; "
        "rocminfo",
        timeout=60,
    )
    wsl_hipconfig = wsl(
        "export PATH=/opt/rocm-wsl/bin:/opt/rocm/core-7.14/bin:/usr/bin; "
        "export HIP_PATH=/opt/rocm/core-7.14; export ROCM_PATH=/opt/rocm/core-7.14; "
        "export LD_LIBRARY_PATH=/opt/rocm-wsl/lib:/opt/rocm/core-7.14/lib; hipconfig --full"
    )
    wsl_sparse = wsl(
        "find /opt/rocm/core-7.14/lib -maxdepth 1 -type f -o "
        "-type l | grep -Ei 'roc(sparse|rand)|hip(sparse|rand)' | sort"
    )
    hip_probe = wsl(
        "hipcc --offload-arch=gfx1201 -x hip "
        "/mnt/d/Github/TES-Programs/Elmer-Projects/scripts/support/hip_device_query.cpp "
        "-o /tmp/tes_hip_device_query && /tmp/tes_hip_device_query"
    )
    sparse_probe = wsl(
        "hipcc --offload-arch=gfx1201 -x hip "
        "/mnt/d/Github/TES-Programs/Elmer-Projects/scripts/support/hip_sparse_query.cpp "
        "-I/opt/rocm/core-7.14/include -L/opt/rocm/core-7.14/lib -lhipsparse "
        "-o /tmp/tes_hip_sparse_query && /tmp/tes_hip_sparse_query"
    )
    runtime = {
        "schema": "tes.amd_gpu_runtime.v1",
        "platform": {
            "os": platform.platform(),
            "python": platform.python_version(),
            "hostname": platform.node(),
        },
        "os_gpu_detection": {
            "raw": gpu_query,
            "parsed": parse_json_stdout(gpu_query),
        },
        "windows_hip_sdk": {
            "root": str(HIP_SDK),
            "exists": HIP_SDK.exists(),
            "hipInfo": hip_info,
            "offload_arch": arch,
            "hipconfig": hipconfig,
            "sparse_runtime_files": sdk_libs,
        },
        "wsl_rocm": {
            "rocminfo": wsl_rocm,
            "hipconfig": wsl_hipconfig,
            "sparse_runtime_files": wsl_sparse,
        },
        "simple_hip_device_query": hip_probe,
        "simple_hipsparse_query": sparse_probe,
        "availability_summary": {
            "os_sees_rx_9070_xt": "AMD Radeon RX 9070 XT" in hip_info.get("stdout", ""),
            "windows_hip_runtime": hip_info.get("returncode") == 0,
            "gfx_architecture": arch.get("stdout", "").strip() or None,
            "windows_sparse_runtime": any("sparse" in path.lower() for path in sdk_libs),
            "wsl_rocm_gpu_agent": "Device Type:             GPU" in wsl_rocm.get("stdout", ""),
        },
    }
    (ARTIFACTS / "amd_gpu_runtime.json").write_text(
        json.dumps(runtime, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    validation = {
        "schema": "tes.hip_device_validation.v1",
        "runtime_artifact": "amd_gpu_runtime.json",
        "device_query_source": "scripts/support/hip_device_query.cpp",
        "sparse_query_source": "scripts/support/hip_sparse_query.cpp",
        "simple_hip_device_query": hip_probe,
        "simple_hipsparse_query": sparse_probe,
        "gpu_execution_proven": (
            hip_probe.get("returncode") == 0
            and "device_memory_probe=PASS" in hip_probe.get("stdout", "")
            and sparse_probe.get("returncode") == 0
        ),
    }
    (ARTIFACTS / "hip_device_validation.json").write_text(
        json.dumps(validation, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(runtime["availability_summary"], ensure_ascii=False))


if __name__ == "__main__":
    main()
