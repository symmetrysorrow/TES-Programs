"""Build and package the TES native tools with CMake.

Usage::

    python -m pip install .

The custom build step invokes CMake for ``tes_cpp`` and bundles the resulting
``posi2pulse`` and ``dump2event`` executables in the ``tes_cpp`` Python
package.  The same setup script works with single- and multi-configuration
generators on Windows, Linux, and macOS.
"""

from pathlib import Path
import os
import platform
import shutil
import subprocess
import sys

from setuptools import find_packages, setup
from setuptools.command.build_ext import build_ext


ROOT = Path(__file__).resolve().parent

# ``python setup.py`` with no command means: build the native tes_cpp tools.
if len(sys.argv) == 1:
    sys.argv.append("build_ext")


def _print_build_help(missing):
    system = platform.system()
    print(f"\nTES native build prerequisites missing: {', '.join(missing)}")
    if system == "Windows":
        print("Recommended: install Visual Studio 2022 with the 'Desktop development with C++' workload and CMake.")
        print("Then run this from 'x64 Native Tools Command Prompt for VS 2022', or ensure cl.exe and cmake.exe are on PATH.")
        print("CMake installer: https://cmake.org/download/")
    elif system == "Darwin":
        print("Recommended: install Xcode Command Line Tools and CMake with Homebrew:")
        print("  xcode-select --install")
        print("  brew install cmake")
    else:
        print("Recommended (Debian/Ubuntu): sudo apt install build-essential cmake")
        print("Recommended (Fedora): sudo dnf install gcc-c++ cmake")
        print("Recommended (Arch): sudo pacman -S base-devel cmake")
    print()


def _raise_build_error(stage, exc):
    """Report a native build error without misdiagnosing it as missing VS."""
    print(f"\nTES native build failed during {stage}.")
    print("See the CMake/MSBuild error above for the actual cause.")
    raise RuntimeError(f"CMake {stage} failed") from exc


class CMakeBuild(build_ext):
    def run(self):
        missing = []
        if shutil.which("cmake") is None:
            missing.append("CMake")
        if missing:
            _print_build_help(missing)
            return
        try:
            subprocess.run(["cmake", "--version"], check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            _print_build_help(["CMake"])
            return

        # Keep all CMake-generated files inside the native project directory.
        # This avoids creating repository-root ``build``/``CMakeFiles`` dirs.
        build_dir = ROOT / "tes_cpp" / "build"
        configure = ["cmake", "-S", str(ROOT / "tes_cpp"), "-B", str(build_dir)]
        # CMake owns compiler discovery for this build.  ``build_ext`` does
        # not initialise ``self.compiler`` until ``super().run()``; querying
        # it here therefore fails when setup.py is invoked directly.
        configure += ["-DCMAKE_BUILD_TYPE=Release"]
        try:
            subprocess.run(configure, check=True)
        except subprocess.CalledProcessError as exc:
            _raise_build_error("configuration", exc)

        build = ["cmake", "--build", str(build_dir), "--config", "Release"]
        jobs = os.environ.get("CMAKE_BUILD_PARALLEL_LEVEL")
        if jobs:
            build += ["--parallel", jobs]
        try:
            subprocess.run(build, check=True)
        except subprocess.CalledProcessError as exc:
            _raise_build_error("compilation", exc)

        super().run()
        package_bin = Path(self.build_lib) / "tes_cpp" / "bin"
        package_bin.mkdir(parents=True, exist_ok=True)
        names = ["posi2pulse.exe", "dump2event.exe"] if os.name == "nt" else ["posi2pulse", "dump2event"]
        # Search only the native CMake output directory.  This also avoids
        # accidentally finding the packaged copy on a subsequent build.
        native_bin = ROOT / "tes_cpp" / "build"
        for name in names:
            candidates = sorted(native_bin.rglob(name))
            if not candidates:
                raise FileNotFoundError(f"CMake did not produce {name}")
            shutil.copy2(candidates[0], package_bin / name)
        print("TES native tools built successfully.")

    def build_extensions(self):
        super().build_extensions()


setup(
    name="tes-programs",
    version="0.1.0",
    description="TES pulse simulation tools",
    package_dir={"": "tes_cpp/python"},
    packages=find_packages("tes_cpp/python"),
    package_data={"tes_cpp": ["bin/*"]},
    include_package_data=True,
    install_requires=["h5py>=3.8", "numpy>=1.23"],
    # Keep setuptools' package build output alongside the CMake output.
    options={"build": {"build_base": str(ROOT / "tes_cpp" / "build" / "python")}},
    cmdclass={"build_ext": CMakeBuild},
)
