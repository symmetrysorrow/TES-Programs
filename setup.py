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
        names = (["posi2pulse.exe", "dump2event.exe", "dump2event_physical_reconstruct.exe",
                  "dumpall_transport_export.exe"] if os.name == "nt" else
                 ["posi2pulse", "dump2event", "dump2event_physical_reconstruct",
                  "dumpall_transport_export"])
        # Search only the native CMake output directory.  This also avoids
        # accidentally finding the packaged copy on a subsequent build.
        native_bin = ROOT / "tes_cpp" / "build"
        for name in names:
            candidates = sorted(native_bin.rglob(name))
            if not candidates:
                raise FileNotFoundError(f"CMake did not produce {name}")
            shutil.copy2(candidates[0], package_bin / name)
        # The opt-in Pb physics converter launches this Python driver from
        # alongside dump2event; keep the script and both helper binaries in
        # the installed tes_cpp/bin directory, not only the source worktree.
        driver = ROOT / "tes_cpp" / "dump2event" / "validation" / "build_physical_hdf5_v2.py"
        if not driver.is_file():
            raise FileNotFoundError(f"Pb HDF5 v2 driver was not found: {driver}")
        shutil.copy2(driver, package_bin / driver.name)
        # Optional source-controlled HATCH coefficients. The exact CSV must
        # be added only after confirming that public redistribution is allowed.
        # Missing table does not block development builds, but conversions
        # without an explicit --stopping-table will fail closed.
        bundled = ROOT / "tes_cpp" / "dump2event" / "data" / "pb_egs5_phits3341.csv"
        if bundled.is_file():
            import hashlib
            expected = "20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9"
            actual = hashlib.sha256(bundled.read_bytes()).hexdigest()
            if actual != expected:
                raise ValueError(f"Bundled Pb HATCH CSV SHA-256 mismatch: {actual}")
            shutil.copy2(bundled, package_bin / bundled.name)
        else:
            print("WARNING: verified Pb HATCH data not bundled; use explicit --stopping-table")
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
