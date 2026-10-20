# HDF5 v2 CI integration

The workflow now installs Python 3.12, `pytest`, `numpy`, and `h5py`; builds the C++ converter and physical-engine test; runs CTest; and runs both `tes_cpp/tests` and `tes_cpp/dump2event/tests` with `PYTHONPATH=tes_cpp/python`. `-rs` prints skip reasons. Input-data physics validation is not performed in CI because PHITS data and private Dumpalls are not supplied to GitHub Actions.

Local development environment: Python 3.14.3, NumPy 2.5.3, h5py 3.16.0, pytest 9.1.1, MSVC toolchain via CMake/MSBuild 18.11.0-1.26427.13+9353559c2. CMake configured and built the `Release` configuration; CTest passed 1/1.

Initial local pytest run before setting the converter executable: 26 passed, 10 skipped. The skips were CLI integration checks requiring `TES_CPP_DUMP2EVENT_EXECUTABLE`. The final local run supplied that variable and passed 41/41 with no skips. CTest passed 1/1.

The pushed commit `02cd2fabbcd9dd00ea867ecb6e162d73948d345b` completed GitHub Actions run [39453919608](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39453919608) successfully. Its `energy-regressions` job passed CMake configure/build, C++ regression tests, and the v1/v2 Python regression tests. Input-data PHITS validation is not part of CI because PHITS data and private Dumpalls are not supplied to GitHub Actions.
