# Build and regression tests

Executed on 2026-10-21 in the research worktree:

- MSVC Release build: `cmake --build build/physdeposit-vs --config Release -j 4` — exit 0.
- CTest: `ctest --test-dir build/physdeposit-vs -C Release --output-on-failure` — 1/1 passed.
- Python: `python -m pytest tes_cpp -q` — 44 passed, 10 skipped.
- The normal CLI was then exercised on the 100- and 1,000-history real PHITS runs; both returned 0 and produced HDF5 v2.
- `scripts/compare_hdf5_to_phits_trace.py` and `scripts/compare_tdeposit_histogram.py` compile with `py_compile` and were executed against both datasets.

No PR #6 files, legacy CLI defaults, PHITS installation files, source inputs, or original Dumpall files were modified. The normal production converter path was not changed in this evidence-only update.
