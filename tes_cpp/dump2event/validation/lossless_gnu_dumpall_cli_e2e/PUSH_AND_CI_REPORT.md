# Push and CI report

## This update

- Target branch: `feat/dump2event-physical-validation-20261020`
- Base at start: `df2e07e7c3517a7b2e836a02412c08a681c85c1f`
- This update records same-run PHITS 3.341 results for the ordinary opt-in `dump2event --physics pb-egs5 --hdf5-version 2` path on 100- and 1,000-history unit-weight Pb/electron runs, and adds reusable raw-to-internal comparison scripts.
- It does not include a change to the default converter, PHITS source/binaries, raw Dumpall, internal trace, private stopping table, or HDF5 outputs.
- The OpenMP merge build attempt failed in GNU assembler (`.tls_common`); threaded serialization remains unvalidated.

## Pre-push checks

- MSVC Release build: passed.
- C++ CTest: 1/1 passed.
- Python tests: 44 passed, 10 skipped.
- `git diff --check`: passed.
- Comparison scripts compile with `python -m py_compile`.
- Only explicitly staged report and script paths are included; pre-existing untracked research directories are excluded.

## Push and CI

- Commit: `43b81439dbce6dab2936ac16a2ca3b105819abd9` (`Document PHITS physical deposit CLI validation`).
- Push to `origin/feat/dump2event-physical-validation-20261020`: succeeded.
- GitHub Actions run [39569564105](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39569564105): **success**; the build, C++ regressions, and Python/v1/v2 regression steps all completed successfully.
- The earlier report-only run `39550926664` also succeeded, but is separate from this commit's result.

PR #6 remains separate and untouched; no merge or PR update is part of this work.
