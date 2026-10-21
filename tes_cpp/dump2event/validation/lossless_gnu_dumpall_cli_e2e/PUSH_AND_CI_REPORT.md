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

To be filled with the commit SHA, push result, and GitHub Actions run for this update after push. The prior report-only commit `5f0c5049c439cdc7580c4442ac2b00525ce8d19f` and its Actions success are historical and do not validate this update.

PR #6 remains separate and untouched; no merge or PR update is part of this work.
