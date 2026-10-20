# Push and CI report

## Push

The research-only implementation was pushed to `feat/dump2event-physical-validation-20261020`; PR #6 and its Draft state were not changed.

| Commit | Change |
|---|---|
| `69780029` | Make HDF5 v2 event appends transactional |
| `f56a9271` | Run HDF5 v2 and dump2event regressions in CI |
| `f8416c83` | Add diagnostic Dumpall-to-HDF5 v2 streaming path |
| `02cd2fabbcd9dd00ea867ecb6e162d73948d345b` | Record real HDF5 v2 physics validation |

At verification, local `HEAD` and `origin/feat/dump2event-physical-validation-20261020` both resolved to `02cd2fabbcd9dd00ea867ecb6e162d73948d345b`. The working tree had no tracked modifications. Two preexisting untracked research directories were left untouched and unstaged: `tes_cpp/dump2event/validation/extended_physics_20261020/` and `tes_cpp/dump2event/validation/weighted_cutoff_positron_photon/`.

PHITS distribution source, private stopping-power tables, raw Dumpalls, generated HDF5 data, and internal per-addition traces were not added to Git.

## Local verification

- MSVC Release CMake build: passed; `dumpall_transport_export` included.
- CTest: 1/1 passed.
- `python -m pytest tes_cpp/tests tes_cpp/dump2event/tests -q -rs`: 41 passed, 0 skipped, with `PYTHONPATH=tes_cpp/python` and `TES_CPP_DUMP2EVENT_EXECUTABLE` set to the Release converter.
- Python bytecode compilation for the new scripts: passed.
- `git diff --check`: passed before the commits.
- Real-data HDF5 readback checks: passed structurally. The 100-history HDF5 ledger matched all 38,077 PHITS internal additions (maximum energy difference `5.7976626971489864e-15 MeV`; addition/event/particle/region mismatches 0; checked positive scoring-position mismatches 0). The 1,000-history histogram matched all 150 bins. The 10,000-history output remains partial: 342 unresolved ledger rows and 13 histogram-bin differences (L1=18).

## GitHub Actions

- Workflow: `dump2event energy regression tests`
- Run: [39453919608](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39453919608), run 33, event `push`
- Head SHA: `02cd2fabbcd9dd00ea867ecb6e162d73948d345b`
- Result: completed, conclusion `success` (2026-10-20 16:03:17 UTC).
- Job `energy-regressions`: success. Checkout, Python setup/dependencies, HDF5 installation, CMake configure/build, C++ regression tests, and v1/v2 Python regression tests all completed successfully.

The successful CI run verifies the pushed code commit. This report-only follow-up may trigger another workflow run when pushed; it does not change the tested implementation. No PR was created, updated, or merged. The physical reconstruction path remains diagnostic and is not enabled by the default converter.
