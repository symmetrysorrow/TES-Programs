# Research branch validation and publication record

## Pre-push checks

- Branch: `feat/dump2event-physical-validation-20261020`, based on `f7f505be30fa7609297db2d116856d8123f9c127`; remote branch was at the same SHA before this change.
- CMake/MSVC Release build succeeded for `dump2event`, `dump2event_physical_reconstruct`, and `test_physical_deposit_engine`.
- CTest: 1/1 passed.
- Python tests: `python -m pytest tes_cpp -q` -> 34 passed, 10 skipped.
- Default converter regression on the 10k sample with `--summary-only --max-histories 1000` processed exactly 1,000 histories and exited successfully.
- Research HDF5 rebuild: 1,954,742/1,954,742 additions, unresolved 0, negatives 0; same-run PHITS internal comparison reports no energy, region or checked positive-position mismatches.
- HDF5 history histogram: 150 bins, count L1=0 vs PHITS 3.341; PHITS printed-value residual L1 `1.554312234475219e-14`.
- `git diff --cached --check` passed. Staged path scan found no PHITS source/build products, Dumpall, HDF5 output, internal trace, local coefficient table, local absolute paths, or credential patterns. The two pre-existing untracked research directories remain unstaged. Private run data is excluded by this directory's `.gitignore`.
- No default `dump2event` energy formula change; PR #6 remains untouched and Draft.

## Push and GitHub Actions

Push was completed to the authorized independent research branch. The implementation commit is `f1f87f24a4d0fe8224f33621c4968f1355aff0f8` (`feat: validate PHITS positron deposit reconstruction`). The push updated `origin/feat/dump2event-physical-validation-20261020` from `f7f505be30fa7609297db2d116856d8123f9c127`.

GitHub Actions run [39539451114](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39539451114), workflow `dump2event energy regression tests`, was triggered by `push` for that exact commit. It completed successfully (`conclusion=success`, run attempt 1; created 2026-10-21 03:35:55 UTC, updated 03:36:26 UTC). No workflow failure remains to investigate.

At final review, the research worktree is on `feat/dump2event-physical-validation-20261020` at `f1f87f24a4d0fe8224f33621c4968f1355aff0f8`. The only untracked paths are two pre-existing research directories, `tes_cpp/dump2event/validation/extended_physics_20261020/` and `tes_cpp/dump2event/validation/weighted_cutoff_positron_photon/`; neither was staged or modified for this push.

Research status remains experimental and scoped. This push is not a claim of general PHITS physics completeness or readiness to integrate into PR #6.
