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

Push is explicitly authorized for this independent research branch after the above gates. The first implementation commit SHA, branch push result, and Actions run URL/status will be recorded after the remote operation. Any Action failure will be investigated before this report is finalized.

Research status remains experimental and scoped. This push is not a claim of general PHITS physics completeness or readiness to integrate into PR #6.
