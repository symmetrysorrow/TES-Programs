# Research branch push and CI report

Push authorization is from the task request. Target branch: `feat/dump2event-physical-validation-20261020`; PR #6 must remain Draft and untouched. Before commit, the worktree is based on `96873b511ce4a68aa049fcfa737867b3221488e4`, equal to the previously observed `origin` tip.

Preliminary checks: MSVC Release build succeeded; CTest passed 1/1; `python -m pytest tes_cpp -q` passed 51 tests and skipped 10 existing environment-gated tests. `git diff --check` is clean. The mono-gamma raw datasets, PHITS source/install, PHITS executable, private Pb stopping table, original inputs, and two protected untracked research directories are excluded from staging.

The PHITS/CLI scope remains experimental. The 2 MeV/1,000 run leaves one below-threshold candidate unresolved and therefore has no published HDF5 file. The independent-seed Co-60 histogram difference and the historical 10,000-history malformed NCOL=4 CLI failure remain limitations.

Commit: pending final review. Push: pending. GitHub Actions: pending.
