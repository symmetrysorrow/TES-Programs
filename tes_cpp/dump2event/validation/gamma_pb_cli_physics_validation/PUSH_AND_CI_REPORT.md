# Research branch push and CI report

Push authorization is from the task request. Target branch: `feat/dump2event-physical-validation-20261020`; PR #6 was left untouched. Before commit, the worktree was based on `96873b511ce4a68aa049fcfa737867b3221488e4`, equal to the observed `origin` tip.

Commit `67d043708932443c07bd48d10e2752f491d994aa` (`Add gamma validation and HDF5 v2 provenance checks`) was pushed to the research branch. GitHub confirmed that the branch points to this SHA. It contains 30 files. The PHITS source/install, executable, private Pb stopping table, raw Dumpall and HDF5 datasets, original inputs, and two protected untracked research directories were excluded from staging.

Checks before push: MSVC Release build succeeded; CTest passed 1/1; `python -m pytest tes_cpp -q` passed 51 tests and skipped 10 existing environment-gated tests. `git diff --cached --check` was clean. The source scan found no hard-coded local absolute paths or secret-like values in the changed implementation files.

GitHub Actions run [#43](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39590935167), for commit `67d043708932443c07bd48d10e2752f491d994aa`, completed with conclusion `success` on 2026-10-21. The GitHub API confirmed PR #6 remains open and Draft, with head branch `fix/dump2event-energy-audit-20261016` at `36dec9f139918d233faf0876736f43d1bbdb965c`; the PR head was not changed.

The PHITS/CLI scope remains experimental. The 2 MeV/1,000 run leaves one below-threshold candidate unresolved and therefore has no published HDF5 file. The independent-seed Co-60 histogram difference and the historical 10,000-history malformed NCOL=4 CLI failure remain limitations. The run #43 success is a code CI result and does not remove those physics limitations.

This report update is a documentation-only follow-up commit; Actions may run again for that commit. The successful run above is for the code commit identified explicitly above.
