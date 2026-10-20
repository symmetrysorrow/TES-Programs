# Commit and CI report

- Branch: `feat/dump2event-physical-validation-20261020`.
- Baseline and remote head before publication: `589f84055f7f1c457487a223d7e0a9ea88196fb0`.
- Parser implementation commit: `abd54f3bbaa389269ce6ac5ec3a792a5e2635071` (`Harden structured Dumpall parsing for truncated records`).
- Evidence and report commit: `142e00f130ef766d232eaec70baef4f394ddc326` (`Document GNU Dumpall truncation and mapping limits`).
- Push command: `git push origin feat/dump2event-physical-validation-20261020`; completed successfully on 2026-10-20 JST.
- GitHub Actions run for `142e00f130ef766d232eaec70baef4f394ddc326`: [dump2event energy regression tests](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39417687378), `energy-regressions` completed with `success`.
- Local checks: MSVC 19.51 targets `dump2event`, `dumpall_record_inspect`, and `test_physical_deposit_engine` built; CTest 1/1 passed; Python 23/23 passed; `git diff --check` clean.
- Scope audit: the pushed diff contains parser/interface/call-site/test files and this validation report folder only. No PHITS distribution source, Dumpall data, stopping-power table, diagnostic executable or install file was included. No energy-deposit equation or physics mode was added. PR #6 and its Draft state were not changed.

The branch still has preexisting untracked validation directories in the worktree (`extended_physics_20261020/` and `weighted_cutoff_positron_photon/`); they were not staged or pushed by this task.
