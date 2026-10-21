# 研究ブランチpushとCI

## Push

- Branch: `feat/dump2event-physical-validation-20261020`
- Pushed HEAD at report update: `60cde10f2b35b1b00c50803deb657ef2392cfdf6`
- Included commits:
  - `83b7ba56` — PHITS 3.341 dead-secondary score gate
  - `fd843124` — 2 MeV HDF5 addition join
  - `98cd6c21` — interrupted 10k run status and fresh 1k HDF5 readback record
  - `60cde10f` — push and CI report
- The working tree still contains uncommitted 10k helper/report files and protected research folders. They were not staged or pushed.

## Local verification

- MSVC Release CMake build: exit 0.
- CTest: 1/1 passed.
- `python -m pytest tes_cpp\dump2event\tests -q`: 34 passed, 10 skipped.
- 2 MeV/1,000-history HDF5 v2 was produced earlier by the normal CLI (1,790,020,441 bytes; recorded CLI elapsed 2,349.8716 s). To avoid another ~39-minute conversion and duplicate 1.79 GB output, this task re-read that existing artifact and freshly re-ran the full 5,978,912-row same-run join and HDF5 completeness audit.
- Readback: 5,978,912/5,978,912 IDs matched; identity, region, and position mismatches 0; max energy error `4.680023729664029e-14 MeV`; no unresolved, invalid, or negative deposits. Machine-readable results are `TWO_MEV_1000_REREAD_JOIN_20261022.json` and `TWO_MEV_1000_REREAD_HDF5_AUDIT_20261022.json`.
- NOCAS 398 child index 0, Q7 `7.5933446597176157e-7 MeV`, is below runtime-measured `eletb(1)=9.9999999999999995e-7 MeV`; the C++ research ledger correctly emits one `NOT_SCORED_BY_PHITS` audit candidate and no child deposit row.

## GitHub and PR #6

- GitHub Actions run [39638078534](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39638078534) for commit `98cd6c21` completed successfully; job `energy-regressions` succeeded.
- GitHub Actions run [39638346141](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39638346141) for the following report commit `60cde10f` also completed successfully; job `energy-regressions` succeeded.
- PR #6 was checked after push: open, Draft, not merged, head remains `36dec9f139918d233faf0876736f43d1bbdb965c`. No change was made to that PR.
- The independent 10,000-history PHITS run was stopped at the user's request and exited `-1`; no successful 10k replay/conversion/join is claimed. Its incomplete output is not a validation run.

The pushed code remains a bounded research implementation. The ordinary converter's default route is unchanged, and the physical reconstruction remains explicitly opt-in. Unverified weights, physical processes, and spatial cases remain outside the validated claim.
