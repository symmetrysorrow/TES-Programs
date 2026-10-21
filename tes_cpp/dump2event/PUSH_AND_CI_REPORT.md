# Research branch push and CI report

## Pre-push review

- Target branch: `feat/dump2event-physical-validation-20261020`.
- The branch started at remote commit `5ec9264ca090324a6562b358197b4930d1b75c07`.
- Existing implementation work was preserved in this isolated worktree; a
  local backup ref was also created before edits.
- The default `dump2event` mode and existing signed `E_deposit` formula are
  unchanged. The experimental algorithm is a separate executable.
- No PHITS source, PHITS executable, Dumpall, generated build output, or
  stopping-power table is included in the intended commit. The table SHA is
  documented, but the nonpublic table itself stays local.
- No local absolute paths are embedded in source. Runtime inputs are explicit
  CLI paths and manifests hash the supplied files/configuration.
- The staged allowlist contains 24 source, documentation, script and small
  validation-report files. A staged-name scan found no Dumpall, executable,
  DLL/library, PHITS Fortran source, PHITS input/output or stopping table.
  `git diff --cached --check` passed; a credential-pattern scan found no
  credential-like strings. The three small comparison CSVs contain only the
  eight boundary rows, 1,000 per-history totals, and 150 histogram bins.
- MSVC 19.51 build passed for all requested targets; CTest passed 1/1; Python
  tests passed 23/23; default converter regression processed 1,000 histories.
- PHITS same-seed internal comparison matched 201,971/201,971 additions and
  all 150 histogram bins after the NCOL=10 boundary branch was added.
- PHITS source/data exclusions, secrets scan and final staged-file allowlist
  are to be recorded after staging.

## Push and Actions

- Commit: `7ecb2e830ca09dc26e3d96519fa6edbe84478abb` (`Add experimental
  PHITS deposit reconstruction ledger`).
- Pushed to `origin/feat/dump2event-physical-validation-20261020`; remote ref
  was read back and matched the commit SHA.
- GitHub Actions run [39393118635](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39393118635)
  completed successfully on `ubuntu-latest` at 2026-10-20 08:25:19 UTC.
  HDF5 setup, CMake configure, build, and energy accounting tests all passed.
- PR #6 was queried after the push and is still open, `draft=true`, and based
  on its original head branch `fix/dump2event-energy-audit-20261016`. No PR #6
  operation or merge was performed.

The research branch is pushed, but PR #6 integration is blocked by the
unsupported physics listed in `UNSUPPORTED_PHYSICS.md`; this is not a complete
physical deposit converter.

## Independent-seed follow-up review (2026-10-20)

### Starting state and preservation

- Isolated worktree: `D:\github\TES-Programs-physdeposit-20261020`.
- Branch: `feat/dump2event-physical-validation-20261020`.
- Starting local and remote HEAD: `540fdd121fee7180d4c74bde460be1a8c5e9be57`.
- No tracked source changes existed at start. Two pre-existing untracked
  physics research directories were left intact and are excluded from staging.
- The experimental C++ implementation was already present on the remote from
  commit `7ecb2e830ca09dc26e3d96519fa6edbe84478abb`; this follow-up changes
  documentation and adds a reproducible bin-edge audit script only.

### Root-cause and implementation review

The exact PHITS 3.341 `tdepstregEVENT` source path and the three serialized-bin
crossings are documented in `INDEPENDENT_SEED_ROOT_CAUSE.md`. The prior L1=6
was reproduced as three histories moving one bin each after eight omitted
NCOL=10 EGS5 boundary contributions totaling `0.018778564270171039 MeV` were
added. After the correction the internal per-history comparison has maximum
absolute difference `7.815970093361102e-14 MeV`, and both PHITS and corrected
Dumpall histograms have L1=0. The user-quoted `2.22e-15 MeV` is not the
PHITS-internal error; it belongs to raw-only code-to-code agreement.

The public implementation audit confirmed explicit PHITS 3.341 / Pb region
1000 / material 2 / EGS5 medium 2 / 1 keV / unit-weight gating, caller-provided
Dumpall and stopping table, wrapper-manifest SHA-256 values and table-source
provenance, wrong-version rejection, unresolved unsupported rows, and the
`EXPERIMENTAL_PARTIAL` status. The new executable is separate from the default
CLI. No negative clipping or normalization was added to the experiment. The
stopping-table documentation was corrected to the actual 11-column parser
header, including the two checked-but-unused positron columns.

### Fresh verification

- Microsoft C++ compiler: MSVC `19.51.36260.0` (Visual Studio 2026 Insider).
- Configure/build: CMake `-S tes_cpp/dump2event -B %TEMP%\\tes_physdeposit_review_build -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release`, then `cmake --build ... --parallel`. Build exit code 0; NMake ignored parallelism and built serially. HDF5 1.14.6 was fetched because no system HDF5 package was found; only third-party C4819 code-page warnings were emitted.
- CTest: `ctest --test-dir %TEMP%\\tes_physdeposit_review_build --output-on-failure` — 1/1 passed.
- Python: with `TES_CPP_DUMP2EVENT_EXECUTABLE` set to the just-built binary, `python -m unittest discover -s tes_cpp/dump2event/tests -v` — 23/23 passed. An earlier pytest run without the executable had 13 passed and 10 skipped; the executable-backed run supersedes that result.
- Default-converter regression: built `dump2event.exe` processed the existing 1,000-history Dumpall with `--summary-only --max-histories 1000 --save-all --history-summary ...`; exit code 0, `Processed 1000 histories`.
- 10,000-history structural parser: the fresh binary read the existing 4,464,636,822-byte `fresh_run_validation/run_10000/dumpall_fresh.dat`, reported 2,036,435 records, 10,000 histories, and terminal NCOL=2. The previously retained `1,983,461`-record regression refers to a different validation artifact and remains documented separately.
- `audit_independent_seed_bins.py` plus `py_compile` passed; it returns the exact three crossing histories and L1 values above from committed derived CSVs.
- WSL was available, but Ubuntu-24.04 has no `cmake` or `g++` (`bash: cmake: command not found`, `bash: g++: command not found`); no WSL build was performed.
- `git diff --check` passed. Source scan found no local absolute paths in the experimental C++ implementation or scripts. Staging must remain limited to the four updated reports and the new audit script; the PHITS distribution, PHITS source, Dumpalls, and nonpublic stopping table are not staged.

### Follow-up push and Actions

The follow-up documentation/script changes were pushed to the research branch
as commit `7d305030de49a54fb26ceb807190ad2cb9c279b5`. GitHub Actions run
[39428338000](https://github.com/symmetrysorrow/TES-Programs/actions/runs/39428338000)
completed successfully: `energy-regressions` passed setup, checkout, HDF5
installation, CMake configure, build, and Python energy-accounting tests.
PR #6 was queried after the push: it remains open, Draft, and at head
`36dec9f139918d233faf0876736f43d1bbdb965c`; no PR #6 change or merge was
performed.
