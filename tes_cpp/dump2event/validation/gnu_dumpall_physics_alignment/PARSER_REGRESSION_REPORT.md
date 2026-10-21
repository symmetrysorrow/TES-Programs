# Parser regression report

Date: 2026-10-20. Worktree: `D:\github\TES-Programs-physdeposit-20261020`, branch `feat/dump2event-physical-validation-20261020`.

The parser’s data schema is based on PHITS 3.341 `ovly15.f:dumpdat`. The preexisting worktree change supports numeric values that continue onto a following physical line. This turn adds an input-source parameter and a fail-closed diagnostic that reports the source, row, record start, NCOL, expected field and observed token count for a malformed particle tuple. The converter now routes files beginning with standalone NCOL=1 through the structural parser without falling back to the permissive legacy line counter after a structural error. It does not fabricate values or make the persisted GNU dumps parseable.

Verified commands/results:

- `cmake --build build/physdeposit-vs --config Release --target test_physical_deposit_engine` — exit 0.
- `ctest --test-dir build/physdeposit-vs -C Release --output-on-failure` — 1/1 passed.
- `$env:TES_CPP_DUMP2EVENT_EXECUTABLE=(Resolve-Path build/physdeposit-vs/Release/dump2event.exe).Path; python -m pytest tes_cpp/dump2event/tests -q` — 23 passed.
- Regression fixtures cover the normal schema, a valid split numeric tuple, an internal truncated field, and EOF truncation.
- A converter-level regression writes a malformed modern NCOL=1 file to a temporary file and verifies `read_dump` rejects it with the input path, rather than accepting a partial history through the legacy parser.
- A read-only, 3,218-byte single-history file at `H:\hata\1332_142_136\1332keV_6\dumpall.dat_1_178256` was wrapped with temporary NCOL=1/2 markers outside the source location to exercise its actual wrapped records: parser exit 0, 4 records, 1 history, NCOL=2 terminal. The original file was not modified.
- Current raw audit executable, linked against the rebuilt parser, stops on the real GNU sample at line 5, record start line 2, NCOL=4, particle header: 11 combined values where 8 are required. It returns exit code 1 and produces no particle-class rows.
- Before the routing change, `dump2event --summary-only` silently returned exit 0 with one partial history for the same malformed file. Afterward it returns exit 1 with the detailed NCOL=4 error. This behavior change is limited to files carrying the PHITS 3.341 standalone NCOL=1 signature.
- The same exact structured failure was reproduced for all three GNU diagnostic files: 2 MeV thick Pb (123,789 internal calls), 511 keV thick Pb (52,897) and 511 keV thin Pb/vacuum (12,440).
- `git diff --check` — clean.

The earlier validated Windows Dumpall regression evidence (100, 1,000, and 10,000 histories) is preserved in the existing reports, but those raw files were not available at the expected PHITS data directory in this worktree session and were not re-read in this turn. The new GNU files are unparseable after the first history prefix, so they cannot support record-level physics regression. No physics formula changed.
