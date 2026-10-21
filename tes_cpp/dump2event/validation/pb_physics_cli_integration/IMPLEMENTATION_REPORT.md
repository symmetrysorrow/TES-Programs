# Pb opt-in physical deposit path: implementation and status

## Implemented

- The ordinary `dump2event` executable now recognizes the explicit `--physics pb-egs5 --hdf5-version 2 --physics-config ... --stopping-table ...` form and dispatches to the existing PHITS-scored `PhysicalDepositLedger` plus the streaming HDF5 v2 writer. The ordinary path and its legacy energy-difference logic remain unchanged.
- `tes_cpp.python.tes_cpp.dump2event.dump2event()` exposes the same opt-in. Default calls still require/use `input_energy`; the Pb mode requires the JSON contract and stopping table and rejects legacy normalize/filter/audit switches.
- Config validation checks PHITS 3.341, region 1000/Pb/material 2/medium 2, 1 keV electron and positron cutoffs, unit-weight scope, MeV and the validated `[t-deposit]` mode. The named PHITS input is checked for explicit seed, source lines, material/cell definition, cutoffs, `icntl`, `negs`, `e-mode`, and tally region/settings; the seed is stored with the input hash in provenance. The parser also checks every Dumpall WT/OLDWT and secondary weight and rejects non-unit values.
- Input Dumpall, PHITS input, JSON contract, and stopping table SHA-256 values are stored in HDF5 provenance. Invalid, unresolved, negative, or failed conversions abort the writer's temporary HDF5; an existing destination is never replaced.
- HDF5 provenance distinguishes terminal-NCOL structural validation, supported-score completeness, undetermined transport completeness, unassessed energy balance, and `all_absorption=undetermined`.

## Validation performed this turn

- MSVC Release build: succeeded.
- CTest: 1/1 passed.
- `python -m pytest tes_cpp -q`: 43 passed, 10 skipped.
- Existing default converter: processed 1,000 histories from the 10,000-history source with exit 0; this confirms the default CLI path still operates. It is not a byte-for-byte event output regression.
- CLI help, Python API option assembly, and the C++→Python driver dispatch were exercised. The integrated CLI was pointed at the GNU diagnostic Dumpall; it failed closed at the first malformed NCOL=4 row and left no output HDF5. A successful opt-in HDF5 conversion was **not** run: the local HATCH stopping table used in the earlier private validation is not present in this worktree or its private output folder.

## Not validated / integration decision

The normal CLI is wired to the research engine, but production conversion remains experimental and must not be treated as complete. The historical 1,954,742-row match used raw/diagnostic artifacts whose run identity is not yet closed: the earlier Dumpall and later PHITS internal-trace Dumpall have different SHA-256 values, and the GNU diagnostic Dumpall currently fails the structured parser at its first NCOL=4 particle header. The old run and later diagnostic are not being assumed identical. A new independent-seed paired run, same-run repaired Dumpall/internal-trace join, successful end-to-end CLI conversion, and gamma-process coverage remain outstanding.

See [provenance audit](VALIDATION_PROVENANCE_AUDIT.md), [unsupported conditions](UNSUPPORTED_CONDITIONS.md), and [final decision](FINAL_DECISION.md).
