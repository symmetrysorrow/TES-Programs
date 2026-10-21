# Push and CI

Push is limited to experimental opt-in wiring and evidence reports, not a validated production physics switch.

## Pre-push checks

- MSVC Release build: succeeded for the standard CLI, physical ledger, transport exporter, and C++ test target.
- CTest: 1/1 passed.
- `python -m pytest tes_cpp -q`: 44 passed, 10 skipped.
- Default converter regression: 1,000 histories from the existing 10k raw input processed with exit 0; this did not compare produced event bytes to a prior file.
- Exclusive 10k ledger-class audit: 1,954,742 rows, zero negative, zero unresolved; counts match the checked expected classes.
- C++→Python CLI failure smoke: GNU diagnostic Dumpall rejected at NCOL=4 physical line 5; exit nonzero and output HDF5 absent.
- `git diff --check` and staged diff check passed. Private run data, PHITS source, stopping table, raw Dumpall, and HDF5 outputs are excluded. The two pre-existing untracked research folders remain untouched and unstaged.
- PR #6 was checked separately and remains open, Draft, unmerged, with head `36dec9f139918d233faf0876736f43d1bbdb965c`.

## Push and CI

Pending the authorized push and GitHub Actions result for the final branch head. Record the SHA, remote head, and workflow URL/status after verification. This report does not claim that the independent-seed or same-run raw/internal provenance gates passed.
