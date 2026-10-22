# Pb HATCH resource and installation status

Date: 2026-10-22 (Asia/Tokyo)

## Data provenance and scope

The table was copied byte-for-byte from the previously validated local source:

`dump2event_physics_validation/egs5_stopping_power_reconstruction/runs/n100_stopping_probe_v2/egs5_stopping_power_coefficients.csv`

- Size: 130,754 bytes
- SHA-256: `20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9`
- Rows: 600; medium 2/PB: 150; `lelke`: 1–150
- PHITS/EGS5 density: 9.78 g/cm³
- Destination in this branch: `tes_cpp/dump2event/data/pb_egs5_phits3341.csv`

The original and destination were checked with the existing importer's
`validate_original` routine. The GitHub connector confirmed the repository is
Private; an authenticated GitHub collaborators query returned only
`symmetrysorrow`. The repository owner reports receiving permission from
PHITS administrator Nobuyuki Shikko for this table within the individual-use
scope. This does not authorize public distribution, mirroring, commercial use,
or broader repository access. The official EGS5 notice limits the manual and
included software to non-commercial use; the bundle remains restricted to
the stated private individual-use scope. See the [PHITS application and
pledge](https://phits.jaea.go.jp/image/application_format.pdf) and [EGS5
copyright notice](https://rcwww.kek.jp/research/egs/slac730-page2.pdf).

## Installation behavior

The normal invocation is:

```powershell
dump2event input.dumpall output.h5
```

The standard command needs only the Dumpall and new HDF5 paths. The native driver
uniquely associates one sibling PHITS input through `file(15)`, checks its
`file(6)` run log reports PHITS 3.341, and validates run-specific seed/source,
cutoff, material, geometry, tally and weight settings. Ambiguous or missing
associations fail closed; no physics JSON is generated on disk. The exact
generated effective contract and SHA-256 are saved in HDF5 v2 provenance.
The installed package's console entry point invokes the bundled native CLI.
The CLI locates the HDF5 v2 driver; the driver resolves the exact bundled CSV
and verifies its SHA-256. `TES_CPP_PB_HATCH_TABLE`, `--stopping-table`, and the
Python API's `stopping_table` argument are optional compatibility/diagnostic
overrides, not requirements. PoST skips the CSV picker when the checked
resource is available. Unsupported material or modified/missing table bytes
remain fail-closed. No physics formula or HDF5 schema changed.

## Verification

Branch: `feat/pb-hatch-private-bundle-20261022`
Base: `origin/main` at `5f129bfa2c4dce6c31ef7230f486c45d87655b89`

- MSVC 19.51.36260 Release build: passed.
- CTest: 1/1 passed.
- Focused dump2event Python/UI tests before the automatic-association follow-up: 81 passed, 7 skipped, 5 subtests passed.
- `pip install .` in a clean Python 3.14 virtual environment: passed. The
  Windows wheel (`tes_cpp_dump2event-0.1.0-cp314-cp314-win_amd64.whl`,
  352,834 bytes, SHA-256 `ad9094854a2d76c559bdcdb85477f91de94ffdd6cb2c11c31bfd1541448175d7`)
  contains the `dump2event.exe` console entry point, native CLI, physical
  reconstructor, transport exporter, HDF5 v2 driver/writer, and CSV.
- Installed-wheel resource resolution with `TES_CPP_PB_HATCH_TABLE` removed
  selected the packaged file and returned the expected SHA-256. CMake install
  resolution without the variable also selected the installed CSV.
- Source-tree PoST test confirmed the real packaged CSV suppresses its
  stopping-table picker.
- The archived 2 MeV/1,000-history validation report records 5,978,912 of
  5,978,912 PHITS additions matched and matching provenance SHA-256. The HDF5
  itself was not present in the checked repo/worktree paths, so it was not
  reread or regenerated in this task.
- Broad Python suite: 565 passed, 10 skipped, 5 subtests passed; 3 unrelated
  `PoST_Simulations` diagnostic tests failed: a high-frequency slope of
  `-39.436646 dB/decade` outside its `-40.0 ±0.5` tolerance, a missing
  `rms_residual_dB` fixture key, and a missing `SimpleNamespace` import. All
  `tes_cpp/dump2event` tests pass.
- The importer's source-validation routine was used. Its CLI privacy check
  requires GitHub CLI, which is not installed here; privacy and collaborator
  scope were independently checked through authenticated GitHub services.

## GitHub status

The feature commit `b5fbdfe3a973f7beafe70f9c73838b4a42cdc6a7` was pushed to the
Private repository on branch `feat/pb-hatch-private-bundle-20261022`. PR #14
is open at https://github.com/symmetrysorrow/TES-Programs/pull/14.

GitHub Actions run #119 (`39730893431`) failed. The job metadata has no steps
or logs, matching the repository's recent pre-job runner failures; this does
not establish a source-code test failure. The Actions service/runner cause is
unresolved, so the PR remains unmerged pending a runnable CI result. PR #6
remains open and Draft at head `36dec9f139918d233faf0876736f43d1bbdb965c`.


## Automatic Pb settings follow-up (2026-10-22)

The default CLI and Python API now take only Dumpall and output paths. The
converter scans sibling `*.inp` files and accepts exactly one whose `[Parameters]
file(15)` resolves to the selected Dumpall. The referenced `[Parameters] file(6)`
log must exist and contain exactly `Version = 3.341`. The selected input then
passes the existing Pb source, region/material/density, cutoff, tally, `esmin`,
seed, `maxcas`/`maxbch`, and unit-weight gates. Explicit physics JSON remains
available; when used, it must point to the unique input associated with that
Dumpall and the run-log version must match. No missing value is guessed.

The internally generated contract JSON and its SHA-256 are added to HDF5 v2
provenance alongside resolved-physics and source/table hashes. The normal Python
API is `tes_cpp.dump2event.dump2event(dump_path, output_path)`; the PoST dialog
no longer asks for a physics JSON and only offers a private CSV picker when the
verified bundled table is unavailable. A wheel-layout test found and fixed the
CMake/scikit-build table path so both the installed driver and UI locate the
resource under `<prefix>/tes_cpp/bin`.

Verification for this follow-up:

- Windows MSVC 19.51 Release build: passed; CTest 1/1 passed.
- Focused dump2event pytest: 84 passed, 7 skipped, 5 subtests passed.
- `pip install . --no-deps --target ...` from `tes_cpp/`: passed; wheel size
  519,284 bytes, SHA-256 `10da46bb9a7c54fccd2b0d54a57f780bc077b1357620ca3d98981975bc6f0480`.
- Installed native CLI invoked with exactly two paths on a software-only fixture.
  It resolved the sibling input, PHITS log, private table, and generated
  provenance, then correctly rejected the deliberately invalid fake Dumpall
  at transport-record parsing. This is an auto-configuration/dispatch smoke
  test, **not** a successful conversion of physical data.
- No real Dumpall was converted in this follow-up. The synthetic CLI smoke exit
  code was 1 with `Dumpall line 1 ... NCOL=0: invalid or non-finite numeric token`.
  Full physical conversion remains unverified by this change.
