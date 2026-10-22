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
dump2event input.dumpall output.h5 --physics-config pb_config.json
```

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
- Focused dump2event Python/UI tests: 81 passed, 7 skipped, 5 subtests passed.
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

No changes were made to PR #6; GitHub reports it open. The CSV is only in
this local feature branch until the branch is pushed to the Private repository.
Record the feature commit, PR, and Actions result here after publication. Do
not merge if the private-repository Actions runner fails before executing
steps; report local tests separately in that case.
