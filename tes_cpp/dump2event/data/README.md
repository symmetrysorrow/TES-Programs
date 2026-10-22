# Frozen Pb HATCH coefficients for dump2event

`pb_egs5_phits3341.csv` is the exact, checksum-pinned 600-row CSV from the
PHITS 3.341 / EGS5 stopping-power probe. Medium 2 is Pb (`PB`) with 150 rows,
indices `lelke=1..150`, at PHITS and EGS5 densities of 9.78 g/cm³.

Expected file size: **130,754 bytes**
Expected SHA-256:

```text
20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9
```

The original local copy was:

```text
dump2event_physics_validation/egs5_stopping_power_reconstruction/runs/
n100_stopping_probe_v2/egs5_stopping_power_coefficients.csv
```

The repository owner reports permission from a PHITS administrator for this
table to be used within the PHITS individual-use scope. This repository is
private; this is not permission for public posting, mirroring, commercial use,
or access beyond that individual-use scope. Do not publish the CSV or add it to
a public repository. The EGS5 notice restricts the manual and included
software to non-commercial use; retain the notice with any applicable copy or
derivative. See the [PHITS application and pledge](https://phits.jaea.go.jp/image/application_format.pdf)
and [EGS5 copyright notice](https://rcwww.kek.jp/research/egs/slac730-page2.pdf).

## Runtime resolution

The Python HDF5 v2 driver and C++ CLI use the bundled table automatically.
The installed CMake/Python package places it next to the physical CLI and the
driver resolves that managed resource. Regular use does not need
`TES_CPP_PB_HATCH_TABLE`, `--stopping-table`, or a Python `stopping_table`
argument. Those explicit paths remain available for diagnostics and must still
match the pinned SHA-256.

The resource resolver fails closed if the bundled file is missing, has a
different hash, or describes an unsupported material. PoST's Pb EGS5 menu
skips the CSV picker when this verified resource is installed. HDF5 v2
provenance records the actual table SHA-256.

## Re-importing the historical table

Use the import tool only with the original file. It validates the source
bytes, size, CSV header, row counts, medium, indices, and densities before
copying. In a GitHub checkout, first confirm that the destination repository is
Private and that access is limited to the permitted individual-use scope:

```powershell
gh repo view --json isPrivate --jq .isPrivate
python tes_cpp/dump2event/tools/import_verified_pb_table.py `
  "D:\path\to\egs5_stopping_power_coefficients.csv" `
  --confirm-private-repository
```

The tool copies only; it does not stage, commit, or push. Never replace the
historical bytes with interpolated or third-party coefficients. This data
resource does not change the physical reconstruction equations or HDF5 v2
schema.
