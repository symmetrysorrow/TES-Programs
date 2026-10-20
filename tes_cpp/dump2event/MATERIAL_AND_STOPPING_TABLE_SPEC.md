# Material and EGS5 stopping-table input specification

## Supported mapping

The current experiment accepts only PHITS region 1000 with Dumpall material
ID 2 mapped to EGS5 medium 2 and material label `PB`. Every record in scope
must carry the expected material ID. Other region/material combinations fail
closed. This mapping is part of the validation configuration, not a general
PHITS material database.

## Table schema

`Egs5StoppingTable::read_csv` consumes the HATCH-derived local CSV with one
row per EGS5 log-energy index and fields:

`medium,material,log_energy_index,eke0,eke1,ededx0,ededx1,phits_density,egs_density`.

The table is not committed because it is a locally generated/nonpublic PHITS
artifact. Reproducible runs should record its SHA-256 and explain the HATCH
source and material setup. The validated local table had SHA-256:

```text
20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9
```

The wrapper manifest additionally hashes the input Dumpall and a canonical
JSON physics configuration. It rejects PHITS version values other than
`3.341`. The C++ executable itself takes the table and Dumpall as explicit
arguments and does not contain local absolute paths.

## Numerical expression

For kinetic energy `E` and EGS5 medium `m`, the experiment follows the
PHITS 3.341 EGS5 helper `egs5edxde`:

```text
lelke = trunc(eke1 * ln(max(E, 0.001 MeV)) + eke0)
dedx(E,m) = ededx1(lelke,m) * ln(E) + ededx0(lelke,m)
E_loss = dedx(E,m) * step
```

The interpolation/indexing details are implemented from the matching PHITS
3.341 source and table layout. The table must be generated for the same HATCH
material, density and EGS5 medium; a table with a different setup is invalid.
No fallback table, density correction, or material guess is applied.

The independent-seed NCOL=10 term uses the raw incoming energy and raw
Dumpall endpoint chord, matching `ein-ecc` for the inspected PHITS branch.
The formula is supported by source and the exact same-seed comparison in
`INDEPENDENT_SEED_ROOT_CAUSE.md`; it is not claimed for other particles,
materials or PHITS versions.

## Provenance boundary

PHITS Fortran source, binary, the HATCH input, and the table itself must stay
outside Git. The source files inspected locally were `src/talls08.f` and
`src/egs5edxde.f`; source hashes from the local installation were recorded in
the validation notes:

```text
talls08.f    8773d0a2aee5a60236dc6c353dea5fb5401d1a8558f5edf71b9383b7bc5448bf
egs5edxde.f  63fb27b4e75babd5e5d14915c63f23590ea2d03aa64b531caf49456aec2f308d
```
