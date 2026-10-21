# Pb HATCH stopping table provenance

The expected historical table was recovered in the main repository's private validation data and its bytes match the recorded expected SHA-256 exactly:

`20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9`

- File: `egs5_stopping_power_reconstruction/runs/n100_stopping_probe_v2/egs5_stopping_power_coefficients.csv`
- Size: 130,754 bytes; 600 rows total.
- Medium 2 / `PB`: 150 rows, `lelke` index 1 through 150.
- Columns: `medium, medname, lelke, eke0, eke1, ededx0, ededx1, pdedx0, pdedx1, phits_density, egs_density`.
- Medium-2 PHITS and EGS5 density fields are both 9.78 g/cm³.
- HDF5 provenance records material Pb, Dumpall material 2, EGS5 medium 2, 1 keV cutoff, and unit-weight-only scope.

This is the exact table required by the existing experimental engine; no substitute was used. The table remains in private local validation data and was not copied or staged into Git. Its originating probe is historical; the current 1,000-history run reuses the byte-verified coefficients and does not regenerate HATCH.
