# Unsupported conditions

- Any PHITS version other than 3.341; regions/materials/media other than region 1000 / Pb / material 2 / EGS5 medium 2.
- Non-unit primary/secondary weights, including dynamic changes.
- Cutoffs other than 1 keV, alternate tally units/modes, or unverified HATCH tables.
- GNU Dumpall serialization with truncated 80-character NCOL=4 fields; no parser-side repair is attempted.
- Photon-local energy loss inferred from parent/child energy differences, positron rest-mass local deposition, and unverified annihilation/relaxation branches.
- Gamma process coverage beyond the measured historical Co-60 tally slice.
- Complete event HDF5 generation when any candidate branch is unresolved.

The current opt-in mode is an experimental path to the existing limited engine; it is not ready as the default converter or for scientific production output.
