# HDF5 v2 integration validation

The opt-in path reuses `EventV2Writer` and the existing POINT/SEGMENT/UNRESOLVED schema. Continuous charged-particle loss retains raw start/end chord, energy, event/particle, region, process, source record and addition ID. Dead-secondary local points use the raw QCLUSTS position. A continuous interval is not collapsed to an invented midpoint; its microscopic distribution and PHITS representative point remain unspecified unless separately validated.

The writer creates a temporary file in the destination directory, validates and atomically publishes only on success. Existing output paths are refused. Negative energies and any unresolved scoring candidates abort and remove the temporary output. `all_absorption` remains `undetermined`; the event `complete` flag means no unresolved candidate score rows in that event, not transport completeness or energy balance.

The existing 10k research HDF5 has 1,954,742 rows and 0 unresolved rows, but the cross-run provenance concern is recorded separately. This turn's integrated writer has not been run on an actual Dumpall because its required local stopping table was unavailable. Transactional failure behavior is covered by the existing writer tests; a complete CLI round-trip test is pending.
