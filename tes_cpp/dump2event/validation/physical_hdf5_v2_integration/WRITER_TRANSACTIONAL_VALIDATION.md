# HDF5 v2 event transaction validation

`append_event` now has a validate/stage phase followed by a buffer commit phase. Validation errors occur before changes to event IDs, addition IDs, offsets, or dataset buffers. If list extension itself raises, buffer lengths are restored. If an HDF5 flush fails, the writer is marked unusable; it cannot publish the temporary file or accept further events.

Checks include duplicate event/particle/addition IDs, absent parent or deposit particle IDs, finite and nonnegative known energies, finite coordinates/time, required POINT and SEGMENT positions, nonempty source-record IDs/processes, HDF5 integer bounds, unresolved/complete consistency, and absorption-state consistency. Existing output is replaced only after successful close and structural validation.

Tests reject negative and NaN energy, duplicate addition IDs, duplicate event IDs, absent parent IDs, missing segment endpoint, and complete events with unresolved deposits. Each caught failure leaves the buffers and global ID sets unchanged; a valid event can then be appended and the resulting file contains only valid events. A separate test preserves an existing output file after failed conversion.

Local result: `tes_cpp/tests/test_event_hdf5_v2.py` and the v1 compatibility tests pass. The complete local Python suite result is recorded in `HDF5_V2_CI_REPORT.md`.
