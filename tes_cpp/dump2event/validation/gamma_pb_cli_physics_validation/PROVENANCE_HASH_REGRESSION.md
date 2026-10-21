# Provenance hash regression

Revision-2 hash unit tests verify the distinction between raw contract bytes and semantic settings, the expected sensitivity of input/table hashes, and rejection of missing or malformed values. The HDF5 writer additionally checks that canonical resolved JSON is valid, compact, key-sorted UTF-8 and hashes exactly to `resolved_physics_sha256`; `physics_config_sha256` must equal that digest. Runtime SHA fields are validated as 64 hexadecimal characters.

Observed same-run 100-history HDF5 readback used CLI JSON and HDF5 `/provenance`. For both mono-gamma runs, the revision-2 alias and resolved hashes matched, and the JSON/result hashes matched the HDF5 attributes. PHITS executable SHA was supplied in the run contract; converter, ledger, exporter, and driver SHA were computed from the actual files. No executable identity was inferred from a version string.

Test command: `python -m pytest tes_cpp -q` — 51 passed, 10 skipped. The skips are existing test-environment skips; they are not counted as passing tests.
