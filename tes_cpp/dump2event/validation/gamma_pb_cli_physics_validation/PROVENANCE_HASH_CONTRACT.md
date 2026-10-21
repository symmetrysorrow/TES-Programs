# HDF5 v2 provenance hash contract

Revision 2 separates byte identity from effective settings. `physics_contract_sha256` is the digest of exact input JSON bytes. `resolved_physics_sha256` hashes sorted, compact UTF-8 canonical JSON of the recognized settings. `physics_config_sha256` remains as a compatibility alias and must equal the resolved hash. `phits_input_sha256` and `source_input_sha256` hash exact input bytes; `dumpall_sha256` and `stopping_table_sha256` hash the exact converted data and table.

For reproducibility, the v2 result also records the hashes of the PHITS executable when explicitly supplied by the contract, the launching `dump2event` binary, the physical reconstruction helper, the transport exporter, and the Python driver. No executable hash is synthesized when its file or explicit contract digest is unavailable.

The resolved object records PHITS version, region/material and EGS5 medium, electron/positron cutoff, unit-weight assumption, energy and deposit semantics, algorithm revision, and the recognized Co-60 or mono-gamma source energy. It is deliberately scoped: it does not canonicalize every PHITS geometry, material card, random seed, or transport card. The exact PHITS input byte hash remains necessary evidence.

`EventV2Writer` now rejects revision-2 provenance if the alias differs, the resolved JSON is invalid/noncanonical, or its SHA-256 does not match `resolved_physics_sha256`. A v2 file without a revision attribute retains the legacy revision-1 interpretation; the dataset layout and HDF5 format version are unchanged.

Regression coverage: JSON whitespace/key-order changes alter the raw contract hash but not the resolved hash; semantic settings alter the resolved hash; PHITS input and stopping table byte changes alter their respective byte hashes; missing files and invalid hashes fail closed. `verify_hdf5_provenance.py` compares CLI JSON to HDF5 readback.
