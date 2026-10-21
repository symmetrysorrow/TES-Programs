# dump2event HDF5 v2 schema

Status: experimental storage/API layer. This schema preserves physics evidence;
it does not claim that every PHITS process is reconstructed or ready for
`Pulse_Ms`.

## File identity and compatibility

Root attributes are `format="tes-dump2event"`, `format_version=2`,
`schema="flat-offset-v2"`, `energy_unit="MeV"`, and `length_unit="cm"`.
Readers must inspect these attributes before interpreting `/events`.

Supported readers distinguish three on-disk forms:

| Form | Identity | Reader behavior |
|---|---|---|
| Legacy C++ hierarchical | format attribute, no version, `/events/<event>/<particle>` groups | Read through the legacy `iter_events()` adapter |
| Python flat v1 | same format, version 1, root flat datasets | Read through the v1 `iter_events()` adapter |
| Flat v2 | same format, version 2, `flat-offset-v2` | Read event summaries with `iter_events()` and deposits with `iter_deposits()` |

Unknown versions and malformed declared schemas fail closed. A v2 event summary
is not the old `event.json` particle mapping. `SEGMENT` rows remain segments;
the legacy API never turns them into a point deposit.

## Groups and columns

All tables are one-dimensional extendible datasets. Coordinate tables have
shape `(N, 3)`. The offset arrays start with zero and have one more element than
their parent table.

### `/events`

`event_id:int64`, `particle_offsets:int64`, `deposit_offsets:int64`,
`deposit_total_MeV:float64`, `complete:bool`,
`all_absorption:utf8` (`confirmed|rejected|undetermined`), and
`unresolved_count:uint64`. The total is the sum of valid deposit rows only; an
event with unresolved rows must be incomplete.

### `/particles`

`event_id:int64`, `particle_id:int64`, `parent_particle_id:int64`,
`parent_valid:bool`, `ityp:int32`, and `track_offsets:int64`. A parent ID of -1
is meaningful only when `parent_valid` is false.

### `/tracks`

`source_record_id:utf8`, `ncol:int16`, `particle_id:int64`, `x_cm/y_cm/z_cm`,
`energy_MeV`, `time`, `time_valid:bool`, and `region_id:int32`. Unknown time is
represented by `time_valid=false`; the stored numeric placeholder is not a
measurement.

### `/deposits`

`addition_id:utf8`, `event_id:int64`, `particle_id:int64` (-1 if unavailable),
`process:utf8`, `region_id:int32`, `energy_MeV:float64`, `energy_valid:bool`,
`kind:uint8` (`1=POINT`, `2=SEGMENT`, `3=UNRESOLVED`),
`start_xyz_cm/end_xyz_cm/score_xyz_cm:float64[3]`,
`start_position_valid/end_position_valid/score_position_valid:bool`,
`source_record_id:utf8`, and
`position_semantics:utf8`.

POINT and SEGMENT energy must be finite and nonnegative. UNRESOLVED has
`energy_valid=false` and is returned as Python `None`, never as a measured zero.
SEGMENT requires both endpoints. The PHITS representative score coordinate is
stored separately from those endpoints. Missing start, end, or score
coordinates are masked by their corresponding `*_position_valid=false` field;
their stored numeric placeholders are not positions.

### `/provenance`

Required attributes: `phits_version`, `dumpall_sha256`,
`physics_config_sha256`, `stopping_table_sha256`, `material`, `medium`,
`cutoff_MeV`, `algorithm_version`, and `position_semantics`. The writer rejects
missing fields. The workflow that creates these values must record their
actual sources; placeholder hashes are only suitable for tests.

Additional scalar attributes are preserved when provided, including
`source_input_sha256`, region/material IDs, weight assumptions, and diagnostic
mode. These are experiment metadata and do not change the v2 dataset layout.
Provenance revision 2 distinguishes byte identity from semantic settings:

| Attribute | Meaning |
|---|---|
| `physics_contract_sha256` | SHA-256 of the exact user-supplied JSON contract bytes, including whitespace and key order. |
| `resolved_physics_sha256` | SHA-256 of the canonical UTF-8 JSON for the resolved, supported semantic settings. |
| `physics_config_sha256` | Compatibility alias for `resolved_physics_sha256`; revision 2 requires equality. |
| `phits_input_sha256` / `source_input_sha256` | SHA-256 of the exact PHITS input bytes used for the run. |
| `dumpall_sha256` | SHA-256 of the exact converted Dumpall bytes. |
| `stopping_table_sha256` | SHA-256 of the exact stopping coefficient table bytes. |
| `phits_executable_sha256` | Optional SHA-256 supplied by the run contract for the PHITS executable. |
| `converter_executable_sha256` | SHA-256 of the `dump2event` executable that launched the v2 driver. |
| `reconstructor_executable_sha256` / `transport_exporter_sha256` | SHA-256 of the C++ physical-ledger and raw-transport helper executables. |
| `driver_sha256` | SHA-256 of the Python HDF5 v2 builder script. |

Revision 2 also stores `resolved_physics_json`, the canonical semantic object
whose digest is `resolved_physics_sha256`. It deliberately excludes file paths,
random seed text, and byte-level hashes. The PHITS input byte digest separately
identifies those exact input bytes. Readers treat a v2 file with no
`provenance_schema_revision` as revision 1 and preserve its existing meaning;
the on-disk datasets and version number do not change.

The resolved object is limited to the settings explicitly parsed by the
experimental Pb workflow. It is not a complete canonicalization of every PHITS
geometry, source, or transport card; consumers must retain the PHITS input hash
and must not interpret this field as proof that all input semantics were
captured.

## Writer and reader API

`tes_cpp.event_hdf5_v2.EventV2Writer` appends event rows into a temporary file.
On successful close it validates offset lengths and atomically replaces the
target. On failure it removes the temporary file and preserves any existing
target. `append_event(event_id, particles, deposits, complete=..., ...)` takes
flat mappings; particle track records are nested under their particle solely as
an input convenience, then stored in `/tracks`.

`tes_cpp.event_hdf5.iter_events()` keeps its legacy behavior for v1 and the old
C++ layout. For v2 it yields only event metadata and offset ranges.
`tes_cpp.event_hdf5.iter_deposits()` yields typed deposit mappings, including
both segment endpoints and score-coordinate validity.

## Scope and migration

The native `dump2event` CLI and both public Python wrappers now use the
PHITS 3.341 Pb/unit-weight EGS5 `PhysicalDepositLedger` plus raw-transport
exporter and streaming HDF5 v2 writer as their **only output converter**.
The obsolete signed reaction-difference converter, energy normalizer, and
legacy JSON/HDF5 v1 converter output path have been removed. The physical
contract JSON, matching Pb HATCH coefficients, and valid Dumpall are mandatory;
unsupported conditions fail closed. Explicit `--physics pb-egs5` remains an
accepted alias but no longer selects a different conversion algorithm.

The legacy v1/hierarchical **readers** remain temporarily available for old
data and for the preexisting `Pulse_Ms` workflow; they are not alternate
dump2event output modes. Do not interpret historical legacy files as v2.

`Pulse_Ms()` currently consumes point-deposit legacy files. It preflights the
format and rejects v2 before opening its pulse output. Segment-to-block
allocation needs a separately validated `event2pulse` model.
