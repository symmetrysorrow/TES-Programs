# 10,000-history normal CLI attempt

**Result: failed closed; no HDF5 was published (CLI exit code 1).** The CLI was executed against the completed PHITS 3.341 Co-60, 10,000-history baseline Dumpall in a separate output directory. The Dumpall is 2,736,869,598 bytes with SHA-256 `c11c51299c6c77a7e8a531a640ed72b31c8455fb2d0a17f2b03d3716a2815060`; its PHITS input SHA-256 is `7401c49903c41317825af7eed2a8c239033ea7505514cecebe0f1524459f007e`. The source run completed with PHITS exit code 0 at 2026-10-21 11:48:11 JST. It failed while the raw transport exporter parsed the first history. The exact error was:

```text
dumpall_transport_export: line 5 (record at line 2, NCOL=4,
field=particle header (particle_id, material_id, ITYP, KTYP, JTYP, MTYP,
RTYP, OLDWT)): combined physical rows contain 11 numeric values; expected 8.
Possible truncation or malformed logical row; values were not inferred.
dump2event: experimental Pb reconstruction failed; subprocess status=2
```

The input raw file and PHITS input were not changed. The CLI refused to publish HDF5 and the PHITS tally from that historical Co-60 run is not used as a substitute for a same-run mono-gamma comparison. The failure occurs before any addition-level join. A new 10,000-history mono-gamma PHITS run was not attempted in this validation window; the 511 keV and 2 MeV same-run experimental matrix remains at 100 histories until the 1,000-history runs complete.

This is a parser/serialization compatibility failure, not a measured physics mismatch. No fields were guessed, and the raw file was not repaired in place. The completed mono-gamma matrix contains full 100-history comparisons for both energies and a full 1,000-history same-run comparison for 0.511 MeV. The 2 MeV / 1,000 CLI attempt separately failed closed on one unsupported ledger branch; it is not counted as a complete conversion. No new 10,000-history mono-gamma PHITS run was attempted.
