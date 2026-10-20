# Real Dumpall to HDF5 v2 execution

All runs used read-only PHITS 3.341 inputs in `dump2event_physics_validation/overnight_physics_reconstruction/runs/`, the explicit Pb/medium-2 stopping table from `egs5_stopping_power_reconstruction/runs/n100_stopping_probe_v2/`, and fresh HDF5 destinations under the ignored `build/physdeposit-vs/` work folder. Commands used `build_physical_hdf5_v2.py` with the built `dump2event_physical_reconstruct.exe` and `dumpall_transport_export.exe`. The script hashes the actual Dumpall, PHITS input, stopping table, and canonical configuration. The table's SHA-256 is `20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9`.

Command form used (the actual per-run file paths were passed as arguments):

```powershell
python tes_cpp/dump2event/validation/build_physical_hdf5_v2.py `
  --dumpall <run>/dumpall.dat --source-input <run>/trace.inp `
  --stopping-table <validated-Pb-table> `
  --reconstructor build/physdeposit-vs/Release/dump2event_physical_reconstruct.exe `
  --track-exporter build/physdeposit-vs/Release/dumpall_transport_export.exe `
  --output build/physdeposit-vs/physical_hdf5_v2_verified_<n>.h5 `
  --work-dir build/physdeposit-vs/hdf5v2_verified_<n>
```

The output was then read with `verify_physical_hdf5_v2.py` in either `--internal-additions` or `--phits-histogram` mode. Every input reached parser `NCOL=2`; PHITS run logs report `PHITS_EXIT=0`.

| Dataset | Dumpall bytes / SHA-256 | Parser records / histories | HDF events / particles / tracks / deposits | unresolved / negative | HDF bytes | seconds |
|---|---|---:|---:|---:|---:|---:|
| `n100_lossless_v2` | 54,036,774 / `d6f1dad78ed7d768a65a14e6b17642ec7ec14cb14e77d157a85c62ded4a34164` | 38,345 / 100 | 100 / 849 / 38,245 / 38,077 | 0 / 0 | 11,523,228 | 5.070 |
| `n1000_independent_seed` | 287,559,946 / `bbfb0f548dd29c105f1ca19c8198f61c285bdbce9bd198bde42ce98d5ace2930` | 204,835 / 1,000 | 1,000 / 5,166 / 203,835 / 201,971 | 0 / 0 | 60,886,302 | about 23 |
| `n10000_independent_seed` | 2,783,824,000 / `d254b3afd13230be84ad8f733c8aaa89811724f832a506a0d5b421d051c06b1c` | 1,983,461 / 10,000 | 10,000 / 49,327 / 1,973,461 / 1,954,742 | 342 / 0 | 587,654,169 | 264.013 |

PHITS input and canonical configuration hashes, in table order:

| Dataset | Input `trace.inp` SHA-256 | Canonical physics configuration SHA-256 |
|---|---|---|
| 100 | `bba748a2247c00635e547ff7159b199a4a0e1b3bc4aebb3059c7b9f6fc7e851e` | `3962a74ecee9e21bed620b1ec14484557e75208a960c067994d186700cee7765` |
| 1,000 | `09392460b09f31b2186dbd8e53ef995435d52a36ac988ec024e075cd9dff5003` | `674febde19e70029a594834ba51be609ed1ee191476cb70135b042263000b21f` |
| 10,000 | `0719b4155c2b0c571e109974bb29558efa9b79f11bd895a11f57f6fd3e2e577e` | `c5cec6ea60f2e1431ac99ef0f52db3e5bd2665e49f937c7495cd68645262f1c5` |

Every input reached the normal `NCOL=2` terminator. The 100-history raw is byte-identical to the separate `n100_stopping_probe_v2` raw used by the preexisting PHITS addition oracle (SHA-256 above). The 1,000- and 10,000-history runs have PHITS normal-end evidence (`PHITS_EXIT=0`), but only the 100-history input has a same-run per-addition internal trace available. The other two are compared to their `[t-deposit]` output at histogram level.

The 10,000-history run exits from the physical ledger with its documented partial status (process exit 1 because 342 rows are unsupported). HDF5 is intentionally emitted as a diagnostic file with those rows explicitly unresolved and 9,991 events complete for the supported tally subset; no event is marked fully absorbed. The measured known-energy sum is `332.05456333632753 MeV`; unresolved energy is unknown and excluded, not set to zero.
