# Runtime and output sizes

| Job | Histories | Wall time | Main artifact size |
|---|---:|---:|---:|
| PHITS 0.511 MeV diagnostic | 100 | 74.1 s | Dumpall 354,930,838 bytes |
| HDF5 v2 CLI 0.511 MeV | 100 | 38.9 s | HDF5 75,482,132 bytes |
| PHITS 2 MeV diagnostic | 100 | 160.8 s | Dumpall 777,546,730 bytes |
| HDF5 v2 CLI 2 MeV | 100 | 107.3 s | HDF5 164,283,437 bytes |
| Default JSON CLI smoke test | 100 input histories | completed | JSON 34,470,912 bytes; 60 retained events |

The 1,000-history runs and CLI conversions are still being measured. The 10,000-history physical CLI attempt failed before publishing output, so no 10,000-history conversion performance number is reported.
