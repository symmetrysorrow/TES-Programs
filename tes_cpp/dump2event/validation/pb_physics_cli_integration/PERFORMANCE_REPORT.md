# Performance report

| Measurement | Result |
|---|---|
| Existing research 10k generation | ~188.02 s for the prior HDF5 rebuild, as previously recorded; this was not the integrated CLI |
| Existing research HDF5 size | 587,741,058 bytes, previous artifact |
| Integrated CLI 10k elapsed time | Not measured |
| Integrated CLI HDF5 size | Not generated |
| Integrated CLI peak memory | Not measured |
| Temporary CSV volume | Not measured for the integrated CLI |

The new path streams C++ CSV exports and grouped event rows into `EventV2Writer`; it does not create a single giant JSON object. Performance must be measured after the exact stopping table and a lossless same-run Dumpall are available. No performance improvement is claimed.
