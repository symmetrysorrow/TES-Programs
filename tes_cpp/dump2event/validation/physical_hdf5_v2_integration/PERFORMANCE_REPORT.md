# HDF5 v2 streaming performance

| Histories | Input bytes | Records | Deposits | HDF5 bytes | End-to-end conversion seconds | Peak observed Python working set |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | 54,036,774 | 38,345 | 38,077 | 11,523,228 | 5.070 | not sampled |
| 1,000 | 287,559,946 | 204,835 | 201,971 | 60,886,302 | about 23 | not sampled |
| 10,000 | 2,783,824,000 | 1,983,461 | 1,954,742 | 587,654,169 | 264.013 | 317,542,400 bytes (~303 MiB) |

The 10,000-history memory figure is the maximum sampled working set of the Python process. C++ child-process memory was not included, so this is not a total process-tree peak. The time includes both raw parser passes and HDF5 writing. Earlier synthetic 10,000-event storage performance is not used as a real-Dumpall benchmark.
