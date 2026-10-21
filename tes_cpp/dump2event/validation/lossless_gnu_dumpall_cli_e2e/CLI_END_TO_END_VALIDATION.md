# Normal CLI to HDF5 v2 end-to-end validation

**Not passed.** The required normal `dump2event` entry point was not run successfully on repaired PHITS data. Two required inputs were unavailable for this run: a lossless same-run Dumpall and the private PHITS 3.341 Pb medium-2 stopping table (expected SHA-256 `20fd86eb0b57483bd19f20a8388d2c1ecb0431a02fa99f0933e2c3c4748514e9`). The historical GNU raw fails the structured parser at the first NCOL=4 tuple. No HDF5 v2 file is presented as a successful result.

The prior code has an explicit opt-in `--physics pb-egs5 --hdf5-version 2` path and previously passed a smoke rejection test on malformed input. That verifies fail-closed behavior, not the requested success path. Python API success with real data was not tested in this task.
