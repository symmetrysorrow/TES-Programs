# CLI and Python API

## CLI

```text
dump2event input.dumpall event.h5 --physics pb-egs5 --hdf5-version 2 \
  --physics-config pb_config.json --stopping-table pb_egs5.csv
```

The option is opt-in. It delegates to the existing ledger and streaming HDF5 v2 writer, requires the physical ledger and transport-export helpers alongside the executable, and locates the installed/repository Python driver or `TES_CPP_PB_PHYSICS_DRIVER`. Python is selected by `TES_CPP_PYTHON`, then `python` on Windows or `python3` elsewhere. Nonzero subprocess status fails the conversion.

## Python API

```python
dump2event(dump_path, output_path,
           physics="pb-egs5",
           physics_config="pb_config.json",
           stopping_table="pb_egs5.csv")
```

Existing calls continue to use `input_energy` and the legacy converter. Pb mode requires `.h5`/`.hdf5`, config and stopping table; it rejects legacy normalize/filter/audit options. The mode is experimental and version/material/weight scoped.

`pb_config.example.json` documents the validated contract. No raw data, local table, PHITS source, or local absolute path is embedded in the sample. Missing binaries or files fail before publishing output. The C++ CLI subprocess integration was compiled and its help/dispatch path inspected; a successful data conversion is not yet verified because the stopping table is unavailable locally.
