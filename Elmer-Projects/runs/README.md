# Runs Directory

This directory is reserved for reproducible frozen run bundles.

## Policy

- Source of truth stays in the repository root as `elmer_project.json`.
- Generated files are copied into each run bundle.
- Mesh and solver outputs are frozen per run.
- A run bundle must contain:
  - `manifest.json`
  - `README.md`
  - copied input files
  - copied generated files
  - copied mesh files
  - copied solver results

## Create a frozen run

```bat
python freeze_repro_run.py current_reference
```

This creates:

- `runs/current_reference/`
- `runs/current_reference.zip`
