# Current Repro Archive

This archive freezes the exact inputs, mesh, and output files used to obtain the current steady-state TES result on 2026-07-11.

## Frozen result

- TES nodal average temperature: `0.168581911088883 K`
- TES volume average temperature: `0.168581952265808 K`
- TES min temperature: `0.150000002 K`
- TES max temperature: `0.168584140 K`
- Solver case: `case_constant_power.sif`
- Mesh directory: `mesh_shifted_merged`

## Included files

- `case_constant_power.sif`
- `generated/tes_shared_variables.sif`
- `generated/tes_case_constant_power.sif`
- `generated/tes_materials.sif`
- `elmer_project.json`
- `elmer_geometry.json`
- `mesh_shifted_merged/mesh.boundary`
- `mesh_shifted_merged/mesh.elements`
- `mesh_shifted_merged/mesh.header`
- `mesh_shifted_merged/mesh.names`
- `mesh_shifted_merged/mesh.nodes`
- `mesh_shifted_merged/entities.sif`
- `mesh_shifted_merged/case_constant_power.ep`
- `mesh_shifted_merged/case_constant_power_t0001.vtu`
- `tmp_extract_tes_avg_ep.py`
- `tmp_extract_tes_volume_avg.py`

## How to rerun

From this archive root:

```bat
ElmerSolver case_constant_power.sif
python tmp_extract_tes_volume_avg.py mesh_shifted_merged mesh_shifted_merged\case_constant_power_t0001.vtu
```

Expected extraction output:

```text
tes_nodal_average=0.168581911088883
tes_volume_average=0.168581952265808
```

## Important note

This archive preserves the current reproducible state exactly as used for the confirmed result. If later project files or generated SIF files change, this archive should still reproduce the frozen result as long as the same Elmer version and mesh files are used.
