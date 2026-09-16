# Phase24 first pre-pulse backend parity

Generated: `2026-09-16T03:22:51.528279+00:00`

## Fixed conditions

- Same mesh, restart, pulse definition, inner circuit, Phase24 assembly, BDF2, and first `18 us` timestep
- One nonlinear iteration per case; this is a matrix/RHS diagnostic, not a convergence qualification
- Only linear backend differs: native HYPRE versus direct MUMPS

## Cases

| backend | exit | ALL DONE | matrix dump | first TES current [µA] |
|---|---:|---:|---|---:|
| Phase24 / HYPRE FlexGMRES + BoomerAMG | 0 | True | `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_baseline_backend_parity_18us\dumps\hypre_firstsolve` | n/a |
| Phase24 / direct MUMPS | 0 | True | `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_baseline_backend_parity_18us\dumps\mumps_firstsolve` | n/a |

## Matrix/RHS comparison

- Matrix (`_a.dat`): `{"available": true, "left": "D:\\Github\\TES-Programs\\Elmer-Projects\\artifacts\\phase24_baseline_backend_parity_18us\\dumps\\hypre_firstsolve_a.dat", "right": "D:\\Github\\TES-Programs\\Elmer-Projects\\artifacts\\phase24_baseline_backend_parity_18us\\dumps\\mumps_firstsolve_a.dat", "left_records": 1447196, "right_records": 1447196, "union_records": 1447196, "nonzero_difference_records": 0, "max_absolute_difference": 0.0, "max_relative_difference": 0.0}`
- RHS (`_b.dat`): `{"available": true, "left": "D:\\Github\\TES-Programs\\Elmer-Projects\\artifacts\\phase24_baseline_backend_parity_18us\\dumps\\hypre_firstsolve_b.dat", "right": "D:\\Github\\TES-Programs\\Elmer-Projects\\artifacts\\phase24_baseline_backend_parity_18us\\dumps\\mumps_firstsolve_b.dat", "left_records": 87534, "right_records": 87534, "union_records": 87534, "nonzero_difference_records": 0, "max_absolute_difference": 0.0, "max_relative_difference": 0.0}`

Interpretation: a material matrix/RHS difference means the Phase24 assembly or constraint path differs. If they match but the first solved current differs, inspect backend solve/scaling and RHS-to-solution handling.

Generated project: `D:\Github\TES-Programs\Elmer-Projects\artifacts\phase24_baseline_backend_parity_18us\phase24_baseline_backend_parity.json`
