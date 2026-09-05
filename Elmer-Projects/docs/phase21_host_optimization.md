# Phase21 host-side optimization and replacement benchmark

Phase21 carries the Phase20 conformal/shared-node HYPRE route forward with a
different objective: reduce host overhead and measure end-to-end replacement
time. It does not enable unconditional AMG/preconditioner reuse. The earlier
`Linear System Refactorize = False` probe remains rejected because both CPU
and GPU failed to converge at timestep 2.

## Implemented low-risk cache

The TES UDF now caches immutable run data:

- TES circuit constants and file paths per circuit instance;
- TES body element connectivity and node-index ranges in
  `tes_parallel_circuit.f90`;
- pulse spatial geometry, discrete normalization, and the temporal interval
  factor once per timestep in `tes_transient_heat_source.f90`.

The parallel circuit UDF emits `TESParallelCircuitProfile` markers containing
the cached element/node counts, CPU-time buckets, and (after rebuilding the
UDF) `SYSTEM_CLOCK` wall buckets. The implementation does not change the
electrical equations, nonlinear relaxation, pulse energy, or temperature
field. Its classification is `SAFE_OPTIMIZATION_WITH_NO_MEASURED_SPEEDUP`.

## Measurement results

The CPU 50-step cache regression was exact on the production mesh:

- TES temperature/current/resistance/Joule power: PASS;
- absorber temperature: PASS;
- all 90,872 temperature nodes: PASS, max difference `0 K`.

The measured CPU wall time was `205.50 s` for the Phase20 baseline and
`210.37 s` for the single Phase21 run (`0.977x`). This is not an optimization
claim; it shows that the cache alone does not overcome the host-side
assembly/solver architecture and run-to-run noise. The new report therefore
leaves the backend verdict as `GPU_VERDICT_PENDING_RUNTIME`; the cache is not
reported as a wall-speedup.

CPU I/O matrix, production 7-step:

- full validation output: `42.90 s`;
- VTU off: `39.80 s`;
- restart result off: `36.00 s`;
- series and iteration CSV off: `40.97 s`.

Only the full-output and VTU-off rows are suitable for production timing
claims. Output-disabled rows are diagnostic overhead measurements.

Fine transient CPU 7-step completed in `158.99 s`. The GPU variants were
prepared but were not rerun because the current WSL session reports
`no CUDA-capable device is detected`; Phase20 GPU results remain the last
validated GPU measurements.

The common replacement artifact now retains per-case states. Mortar steady,
1-step, and 7-step are `DONE`; conformal CPU steady, 1-step, and 7-step are
`DONE` with the provisional `1e-7` HYPRE tolerance; all conformal GPU cases are
`BLOCKED_NO_CUDA`. It reports `S_algorithmic` for completed CPU windows and
leaves `S_gpu`/`S_total` null until a GPU runtime is available.

The Phase22 tolerance sweep completed conformal CPU steady and 7-step cases at
`1e-6`, `1e-7`, and `1e-8`. The common Mortar field offset is about
`1.91e-2 K`; this is much larger than the within-conformal tolerance error and
is therefore retained as a separate route/mesh/reference gate. Electrical and
pulse-series comparison is explicitly `REFERENCE_UNAVAILABLE` because the
retained Mortar run emitted no machine-readable TES series. Thus `1e-7` is a
provisional common-replacement candidate, not a final physical-accuracy
verdict.

The wall report separates total solver wall, output/I/O wall, HYPRE setup,
HYPRE Krylov, UDF wall (when the rebuilt UDF marker is present), and an
unclassified residual. HeatSolve native instrumentation is staged in the
Phase22 native integration source; the installed solver must be rebuilt before
its `HeatSolveWallBreakdown` and local-FEM/global-insertion markers populate
production artifacts. CPU_TIME values remain a separate CPU profile and are
never added to the wall breakdown.

## Artifacts

See [transient_wall_breakdown.json](../artifacts/phase21_host_optimization/transient_wall_breakdown.json),
[circuit_profile.json](../artifacts/phase21_host_optimization/circuit_profile.json),
[static_dynamic_matrix_analysis.json](../artifacts/phase21_host_optimization/static_dynamic_matrix_analysis.json),
[udf_cache_benchmark.json](../artifacts/phase21_host_optimization/udf_cache_benchmark.json),
[fine_transient_cpu_gpu.json](../artifacts/phase21_host_optimization/fine_transient_cpu_gpu.json),
[io_overhead_benchmark.json](../artifacts/phase21_host_optimization/io_overhead_benchmark.json),
[mortar_conformal_replacement_benchmark.json](../artifacts/phase21_host_optimization/mortar_conformal_replacement_benchmark.json),
and [production_backend_recommendation.json](../artifacts/phase21_host_optimization/production_backend_recommendation.json).
The tolerance study is [hypre_physical_tolerance_study.json](../artifacts/phase22_physical_tolerance/hypre_physical_tolerance_study.json).

The canonical preparation commands are:

```powershell
python scripts/prep/prepare_phase21_host_optimization.py
python scripts/prep/prepare_phase21_mortar_conformal_replacement.py
python scripts/analysis/assemble_phase21_host_optimization_reports.py
python scripts/analysis/assemble_phase22_physical_tolerance.py
```
