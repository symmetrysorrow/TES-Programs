# Phase24 HYPRE native residual trace

Generated: `2026-09-16T07:48:01.082111+00:00`

- Case: `1e-8`, `Linear System Max Iterations=10000`
- Diagnostic environment: `PHASE24_HYPRE_RESIDUAL_TRACE=1`
- Exit code: `1`
- ALL DONE: `False`
- Trace enabled marker: `True`
- Candidate trace lines: `20005`
- Native iteration records: `20000`
- Native minimum residual: `1.267862e-07`
- Native final residual: `1.267862e-07`

## Interpretation

The complete native trace is in `solver.log`. The head and tail below are only a compact index. Inspect whether residuals decrease steadily, stagnate, or become unstable before the 10000-iteration stop.

### Trace head

```text
PHASE24_HYPRE_RESIDUAL_TRACE_ENABLED method=901 max_iterations=10000 tolerance=1.00000000000000002e-08
Initial L2 norm of residual: 4.042852e-05
    1    4.030830e-05    0.997026   4.030830e-05
    2    2.613609e-05    0.648405   2.613609e-05
    3    2.154738e-05    0.824430   2.154738e-05
    4    2.136904e-05    0.991723   2.136904e-05
    5    1.638112e-05    0.766582   1.638112e-05
    6    1.156859e-05    0.706215   1.156859e-05
    7    1.073996e-05    0.928372   1.073996e-05
    8    1.055812e-05    0.983068   1.055812e-05
    9    8.886301e-06    0.841656   8.886301e-06
   10    7.371882e-06    0.829578   7.371882e-06
   11    7.049430e-06    0.956259   7.049430e-06
   12    6.936624e-06    0.983998   6.936624e-06
   13    6.278445e-06    0.905115   6.278445e-06
   14    5.636906e-06    0.897819   5.636906e-06
   15    5.550144e-06    0.984608   5.550144e-06
   16    5.471149e-06    0.985767   5.471149e-06
   17    5.038566e-06    0.920934   5.038566e-06
   18    4.779655e-06    0.948614   4.779655e-06
```

### Trace tail

```text
 9982    1.268148e-07    1.000000   1.268148e-07
 9983    1.268142e-07    0.999995   1.268142e-07
 9984    1.268141e-07    0.999999   1.268141e-07
 9985    1.268135e-07    0.999995   1.268135e-07
 9986    1.268133e-07    0.999999   1.268133e-07
 9987    1.268132e-07    0.999999   1.268132e-07
 9988    1.268122e-07    0.999992   1.268122e-07
 9989    1.268093e-07    0.999976   1.268093e-07
 9990    1.268089e-07    0.999998   1.268089e-07
 9991    1.268075e-07    0.999989   1.268075e-07
 9992    1.268047e-07    0.999978   1.268047e-07
 9993    1.268043e-07    0.999997   1.268043e-07
 9994    1.268012e-07    0.999975   1.268012e-07
 9995    1.267935e-07    0.999940   1.267935e-07
 9996    1.267935e-07    1.000000   1.267935e-07
 9997    1.267912e-07    0.999982   1.267912e-07
 9998    1.267872e-07    0.999968   1.267872e-07
 9999    1.267863e-07    0.999993   1.267863e-07
 10000    1.267862e-07    0.999999   1.267862e-07
PHASE24_HYPRE_SOLVE_FAILURE phase=backend_status method=901 solve_status=256 iterations=10000 final_relative_residual=1.26786240530270873e-07 requested_tolerance=1.00000000000000002e-08 max_iterations=10000 matrix_epoch=2 preconditioner_epoch=2 case=3 hypre_error=256 hypre_global_error=256
```

Solver log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_hypre_tol1e8_max10000_1us\solver.log`
Launcher log: `D:\Github\TES-Programs\Elmer-Projects\results\case_phase24_short_hypre_tol1e8_max10000_1us\residual_trace_launcher.log`
