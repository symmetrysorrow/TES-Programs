# Independent analysis prompt

You are reviewing an evidence bundle for Elmer FEM performance. Read `README.md` first, then inspect every included solver log, SIF, manifest, CSV, and source note. The user wants a DualTES transient runtime close to 2x the single-pixel transient runtime. A full 75 ms rerun is not needed for the next analysis; two four-step diagnostics already exist.

Please independently:
- reconstruct comparable timings, distinguishing CPUTime, RealTime, rank-local, and global timings;
- trace the code paths represented by each timer and identify intervals omitted by the current instrumentation;
- explain the outer-vs-inner HYPRE solve timing discrepancy without assuming its cause;
- compare the short DualTES and single-pixel profiles, accounting for different mesh families and restart convergence histories;
- recommend the smallest next measurement and likely optimization candidates, with numerical-accuracy checks.

Avoid conclusions based on subtracting the 103.80 s Assembly CPUTime from the 37.204 s Phase24 wall-event total. Source-to-runtime binary provenance is not yet established. State uncertainties explicitly and cite filenames/line numbers or log entries for each conclusion.
