# Single-pixel parameter search

This search is the Elmer-versus-COMSOL single-pixel validation and one-factor
model-sensitivity stage.  There is no experimental single-pixel transient under
the same conditions as the PoST measurement, so it must not be used to fit the
PoST experimental waveform.  The experimental multivariate fit is implemented
separately in `scripts/search/post_multivariate_search.py`.

## First-stage variables

The initial one-factor-at-a-time scan contains the baseline plus ±20% changes
to five quantities:

- Pb heat capacity
- TES heat capacity
- membrane-stack effective heat capacity
- Stycast thermal conductivity (used as the first contact-path proxy)
- Pb thermal conductivity

Density is fixed.  `R_0`, `alpha`, `beta`, and `L_tes` remain fixed in this
first scan; they should be released only if the thermal variables cannot
explain the residual waveform shape.

For every candidate, the runner first solves the steady case and adjusts `G0`
with a bracketed log-space secant/bisection until the converged TES temperature
is within 0.01 mK of `T_0`.  This is required even for the nominal material set
because the coarse screening mesh has a slightly different effective thermal
balance.  Heat-capacity-only candidates should converge to the same calibrated
`G0`, while conductivity candidates are allowed to move it.
After the baseline calibration is complete, heat-capacity-only candidates reuse
that `G0` and verify it with one steady run.  Conductivity candidates start
their bracketed solve from the calibrated baseline value.

The screening pulse uses the coarse `mesh_shifted_merged` mesh and ends 0.59 ms
after the pulse.  It writes only the TES series CSV, not VTU files.  This stage
is for ranking sensitivities; promising candidates must later be rerun on the
3x mesh and over the full experimental decay window.

The generated search cases explicitly route result and restart files back to
`work/meshes/mesh_shifted_merged`.  This avoids the legacy Elmer behavior that
otherwise prepends a repository-root `mesh_shifted_merged` directory.

## Commands

Run from `Elmer-Projects`:

```powershell
python scripts/search/single_pixel_search.py prepare
python scripts/search/single_pixel_search.py score-existing
python scripts/search/single_pixel_search.py dry-run baseline
python scripts/search/single_pixel_search.py run baseline
python scripts/search/single_pixel_search.py run-all
```

`run-all --limit 3` executes only the first three candidates.  Each run is
cached under `results/spsearch_*`; search metadata, scores, and the leaderboard
are written under `artifacts/search/single_pixel`.

## Waveform objective

The reference and Elmer signals are treated as different observables.  Each is
baseline corrected and divided by its own post-pulse peak, so the score compares
shape rather than volts against amperes.  The baseline statistic defaults to
the pre-pulse median.  Signal column, unit scale, and pulse polarity are
configured independently in `reference_data` and `simulation_data` in
`single_pixel_search_config.json`.

The aligned waveform is separated into rise, peak, decay, and tail regions.
Each region receives one RMSE before the weighted mean is formed; this prevents
a densely sampled tail from overwhelming the rise and peak.  A bounded time
shift is optimized for every candidate, but a weak quadratic penalty prevents
the shift from freely absorbing a physical rise-time mismatch.  Both shifted
and unshifted scores are retained.

Every scored candidate writes `aligned_waveform.csv` beside `score.json`.  It
contains the common time grid, normalized reference, normalized simulation,
residual, and region label.

For an experimental voltage trace, change the reference configuration, for
example:

```json
"reference_data": {
  "delimiter": ",",
  "comments": "#",
  "time_column": 0,
  "signal_column": 1,
  "time_scale_to_ms": 1000.0,
  "signal_scale": 1.0,
  "response_direction": "auto"
}
```

`response_direction` may be `rise`, `drop`, or `auto`.  Filtering in the
experimental readout still changes shape; either apply the same transfer
function to the Elmer current before scoring or add it explicitly to the
readout model.

## Recommended execution order

Run the nominal case first because all later candidates need its calibrated
coarse-mesh `G0` as a starting point:

```powershell
python scripts/search/single_pixel_search.py prepare
python scripts/search/single_pixel_search.py run baseline
python scripts/search/single_pixel_search.py run-all
```

After all low/high one-factor candidates have scores, produce the sensitivity
ranking and the correlation matrix of sensitivity waveforms:

```powershell
python scripts/search/single_pixel_search.py analyze-sensitivity
```

The outputs are:

- `artifacts/search/single_pixel/sensitivity_summary.csv`
- `artifacts/search/single_pixel/sensitivity_correlation.csv`

Large absolute correlations identify parameter pairs that the single waveform
cannot distinguish reliably.  Keep only the influential, weakly redundant
variables for the next stage, or narrow their ranges before continuing.

Do not run a single-pixel multivariate waveform fit against PoST data.  After
the one-factor sensitivity results have selected the influential, weakly
redundant shared variables, run the PoST-specific multivariate workflow:

```powershell
python scripts/search/post_multivariate_search.py prepare-references
python scripts/search/post_multivariate_search.py prepare
python scripts/search/post_multivariate_search.py run-all --limit 3
python scripts/search/post_multivariate_search.py run-all --skip-scored
```

That runner uses the dual-TES 20 mm geometry, both experimental channels, five
ordered pulse-position groups, and the measured readout filter.  Its default is
48 reproducible log-space Latin-hypercube samples.

A weak variable-specific log prior is added to the phase-balanced waveform
objective.  The initial search ranges are intentionally moderate.  Expand a
range only when its best candidates sit near a bound and the sensitivity shape
is physically plausible.

## Releasing electrical parameters

Do not release `R_0`, `alpha`, `beta`, and `L_tes` together with every material
property in the first multivariate run.  First fit the thermal subset above.
If structured residuals remain, add electrical variables in this order:

1. `L_tes` or the measured readout transfer function for a global rise-time
   mismatch.
2. `alpha` and `beta` for peak-region and pulse-height-dependent mismatch.
3. `R_0` only with a steady-current constraint in addition to the normalized
   waveform objective.

Because this single-pixel stage has no matching experimental transient, it does
not produce an experimental posterior.  It produces sensitivity information
and steady-state constraints.  Those results define variable selection and
bounded priors for the PoST fit; PoST-specific contact area, contact thickness,
absorber-to-TES coupling, and left/right variation remain geometry-stage
parameters.
