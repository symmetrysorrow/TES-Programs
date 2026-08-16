# PoST multivariate search

The multivariate fit is now performed on the PoST/dual-TES geometry, because
there is no experimental transient waveform from a single-pixel device under
the same operating conditions.  The single-pixel runner remains useful only
for Elmer-versus-COMSOL validation and one-factor model sensitivity.

## What is fitted

Each candidate uses `mesh_dual_20mm_localrefine` and runs one shared steady
case plus a single pulse position, `case_dual20_lr_pos30` (a hit near one
edge of the absorber).  The experimental side matches this with two
single-event references, chosen the same way `PoST_Simulations/Fitting_LSE.py`
picks its high/low pulses: among the selected events, the one with the
largest raw CH0 peak (`high`) and the one with the smallest raw CH0 peak
(`low`), with no averaging and no left/right fraction involved in the
selection itself.

`high` is compared directly against the `pos30` simulation.  `low` reuses the
*same* `pos30` run with `swap_simulation_channels: true` — by the mesh's
left/right mirror symmetry, swapping which simulated channel is called "left"
turns a near-left hit into the equivalent of a near-right hit, so the pipeline
never has to solve a second, mirrored position in Elmer.  A five-target
version of this search (`left_edge`, `left_mid`, `center`, `right_mid`,
`right_edge`, selected by CH0-fraction percentile bands of 16 events each)
was used earlier; it was replaced because two of its five pulse cases
(`right_mid`, `right_edge`) turned out to be byte-for-byte duplicates of
`left_mid`/`left_edge` — same physics, only the output filenames differed —
and because a single extreme event, while noisier than a 16-event average,
was found (on real data) to reach further into the position range than the
90th/10th percentile bands did.  The selected keys and their actual mean
channel fractions are recorded in
`artifacts/search/post_multivariate/reference_manifest.json`.

Both experimental channels are baseline corrected and filtered with the
measured second-order 2 kHz Bessel readout.  Before forming the position proxy,
each channel is divided by its median selected-event peak.  This removes the
large CH0/CH1 relative readout-gain offset without fitting it as a thermal
parameter; the applied factors are recorded in the reference manifest.  The
two Elmer TES-current drops are resampled at 100 kHz and receive the same
filter.  A pair is then divided by the sum of its two channel peaks.  This
removes the remaining unknown global current-to-voltage gain but preserves the
corrected left/right peak fraction, which is included explicitly in the
objective.  One common time shift is fitted across all positions and both
channels; separate per-channel or per-position shifts are not allowed.

The initial variables remain the shared thermal subset selected by the
single-pixel sensitivity study:

- Pb heat capacity;
- TES heat capacity;
- membrane effective heat capacity;
- Stycast thermal conductivity;
- Pb thermal conductivity.

`G0` is not treated as a freely waveform-fitted variable.  For every candidate
it is calibrated against the mean steady temperature of the left and right TES
stacks, keeping the operating point separate from the transient shape fit.

## Current CH0-only scoring

The current configured method uses CH0 only, because CH1 from this
measurement set is not considered sufficiently reliable for fitting.  The
CH0 waveform from the CH0-largest event is compared with the larger Elmer
response.  The CH0 waveform from the CH0-smallest event is compared with the
smaller Elmer response, obtained by the existing mirrored-channel mapping.

Each retained waveform is normalized by its own peak, so the waveform term
tests shape.  The high/low CH0 peak ratio is included separately as a squared
log-ratio error; since both measurements are CH0, it does not depend on the
unknown CH0-to-CH1 gain.  CH1 remains available in the reference CSV for
diagnostics, but does not contribute to the objective.  One common time shift
is fitted across high and low.

This mode writes separate score and trace files, preserving the older paired
CH0/CH1 results.  Use the rescore-all command with skip-scored to score
already-computed Elmer traces under the CH0-only method without rerunning
Elmer.

## Commands

Run from `Elmer-Projects`:

```powershell
python scripts/search/post_multivariate_search.py prepare-references
python scripts/search/post_multivariate_search.py prepare
python scripts/search/post_multivariate_search.py run-all --limit 3
python scripts/search/post_multivariate_search.py run-all --skip-scored
```

`prepare-references` reads the experiment location configured in
`post_multivariate_search_config.json`, selects representative event groups,
and writes reproducible averaged reference CSV files.  `prepare` creates the
48 deterministic Latin-hypercube candidates without running Elmer.  Use
`dry-run <candidate-id>` to inspect the generated dependency chain before the
first solver run.

Candidate scores and aligned traces are written below
`artifacts/search/post_multivariate/candidates`.  `leaderboard.csv` ranks the
combined high/low objective.

## Interpretation limits

`high` and `low` are single events, not averages, so each carries more shot
noise than the old 16-event bands and is not perfectly reproducible run to
run (a different measurement session could select a different pair of
extreme events). Their CH0 fraction is also a byproduct of ranking by raw
CH0 peak, not a controlled position proxy, so it constrains positional
ordering only loosely. Consequently this stage can constrain shared thermal
parameters and gross absorber-to-TES transport, but it cannot uniquely
separate material conductivity from exact contact geometry.  Contact area,
contact thickness, and left/right fabrication differences should be added only
after the shared-material search is stable, and should retain priors rather
than replacing the shared parameters with fixed values.
