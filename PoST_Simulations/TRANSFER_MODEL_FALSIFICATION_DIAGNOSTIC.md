# Low-order transfer-model falsification diagnostic

## Purpose

The detector-side diagnostics have not produced a stable parameter change that
raises the 5--15 kHz deficit while lowering the 40--100 kHz excess. The
existing required-transfer diagnostic already reports

`H_required(f) = measured normalized ASD / model normalized ASD`.

This diagnostic asks a narrower question: can that residual shape be
represented by a smooth low-order single-path transfer magnitude built only
from stable real poles and minimum-phase real zeros?

It does not alter the production optimizer, detector model, analog Bessel,
alias fold, or digital analysis filter.

## Tested transfer families

The script fits the following normalized magnitude families:

- unity;
- 1 pole;
- 1 pole / 1 zero;
- 2 poles / 1 zero;
- 2 poles / 2 zeros;
- 3 poles / 2 zeros;
- 3 poles / 3 zeros.

Every family is proper and normalized to unity at 1 kHz. Pole and zero corner
frequencies are positive, so the tested sections are stable and minimum phase.

The default corner search is 300 Hz to 2 MHz. The search uses differential
evolution followed by a least-squares polish.

## Run

From the repository root in PowerShell:

```powershell
python PoST_Simulations/subScript/transfer_model_falsification_diagnostic.py --summary PoST_Simulations/.noise_optimization_work_rsh_sweep/summary.json
```

Default output:

```text
PoST_Simulations/.noise_optimization_work_rsh_sweep/transfer_model_falsification_diagnostic.json
```

Useful optional controls are:

```text
--corner-min-hz
--corner-max-hz
--de-maxiter
--seed
--rms-threshold-db
--max-threshold-db
```

The default screen tolerance is 1 dB RMS and 3 dB maximum absolute residual.
Those values are explicitly a configurable diagnostic tolerance, not a
physical prior.

## Output interpretation

For every family the JSON records:

- fitted pole and zero corners;
- RMS, mean, and maximum absolute residual in dB;
- residuals by the same fit bands used by the noise optimizer;
- the production shape score after multiplying the frozen detector spectrum
  by the fitted residual transfer;
- whether the configured screen tolerance is met.

The output also reports the smallest passing family, the family with minimum
RMS dB error, and the family with the best corrected production shape score.

A positive result means only that a low-order smooth residual magnitude has
enough shape leverage. Magnitude-only fitting cannot identify where that
response lives in the acquisition chain, and it does not determine the actual
phase or distinguish minimum-phase behavior from a system with additional
all-pass/non-minimum-phase structure.

A negative result is stronger: within the stated corner-frequency box and
numerical tolerance, the tested low-order single-path pole/zero magnitude
models cannot explain the required residual transfer. That would argue
against interpreting the current discrepancy as one missing simple readout
filter.
