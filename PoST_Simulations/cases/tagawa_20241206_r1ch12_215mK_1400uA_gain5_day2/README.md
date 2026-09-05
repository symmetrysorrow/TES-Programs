# Target case: r1ch12, 215 mK, 1400 uA, gain 5

This is the versioned case record for the acquisition requested in the noise-mismatch investigation. It is deliberately **not runnable yet**: the exact target TES physics and readout calibration are not present in the repository or in the target folder inspected on 2026-09-05. Unknown values are `null` in [`input.json`](input.json); they must not be replaced with values from the generic [`PoST_Simulations/input.json`](../../input.json).

The folder name supplies useful labels (`215mK`, `1400uA1400uA`, `gain5`), but those labels are recorded as provenance only. In particular, the 1400 uA label is not promoted to calibrated TES current, and the 215 mK label is not promoted to `T_c` or to a measured TES temperature. See [`provenance.json`](provenance.json) for the source and confidence of every value.

Known acquisition facts are sufficient to fix the FFT grid: `rate=500000 Hz`, `samples=100000`, `cutoff=10000 Hz`, and `df=5 Hz`. The target experimental ASD is summarized at the requested frequencies in [`comparison_summary.json`](comparison_summary.json); the current readiness decision is in [`audit_status.json`](audit_status.json). The simulation columns are intentionally `null`, because an absolute no-added-noise comparison would otherwise manufacture a physical operating point.

Regenerate the versioned summary from the raw target records with:

```powershell
python PoST_Simulations/subScript/target_case_report.py `
  --experiment-path G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2 `
  --case-dir PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2 `
  --output PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/comparison_summary.json
```

To unblock the comparison, add an independently sourced target calibration/operating-point record covering at least `T_c`, `T_bath`, TES `I/R`, `alpha`, `beta`, `L`, `n`, `C_abs`, `C_tes`, `G_abs-abs`, `G_abs-tes`, `G_tes-bath`, `n_abs`, and the CH0 output calibration. Then copy those values into `input.json`, rerun the operating-point/stability and Python/C++ parity audits, and only then populate the simulation ASD and ratio columns.
