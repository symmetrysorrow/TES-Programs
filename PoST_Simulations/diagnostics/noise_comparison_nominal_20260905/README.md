# Nominal noise comparison (2026-09-05)

This directory contains a local generated diagnostic copy. It is not currently
versioned as a complete artifact: the repository `.gitignore` excludes the
CSV (`*.csv`) and the PoST-specific ignore excludes the PNG (`*.png`). The
small, versionable target-case summary is under
`PoST_Simulations/cases/tagawa_20241206_r1ch12_215mK_1400uA_gain5_day2/`.

- `comparison.csv`: generated comparison data (local, ignored) from
  `subScript/midband_residual.py`.
- `comparison.png`: comparison plot (local, ignored).
- `model_audit.json`: operating-point, Jacobian, and stability diagnostic
  (versionable JSON, if explicitly added).

The comparison used:

- `PoST_Simulations/input.json` as the nominal simulation input;
- the experiment directory configured by `--experiment-path`;
- the experiment's 500 kS/s, 100,000-sample, 10 kHz analysis settings;
- `white/readout ASD = 0`;
- no fitted parameter override and no empirical residual noise.

Reproduce from the repository root:

```powershell
python PoST_Simulations/subScript/midband_residual.py `
  --work-dir D:\desktop\post_midband_residual_audit `
  --csv D:\desktop\post_midband_residual_audit.csv `
  --plot D:\desktop\post_midband_residual_audit.png `
  --experiment-path G:\tagawa\20241206\r1ch12_215mK_1400uA1400uA_difftrig5e-5_rate500k_samples100k_gain5_day2 `
  --base-input-path PoST_Simulations/input.json

python PoST_Simulations/subScript/noise_model_audit.py `
  --input PoST_Simulations/input.json `
  --json D:\desktop\post_midband_model_audit.json
```
