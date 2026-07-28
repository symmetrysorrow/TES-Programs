from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
comsol = np.loadtxt(ROOT / "reference" / "SignglePixel.txt", comments="%", encoding="utf-8")
elmer = np.genfromtxt(ROOT / "tes_mpi_comsol_grid_series.csv", delimiter=",", names=True)
tc, ic = comsol[:, 0], comsol[:, 3] * 1e-6
te, ie = elmer["time_s"] * 1e3, elmer["tes_current_A"]
end = min(tc.max(), te.max())

def calc(t, raw, name):
    base = float(np.mean(raw[t < 20.0])); y = raw - base; m = t <= end
    return {"series": name, "baseline_A": base, "min_A": float(y[m].min()), "max_A": float(y[m].max()), "min_max_diff_A": float(np.ptp(y[m])), "window_end_ms": float(end), "baseline_window": "t < 20 ms"}

rows = [calc(tc, ic, "COMSOL single-pixel"), calc(te, ie, "Elmer 4rank all-tet (completed series)")]
out = ROOT / "artifacts" / "single_pixel_comparison_now"; out.mkdir(parents=True, exist_ok=True)
with (out / "metrics.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
grid = np.unique(np.concatenate([tc[tc <= end], te[te <= end]]))
bc, be = ic - rows[0]["baseline_A"], ie - rows[1]["baseline_A"]
with (out / "timeseries.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["time_ms", "COMSOL_current_A", "Elmer_current_A"])
    for t in grid: w.writerow([t, np.interp(t, tc, bc), np.interp(t, te, be)])
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(tc[tc <= end], bc[tc <= end]*1e6, label="COMSOL single-pixel")
ax.plot(te[te <= end], be[te <= end]*1e6, label="Elmer 4rank all-tet")
ax.axvline(20.02, color="k", ls="--", lw=.8, label="pulse start")
ax.set(xlabel="Time [ms]", ylabel="Baseline-corrected current [µA]", title="Single-pixel current comparison (available interval)")
ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(out / "timeseries.png", dpi=180); plt.close(fig)
print(out / "metrics.csv")
