from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
COMSOL = ROOT.parent / "PoST_Simulations" / "ComsolResults" / "PoST_130.txt"
ELMER = ROOT / "results" / "case_dual20_pos30"
OUT = ROOT / "artifacts" / "dual_current_comparison"
OUT.mkdir(parents=True, exist_ok=True)

def crossing(t, y, level, rising):
    idx = np.where(y >= level if rising else y <= level)[0]
    if len(idx) == 0: return np.nan
    i = int(idx[0])
    if i == 0: return float(t[0])
    return float(t[i-1] + (level-y[i-1])*(t[i]-t[i-1])/(y[i]-y[i-1]))

def metrics(t_ms, raw, name, end_ms):
    base = float(np.mean(raw[t_ms < 20.0]))
    y = raw - base
    mask = t_ms <= end_ms
    tt, yy = t_ms[mask], y[mask]
    peak_i = int(np.argmax(np.abs(yy)))
    peak = float(yy[peak_i])
    lo, hi = 0.1*peak, 0.9*peak
    if peak < 0:
        rise = crossing(tt, yy, lo, False)
        fall = crossing(tt[peak_i:], yy[peak_i:], hi, True)
    else:
        rise = crossing(tt, yy, lo, True)
        fall = crossing(tt[peak_i:], yy[peak_i:], hi, False)
    return {"series": name, "baseline_A": base, "min_A": float(np.min(yy)),
            "max_A": float(np.max(yy)), "min_max_diff_A": float(np.ptp(yy)),
            "peak_corrected_A": peak, "rise_time_10_90_ms": abs(float(rise-tt[0])) if np.isfinite(rise) else np.nan,
            "fall_time_90_10_ms": abs(float(fall-rise)) if np.isfinite(fall) and np.isfinite(rise) else np.nan,
            "window_end_ms": float(end_ms), "baseline_window": "t < 20 ms"}, y

comsol = np.loadtxt(COMSOL, comments="%", encoding="utf-8")
tc = comsol[:, 0]
elmer_l = np.genfromtxt(ELMER / "tes_dual20_pos30_L_series.csv", delimiter=",", names=True)
te = elmer_l["time_s"] * 1e3
rows, series = [], []
for name, col in (("COMSOL_L", 3), ("COMSOL_R", 4)):
    row, y = metrics(tc, comsol[:, col], name, 195.0)
    rows.append(row); series.append((tc, y, name))
for side in ("L", "R"):
    dat = np.genfromtxt(ELMER / f"tes_dual20_pos30_{side}_series.csv", delimiter=",", names=True)
    row, y = metrics(dat["time_s"]*1e3, dat["tes_current_A"], f"Elmer_{side}", 100.0)
    rows.append(row); series.append((dat["time_s"]*1e3, y, f"Elmer_{side} (completed 100 ms)"))

with (OUT / "metrics.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
with (OUT / "comparison_timeseries.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["time_ms", "COMSOL_L_A", "COMSOL_R_A", "Elmer_L_A", "Elmer_R_A"])
    grid = np.unique(np.concatenate([tc, te]))
    vals = [np.interp(grid, t, y, left=np.nan, right=np.nan) for t, y, _ in series]
    for i, t in enumerate(grid): w.writerow([t, *(v[i] for v in vals)])

fig, ax = plt.subplots(figsize=(10, 5.5))
for t, y, label in series: ax.plot(t, y*1e6, label=label, lw=1.1)
ax.axvline(20.02, color="k", ls="--", lw=.8, label="pulse start")
ax.set(xlabel="Time [ms]", ylabel="Baseline-corrected current [µA]", title="Dual-TES current comparison")
ax.grid(True, alpha=.25); ax.legend(ncol=2); fig.tight_layout()
fig.savefig(OUT / "comparison_timeseries.png", dpi=180); plt.close(fig)
print(OUT / "metrics.csv")
print(OUT / "comparison_timeseries.png")
