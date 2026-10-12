#!/usr/bin/env python3
"""Recover an unflushed final TES row from the accepted iteration trace."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RUN = ROOT / "runs" / "attempt6" / "transient"
audit = {"cases": {}, "source": "final accepted iteration CSV row", "writes_original_series": False}
for circuit in (1, 2):
    series = RUN / f"dualtes_wsl_transient_20261013_newton_{circuit}_series.csv"
    iterations = RUN / f"dualtes_wsl_transient_20261013_newton_iterations_{circuit}.csv"
    out = RUN / f"dualtes_wsl_transient_20261013_newton_{circuit}_series_reconciled.csv"
    with series.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, skipinitialspace=True))
        header = f"time_s,tes_temperature_K,tes_current_A,tes_resistance_ohm,tes_power_W"
    with iterations.open(newline="", encoding="utf-8") as f:
        iteration_rows = list(csv.DictReader(f, skipinitialspace=True))
    final = iteration_rows[-1]
    final_time = float(final["time_s"])
    old_time = float(rows[-1]["time_s"])
    appended = final_time > old_time
    with out.open("w", newline="", encoding="utf-8") as f:
        f.write(header + "\n")
        writer = csv.writer(f, lineterminator="\n")
        for row in rows:
            writer.writerow([row[k] for k in ("time_s", "tes_temperature_K", "tes_current_A", "tes_resistance_ohm", "tes_power_W")])
        if appended:
            writer.writerow([final["time_s"], final["tes_temperature_K"], final["raw_current_A"], final["tes_resistance_ohm"], final["raw_power_W"]])
    audit["cases"][str(circuit)] = {
        "original_series_rows": len(rows), "reconciled_rows": len(rows) + int(appended),
        "original_last_time_s": old_time, "accepted_iteration_time_s": final_time,
        "appended": appended, "final_tes_temperature_K": float(final["tes_temperature_K"]),
        "final_tes_current_A": float(final["raw_current_A"]), "final_tes_resistance_ohm": float(final["tes_resistance_ohm"]),
        "final_tes_power_W": float(final["raw_power_W"]), "output": out.name,
    }
(ROOT / "series_flush_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
print(json.dumps(audit, indent=2))
