"""Generate diagnostic plots from pulse-contamination audit artifacts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case-dir", type=Path, required=True)
    args = parser.parse_args()
    import matplotlib.pyplot as plt
    case, plot = args.case_dir, args.case_dir / "plots"
    plot.mkdir(parents=True, exist_ok=True)
    library = json.loads((case / "pulse_template_library.json").read_text(encoding="utf-8"))
    thresholds = json.loads((case / "pulse_detection_thresholds.json").read_text(encoding="utf-8"))
    classified = json.loads((case / "noise_record_pulse_classification.json").read_text(encoding="utf-8"))
    partition = json.loads((case / "pulse_partitioned_noise_spectra.json").read_text(encoding="utf-8"))
    prediction = json.loads((case / "pulse_contamination_psd_prediction.json").read_text(encoding="utf-8"))
    coherence = json.loads((case / "partitioned_coherence.json").read_text(encoding="utf-8"))
    colors = {"CH0": "tab:blue", "CH1": "tab:orange"}
    for channel in ("CH0", "CH1"):
        templates = library["channels"][channel]["templates"]
        plt.figure(figsize=(8, 4))
        for name, row in templates.items():
            values = np.asarray(row["values"])
            step = int(row.get("sample_step", 1))
            plt.plot(np.arange(len(values)) * step / 500000.0 * 1e3, values, label=name)
        plt.xlabel("time from peak [ms]"); plt.ylabel("normalized amplitude"); plt.title(f"{channel} pulse templates"); plt.legend(fontsize=8); plt.tight_layout(); plt.savefig(plot / "pulse_templates.png" if channel == "CH0" else plot / "pulse_templates_ch1.png", dpi=160); plt.close()
        scores = [row[channel]["matched_score"] for row in classified["records"].values() if row[channel]["classification"] != "not_accepted"]
        plt.figure(figsize=(7, 4)); plt.hist(scores, bins=40, color=colors[channel], alpha=0.75); plt.axvline(thresholds["channels"][channel]["thresholds_by_template"]["full_pulse"]["fpr_0.001"], color="black", label="primary fixed FPR threshold"); plt.xlabel("best matched score"); plt.ylabel("records"); plt.title(f"{channel} pulse score distribution"); plt.legend(); plt.tight_layout(); plt.savefig(plot / "pulse_score_distribution.png" if channel == "CH0" else plot / "pulse_score_distribution_ch1.png", dpi=160); plt.close()
    ch0 = partition["subsets"]["CH0"]
    freq = np.asarray(ch0["all_accepted"]["frequencies_Hz"])
    example_scores = sorted((row["CH0"]["matched_score"] for row in classified["records"].values() if row["CH0"]["classification"] != "not_accepted"), reverse=True)[:40]
    plt.figure(figsize=(8, 4)); plt.plot(np.arange(len(example_scores)), example_scores, "o-"); plt.axhline(thresholds["channels"]["CH0"]["thresholds_by_template"]["full_pulse"]["fpr_0.001"], color="black", label="primary fixed FPR threshold"); plt.xlabel("top detected records [score rank]"); plt.ylabel("matched score"); plt.title("CH0 pulse detection examples"); plt.legend(); plt.tight_layout(); plt.savefig(plot / "pulse_detection_examples.png", dpi=160); plt.close()
    for filename, names in (("clean_vs_all_asd.png", ("pulse_free", "all_accepted")), ("clean_vs_contaminated_asd.png", ("pulse_free", "definitely_contaminated"))):
        plt.figure(figsize=(7, 4))
        for name in names:
            row = ch0[name]
            if row["asd"] is not None: plt.loglog(freq, row["asd"], label=name)
        plt.xlabel("frequency [Hz]"); plt.ylabel("ASD"); plt.grid(True, which="both"); plt.legend(); plt.tight_layout(); plt.savefig(plot / filename, dpi=160); plt.close()
    pred = prediction["channels"]["CH0"]
    plt.figure(figsize=(7, 4)); plt.loglog(freq, pred["pulse_only_asd"], label="pulse-only"); plt.loglog(freq, pred["predicted_all_asd"], label="clean + pulse-only"); plt.loglog(freq, pred["measured_all_asd"], "--", label="measured all"); plt.xlabel("frequency [Hz]"); plt.ylabel("ASD"); plt.grid(True, which="both"); plt.legend(); plt.tight_layout(); plt.savefig(plot / "predicted_vs_measured_all_psd.png", dpi=160); plt.close()
    plt.figure(figsize=(7, 4)); plt.loglog(freq, np.asarray(pred["pulse_only_asd"]) ** 2, label="pulse-only PSD"); plt.xlabel("frequency [Hz]"); plt.ylabel("PSD"); plt.grid(True, which="both"); plt.legend(); plt.tight_layout(); plt.savefig(plot / "pulse_only_psd.png", dpi=160); plt.close()
    simulation = json.loads((case / "pulse_free_simulation_comparison.json").read_text(encoding="utf-8"))
    x = np.asarray(simulation["frequencies_Hz"]); plt.figure(figsize=(7, 4)); plt.fill_between(x, [row["proxy_min"] for row in simulation["rows"]], [row["proxy_max"] for row in simulation["rows"]], alpha=0.25, label="pulse-gated model envelope"); plt.plot(x, [row["proxy_median"] for row in simulation["rows"]], label="model median"); plt.plot(x, [row["pulse_free_experimental_normalized"] for row in simulation["rows"]], "o-", label="pulse-free experiment"); plt.xscale("log"); plt.yscale("log"); plt.xlabel("frequency [Hz]"); plt.ylabel("normalized ASD"); plt.grid(True, which="both"); plt.legend(fontsize=8); plt.tight_layout(); plt.savefig(plot / "clean_vs_simulation.png", dpi=160); plt.close()
    plt.figure(figsize=(7, 4))
    for name, row in coherence["subsets"].items():
        if row: plt.semilogx(row["frequencies_Hz"], row["coherence"], label=name)
    plt.xlabel("frequency [Hz]"); plt.ylabel("magnitude-squared coherence"); plt.ylim(0, 1.05); plt.grid(True, which="both"); plt.legend(); plt.tight_layout(); plt.savefig(plot / "coherence_all_clean_contaminated.png", dpi=160); plt.close()


if __name__ == "__main__":
    main()
