"""Create an auditable next-run plan from completed SinglePixel convergence pilots.

This does not claim that Elmer's BDF solver is internally adaptive.  Instead it
implements the safer outer-loop equivalent: compare one refinement at a time,
accept an axis only below a specified waveform-change budget, and propose the
next spatial or time-step pilot.  The generated JSON is intended to be read by
the case-preparation script (or by a user) before launching the next run.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "singlepixel_resolution_optimization.json"
DEFAULT_METRICS = ROOT / "artifacts" / "comparison" / "stycast_z16_resolution" / "metrics.json"
DEFAULT_OUT = ROOT / "artifacts" / "optimization" / "singlepixel_resolution_plan.json"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--metrics", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return parser.parse_args()


def pair_result(pairs: dict[str, dict[str, float]], name: str, limit: float) -> dict[str, object]:
    item = pairs.get(name)
    if item is None:
        return {"status": "not_measured", "pair": name}
    value = float(item["pct_comsol_peak"])
    return {
        "status": "accepted" if value <= limit else "refine_further",
        "pair": name,
        "waveform_change_pct_of_peak": value,
        "limit_pct_of_peak": limit,
        "max_difference_time_us": float(item["time_us"]),
    }


def main() -> None:
    args = parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    metrics = json.loads(args.metrics.read_text(encoding="utf-8"))
    limit = float(config["waveform_change_limit_pct_of_peak"])
    pairs = metrics["pairwise"]

    stycast = pair_result(pairs, "z16_1p25us_minus_z32_1p25us", limit)
    time = pair_result(pairs, "z32_1p25us_minus_z32_0p625us", limit)
    spatial = {
        "stycast_layers": {
            **stycast,
            "selected_layers": 32 if stycast["status"] == "accepted" else 64,
            "next_case": None if stycast["status"] == "accepted" else "stycast 64 layers, 0.625 us",
        },
        "other_axes": {
            name: {
                "status": "not_measured",
                "recommended_first_comparison": [values[0], values[1]],
                "note": "Keep all already accepted axes fixed; use the selected 0.625 us time step.",
            }
            for name, values in config["spatial_axes"].items()
            if name != "stycast_layers"
        },
    }
    time_plan = {
        **time,
        "selected_early_step_us": 0.625 if time["status"] == "accepted" else 0.3125,
        "schedule_to_validate": [
            {"start_us": 0.0, "end_us": 105.0, "step_us": 0.625, "purpose": "initial response"},
            {"start_us": 105.0, "end_us": 225.0, "step_us": 1.25, "purpose": "first tail coarsening pilot"},
            {"start_us": 225.0, "end_us": None, "step_us": 2.5, "purpose": "provisional late-tail step; accept only after comparison"},
        ],
        "acceptance_rule": (
            f"At each switch, compare against a run retaining the preceding step; "
            f"accept only when the maximum waveform change is <= {limit:.3g}% of the COMSOL peak."
        ),
    }
    result = {
        "input": {"config": str(args.config.relative_to(ROOT)), "metrics": str(args.metrics.relative_to(ROOT))},
        "criterion": {"waveform_change_limit_pct_of_peak": limit, "window_us": config["comparison_window_us"]},
        "current_recommendation": {
            "mesh": {
                "stycast_layers": spatial["stycast_layers"]["selected_layers"],
                "tes_layers": 1,
                "note": "TES is intentionally held at one layer until its own one-axis pilot is measured; its 0.16 um thermal time scale is far shorter than the present microsecond response.",
            },
            "time": {"initial_step_us": time_plan["selected_early_step_us"]},
        },
        "spatial": spatial,
        "time": time_plan,
        "execution_order": [
            "Use the accepted Stycast=32 / early dt=0.625 us baseline.",
            "Run the 105--225 us tail-coarsening pair before using 1.25 us in production.",
            "For each remaining spatial axis, run only its first two candidates; continue only when that pair exceeds the budget.",
            "Repeat lateral stack and absorber local-size pilots last, because they change the total element count most strongly.",
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
