"""Record the production pulse gate when the production mesh is unavailable."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    benchmark = json.loads(args.benchmark.read_text(encoding="utf-8"))
    report = {
        "status": "NOT_RUN",
        "dependency": str(args.benchmark.resolve()),
        "reason": "Production pulse waveform is gated on a valid current-branch production conformal mesh and matching CPU/GPU benchmark.",
        "benchmark_status": benchmark.get("status"),
        "small_mesh_reference": "artifacts/phase20_conformal/pulse_waveform_parity.json",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
