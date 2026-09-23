#!/usr/bin/env python3
"""Run the three-point frozen-power MUMPS diagnostic for one density probe."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.support import run_phase24_thermal_network_localization as base  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--key", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--template",
        type=Path,
        default=ROOT / "artifacts/phase24_thermal_network_localization/phase24_membrane_stycast_variant_v2_1p00P.sif",
    )
    parser.add_argument("--repeat", action="store_true")
    args = parser.parse_args()

    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    base.OUT = args.output
    case = base.MeshCase(
        args.key,
        f"Phase24 Stycast density probe {args.key}",
        args.mesh,
        args.template,
        101,
        1804,
        (
            base.Interface("TES_to_membrane", "TES->membrane", 1104, 1305, 101, 103),
            base.Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, 101, 102),
            base.Interface("Stycast_to_substrate", "Stycast->substrate", 1205, 1004, 102, 108),
        ),
    )
    manifest = [base.run_case(case, fraction, repeat=args.repeat) for fraction in (0.95, 1.0, 1.05)]
    (args.output / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
