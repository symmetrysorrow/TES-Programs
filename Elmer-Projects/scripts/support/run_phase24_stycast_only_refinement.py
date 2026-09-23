#!/usr/bin/env python3
"""Run the fixed-power three-point TES-frozen Stycast refinement control."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.support import run_phase24_thermal_network_localization as base  # noqa: E402


OUT = ROOT / "artifacts/phase24_stycast_only_refinement"
MESH = "mesh_phase24_stycast_only_refined"
TEMPLATE = ROOT / "artifacts/phase24_thermal_network_localization/phase24_membrane_stycast_variant_v2_1p00P.sif"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    base.OUT = OUT
    case = base.MeshCase(
        "phase24_stycast_only_refined",
        "Phase24 Stycast-side-only refined",
        MESH,
        TEMPLATE,
        101,
        1804,
        (
            base.Interface("TES_to_membrane", "TES->membrane", 1104, 1305, 101, 103),
            base.Interface("TES_to_Stycast", "TES->Stycast", 1105, 1204, 101, 102),
            base.Interface("Stycast_to_substrate", "Stycast->substrate", 1205, 1004, 102, 108),
        ),
    )
    manifest = [base.run_case(case, fraction, repeat=False) for fraction in (0.95, 1.0, 1.05)]
    (OUT / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
