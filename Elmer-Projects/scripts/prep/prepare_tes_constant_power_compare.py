from __future__ import annotations

import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TES_JSON = ROOT / "reference" / "tes.json"
ELMER_PROJECT_JSON = ROOT / "elmer_project.json"
GENERATED_DIR = ROOT / "generated"


def _apply_y_shift_only(base: dict, source: dict) -> None:
    source_children = {child["name"]: child for child in source["geometry"]["children"]}
    for child in base["geometry"]["children"]:
        name = child.get("name")
        if name not in {"SiO2_1", "Si_1", "SiNx", "Si_2", "SiO2_2"}:
            continue
        source_group = source_children[name]
        for base_part, source_part in zip(child["children"], source_group["children"]):
            if source_part.get("body_mode") != "add":
                continue
            base_part["y"] = source_part["y"]
            base_part["y_expr"] = source_part["y_expr"]


def main() -> None:
    variant = sys.argv[1] if len(sys.argv) > 1 else "base"
    base = json.loads(TES_JSON.read_text(encoding="utf-8"))
    source = json.loads(ELMER_PROJECT_JSON.read_text(encoding="utf-8"))

    base["elmer_cases"] = source.get("elmer_cases", {})
    stycast = next(
        child for child in base["geometry"]["children"] if child.get("name") == "Stycast"
    )
    base["elmer_overrides"] = {
        "stycast_shape": stycast["shape"],
        "stycast_x": stycast["x"],
        "stycast_y": stycast["y"],
        "stycast_z": stycast["z"],
        "stycast_dx": stycast["dx"],
        "stycast_dy": stycast["dy"],
        "stycast_dz": stycast["dz"],
        "fragment_mortar_interfaces": False,
    }

    if variant == "yshift":
        _apply_y_shift_only(base, source)
        out_json = GENERATED_DIR / "tes_constant_power_compare_yshift.json"
    elif variant == "fragment":
        base["elmer_overrides"]["fragment_mortar_interfaces"] = True
        out_json = GENERATED_DIR / "tes_constant_power_compare_fragment.json"
    else:
        out_json = GENERATED_DIR / "tes_constant_power_compare.json"

    GENERATED_DIR.mkdir(exist_ok=True)
    out_json.write_text(json.dumps(base, indent=2), encoding="utf-8")
    print(str(out_json))


if __name__ == "__main__":
    main()
