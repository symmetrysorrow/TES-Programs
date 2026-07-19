from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parent


def copy_many(paths: list[Path], destination: Path) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    for path in paths:
        shutil.copy2(path, destination / path.name)


def write_manifest(run_dir: Path, run_name: str) -> None:
    manifest = {
        "run_name": run_name,
        "source_of_truth": "elmer_project.json",
        "source_files": [
            "elmer_project.json",
        ],
        "generated_files": [
            "generated/tes_shared_variables.sif",
            "generated/tes_case_constant_power.sif",
            "generated/tes_materials.sif",
            "case_constant_power.sif",
        ],
        "mesh_dir": "mesh_shifted_merged",
        "result_files": [
            "mesh_shifted_merged/case_constant_power.ep",
            "mesh_shifted_merged/case_constant_power_t0001.vtu",
        ],
        "expected_values": {
            "tes_nodal_average_k": 0.168581911088883,
            "tes_volume_average_k": 0.168581952265808,
        },
        "rerun_commands": [
            "ElmerSolver case_constant_power.sif",
            "python scripts\\analysis\\extract_tes_volume_avg.py mesh_shifted_merged mesh_shifted_merged\\case_constant_power_t0001.vtu",
        ],
        "notes": [
            "Do not hand-edit files under generated/ inside a frozen run.",
            "This run is intended as a reproducible snapshot of the confirmed 2026-07-11 result.",
        ],
    }
    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def write_readme(run_dir: Path, run_name: str) -> None:
    text = f"""# Repro Run: {run_name}

This directory is a frozen reproducible run bundle.

## Policy

- Source of truth: `elmer_project.json`
- Generated files: `generated/*.sif` and `case_constant_power.sif`
- Mesh and results are frozen under `mesh_shifted_merged/`

## Expected result

- TES nodal average temperature: `0.168581911088883 K`
- TES volume average temperature: `0.168581952265808 K`

## Rerun

```bat
ElmerSolver case_constant_power.sif
python scripts\\analysis\\extract_tes_volume_avg.py mesh_shifted_merged mesh_shifted_merged\\case_constant_power_t0001.vtu
```

## Files

- `manifest.json`: machine-readable run definition
- `generated/`: generated solver inputs copied into the run
- `mesh_shifted_merged/`: frozen mesh and result files
"""
    (run_dir / "README.md").write_text(text, encoding="utf-8")


def make_zip(run_dir: Path) -> Path:
    zip_path = run_dir.with_suffix(".zip")
    if zip_path.exists():
        zip_path.unlink()
    with ZipFile(zip_path, "w", compression=ZIP_DEFLATED) as zf:
        for path in sorted(run_dir.rglob("*")):
            if path.is_file():
                zf.write(path, path.relative_to(run_dir.parent))
    return zip_path


def freeze_run(run_name: str) -> tuple[Path, Path]:
    run_dir = ROOT / "runs" / run_name
    generated_dir = run_dir / "generated"
    mesh_dir = run_dir / "mesh_shifted_merged"

    generated_sources = [
        ROOT / "generated" / "tes_shared_variables.sif",
        ROOT / "generated" / "tes_case_constant_power.sif",
        ROOT / "generated" / "tes_materials.sif",
    ]
    root_sources = [
        ROOT / "case_constant_power.sif",
        ROOT / "elmer_project.json",
        ROOT / "scripts" / "analysis" / "extract_tes_avg_ep.py",
        ROOT / "scripts" / "analysis" / "extract_tes_volume_avg.py",
    ]
    mesh_sources = [
        ROOT / "mesh_shifted_merged" / "mesh.boundary",
        ROOT / "mesh_shifted_merged" / "mesh.elements",
        ROOT / "mesh_shifted_merged" / "mesh.header",
        ROOT / "mesh_shifted_merged" / "mesh.names",
        ROOT / "mesh_shifted_merged" / "mesh.nodes",
        ROOT / "mesh_shifted_merged" / "entities.sif",
        ROOT / "mesh_shifted_merged" / "case_constant_power.ep",
        ROOT / "mesh_shifted_merged" / "case_constant_power_t0001.vtu",
    ]

    if run_dir.exists():
        shutil.rmtree(run_dir)

    copy_many(root_sources, run_dir)
    copy_many(generated_sources, generated_dir)
    copy_many(mesh_sources, mesh_dir)
    write_manifest(run_dir, run_name)
    write_readme(run_dir, run_name)
    zip_path = make_zip(run_dir)
    return run_dir, zip_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Freeze the current reproducible Elmer run.")
    parser.add_argument("run_name", nargs="?", default="current_reference")
    args = parser.parse_args()

    run_dir, zip_path = freeze_run(args.run_name)
    print(run_dir)
    print(zip_path)


if __name__ == "__main__":
    main()
