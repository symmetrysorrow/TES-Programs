from __future__ import annotations

from pathlib import Path

import pyvista as pv


REPO_ROOT = Path(__file__).resolve().parents[2]
MESH_DIR = REPO_ROOT / "mesh_shifted_merged"
OUT_DIR = REPO_ROOT / "generated"


def pick_mesh_file() -> Path:
    preferred = MESH_DIR / "case_tes_shunt_transient_t0001.vtu"
    if preferred.exists():
        return preferred

    candidates = sorted(MESH_DIR.glob("*_t*.vtu"))
    if candidates:
        return candidates[0]

    raise FileNotFoundError(f"No VTU file found in {MESH_DIR}")


def main() -> Path:
    mesh_path = pick_mesh_file()
    mesh = pv.read(mesh_path)
    surface = mesh.extract_surface().triangulate()

    OUT_DIR.mkdir(exist_ok=True)
    out_path = OUT_DIR / "elmer_mesh_3d_preview.png"

    pv.OFF_SCREEN = True
    plotter = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    plotter.set_background("#fbfaf7")
    plotter.add_mesh(
        surface,
        color="#4e79a7",
        opacity=0.92,
        smooth_shading=True,
        show_edges=True,
        edge_color="#233042",
        line_width=1,
    )
    plotter.add_text(
        f"Elmer mesh 3D preview\n{mesh_path.relative_to(REPO_ROOT)}",
        position="upper_left",
        font_size=16,
        color="black",
    )
    plotter.show_axes()
    plotter.camera_position = "iso"
    plotter.camera.zoom(1.25)
    plotter.show(screenshot=str(out_path))
    plotter.close()
    return out_path


if __name__ == "__main__":
    print(main())
