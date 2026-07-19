from __future__ import annotations

from pathlib import Path

import meshio
import pyvista as pv
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
MESH_PATH = REPO_ROOT / "gmsh" / "project.msh"
OUT_DIR = REPO_ROOT / "generated"

BODY_INFO = {
    100: {"name": "abs", "color": "#c9c9c9", "opacity": 0.14},
    101: {"name": "TES", "color": "#d94841", "opacity": 1.0},
    102: {"name": "Stycast", "color": "#f2a93b", "opacity": 0.95},
    103: {"name": "Membrane_SiNx", "color": "#4e79a7", "opacity": 0.95},
    104: {"name": "SiO2_1", "color": "#9fd3c7", "opacity": 0.20},
    105: {"name": "Si_1", "color": "#59a14f", "opacity": 0.20},
    106: {"name": "SiNx", "color": "#2f5597", "opacity": 0.28},
    107: {"name": "Si_2", "color": "#8cd17d", "opacity": 0.10},
    108: {"name": "SiO2_2", "color": "#bde0fe", "opacity": 0.08},
    109: {"name": "Membrane_Si1", "color": "#2ca02c", "opacity": 0.95},
}
SECTION_Z_SCALE = 5000.0


def tetra_grid(points: np.ndarray, tetra_cells: np.ndarray) -> pv.UnstructuredGrid:
    n = tetra_cells.shape[0]
    cells = np.hstack([np.full((n, 1), 4, dtype=np.int64), tetra_cells.astype(np.int64)]).ravel()
    celltypes = np.full(n, pv.CellType.TETRA, dtype=np.uint8)
    return pv.UnstructuredGrid(cells, celltypes, points)


def load_body_surfaces() -> dict[int, pv.PolyData]:
    mesh = meshio.read(MESH_PATH)
    tetra = mesh.cells_dict["tetra"]
    physical = np.asarray(mesh.cell_data_dict["gmsh:physical"]["tetra"])
    surfaces: dict[int, pv.PolyData] = {}

    for body_id in sorted(np.unique(physical)):
        mask = physical == body_id
        body_grid = tetra_grid(mesh.points, tetra[mask])
        surfaces[int(body_id)] = body_grid.extract_surface().triangulate()

    return surfaces


def add_bodies(plotter: pv.Plotter, surfaces: dict[int, pv.PolyData], zoom: bool) -> None:
    for body_id, surface in surfaces.items():
        info = BODY_INFO.get(body_id, {"name": str(body_id), "color": "white", "opacity": 1.0})
        opacity = info["opacity"]
        if zoom and body_id in {100, 107, 108}:
            opacity *= 0.35
        plotter.add_mesh(
            surface,
            color=info["color"],
            opacity=opacity,
            smooth_shading=True,
            show_edges=False,
            name=info["name"],
            label=info["name"],
        )


def render_views() -> tuple[Path, Path]:
    pv.OFF_SCREEN = True
    OUT_DIR.mkdir(exist_ok=True)
    surfaces = load_body_surfaces()

    full_path = OUT_DIR / "geometry_preview_full.png"
    zoom_path = OUT_DIR / "geometry_preview_zoom.png"

    full = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    full.set_background("#fbfaf7")
    add_bodies(full, surfaces, zoom=False)
    full.add_legend(face="circle", bcolor="white", border=True, size=(0.2, 0.35))
    full.camera_position = [
        (0.0038, -0.0048, 0.0030),
        (0.0, 0.0, 0.00018),
        (0.0, 0.0, 1.0),
    ]
    full.show(screenshot=str(full_path))
    full.close()

    zoom = pv.Plotter(off_screen=True, window_size=(1800, 1200))
    zoom.set_background("#fbfaf7")
    add_bodies(zoom, surfaces, zoom=True)
    zoom.add_legend(face="circle", bcolor="white", border=True, size=(0.2, 0.35))
    zoom.camera_position = [
        (0.00095, -0.00105, 0.00085),
        (0.0, 0.0, 0.000365),
        (0.0, 0.0, 1.0),
    ]
    zoom.camera.zoom(1.35)
    zoom.show(screenshot=str(zoom_path))
    zoom.close()

    return full_path, zoom_path


def render_sections() -> tuple[Path, Path]:
    pv.OFF_SCREEN = True
    OUT_DIR.mkdir(exist_ok=True)
    surfaces = load_body_surfaces()

    xz_path = OUT_DIR / "geometry_section_xz.png"
    yz_path = OUT_DIR / "geometry_section_yz.png"
    xz_scaled_path = OUT_DIR / "geometry_section_xz_scaled.png"
    yz_scaled_path = OUT_DIR / "geometry_section_yz_scaled.png"

    def setup_section_plot() -> pv.Plotter:
        plotter = pv.Plotter(off_screen=True, window_size=(1800, 1200))
        plotter.set_background("#fbfaf7")
        return plotter

    def add_section(plotter: pv.Plotter, normal: tuple[float, float, float], origin: tuple[float, float, float]) -> None:
        for body_id, surface in surfaces.items():
            info = BODY_INFO.get(body_id, {"name": str(body_id), "color": "white"})
            section = surface.slice(normal=normal, origin=origin)
            if section.n_points == 0:
                continue
            width = 8 if body_id in {101, 102, 103, 109} else 5
            plotter.add_mesh(
                section,
                color=info["color"],
                line_width=width,
                name=info["name"],
                label=info["name"],
            )

    def add_scaled_section(
        plotter: pv.Plotter,
        normal: tuple[float, float, float],
        origin: tuple[float, float, float],
    ) -> None:
        for body_id, surface in surfaces.items():
            info = BODY_INFO.get(body_id, {"name": str(body_id), "color": "white"})
            section = surface.slice(normal=normal, origin=origin)
            if section.n_points == 0:
                continue
            scaled = section.copy()
            pts = scaled.points.copy()
            pts[:, 2] *= SECTION_Z_SCALE
            scaled.points = pts
            width = 8 if body_id in {101, 102, 103, 109} else 5
            plotter.add_mesh(
                scaled,
                color=info["color"],
                line_width=width,
                name=info["name"],
                label=info["name"],
            )

    xz = setup_section_plot()
    add_section(xz, normal=(0.0, 1.0, 0.0), origin=(0.0, 0.0, 0.0))
    xz.add_text("Center X-Z Section (y = 0)", position="upper_left", font_size=16, color="black")
    xz.add_legend(face="circle", bcolor="white", border=True, size=(0.2, 0.35))
    xz.camera_position = "xy"
    xz.camera.zoom(1.6)
    xz.show_bounds(
        xtitle="X [m]",
        ytitle="Z [m]",
        show_zaxis=False,
        location="outer",
        font_size=14,
    )
    xz.show(screenshot=str(xz_path))
    xz.close()

    yz = setup_section_plot()
    add_section(yz, normal=(1.0, 0.0, 0.0), origin=(0.0, 0.0, 0.0))
    yz.add_text("Center Y-Z Section (x = 0)", position="upper_left", font_size=16, color="black")
    yz.add_legend(face="circle", bcolor="white", border=True, size=(0.2, 0.35))
    yz.camera_position = "yz"
    yz.camera.zoom(1.6)
    yz.show_bounds(
        xtitle="Y [m]",
        ytitle="Z [m]",
        show_zaxis=False,
        location="outer",
        font_size=14,
    )
    yz.show(screenshot=str(yz_path))
    yz.close()

    xz_scaled = setup_section_plot()
    add_scaled_section(xz_scaled, normal=(0.0, 1.0, 0.0), origin=(0.0, 0.0, 0.0))
    xz_scaled.add_text(
        f"Center X-Z Section (y = 0, z x{int(SECTION_Z_SCALE)})",
        position="upper_left",
        font_size=16,
        color="black",
    )
    xz_scaled.add_legend(face="circle", bcolor="white", border=True, size=(0.2, 0.35))
    xz_scaled.camera_position = "xy"
    xz_scaled.camera.zoom(1.2)
    xz_scaled.show_bounds(
        xtitle="X [m]",
        ytitle=f"Z x{int(SECTION_Z_SCALE)}",
        show_zaxis=False,
        location="outer",
        font_size=14,
    )
    xz_scaled.show(screenshot=str(xz_scaled_path))
    xz_scaled.close()

    yz_scaled = setup_section_plot()
    add_scaled_section(yz_scaled, normal=(1.0, 0.0, 0.0), origin=(0.0, 0.0, 0.0))
    yz_scaled.add_text(
        f"Center Y-Z Section (x = 0, z x{int(SECTION_Z_SCALE)})",
        position="upper_left",
        font_size=16,
        color="black",
    )
    yz_scaled.add_legend(face="circle", bcolor="white", border=True, size=(0.2, 0.35))
    yz_scaled.camera_position = "yz"
    yz_scaled.camera.zoom(1.2)
    yz_scaled.show_bounds(
        xtitle="Y [m]",
        ytitle=f"Z x{int(SECTION_Z_SCALE)}",
        show_zaxis=False,
        location="outer",
        font_size=14,
    )
    yz_scaled.show(screenshot=str(yz_scaled_path))
    yz_scaled.close()

    return xz_path, yz_path, xz_scaled_path, yz_scaled_path


if __name__ == "__main__":
    full_path, zoom_path = render_views()
    xz_path, yz_path, xz_scaled_path, yz_scaled_path = render_sections()
    print(full_path)
    print(zoom_path)
    print(xz_path)
    print(yz_path)
    print(xz_scaled_path)
    print(yz_scaled_path)
