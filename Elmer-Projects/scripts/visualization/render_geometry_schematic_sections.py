from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_JSON = REPO_ROOT / "elmer_project.json"
OUT_DIR = REPO_ROOT / "generated"

COLORS = {
    "abs": "#c9c9c9",
    "TES": "#d94841",
    "Stycast": "#f2a93b",
    "Membrane_SiNx": "#4e79a7",
    "SiO2_1": "#9fd3c7",
    "Si_1": "#59a14f",
    "SiNx": "#2f5597",
    "Si_2": "#8cd17d",
    "SiO2_2": "#bde0fe",
    "Membrane_Si1": "#2ca02c",
}


def load_params() -> dict:
    return json.loads(PROJECT_JSON.read_text(encoding="utf-8"))["parameters"]


def layer_bounds(params: dict) -> dict[str, tuple[float, float]]:
    si2_top = params["Si_2_dz"] / 2.0
    sio2_1_top = si2_top + params["SiO2_1_dz"]
    si1_top = sio2_1_top + params["Si_1_dz"]
    sinx_top = si1_top + params["SiNx_dz"]
    tes_top = sinx_top + params["TES_dz"]
    stycast_top = tes_top + params["Stycast_dz"]
    abs_top = stycast_top + params["abs_dz"]
    sio2_2_bottom = -params["Si_2_dz"] / 2.0 - params["SiO2_2_dz"]

    return {
        "SiO2_2": (sio2_2_bottom, -params["Si_2_dz"] / 2.0),
        "Si_2": (-params["Si_2_dz"] / 2.0, params["Si_2_dz"] / 2.0),
        "SiO2_1": (si2_top, sio2_1_top),
        "Si_1": (sio2_1_top, si1_top),
        "SiNx": (si1_top, sinx_top),
        "Membrane_Si1": (sio2_1_top, si1_top),
        "Membrane_SiNx": (si1_top, sinx_top),
        "TES": (sinx_top, tes_top),
        "Stycast": (tes_top, stycast_top),
        "abs": (stycast_top, abs_top),
    }


def add_rect(ax, x0, x1, z0, z1, color, label=None, alpha=1.0, lw=1.5):
    ax.add_patch(
        Rectangle(
            (x0 * 1e6, z0 * 1e6),
            (x1 - x0) * 1e6,
            (z1 - z0) * 1e6,
            facecolor=color,
            edgecolor="black",
            linewidth=lw,
            alpha=alpha,
            label=label,
        )
    )


def render_xz(params: dict, bounds: dict) -> Path:
    fig, ax = plt.subplots(figsize=(14, 8))

    half_sub = params["Si_dx"] / 2.0
    half_mem = params["membrane_dx"] / 2.0
    half_tes = params["TES_Au_dx"] / 2.0
    half_sty = params["Stycast_dx"] / 2.0
    half_abs = params["abs_dx"] / 2.0

    add_rect(ax, -half_sub, half_sub, *bounds["SiO2_2"], COLORS["SiO2_2"], "SiO2_2")
    add_rect(ax, -half_sub, half_sub, *bounds["Si_2"], COLORS["Si_2"], "Si_2")
    add_rect(ax, -half_sub, half_sub, *bounds["SiO2_1"], COLORS["SiO2_1"], "SiO2_1")

    add_rect(ax, -half_sub, -half_mem, *bounds["Si_1"], COLORS["Si_1"], "Si_1")
    add_rect(ax, half_mem, half_sub, *bounds["Si_1"], COLORS["Si_1"])
    add_rect(ax, -half_sub, -half_mem, *bounds["SiNx"], COLORS["SiNx"], "SiNx")
    add_rect(ax, half_mem, half_sub, *bounds["SiNx"], COLORS["SiNx"])

    add_rect(ax, -half_mem, half_mem, *bounds["Membrane_Si1"], COLORS["Membrane_Si1"], "Membrane_Si1")
    add_rect(ax, -half_mem, half_mem, *bounds["Membrane_SiNx"], COLORS["Membrane_SiNx"], "Membrane_SiNx")
    add_rect(ax, -half_tes, half_tes, *bounds["TES"], COLORS["TES"], "TES")
    add_rect(ax, -half_sty, half_sty, *bounds["Stycast"], COLORS["Stycast"], "Stycast")
    add_rect(ax, -half_abs, half_abs, *bounds["abs"], COLORS["abs"], "abs", alpha=0.35)

    ax.set_title("Center X-Z Section (schematic from active geometry)")
    ax.set_xlabel("X [um]")
    ax.set_ylabel("Z [um]")
    ax.set_xlim(-1700, 1700)
    ax.set_ylim(bounds["SiO2_2"][0] * 1e6 - 15, bounds["abs"][1] * 1e6 + 25)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.legend(loc="upper right", ncol=2, framealpha=0.95)
    fig.tight_layout()

    out = OUT_DIR / "geometry_section_xz_schematic.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def render_yz(params: dict, bounds: dict) -> Path:
    fig, ax = plt.subplots(figsize=(14, 8))

    y_sub_min = -0.001 - params["Si_dy"] / 2.0
    y_sub_max = -0.001 + params["Si_dy"] / 2.0
    half_mem = params["membrane_dy"] / 2.0
    half_tes = params["TES_Au_dy"] / 2.0
    half_sty = params["Stycast_dy"] / 2.0
    half_abs = params["abs_dy"] / 2.0

    add_rect(ax, y_sub_min, y_sub_max, *bounds["SiO2_2"], COLORS["SiO2_2"], "SiO2_2")
    add_rect(ax, y_sub_min, y_sub_max, *bounds["Si_2"], COLORS["Si_2"], "Si_2")
    add_rect(ax, y_sub_min, y_sub_max, *bounds["SiO2_1"], COLORS["SiO2_1"], "SiO2_1")

    add_rect(ax, y_sub_min, -half_mem, *bounds["Si_1"], COLORS["Si_1"], "Si_1")
    add_rect(ax, half_mem, y_sub_max, *bounds["Si_1"], COLORS["Si_1"])
    add_rect(ax, y_sub_min, -half_mem, *bounds["SiNx"], COLORS["SiNx"], "SiNx")
    add_rect(ax, half_mem, y_sub_max, *bounds["SiNx"], COLORS["SiNx"])

    add_rect(ax, -half_mem, half_mem, *bounds["Membrane_Si1"], COLORS["Membrane_Si1"], "Membrane_Si1")
    add_rect(ax, -half_mem, half_mem, *bounds["Membrane_SiNx"], COLORS["Membrane_SiNx"], "Membrane_SiNx")
    add_rect(ax, -half_tes, half_tes, *bounds["TES"], COLORS["TES"], "TES")
    add_rect(ax, -half_sty, half_sty, *bounds["Stycast"], COLORS["Stycast"], "Stycast")
    add_rect(ax, -half_abs, half_abs, *bounds["abs"], COLORS["abs"], "abs", alpha=0.35)

    ax.set_title("Center Y-Z Section (schematic from active geometry)")
    ax.set_xlabel("Y [um]")
    ax.set_ylabel("Z [um]")
    ax.set_xlim(y_sub_min * 1e6 - 100, y_sub_max * 1e6 + 100)
    ax.set_ylim(bounds["SiO2_2"][0] * 1e6 - 15, bounds["abs"][1] * 1e6 + 25)
    ax.grid(True, alpha=0.25, linestyle=":")
    ax.legend(loc="upper right", ncol=2, framealpha=0.95)
    fig.tight_layout()

    out = OUT_DIR / "geometry_section_yz_schematic.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


if __name__ == "__main__":
    OUT_DIR.mkdir(exist_ok=True)
    params = load_params()
    bounds = layer_bounds(params)
    print(render_xz(params, bounds))
    print(render_yz(params, bounds))
