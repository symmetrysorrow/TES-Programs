"""Independent cross-check of the dual-TES series CSVs against the VTU field
data, plus a coarse pulse energy-balance test (docs/dual_tes_plan.md gate
D-3.4 amplitude verification).

For each requested snapshot it recomputes the TES_L (body 101) and TES_R
(body 110) temperatures directly from the VTU temperature field -- both the
nodal arithmetic mean (which is what the UDF's own AverageTemperature is, so
this is the apples-to-apples check against the series CSV) and the
volume-weighted average -- and compares to the series CSV value at the same
solver time. Snapshot -> solver-time mapping is recovered from the solver
log (each 'Saving in unstructured VTK XML' line is attributed to the most
recent 'MAIN: Time: step/total: time' line), because the Elmer VTU files
carry no embedded TIME field and the staged Output Intervals make the
snapshot index -> step map non-obvious.

The energy check integrates rho*cp*(T_post - T_pre) over the whole mesh
(materials + body->material map read from the case SIF) and compares to the
1332 keV pulse energy.

Read-only. Usage: python scripts/analysis/dual_vtu_verify.py
"""
from __future__ import annotations

import re
from pathlib import Path

import meshio
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MESH = ROOT / "mesh_dual_base"

BODY = {"L": 101, "R": 110}
T_STEADY = {"L": 0.1672581, "R": 0.1673430}
KEV_1332 = 1332 * 1.602176634e-16  # J

# Material id -> (density kg/m^3, heat capacity J/(kg K)); target body id ->
# material id. Both are read straight from the case SIF at runtime below so
# they cannot drift from what actually ran; these constants are only the
# fallback / documentation of the expected values.


def load_nodes() -> np.ndarray:
    ids, xyz = [], []
    with (MESH / "mesh.nodes").open() as f:
        for line in f:
            p = line.split()
            if not p:
                continue
            ids.append(int(p[0]))
            xyz.append((float(p[2]), float(p[3]), float(p[4])))
    return np.array(xyz)[np.argsort(ids)]  # row (node_id-1) -> xyz


def load_elements() -> tuple[np.ndarray, np.ndarray]:
    body, conn = [], []
    with (MESH / "mesh.elements").open() as f:
        for line in f:
            p = line.split()
            if not p or int(p[2]) != 504:
                continue
            body.append(int(p[1]))
            conn.append([int(x) for x in p[3:]])
    return np.array(body), np.array(conn, dtype=int)


def tet_volumes(conn: np.ndarray, nodes: np.ndarray) -> np.ndarray:
    p0, p1, p2, p3 = (nodes[conn[:, i] - 1] for i in range(4))
    return np.abs(np.einsum("ij,ij->i", np.cross(p1 - p0, p2 - p0), p3 - p0)) / 6.0


def snapshot_times(case: str) -> list[tuple[int, float]]:
    time_re = re.compile(r"MAIN: Time:\s+(\d+)/\d+:\s+([\dExe.+-]+)")
    save_re = re.compile(r"Saving in unstructured VTK XML")
    cur = (None, None)
    out = []
    for ln in (ROOT / "results" / case / "solver.log").read_text(errors="replace").splitlines():
        m = time_re.search(ln)
        if m:
            cur = (int(m.group(1)), float(m.group(2)))
        elif save_re.search(ln):
            out.append(cur)
    return out


def read_temperature(case: str, snap: int) -> np.ndarray:
    m = meshio.read(ROOT / "results" / case / f"{case}_t{snap:04d}.vtu")
    return np.asarray(m.point_data["temperature"]).reshape(-1)


def load_csv(case: str, base: str, side: str) -> pd.DataFrame:
    for parent in (ROOT / "results" / case, ROOT):
        for tag in (side, side.lower()):
            p = parent / f"{base}_{tag}_series.csv"
            if p.exists():
                return pd.read_csv(p, names=["time_s", "T", "I", "R", "P"], header=0)
    raise FileNotFoundError(f"series CSV for {case} side {side}")


def parse_materials(sif: Path) -> tuple[dict[int, tuple[float, float]], dict[int, int]]:
    text = sif.read_text(encoding="utf-8")
    mats: dict[int, tuple[float, float]] = {}
    for m in re.finditer(r"Material (\d+)\n(.*?)\nEnd", text, re.S):
        mid = int(m.group(1))
        blk = m.group(2)
        rho = float(re.search(r"Density = ([\dEe.+-]+)", blk).group(1))
        cp = float(re.search(r"Heat Capacity = ([\dEe.+-]+)", blk).group(1))
        mats[mid] = (rho, cp)
    body_mat: dict[int, int] = {}
    for m in re.finditer(r"Body \d+\n.*?Target Bodies\(1\) = (\d+).*?Material = (\d+)", text, re.S):
        body_mat[int(m.group(1))] = int(m.group(2))
    return mats, body_mat


def main() -> None:
    nodes = load_nodes()
    ebody, econn = load_elements()
    evol = tet_volumes(econn, nodes)
    tet = {s: econn[ebody == BODY[s]] for s in BODY}
    vol = {s: evol[ebody == BODY[s]] for s in BODY}

    # Coordinate alignment check (VTU point i == mesh node i+1).
    t1 = meshio.read(ROOT / "results" / "case_dual_pulse_offset" / "case_dual_pulse_offset_t0001.vtu")
    print(f"VTU<->mesh.nodes max coord mismatch: {np.abs(t1.points - nodes).max():.2e} m (expect ~0)")
    print(f"TES_L mesh volume = {vol['L'].sum():.4e} m^3 ({len(tet['L'])} tets)")
    print(f"TES_R mesh volume = {vol['R'].sum():.4e} m^3 ({len(tet['R'])} tets)")
    print()

    def volavg(temp, s):
        return float((temp[tet[s] - 1].mean(axis=1) * vol[s]).sum() / vol[s].sum())

    def nodalmean(temp, s):
        return float(temp[np.unique(tet[s].reshape(-1)) - 1].mean())

    def verify(case, base, snaps, labels):
        print("=" * 84)
        print(f"{case}: TES temperature, VTU field vs series CSV")
        print("=" * 84)
        sv = snapshot_times(case)
        csv = {s: load_csv(case, base, s) for s in BODY}
        rows = []
        for snap, lab in zip(snaps, labels):
            step, vt = sv[snap - 1]
            temp = read_temperature(case, snap)
            for s in BODY:
                va, nm = volavg(temp, s), nodalmean(temp, s)
                ci = (csv[s]["time_s"] - vt).abs().idxmin()
                cT = float(csv[s].loc[ci, "T"])
                rows.append({
                    "region": lab, "snap": f"t{snap:04d}", "t_ms": round(vt * 1e3, 4),
                    "side": f"TES_{s}",
                    "dT_vtu_vol_uK": round((va - T_STEADY[s]) * 1e6, 2),
                    "dT_vtu_nodal_uK": round((nm - T_STEADY[s]) * 1e6, 2),
                    "dT_csv_uK": round((cT - T_STEADY[s]) * 1e6, 2),
                    "nodal-csv_uK": round((nm - cT) * 1e6, 2),
                })
        df = pd.DataFrame(rows)
        pd.set_option("display.width", 200)
        print(df.to_string(index=False))
        print()
        return df

    off_snaps = [2, 15, 16, 17, 18, 21]
    off_labels = ["pre-pulse", "R-peak~20.1ms", "L-rise~20.5ms", "L-peak~20.9ms", "plateau~32ms", "late~92ms"]
    verify("case_dual_pulse_offset", "tes_dual_pulse_offset", off_snaps, off_labels)
    verify("case_dual_pulse_center", "tes_dual_pulse_center", [2, 11, 17, 21],
           ["pre-pulse", "post-pulse", "peak~20.9ms", "late~92ms"])

    # Energy balance (offset case).
    print("=" * 84)
    print("Energy balance: integral rho*cp*(T - T_pre) over whole mesh vs 1332 keV")
    print("=" * 84)
    mats, body_mat = parse_materials(ROOT / "generated" / "cases" / "case_dual_pulse_offset.sif")
    rho_cp = np.array([mats[body_mat[b]][0] * mats[body_mat[b]][1] for b in ebody])
    T_pre = read_temperature("case_dual_pulse_offset", 2)  # t=0.02, just before pulse
    sv = snapshot_times("case_dual_pulse_offset")
    print(f"E_pulse (1332 keV) = {KEV_1332:.4e} J")
    for snap in (6, 7, 11, 15, 17, 18, 21):
        step, vt = sv[snap - 1]
        T = read_temperature("case_dual_pulse_offset", snap)
        dT = (T[econn - 1] - T_pre[econn - 1]).mean(axis=1)
        E = rho_cp * evol * dT
        E_abs = E[ebody == 100].sum()
        print(f"  t{snap:04d} t={vt * 1e3:8.4f}ms: E_total={E.sum():.4e} J "
              f"(abs={E_abs:.4e})  E/E_pulse={E.sum() / KEV_1332:.3f}")


if __name__ == "__main__":
    main()
