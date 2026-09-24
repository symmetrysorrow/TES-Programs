"""Materialize the Phase24 native branch/operator audit deliverables.

Reads ``artifacts/phase24_native_branch_operator_audit/raw/*.json`` written by
``analyze_phase24_native_branch_operator_audit.py`` and writes the CSV/JSON
tables and ``summary.md`` / ``next_controlled_fix.md``.
"""
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "artifacts/phase24_native_branch_operator_audit"
RAW = OUT / "raw"
P0 = 3.203004762115138e-10
TBATH = 0.15
TAGS = ("0p95P", "1p00P", "1p05P")
FRACTION = {"0p95P": 0.95, "1p00P": 1.0, "1p05P": 1.05}
CASE_KEYS = ("historical", "best_phase24", "refined_parent")
SHEETS = ("Membrane_SiNx", "Membrane_Si1", "SiO2_1")
SHEET_BRANCH = {
    "Membrane_SiNx": ("A", "Membrane_SiNx window sheet (z191-192um) -> SiNx/Si_1 frame at window perimeter"),
    "Membrane_Si1": ("B", "Membrane_Si1 window sheet (z176-191um) -> Si_1 frame at window perimeter"),
    "SiO2_1": ("C", "SiO2_1 window sheet (z175-176um, back-etched) -> SiO2_1 frame at window perimeter"),
}


def load() -> dict[str, dict[str, Any]]:
    return {k: json.loads((RAW / f"{k}.json").read_text(encoding="utf-8")) for k in CASE_KEYS}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields: list[str] = []
    for row in rows:
        fields.extend(k for k in row if k not in fields)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def dump(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def dT_src(point: dict[str, Any]) -> float:
    return point["traces_K"]["TES_source_weighted"] - TBATH


def g_dissipation(point: dict[str, Any], D: float) -> float:
    return D / dT_src(point) ** 2


def sym(case: dict[str, Any], getter) -> float:
    """Symmetric derivative d(.)/dP from the 0.95/1.05 P0 points."""
    return (getter(case["points"]["1p05P"]) - getter(case["points"]["0p95P"])) / (0.10 * P0)


# --------------------------------------------------------------------------


def branch_graph(R) -> dict[str, Any]:
    out: dict[str, Any] = {"cases": {}}
    for key, case in R.items():
        p = case["points"]["1p00P"]
        out["cases"][key] = {
            "label": case["label"],
            "mesh": case["provenance"]["mesh"],
            "body_adjacency": case["graph"]["edges"],
            "bath": case["graph"]["bath"],
            "evidence": {
                "Q_TES_to_Membrane_over_P": p["Q_TES_to_Membrane_W"] / p["P_source_W"],
                "Q_Stycast_abs_dead_end_W": p["Q_into_Stycast_abs_W"],
                "frame_top_SiNx_T_K": p["traces_K"]["frame_top_SiNx"],
                "window_perimeter_sidewall_T_K": [p["traces_K"]["MembraneSiNx_sidewall_to_SiNx"], p["traces_K"]["MembraneSi1_sidewall_to_Si_1"]],
                "TES_T_K": p["traces_K"]["TES_average_harness"],
                "window_perimeter_to_bath_conductance_W_per_K": case["operators"]["Lambda_window_perimeter_to_bath"]["constant_mode_total_reaction_W_per_K"],
            },
        }
    out["thermal_graph"] = {
        "nodes": ["TES (source, k=68, isothermal within 6uK)", "Stycast+abs (dead end, net Q=0)",
                  "Membrane-side TES trace (z=192um, TES footprint)",
                  "branch A: Membrane_SiNx window sheet", "branch B: Membrane_Si1 window sheet",
                  "branch C: SiO2_1 window sheet (back-etched, bottom free)",
                  "frame (SiNx, Si_1, SiO2_1 frame, Si_2 frame) ~ bath", "SiO2_2", "bath 0.15 K (SiO2_2 zmin)"],
        "edges": [
            ["TES", "Stycast+abs", "mortar (dead end)"],
            ["TES", "Membrane-side trace", "historical: mortar (2218 rows); Phase24: conformal (+18 hanging mortar rows in best)"],
            ["Membrane-side trace", "branch A", "conformal"],
            ["branch A", "branch B", "conformal face z=191 (vertical exchange, ~97-99% of P)"],
            ["branch B", "branch C", "conformal face z=176 (vertical exchange, ~52-53% of P)"],
            ["branch A", "frame", "sidewall z191-192 at s=350um"],
            ["branch B", "frame", "sidewall z176-191 at s=350um"],
            ["branch C", "frame", "continuous SiO2_1 across s=350um"],
            ["frame", "SiO2_2", "series, conductance ~7.5e-4 W/K (negligible resistance)"],
            ["SiO2_2", "bath", "Dirichlet"],
        ],
        "parallel_structure": "A, B, C are parallel lateral sheets between the TES footprint and the window perimeter, strongly coupled vertically (not independent).",
        "correction_to_previous_labels": "The previous 'Si1/SiNx' and 'Si2/SiO2_2' parallel branches do not exist in the actual adjacency: Si_2/SiO2_2 carry heat only in the frame, in series, below the window perimeter; the window region of Si_2 is etched (SiO2_1 window nodes are not shared with Si_2).",
    }
    return out


def native_flux_rows(R) -> list[dict[str, Any]]:
    rows = []
    for key, case in R.items():
        for t in TAGS:
            p = case["points"][t]
            P = p["P_source_W"]
            coupling = "mortar" if key == "historical" else ("conformal + 18 hanging-node mortar rows" if p["Q_TES_to_Membrane_mortar_part_W"] else "conformal")
            base = {"case": key, "P_fraction": FRACTION[t], "P_W": P}
            rows.append({**base, "branch_id": "B0", "body_pair": "TES->Membrane_SiNx", "interface_ids": "TES__zmin / Membrane_SiNx__zmax(under TES)",
                         "coupling_type": coupling, "Q_W": p["Q_TES_to_Membrane_W"], "Q_over_P": p["Q_TES_to_Membrane_W"] / P,
                         "Q_mortar_part_W": p["Q_TES_to_Membrane_mortar_part_W"], "Q_conformal_part_W": p["Q_TES_to_Membrane_conformal_part_W"],
                         "sign": "+ = out of TES", "source": "native assembled edge flow across closed TES/Stycast/abs cut + mortar reaction C^T lambda",
                         "closure_W": p["Q_TES_to_Membrane_W"] + p["Q_into_Stycast_abs_W"] - P})
            rows.append({**base, "branch_id": "B0s", "body_pair": "Stycast+abs->TES", "interface_ids": "TES__zmax / Stycast__zmin (mortar)",
                         "coupling_type": "mortar", "Q_W": p["Q_into_Stycast_abs_W"], "Q_over_P": p["Q_into_Stycast_abs_W"] / P,
                         "sign": "+ = out of Stycast/abs", "source": "native edge flow + mortar reaction; dead end", "closure_W": p["Q_into_Stycast_abs_W"]})
            for sheet in SHEETS:
                b = p["branch_sheets"][sheet]
                bid, desc = SHEET_BRANCH[sheet]
                exch = {k: v for k, v in b.items() if k.startswith("Q_in_from_") and k != "Q_in_from_TES_W"}
                rows.append({**base, "branch_id": bid, "body_pair": desc, "interface_ids": "window perimeter s=350um",
                             "coupling_type": "conformal", "Q_W": b["Q_exit_perimeter_W"], "Q_over_P": b["Q_exit_perimeter_W"] / P,
                             "Q_in_from_TES_W": b["Q_in_from_TES_W"], **exch,
                             "sign": "Q = exit into frame (+ = leaving sheet); Q_in_* = entering sheet",
                             "source": "element reaction r=K_sheet (T-Tb) with element matrices verified against native K",
                             "closure_W": b["balance_residual_W"]})
            rows.append({**base, "branch_id": "M", "body_pair": "frame/SiO2_2->bath", "interface_ids": "SiO2_2__zmin (Dirichlet)",
                         "coupling_type": "Dirichlet", "Q_W": p["Q_bath_W"], "Q_over_P": p["Q_bath_W"] / P,
                         "sign": "+ = into bath", "source": "native edge flow into Dirichlet rows", "closure_W": p["Q_bath_W"] - P})
            sumexit = sum(p["branch_sheets"][s]["Q_exit_perimeter_W"] for s in SHEETS)
            rows.append({**base, "branch_id": "A+B+C", "body_pair": "sum of parallel sheet exits", "Q_W": sumexit, "Q_over_P": sumexit / P,
                         "source": "sum", "closure_W": sumexit - p["Q_TES_to_Membrane_W"]})
            for label, cut in p["lateral_cuts_edge_attribution"].items():
                rows.append({**base, "branch_id": "edge-attribution", "body_pair": label, "Q_W": cut["total"], "Q_over_P": cut["total"] / P,
                             "source": "native edge flow, layer-attributed by edge z", "status": "unavailable as branch Q: obtuse-tet edge circulation",
                             "negative_edge_flow_over_P": cut["negative_edge_flow_sum"] / P})
    return rows


def trace_rows(R) -> list[dict[str, Any]]:
    rows = []
    for key, case in R.items():
        for t in TAGS:
            p = case["points"][t]
            for sheet in SHEETS:
                b = p["branch_sheets"][sheet]
                rows.append({"case": key, "P_fraction": FRACTION[t], "branch_id": SHEET_BRANCH[sheet][0], "sheet": sheet,
                             "hot_trace": "sheet nodes inside TES footprint (mean)", "T_hot_K": b["T_hot_footprint_mean_K"],
                             "cold_trace": "sheet nodes on window perimeter (mean)", "T_cold_K": b["T_cold_perimeter_mean_K"],
                             "dT_K": b["dT_branch_K"], "native_residual_primal_inf": p["native_residual_primal_inf"],
                             "native_balance_rel": p["native_balance_rel"]})
            tr = p["traces_K"]
            rows.append({"case": key, "P_fraction": FRACTION[t], "branch_id": "nodes", "sheet": "",
                         **{f"T_{k}_K": v for k, v in tr.items()},
                         "dT_TES_trace_jump_K": tr["TES_bottom"] - tr["Membrane_top_under_TES"]})
    return rows


def conductance_rows(R) -> tuple[list[dict[str, Any]], dict[str, dict[str, float]]]:
    rows = []
    g: dict[str, dict[str, float]] = {}
    for key, case in R.items():
        g[key] = {}
        p1 = case["points"]["1p00P"]
        for sheet in SHEETS:
            q = lambda p, s=sheet: p["branch_sheets"][s]["Q_exit_perimeter_W"]
            d = lambda p, s=sheet: p["branch_sheets"][s]["dT_branch_K"]
            dq, ddt = sym(case, q), sym(case, d)
            gd = g_dissipation(p1, p1["branch_sheets"][sheet]["dissipation_WK"])
            gq = q(p1) / dT_src(p1)
            g[key][sheet] = gd
            rows.append({"case": key, "branch_id": SHEET_BRANCH[sheet][0], "branch": sheet,
                         "dQ_dP": dq, "d_dT_dP_K_per_W": ddt, "G_own_differential_W_per_K": dq / ddt, "R_own_differential_K_per_W": ddt / dq,
                         "G_exitflow_common_dT_W_per_K": gq, "G_dissipation_W_per_K": gd,
                         "definition": "G_dissipation = sum_e T_e^T K_e T_e / dT_src^2 = dG/dln(k_sheet); sum over all bodies = G_src exactly"})
        for body, D in p1["body_dissipation_WK"].items():
            if body in SHEETS:
                continue
            gd = g_dissipation(p1, D)
            g[key][body] = gd
        # SiO2_1 frame part
        g[key]["SiO2_1_frame"] = g_dissipation(p1, p1["body_dissipation_WK"]["SiO2_1"] - p1["branch_sheets"]["SiO2_1"]["dissipation_WK"])
        g[key]["G_src"] = p1["G_src_W_per_K"]
        g[key]["G_eff"] = case["G_eff_harness_symmetric_W_per_K"]
        g[key]["G_sum_components"] = sum(v for k, v in g[key].items() if k not in ("G_src", "G_eff", "SiO2_1"))  # placeholder, fixed below
    for key in g:
        comps = {k: v for k, v in g[key].items() if k not in ("G_src", "G_eff", "G_sum_components")}
        g[key]["G_sum_components"] = sum(comps.values())
    # comparison rows
    for other in ("best_phase24", "refined_parent"):
        for comp in list(SHEETS) + ["SiO2_1_frame", "TES", "Stycast", "abs", "SiNx", "Si_1", "Si_2", "SiO2_2", "G_src", "G_eff"]:
            h, o = g["historical"].get(comp, float("nan")), g[other].get(comp, float("nan"))
            rows.append({"case": f"{other}_vs_historical", "branch": comp, "G_historical_W_per_K": h, f"G_case_W_per_K": o,
                         "ratio_case_over_historical": o / h if h else float("nan"), "excess_W_per_K": o - h})
    return rows, g


def gap_accounting(R, g) -> dict[str, Any]:
    def split(a: str, b: str) -> dict[str, Any]:
        gap_src = g[b]["G_src"] - g[a]["G_src"]
        gap_eff = g[b]["G_eff"] - g[a]["G_eff"]
        comps = {}
        for comp in ("Membrane_SiNx", "Membrane_Si1", "SiO2_1", "SiO2_1_frame", "TES", "Stycast", "abs", "SiNx", "Si_1", "Si_2", "SiO2_2"):
            d = g[b].get(comp, 0.0) - g[a].get(comp, 0.0)
            comps[comp] = {"dG_W_per_K": d, "fraction_of_gap": d / gap_src}
        explained = sum(v["dG_W_per_K"] for v in comps.values())
        cats = {
            "TES_to_Membrane_operator (TES body + coupling)": comps["TES"]["dG_W_per_K"] + comps["Stycast"]["dG_W_per_K"] + comps["abs"]["dG_W_per_K"],
            "branch_A_Membrane_SiNx_sheet": comps["Membrane_SiNx"]["dG_W_per_K"],
            "branch_B_Membrane_Si1_sheet": comps["Membrane_Si1"]["dG_W_per_K"],
            "branch_C_SiO2_1_window_sheet": comps["SiO2_1"]["dG_W_per_K"],
            "merge_frame_to_bath_network": sum(comps[c]["dG_W_per_K"] for c in ("SiO2_1_frame", "SiNx", "Si_1", "Si_2", "SiO2_2")),
        }
        unexplained = gap_eff - explained
        return {
            "G_reference": g[a]["G_eff"], "G_case": g[b]["G_eff"], "ratio": g[b]["G_eff"] / g[a]["G_eff"],
            "gap_G_eff_W_per_K": gap_eff, "gap_G_src_W_per_K": gap_src,
            "components": comps,
            "categories_W_per_K": cats,
            "categories_fraction_of_gap_eff": {k: v / gap_eff for k, v in cats.items()},
            "explained_gap_W_per_K": explained,
            "explained_fraction": explained / gap_eff,
            "unexplained_W_per_K": unexplained,
            "unexplained_fraction": unexplained / gap_eff,
            "unexplained_breakdown": {
                "G_eff(TES element-node mean) vs G_src(source-weighted) definition": (g[b]["G_eff"] - g[b]["G_src"]) - (g[a]["G_eff"] - g[a]["G_src"]),
                "dissipation closure (MUMPS balance precision)": g[b]["G_src"] - g[b]["G_sum_components"] - (g[a]["G_src"] - g[a]["G_sum_components"]),
            },
        }

    return {
        "method": ("Exact additive decomposition G_src = sum_body D_body/dT_src^2 with D_body = sum_e T_e^T K_e T_e "
                   "(element matrices verified against the native K; D_body/dT^2 = dG/dln k_body, the first-order sensitivity "
                   "of G to that body's discrete operator). System is linear (A identical at 0.95/1.00/1.05 P0), so differential = secant."),
        "best_phase24_minus_historical": split("historical", "best_phase24"),
        "refined_parent_minus_best_phase24_controlled_improvement": split("best_phase24", "refined_parent"),
        "refined_parent_minus_historical": split("historical", "refined_parent"),
    }


def operator_tables(R) -> tuple[dict, list, dict, list]:
    tes_json: dict[str, Any] = {}
    tes_rows: list[dict[str, Any]] = []
    net_json: dict[str, Any] = {}
    net_rows: list[dict[str, Any]] = []
    ref = json.loads((RAW / "annulus_reference.json").read_text(encoding="utf-8"))
    s_ref_ntd = ref["ntd_uniform_footprint"]["h1um"]
    s_ref_dir = ref["dirichlet"]["h1um"]
    for key, case in R.items():
        ops = case["operators"]
        tes_json[key] = {
            "TES_Membrane_local_operator(membrane side, Membrane_SiNx layer, grounded below)": ops["TES_Membrane_local_operator"],
            "TES_side_native_NtD(full native saddle A, TES source modes)": ops["TES_side_native_NtD"],
            "mortar_slave_Dirichlet_wellposedness": ops.get("TES_Membrane_mortar_slave_Dirichlet_wellposedness", "no TES/Membrane mortar rows beyond hanging nodes" if key != "historical" else None),
            "TES_dissipation_G_W_per_K": g_dissipation(case["points"]["1p00P"], case["points"]["1p00P"]["body_dissipation_WK"]["TES"]),
            "TES_bottom_minus_membrane_top_trace_K_at_P0": case["points"]["1p00P"]["traces_K"]["TES_bottom"] - case["points"]["1p00P"]["traces_K"]["Membrane_top_under_TES"],
        }
        net_json[key] = {
            "Lambda_membrane_trace_to_bath(network DtN)": ops["Lambda_membrane_trace_to_bath"],
            "sheet_operators(parallel branches)": {
                s: {**v, "NtD_shape_factor_reference_2D": s_ref_ntd, "NtD_relative_error": v["NtD_shape_factor_mesh"] / s_ref_ntd - 1.0,
                    "Dirichlet_shape_factor_reference_2D": s_ref_dir, "Dirichlet_relative_error": v["shape_factor_mesh"] / s_ref_dir - 1.0}
                for s, v in ops["sheet_annulus_operators"].items()},
            "Lambda_window_perimeter_to_bath": ops["Lambda_window_perimeter_to_bath"],
        }
        for row in case["mode_rows"]:
            target = tes_rows if row["operator"] == "TES_Membrane_local_operator" else net_rows
            target.append({"case": key, **row})
        for op_name, target in (("TES_Membrane_local_operator", tes_rows), ("Lambda_membrane_trace_to_bath", net_rows)):
            o = ops[op_name]
            for mode in ("constant", "linear_x", "linear_y", "radial_quadratic"):
                m = o[f"mode_{mode}"]
                target.append({"case": key, "operator": op_name, "mode": mode, "s_bin_um": "TOTAL",
                               "trace_nodes": o["trace_dof_count"], "reaction_W_per_K": m["total_reaction_W_per_K"],
                               "energy_W_per_K": m["energy_W_per_K"], "normalized_response_W_per_K_m2": m["normalized_response_W_per_K_m2"]})
        ntd = ops["TES_side_native_NtD"]
        for mode in ("constant", "linear_x", "linear_y", "radial_quadratic"):
            m = ntd[f"mode_{mode}"]
            tes_rows.append({"case": key, "operator": "TES_side_native_NtD", "mode": mode, "s_bin_um": "TOTAL",
                             "energy_K_per_W": m["energy_K_per_W"], "TES_mean_response_K_per_W": m["TES_mean_temperature_response_K_per_W"],
                             "net_source_W": m["net_source_W"]})
    # ratios
    for name, js, key_op in (("tes", tes_json, None), ("net", net_json, None)):
        pass
    comp = {}
    for mode in ("constant", "linear_x", "linear_y", "radial_quadratic"):
        h = R["historical"]["operators"]
        b = R["best_phase24"]["operators"]
        comp[mode] = {
            "Lambda_membrane_energy_ratio_best_over_hist": b["Lambda_membrane_trace_to_bath"][f"mode_{mode}"]["energy_W_per_K"] / h["Lambda_membrane_trace_to_bath"][f"mode_{mode}"]["energy_W_per_K"],
            "Lambda_membrane_normalized_ratio_best_over_hist": b["Lambda_membrane_trace_to_bath"][f"mode_{mode}"]["normalized_response_W_per_K_m2"] / h["Lambda_membrane_trace_to_bath"][f"mode_{mode}"]["normalized_response_W_per_K_m2"],
            "TES_Membrane_local_normalized_ratio_best_over_hist": b["TES_Membrane_local_operator"][f"mode_{mode}"]["normalized_response_W_per_K_m2"] / h["TES_Membrane_local_operator"][f"mode_{mode}"]["normalized_response_W_per_K_m2"],
            "TES_side_NtD_energy_ratio_hist_over_best": h["TES_side_native_NtD"][f"mode_{mode}"]["energy_K_per_W"] / b["TES_side_native_NtD"][f"mode_{mode}"]["energy_K_per_W"],
        }
    tes_json["comparison_best_over_historical"] = {m: {k: v for k, v in c.items() if "TES" in k} for m, c in comp.items()}
    net_json["comparison_best_over_historical"] = {m: {k: v for k, v in c.items() if "Lambda" in k} for m, c in comp.items()}
    net_json["reference_2D_annulus"] = ref
    return tes_json, tes_rows, net_json, net_rows


SAVESCALARS = {
    "historical": ROOT / "artifacts/phase24_thermal_network_localization/native_scalars/historical_{t}.dat",
    "best_phase24": ROOT / "artifacts/p24trace/ms/native_scalars/ms_hist_{t}.dat",
    "refined_parent": ROOT / "artifacts/p24d10b/native_scalars/phase24_historical_density_{t}.dat",
}


def savescalars_bath(key: str, t: str) -> float:
    path = Path(str(SAVESCALARS[key]).format(t=t))
    if not path.is_file():
        return float("nan")
    last = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()][-1]
    return -float(last[0])


def operator_decomposition(R) -> dict[str, Any]:
    """Exact multiplicative split: coupled G = G_uniform_flux_stack x edge-loading factor."""
    out: dict[str, Any] = {}
    for key, case in R.items():
        gu = case["operators"]["Stack_NtD_uniform_footprint_flux"]["G_uniform_flux_W_per_K"]
        ge = case["G_eff_harness_symmetric_W_per_K"]
        out[key] = {"G_eff": ge, "G_uniform_flux_stack": gu, "edge_loading_factor": ge / gu,
                    "uniform_flux_mode_energies_K_per_W": {m: case["operators"]["Stack_NtD_uniform_footprint_flux"][f"mode_{m}"]["energy_K_per_W"]
                                                            for m in ("constant", "linear_x", "linear_y", "radial_quadratic")}}
    for a, bkey in (("historical", "best_phase24"), ("best_phase24", "refined_parent"), ("historical", "refined_parent")):
        ratio = out[bkey]["G_eff"] / out[a]["G_eff"]
        bulk = out[bkey]["G_uniform_flux_stack"] / out[a]["G_uniform_flux_stack"]
        edge = out[bkey]["edge_loading_factor"] / out[a]["edge_loading_factor"]
        out[f"{bkey}_over_{a}"] = {
            "G_ratio": ratio, "bulk_stack_operator_ratio(uniform flux)": bulk, "TES_edge_loading_ratio": edge,
            "log_share_bulk": math.log(bulk) / math.log(ratio), "log_share_edge": math.log(edge) / math.log(ratio),
            "mode_energy_ratio_uniform_flux(a_over_b)": {m: out[a]["uniform_flux_mode_energies_K_per_W"][m] / out[bkey]["uniform_flux_mode_energies_K_per_W"][m]
                                                         for m in ("constant", "linear_x", "linear_y", "radial_quadratic")},
        }
    return out


def zone_rows(R) -> list[dict[str, Any]]:
    rows = []
    P = {k: R[k]["points"]["1p00P"]["P_source_W"] for k in R}
    zh = R["historical"]["points"]["1p00P"]["zone_dissipation_WK"]
    for other in ("best_phase24", "refined_parent"):
        zo = R[other]["points"]["1p00P"]["zone_dissipation_WK"]
        dR_tot = sum(zo.values()) / P[other] ** 2 - sum(zh.values()) / P["historical"] ** 2
        for z in zh:
            rh, ro = zh[z] / P["historical"] ** 2, zo[z] / P[other] ** 2
            if max(abs(rh), abs(ro)) < 1e3:
                continue
            body, band = z.split("|")
            rows.append({"comparison": f"{other}_minus_historical", "body": body, "s_band_um": band,
                         "R_historical_K_per_W": rh, "R_case_K_per_W": ro, "dR_K_per_W": ro - rh,
                         "fraction_of_dR_total": (ro - rh) / dR_tot})
    return rows


def fmt(x: float, d: int = 4) -> str:
    return f"{x:.{d}e}"


def main() -> int:
    R = load()
    dump(OUT / "actual_branch_graph.json", branch_graph(R))
    write_csv(OUT / "native_branch_flux.csv", native_flux_rows(R))
    write_csv(OUT / "branch_trace_temperatures.csv", trace_rows(R))
    crow, g = conductance_rows(R)
    write_csv(OUT / "branch_differential_conductance.csv", crow)
    acc = gap_accounting(R, g)
    dump(OUT / "remaining_gap_accounting.json", acc)
    tes_json, tes_rows, net_json, net_rows = operator_tables(R)
    dump(OUT / "tes_membrane_condensed_operator.json", tes_json)
    write_csv(OUT / "tes_membrane_mode_response.csv", tes_rows)
    dump(OUT / "membrane_substrate_operator.json", net_json)
    write_csv(OUT / "membrane_substrate_mode_response.csv", net_rows)
    opdec = operator_decomposition(R)
    dump(OUT / "operator_gap_decomposition.json", opdec)
    write_csv(OUT / "zone_resistance_accounting.csv", zone_rows(R))
    dump(OUT / "tes_edge_band_mesh_statistics.json", {k: R[k]["mesh_stats"] for k in R})
    xc = []
    for key in R:
        for t in TAGS:
            p = R[key]["points"][t]
            xc.append({"case": key, "P_fraction": FRACTION[t], "Q_bath_reaction_W": p["Q_bath_W"],
                       "SaveScalars_bath_diffusive_flux_W": savescalars_bath(key, t),
                       "relative_difference": savescalars_bath(key, t) / p["Q_bath_W"] - 1.0,
                       "note": "SaveScalars flux is element-gradient based (not a reaction); agreement is a cross-check only"})
    write_csv(OUT / "bath_flux_crosscheck.csv", xc)
    print(json.dumps({k: v for k, v in opdec.items() if "_over_" in k}, indent=1))
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dump(OUT / "provenance.json", {"audit_parent_commit": head, "cases": {k: R[k]["provenance"] for k in R},
                                  "element_reassembly_verification": {k: R[k]["element_reassembly_verification"] for k in R},
                                  "solver_reruns": "none (existing opt-in native captures reused)",
                                  "native_source_changes": "none"})
    print(json.dumps({k: acc["best_phase24_minus_historical"][k] for k in ("ratio", "explained_fraction", "unexplained_fraction")}, indent=1))
    print(json.dumps(acc["best_phase24_minus_historical"]["categories_fraction_of_gap_eff"], indent=1))
    print(json.dumps(acc["refined_parent_minus_best_phase24_controlled_improvement"]["categories_fraction_of_gap_eff"], indent=1))
    print(json.dumps(tes_json["comparison_best_over_historical"], indent=1))
    print(json.dumps(net_json["comparison_best_over_historical"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
