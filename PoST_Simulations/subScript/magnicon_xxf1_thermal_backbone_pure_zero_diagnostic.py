"""Thermal-backbone profile with c2=0, capped-c2, and pure-zero branches.

The existing separated-substrate diagnostic held several thermal-backbone
quantities at an inherited snapshot.  This diagnostic relaxes the specific
assumptions identified as weakly constrained:

* one shared effective TES-to-bath DC conductance scale is profiled broadly;
* T_bath is allowed to move modestly below the nominal 215 mK run label;
* Pb absorber thickness is shared and profiled over the measured
  0.5 +/- 0.1 mm range;
* C_abs and G_abs-abs are recomputed from that same thickness using the
  current Elmer Pb rho, cp, and k, rather than inherited fit values.

The separated local substrate node, L=0.1--1.3 nH profile, documented
Magnicon 10 kHz Bessel response, and SRS SIM965 response remain in place.

Three readout branches are compared:
  1. c2 = 0;
  2. historical capped c2 reference;
  3. a one-parameter pure-zero response sqrt(1 + (f/fz)^2).

The pure-zero branch is an interpretable reparameterization of the c2
numerator with no artificial c2 upper-bound censoring over its configured
zero-frequency interval.  It is not, by itself, a component identification.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
WORK_DIR = ROOT / ".noise_optimization_work_rsh_sweep"

DEFAULT_CONFIG = (
    CONFIG_DIR
    / "magnicon_xxf1_thermal_backbone_pure_zero_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR
    / "magnicon_xxf1_thermal_backbone_pure_zero_diagnostic.json"
)
DEFAULT_FIGURE = (
    WORK_DIR
    / "magnicon_xxf1_thermal_backbone_pure_zero_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from lib.tes_noise_model import THERMAL_EXTENSION_TES_BATH_SERIES  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_separated_substrate_joint_c2_diagnostic as separated  # noqa: E402
from subScript import magnicon_xxf1_shared_L_Ctes_joint_c2_diagnostic as shared  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402


BRANCH_ZERO = "c2_zero"
BRANCH_C2 = "c2_reference"
BRANCH_PURE_ZERO = "pure_zero"
BRANCHES = (BRANCH_ZERO, BRANCH_C2, BRANCH_PURE_ZERO)
DAY_PARAMETERS = ("alpha", "beta", "T_bath")


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def pb_absorber_values(
    thickness_m: float,
    *,
    length_m: float,
    width_m: float,
    rho_kg_per_m3: float,
    cp_J_per_kgK: float,
    k_W_per_mK: float,
) -> dict:
    thickness = float(thickness_m)
    length = float(length_m)
    width = float(width_m)
    rho = float(rho_kg_per_m3)
    cp = float(cp_J_per_kgK)
    conductivity = float(k_W_per_mK)
    if min(thickness, length, width, rho, cp, conductivity) <= 0.0:
        raise ValueError("Pb geometry and material values must be positive")
    volume = length * width * thickness
    area = width * thickness
    return {
        "thickness_m": thickness,
        "volume_m3": float(volume),
        "C_abs_J_per_K": float(rho * cp * volume),
        "G_abs_abs_W_per_K": float(conductivity * area / length),
    }


def pure_zero_to_c2(zero_hz: float, scale_hz: float) -> float:
    zero = float(zero_hz)
    scale = float(scale_hz)
    if zero <= 0.0 or scale <= 0.0:
        raise ValueError("zero and reference scale must be positive")
    return float((scale / zero) ** 2)


def c2_to_pure_zero(c2: float, scale_hz: float) -> float | None:
    value = float(c2)
    scale = float(scale_hz)
    if value <= 0.0:
        return None
    return float(scale / np.sqrt(value))


def _shared_bounds(config: dict, base_config: dict, material_C_tes: float) -> dict:
    profile = base_config["shared_physical_profile"]
    thermal = config["shared_thermal_backbone_profile"]

    bounds = {
        "C_tes": (
            material_C_tes
            * float(profile["C_tes_material_multiplier"]["min"]),
            material_C_tes
            * float(profile["C_tes_material_multiplier"]["max"]),
        ),
        "L": (
            float(profile["L_H"]["min"]),
            float(profile["L_H"]["max"]),
        ),
        "C_substrate": (
            float(profile["C_substrate_J_per_K"]["min"]),
            float(profile["C_substrate_J_per_K"]["max"]),
        ),
        "G_ratio": (
            float(
                profile["G_tes_substrate_over_G_substrate_bath"]["min"]
            ),
            float(
                profile["G_tes_substrate_over_G_substrate_bath"]["max"]
            ),
        ),
        "G_tes_bath_scale": (
            float(thermal["G_tes_bath_scale_to_snapshot"]["min"]),
            float(thermal["G_tes_bath_scale_to_snapshot"]["max"]),
        ),
        "Pb_thickness_m": (
            float(thermal["Pb_absorber"]["thickness_m"]["min"]),
            float(thermal["Pb_absorber"]["thickness_m"]["max"]),
        ),
    }
    for name, (lower, upper) in bounds.items():
        if not 0.0 < lower < upper:
            raise ValueError(f"invalid shared bounds for {name}")
    return bounds


def _day_bounds(config: dict, reference_problem: dict, repeat_problem: dict):
    override = config.get("day_specific_bounds", {})
    t_cfg = override.get("T_bath_K")
    result = {}
    for prefix, problem in (
        ("reference", reference_problem),
        ("repeat", repeat_problem),
    ):
        rows = {}
        for name in DAY_PARAMETERS:
            if name == "T_bath" and t_cfg is not None:
                rows[name] = (
                    float(t_cfg["min"]),
                    float(t_cfg["max"]),
                )
            else:
                low, high = problem["detector_bounds"][name]
                rows[name] = (float(low), float(high))
            if not rows[name][0] < rows[name][1]:
                raise ValueError(f"invalid {prefix} bound for {name}")
        result[prefix] = rows
    return result


def build_spec(
    *,
    branch: str,
    config: dict,
    base_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    if branch not in BRANCHES:
        raise ValueError(f"unknown branch: {branch}")

    shared_bounds = _shared_bounds(config, base_config, material_C_tes)
    shared_names = tuple(shared_bounds)
    day_bounds = _day_bounds(config, reference_problem, repeat_problem)

    white_bounds = {}
    for prefix, problem in (
        ("reference", reference_problem),
        ("repeat", repeat_problem),
    ):
        white_cfg = problem["continuum_config"]["white_profile"]
        white_bounds[prefix] = (
            float(white_cfg["min_scale"]),
            float(white_cfg["max_scale"]),
        )

    names = []
    bounds = []
    for name in shared_names:
        low, high = shared_bounds[name]
        names.append(f"shared_log10_{name}")
        bounds.append((np.log10(low), np.log10(high)))

    for prefix in ("reference", "repeat"):
        for name in DAY_PARAMETERS:
            low, high = day_bounds[prefix][name]
            names.append(f"{prefix}_{name}")
            bounds.append(
                (
                    competition.encode_value(name, low),
                    competition.encode_value(name, high),
                )
            )

    for prefix in ("reference", "repeat"):
        low, high = white_bounds[prefix]
        names.append(f"{prefix}_log10_white_scale")
        bounds.append((np.log10(low), np.log10(high)))

    branch_cfg = config["readout_branches"][branch]
    if branch == BRANCH_C2:
        names.append("shared_c2")
        bounds.append((0.0, float(branch_cfg["upper_bound"])))
    elif branch == BRANCH_PURE_ZERO:
        zero_cfg = branch_cfg["zero_Hz"]
        names.append("shared_log10_pure_zero_Hz")
        bounds.append(
            (
                np.log10(float(zero_cfg["min"])),
                np.log10(float(zero_cfg["max"])),
            )
        )

    return {
        "branch": branch,
        "names": tuple(names),
        "bounds": tuple(bounds),
        "shared_names": shared_names,
        "shared_bounds": shared_bounds,
        "day_bounds": day_bounds,
        "white_bounds": white_bounds,
    }


def decode_vector(
    vector,
    *,
    spec: dict,
    config: dict,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
) -> dict:
    vector = np.asarray(vector, dtype=float)
    if vector.size != len(spec["names"]):
        raise ValueError("joint vector length mismatch")

    offset = 0
    shared_values = {}
    for name in spec["shared_names"]:
        shared_values[name] = float(10.0 ** vector[offset])
        offset += 1

    pb_cfg = config["shared_thermal_backbone_profile"]["Pb_absorber"]
    pb = pb_absorber_values(
        shared_values["Pb_thickness_m"],
        length_m=float(pb_cfg["length_m"]),
        width_m=float(pb_cfg["width_m"]),
        rho_kg_per_m3=float(pb_cfg["rho_kg_per_m3"]),
        cp_J_per_kgK=float(pb_cfg["cp_J_per_kgK"]),
        k_W_per_mK=float(pb_cfg["k_W_per_mK"]),
    )

    detectors = {}
    for prefix, problem in (
        ("reference", reference_problem),
        ("repeat", repeat_problem),
    ):
        detector = copy.deepcopy(problem["baseline_detector"])
        detector["C_tes"] = shared_values["C_tes"]
        detector["L"] = shared_values["L"]
        detector["C_substrate"] = shared_values["C_substrate"]
        detector["thermal_extension"] = THERMAL_EXTENSION_TES_BATH_SERIES

        detector["G_tes-bath"] = (
            float(problem["baseline_detector"]["G_tes-bath"])
            * shared_values["G_tes_bath_scale"]
        )
        g1, g2 = separated.split_series_conductance(
            detector["G_tes-bath"],
            shared_values["G_ratio"],
        )
        detector["G_tes-substrate"] = g1
        detector["G_substrate-bath"] = g2

        detector["C_abs"] = pb["C_abs_J_per_K"]
        detector["G_abs-abs"] = pb["G_abs_abs_W_per_K"]

        for name in DAY_PARAMETERS:
            detector[name] = competition.decode_value(
                name,
                vector[offset],
            )
            offset += 1
        detectors[prefix] = detector

    white_scales = {}
    for prefix in ("reference", "repeat"):
        white_scales[prefix] = float(10.0 ** vector[offset])
        offset += 1

    branch = spec["branch"]
    pure_zero_hz = None
    c2 = 0.0
    if branch == BRANCH_C2:
        c2 = float(vector[offset])
        offset += 1
    elif branch == BRANCH_PURE_ZERO:
        pure_zero_hz = float(10.0 ** vector[offset])
        offset += 1
        scale_ref = float(reference_problem["scale_hz"])
        scale_rep = float(repeat_problem["scale_hz"])
        if not np.isclose(scale_ref, scale_rep, rtol=0.0, atol=1e-12):
            raise ValueError("reference/repeat readout scales differ")
        c2 = pure_zero_to_c2(pure_zero_hz, scale_ref)

    if offset != vector.size:
        raise ValueError("joint vector decode mismatch")

    readout = dict(readout_reference)
    readout["c2"] = c2
    return {
        "shared": shared_values,
        "pb": pb,
        "detectors": detectors,
        "white_scales": white_scales,
        "readout": readout,
        "pure_zero_Hz": pure_zero_hz,
    }


def warm_vector(
    *,
    spec: dict,
    config: dict,
    reference_problem: dict,
    repeat_problem: dict,
    prior_solution: dict | None = None,
) -> np.ndarray:
    values = []
    for name in spec["shared_names"]:
        low, high = spec["shared_bounds"][name]
        value = float(np.sqrt(low * high))
        if prior_solution is not None:
            key_map = {
                "C_tes": "shared_C_tes_J_per_K",
                "L": "shared_L_H",
                "C_substrate": "shared_C_substrate_J_per_K",
                "G_ratio": "shared_G_tes_substrate_over_G_substrate_bath",
                "G_tes_bath_scale": "shared_G_tes_bath_scale_to_snapshot",
                "Pb_thickness_m": "shared_Pb_absorber_thickness_m",
            }
            candidate = prior_solution.get(key_map[name])
            if candidate is not None:
                value = float(np.clip(candidate, low, high))
        values.append(np.log10(value))

    for prefix, problem in (
        ("reference", reference_problem),
        ("repeat", repeat_problem),
    ):
        detector = problem["baseline_detector"]
        if prior_solution is not None:
            detector = prior_solution[f"{prefix}_day"]["detector_candidate"]
        for name in DAY_PARAMETERS:
            low, high = spec["day_bounds"][prefix][name]
            value = float(np.clip(detector[name], low, high))
            values.append(competition.encode_value(name, value))

    for prefix in ("reference", "repeat"):
        white = 1.0
        if prior_solution is not None:
            white = float(prior_solution[f"{prefix}_day"]["white_scale"])
        low, high = spec["white_bounds"][prefix]
        values.append(np.log10(float(np.clip(white, low, high))))

    if spec["branch"] == BRANCH_C2:
        guess = 100.0
        if prior_solution is not None:
            guess = float(prior_solution.get("shared_c2", guess))
        upper = config["readout_branches"][BRANCH_C2]["upper_bound"]
        values.append(float(np.clip(guess, 0.0, float(upper))))
    elif spec["branch"] == BRANCH_PURE_ZERO:
        guess = 600.0
        if prior_solution is not None:
            c2 = float(prior_solution.get("shared_c2", 0.0))
            mapped = c2_to_pure_zero(c2, reference_problem["scale_hz"])
            if mapped is not None:
                guess = mapped
        zero_cfg = config["readout_branches"][BRANCH_PURE_ZERO]["zero_Hz"]
        guess = float(
            np.clip(guess, float(zero_cfg["min"]), float(zero_cfg["max"]))
        )
        values.append(np.log10(guess))

    return np.asarray(values, dtype=float)


def fit_branch(
    *,
    branch: str,
    config: dict,
    base_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
    warm_from: dict | None = None,
) -> dict:
    spec = build_spec(
        branch=branch,
        config=config,
        base_config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    lower = np.asarray([row[0] for row in spec["bounds"]], dtype=float)
    upper = np.asarray([row[1] for row in spec["bounds"]], dtype=float)
    penalty = float(reference_problem["optimizer_cfg"]["instability_penalty"])

    template_ref = opt.weighted_residual_vector(
        np.ones_like(reference_problem["target"]),
        reference_problem["target"],
        reference_problem["frequency"],
        reference_problem["fit_args"],
    )
    template_rep = opt.weighted_residual_vector(
        np.ones_like(repeat_problem["target"]),
        repeat_problem["target"],
        repeat_problem["frequency"],
        repeat_problem["fit_args"],
    )
    invalid_residual = np.full(
        template_ref.size + template_rep.size,
        np.sqrt(penalty),
        dtype=float,
    )

    evaluation_count = 0
    invalid_reference = 0
    invalid_repeat = 0

    def evaluate(vector):
        nonlocal evaluation_count, invalid_reference, invalid_repeat
        evaluation_count += 1
        decoded = decode_vector(
            vector,
            spec=spec,
            config=config,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
            readout_reference=readout_reference,
        )
        models = {}
        points = {}
        for prefix, problem in (
            ("reference", reference_problem),
            ("repeat", repeat_problem),
        ):
            white_asd = (
                float(problem["tracked_white"])
                * decoded["white_scales"][prefix]
            )
            try:
                model, point = residual.model_for_candidate(
                    decoded["detectors"][prefix],
                    decoded["readout"],
                    problem["frequency"],
                    problem["scale_hz"],
                    white_asd,
                )
            except Exception:
                model, point = None, {}
            if model is None:
                if prefix == "reference":
                    invalid_reference += 1
                else:
                    invalid_repeat += 1
                return None, decoded, points.get("reference", {}), point
            models[prefix] = np.asarray(model, dtype=float)
            points[prefix] = point
        return models, decoded, points["reference"], points["repeat"]

    def objective(vector):
        models, _, _, _ = evaluate(vector)
        if models is None:
            return penalty
        scores = []
        for prefix, problem in (
            ("reference", reference_problem),
            ("repeat", repeat_problem),
        ):
            scores.append(
                opt.fit_score(
                    models[prefix],
                    problem["target"],
                    problem["frequency"],
                    problem["fit_args"],
                )
            )
        return float(np.mean(scores))

    def ls_residual(vector):
        models, _, _, _ = evaluate(vector)
        if models is None:
            return invalid_residual
        parts = []
        for prefix, problem in (
            ("reference", reference_problem),
            ("repeat", repeat_problem),
        ):
            parts.append(
                opt.weighted_residual_vector(
                    models[prefix],
                    problem["target"],
                    problem["frequency"],
                    problem["fit_args"],
                )
            )
        return np.concatenate(parts) / np.sqrt(2.0)

    optimizer = config["optimizer"]
    seed = (
        int(reference_problem["optimizer_cfg"]["seed"])
        + int(optimizer["seed_offset"])
        + 1000 * BRANCHES.index(branch)
    )
    de = differential_evolution(
        objective,
        spec["bounds"],
        seed=seed,
        maxiter=int(optimizer["DE_maxiter"]),
        popsize=int(optimizer["DE_popsize"]),
        tol=1.0e-8,
        polish=False,
        workers=1,
        updating="immediate",
    )
    candidates = [("de", np.asarray(de.x, dtype=float))]

    prior_solution = warm_from["solution"] if warm_from is not None else None
    warm = warm_vector(
        spec=spec,
        config=config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        prior_solution=prior_solution,
    )
    if np.all(warm >= lower) and np.all(warm <= upper):
        candidates.append(("physical_midpoint_warm", warm))

    for label, vector in list(candidates):
        if objective(vector) >= penalty:
            continue
        try:
            ls = least_squares(
                ls_residual,
                vector,
                bounds=(lower, upper),
                max_nfev=int(optimizer["least_squares_max_nfev"]),
                x_scale="jac",
            )
        except Exception:
            continue
        candidates.append(
            (f"{label}_least_squares", np.asarray(ls.x, dtype=float))
        )

    scored = []
    for label, vector in candidates:
        score = objective(vector)
        if score < penalty:
            scored.append((float(score), label, vector))
    if not scored:
        raise RuntimeError(f"no stable solution for {branch}")

    score, source, vector = min(scored, key=lambda row: row[0])
    models, decoded, point_ref, point_rep = evaluate(vector)
    if models is None:
        raise RuntimeError("selected solution became invalid")

    ref_metrics = shared._day_metrics(reference_problem, models["reference"])
    rep_metrics = shared._day_metrics(repeat_problem, models["repeat"])
    combined_rms = shared._combined_rms(ref_metrics, rep_metrics)

    shared_hits = {
        name: shared._encoded_boundary(
            decoded["shared"][name],
            *spec["shared_bounds"][name],
            log=True,
        )
        for name in spec["shared_names"]
    }
    white_hits = {
        prefix: shared._encoded_boundary(
            decoded["white_scales"][prefix],
            *spec["white_bounds"][prefix],
            log=True,
        )
        for prefix in ("reference", "repeat")
    }

    readout_hit = None
    if branch == BRANCH_C2:
        readout_hit = shared._encoded_boundary(
            decoded["readout"]["c2"],
            0.0,
            float(config["readout_branches"][BRANCH_C2]["upper_bound"]),
            log=False,
        )
    elif branch == BRANCH_PURE_ZERO:
        zcfg = config["readout_branches"][BRANCH_PURE_ZERO]["zero_Hz"]
        readout_hit = shared._encoded_boundary(
            decoded["pure_zero_Hz"],
            float(zcfg["min"]),
            float(zcfg["max"]),
            log=True,
        )

    def day_solution(prefix, problem, point, metrics):
        detector = decoded["detectors"][prefix]
        return {
            "detector_candidate": {
                key: (
                    str(detector[key])
                    if key == "thermal_extension"
                    else float(detector[key])
                )
                for key in (
                    "alpha",
                    "beta",
                    "T_bath",
                    "C_tes",
                    "C_substrate",
                    "C_abs",
                    "L",
                    "R",
                    "G_tes-bath",
                    "G_tes-substrate",
                    "G_substrate-bath",
                    "G_tes-stycast",
                    "G_stycast-abs",
                    "G_abs-abs",
                    "thermal_extension",
                )
                if key in detector
            },
            "day_specific_boundary_hits": {
                name: shared._encoded_boundary(
                    detector[name],
                    *spec["day_bounds"][prefix][name],
                    log=False,
                )
                for name in DAY_PARAMETERS
            },
            "white_scale": float(decoded["white_scales"][prefix]),
            "white_asd_A_rtHz": float(
                problem["tracked_white"] * decoded["white_scales"][prefix]
            ),
            "white_boundary_hit": white_hits[prefix],
            "operating_point": {
                "valid": bool(point.get("valid")),
                "stable": bool(point.get("stable")),
                "current_A": (
                    float(point["current_A"])
                    if point.get("current_A") is not None
                    else None
                ),
                "joule_power_W": (
                    float(point["joule_power_W"])
                    if point.get("joule_power_W") is not None
                    else None
                ),
                "f_substrate_node_Hz": (
                    float(point["f_substrate_node_Hz"])
                    if point.get("f_substrate_node_Hz") is not None
                    else None
                ),
            },
            "metrics": metrics,
        }

    shared_values = decoded["shared"]
    historical = float(
        base_config["shared_physical_profile"]["C_substrate_J_per_K"][
            "historical_PoST_total_reference_J_per_K"
        ]
    )
    local_total = shared_values["C_tes"] + shared_values["C_substrate"]

    solution = {
        "shared_C_tes_J_per_K": float(shared_values["C_tes"]),
        "shared_C_tes_over_material": float(
            shared_values["C_tes"] / material_C_tes
        ),
        "shared_C_substrate_J_per_K": float(shared_values["C_substrate"]),
        "shared_local_C_total_J_per_K": float(local_total),
        "shared_local_C_total_over_historical_PoST": float(
            local_total / historical
        ),
        "shared_L_H": float(shared_values["L"]),
        "shared_G_tes_substrate_over_G_substrate_bath": float(
            shared_values["G_ratio"]
        ),
        "shared_G_tes_bath_scale_to_snapshot": float(
            shared_values["G_tes_bath_scale"]
        ),
        "shared_Pb_absorber_thickness_m": float(
            shared_values["Pb_thickness_m"]
        ),
        "shared_Pb_C_abs_J_per_K": float(decoded["pb"]["C_abs_J_per_K"]),
        "shared_Pb_G_abs_abs_W_per_K": float(
            decoded["pb"]["G_abs_abs_W_per_K"]
        ),
        "shared_c2": float(decoded["readout"]["c2"]),
        "shared_pure_zero_Hz": (
            float(decoded["pure_zero_Hz"])
            if decoded["pure_zero_Hz"] is not None
            else None
        ),
        "shared_physical_boundary_hits": shared_hits,
        "readout_boundary_hit": readout_hit,
        "reference_day": day_solution(
            "reference", reference_problem, point_ref, ref_metrics
        ),
        "repeat_day": day_solution(
            "repeat", repeat_problem, point_rep, rep_metrics
        ),
    }
    return {
        "name": branch,
        "joint_shape_score": float(score),
        "joint_continuum_rms_dB": float(combined_rms),
        "best_candidate_source": source,
        "evaluation_count": int(evaluation_count),
        "invalid_reference_evaluations": int(invalid_reference),
        "invalid_repeat_evaluations": int(invalid_repeat),
        "free_parameter_names": list(spec["names"]),
        "n_free_parameters": int(len(spec["names"])),
        "solution": solution,
        "_models": models,
    }


def compare_branches(zero: dict, c2: dict, pure: dict, screen: dict) -> dict:
    ratio_pure_zero = float(
        pure["joint_shape_score"] / zero["joint_shape_score"]
    )
    rms_gain = float(
        zero["joint_continuum_rms_dB"] - pure["joint_continuum_rms_dB"]
    )
    rms_delta_c2 = float(
        pure["joint_continuum_rms_dB"] - c2["joint_continuum_rms_dB"]
    )
    material = bool(
        ratio_pure_zero <= float(screen["max_score_ratio_to_c2_zero"])
        and rms_gain >= float(screen["min_joint_rms_improvement_dB"])
    )
    near_c2 = bool(
        rms_delta_c2
        <= float(screen["max_rms_delta_from_c2_reference_dB"])
    )

    c2_hit = c2["solution"]["readout_boundary_hit"]
    c2_censored = bool(c2_hit and c2_hit.get("at_upper"))
    pure_hit = pure["solution"]["readout_boundary_hit"]
    pure_censored = bool(
        pure_hit
        and (pure_hit.get("at_lower") or pure_hit.get("at_upper"))
    )

    if material and near_c2 and not pure_censored:
        classification = (
            "pure_zero_replaces_c2_after_thermal_backbone_profile"
        )
    elif material:
        classification = (
            "pure_zero_material_after_thermal_backbone_profile"
        )
    else:
        classification = (
            "thermal_backbone_profile_reduces_zero_like_readout_need"
        )

    return {
        "classification": classification,
        "pure_zero_material_improvement_vs_c2_zero": material,
        "pure_zero_within_c2_reference_rms_screen": near_c2,
        "pure_zero_score_ratio_to_c2_zero": ratio_pure_zero,
        "pure_zero_rms_improvement_vs_c2_zero_dB": rms_gain,
        "pure_zero_rms_delta_from_c2_reference_dB": rms_delta_c2,
        "c2_reference_upper_bound_censored": c2_censored,
        "pure_zero_boundary_censored": pure_censored,
        "screen": screen,
    }


def run(config: dict, config_path: Path) -> dict:
    base_path = resolve_config_path(
        config["base_separated_substrate_config"], config_path
    )
    base_config = json.loads(base_path.read_text(encoding="utf-8"))
    stack = compensation._load_stack(base_config, base_path)
    base_cross = stack["base_cross"]
    magnicon_config = stack["magnicon_config"]

    reference_problem = crossday._problem_for_case(
        magnicon_config,
        stack["magnicon_config_path"],
        str(base_cross["reference_case_label"]),
    )
    repeat_problem = crossday._problem_for_case(
        magnicon_config,
        stack["magnicon_config_path"],
        str(base_cross["repeat_case_label"]),
    )
    if not np.allclose(
        reference_problem["frequency"],
        repeat_problem["frequency"],
        rtol=0.0,
        atol=0.0,
    ):
        raise ValueError("reference/repeat frequency grids differ")

    filter_cfg = base_config["magnicon_filter"]
    canonical = xxf1.second_order_bessel_canonical(
        float(filter_cfg["cutoff_Hz"]),
        str(filter_cfg["normalization"]),
    )
    readout_reference = {
        "pole_Hz": float(canonical["pole_Hz"]),
        "pole_Q": float(canonical["pole_Q"]),
        "c2": 0.0,
        "c4": 0.0,
    }
    material_C_tes = float(opt.C_TES_MATERIAL_J_PER_K)

    zero = fit_branch(
        branch=BRANCH_ZERO,
        config=config,
        base_config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
    )
    c2 = fit_branch(
        branch=BRANCH_C2,
        config=config,
        base_config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
        warm_from=zero,
    )
    pure = fit_branch(
        branch=BRANCH_PURE_ZERO,
        config=config,
        base_config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
        warm_from=c2,
    )
    comparison = compare_branches(
        zero, c2, pure, config["materiality_screen"]
    )

    pb_cfg = config["shared_thermal_backbone_profile"]["Pb_absorber"]
    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "After allowing modestly lower T_bath, broad shared effective "
            "G_tes-bath, and Pb C/G tied to measured 0.4-0.6 mm thickness, "
            "does a shared zero-like readout numerator remain necessary?"
        ),
        "fit_semantics": {
            "shared_across_days": [
                "C_tes",
                "L",
                "C_substrate",
                "G_tes-substrate/G_substrate-bath",
                "G_tes-bath scale to inherited snapshot",
                "Pb absorber thickness",
            ],
            "day_specific": [
                "alpha",
                "beta",
                "T_bath",
                "post-filter white ASD",
            ],
            "Pb_absorber": {
                "length_m": float(pb_cfg["length_m"]),
                "width_m": float(pb_cfg["width_m"]),
                "rho_kg_per_m3": float(pb_cfg["rho_kg_per_m3"]),
                "cp_J_per_kgK": float(pb_cfg["cp_J_per_kgK"]),
                "k_W_per_mK": float(pb_cfg["k_W_per_mK"]),
                "material_source": str(pb_cfg["material_source"]),
                "C_formula": "rho*cp*length*width*thickness",
                "G_formula": "k*width*thickness/length",
            },
            "Magnicon": {
                "normalization": str(filter_cfg["normalization"]),
                "cutoff_Hz": float(filter_cfg["cutoff_Hz"]),
                "pole_Hz": float(canonical["pole_Hz"]),
                "pole_Q": float(canonical["pole_Q"]),
                "c4": 0.0,
            },
            "external_analog_filter": {
                "model": "SRS SIM965",
                "front_panel_cutoff_Hz": float(
                    opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ
                ),
                "mode": "Bessel",
                "slope_dB_per_oct": 24,
                "order": int(opt.TARGET_HARDWARE_BESSEL_ORDER),
                "coupling": "DC",
                "scipy_norm": str(opt.TARGET_HARDWARE_BESSEL_NORM),
                "nominal_minus3dB_Hz": float(
                    0.6604 * opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ
                ),
            },
        },
        "c2_zero": {
            key: value for key, value in zero.items() if not key.startswith("_")
        },
        "c2_reference": {
            key: value for key, value in c2.items() if not key.startswith("_")
        },
        "pure_zero": {
            key: value for key, value in pure.items() if not key.startswith("_")
        },
        "comparison": comparison,
        "interpretation": {
            "classification": comparison["classification"],
            "pure_zero_Hz": float(
                pure["solution"]["shared_pure_zero_Hz"]
            ),
            "pure_zero_equivalent_c2": float(
                pure["solution"]["shared_c2"]
            ),
            "c2_reference": float(c2["solution"]["shared_c2"]),
            "guardrail": str(config["guardrail"]),
        },
        "inputs": {
            "config": str(config_path),
            "base_separated_substrate_config": str(base_path),
            "base_cross_day_config": str(stack["base_cross_path"]),
            "base_magnicon_config": str(stack["magnicon_config_path"]),
        },
        "_plot": {
            "frequency_Hz": reference_problem["frequency"].tolist(),
            "reference_target": reference_problem["target"].tolist(),
            "repeat_target": repeat_problem["target"].tolist(),
            "reference_raw": reference_problem["raw_target"].tolist(),
            "repeat_raw": repeat_problem["raw_target"].tolist(),
            "reference_zero": zero["_models"]["reference"].tolist(),
            "repeat_zero": zero["_models"]["repeat"].tolist(),
            "reference_c2": c2["_models"]["reference"].tolist(),
            "repeat_c2": c2["_models"]["repeat"].tolist(),
            "reference_pure": pure["_models"]["reference"].tolist(),
            "repeat_pure": pure["_models"]["repeat"].tolist(),
        },
    }


def _residual_db(model, target):
    return 20.0 * np.log10(
        np.asarray(model, dtype=float) / np.asarray(target, dtype=float)
    )


def make_plot(result: dict, output: Path, show=False) -> None:
    import matplotlib.pyplot as plt

    p = result["_plot"]
    f = np.asarray(p["frequency_Hz"], dtype=float)
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(13.8, 8.8),
        sharex="col",
        gridspec_kw={"height_ratios": [2.0, 1.0]},
    )
    rows = [
        (
            0,
            "2024-12-06 reference",
            np.asarray(p["reference_target"], dtype=float),
            np.asarray(p["reference_raw"], dtype=float),
            np.asarray(p["reference_zero"], dtype=float),
            np.asarray(p["reference_c2"], dtype=float),
            np.asarray(p["reference_pure"], dtype=float),
        ),
        (
            1,
            "2024-12-05 repeat",
            np.asarray(p["repeat_target"], dtype=float),
            np.asarray(p["repeat_raw"], dtype=float),
            np.asarray(p["repeat_zero"], dtype=float),
            np.asarray(p["repeat_c2"], dtype=float),
            np.asarray(p["repeat_pure"], dtype=float),
        ),
    ]
    for col, label, target, raw, zero, c2, pure in rows:
        top = axes[0, col]
        bottom = axes[1, col]
        top.loglog(f, raw, linewidth=0.7, alpha=0.35, label="Raw ASD")
        top.loglog(f, target, linewidth=1.5, label="Continuum target")
        top.loglog(f, zero, linewidth=1.2, label="thermal profile, c2=0")
        top.loglog(f, c2, linewidth=1.2, label="thermal profile, c2 ref")
        top.loglog(f, pure, linewidth=1.4, label="thermal profile, pure zero")
        top.set_title(label, fontsize=10)
        top.set_ylabel("Normalized ASD")
        top.grid(True, which="both", alpha=0.2)
        top.legend(frameon=False, fontsize=7)

        for model, name in (
            (zero, "c2=0"),
            (c2, "c2 reference"),
            (pure, "pure zero"),
        ):
            bottom.semilogx(
                f,
                _residual_db(model, target),
                linewidth=1.1,
                label=name,
            )
        bottom.axhline(0.0, linewidth=0.9)
        bottom.set_xlabel("Frequency [Hz]")
        bottom.set_ylabel("Model / continuum [dB]")
        bottom.grid(True, which="both", alpha=0.2)
        bottom.legend(frameon=False, fontsize=7)

    z = result["c2_zero"]["solution"]
    q = result["pure_zero"]["solution"]
    fig.suptitle(
        result["interpretation"]["classification"]
        + "\n"
        + (
            f"c2=0: Gbath scale={z['shared_G_tes_bath_scale_to_snapshot']:.3g}, "
            f"Pb t={z['shared_Pb_absorber_thickness_m']*1e3:.3f} mm; "
            f"pure zero={q['shared_pure_zero_Hz']:.3g} Hz, "
            f"Gbath scale={q['shared_G_tes_bath_scale_to_snapshot']:.3g}, "
            f"Pb t={q['shared_Pb_absorber_thickness_m']*1e3:.3f} mm"
        ),
        fontsize=9,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result: dict) -> dict:
    return {key: value for key, value in result.items() if key != "_plot"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure", type=Path, default=DEFAULT_FIGURE)
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(config, args.config)
    make_plot(result, args.figure, show=args.show)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(cleaned_result(result), indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "classification": result["interpretation"]["classification"],
                "c2_zero": {
                    "joint_rms_dB": result["c2_zero"][
                        "joint_continuum_rms_dB"
                    ],
                    "solution": result["c2_zero"]["solution"],
                },
                "c2_reference": {
                    "joint_rms_dB": result["c2_reference"][
                        "joint_continuum_rms_dB"
                    ],
                    "c2": result["c2_reference"]["solution"]["shared_c2"],
                },
                "pure_zero": {
                    "joint_rms_dB": result["pure_zero"][
                        "joint_continuum_rms_dB"
                    ],
                    "zero_Hz": result["pure_zero"]["solution"][
                        "shared_pure_zero_Hz"
                    ],
                    "equivalent_c2": result["pure_zero"]["solution"][
                        "shared_c2"
                    ],
                },
                "comparison": result["comparison"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
