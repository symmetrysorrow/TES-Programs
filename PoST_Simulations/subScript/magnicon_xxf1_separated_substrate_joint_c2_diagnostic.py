"""Joint c2 test with an explicit TES--substrate--bath thermal node.

This diagnostic asks whether the heat capacity historically folded into C_tes
is better represented by a separate local Si/substrate-membrane node.

Topology for each TES side:

    TES -- G_tes-substrate -- C_substrate -- G_substrate-bath -- bath
     |
     +-- existing TES--Stycast--Pb path

The two substrate conductances are parameterized by one shared ratio while
their series-equivalent conductance is constrained to equal the existing
G_tes-bath exactly for each day.  Therefore this diagnostic adds a thermal
time scale without changing the reduced model's DC TES-to-bath conductance.

Shared across the adjacent-day datasets:
  * C_tes
  * L
  * C_substrate
  * G_tes-substrate / G_substrate-bath
  * optional c2 in the c2-free branch

Day specific:
  * alpha, beta, T_bath
  * post-filter white ASD scale

The Magnicon XXF-1 response remains the nominal phase-normalized 10 kHz
second-order Bessel with c4=0.

This is diagnostic-only.  C_substrate is a lumped local thermal body, not a
claim that the geometrically resolved silicon is isothermal or that the
historical PoST C_tes=7.9e-12 J/K was explicitly calculated from silicon.
"""

from __future__ import annotations

import argparse
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
    / "magnicon_xxf1_separated_substrate_joint_c2_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR
    / "magnicon_xxf1_separated_substrate_joint_c2_diagnostic.json"
)
DEFAULT_FIGURE = (
    WORK_DIR
    / "magnicon_xxf1_separated_substrate_joint_c2_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from lib.tes_noise_model import (  # noqa: E402
    THERMAL_EXTENSION_TES_BATH_SERIES,
)
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_shared_L_Ctes_joint_c2_diagnostic as shared  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402


DAY_PARAMETERS = ("alpha", "beta", "T_bath")


def split_series_conductance(
    equivalent: float,
    tes_over_bath_ratio: float,
) -> tuple[float, float]:
    """Return G_tes-substrate, G_substrate-bath preserving series G."""
    equivalent = float(equivalent)
    ratio = float(tes_over_bath_ratio)
    if equivalent <= 0.0 or ratio <= 0.0:
        raise ValueError("equivalent conductance and ratio must be positive")
    g_tes_substrate = equivalent * (1.0 + ratio)
    g_substrate_bath = equivalent * (1.0 + ratio) / ratio
    return float(g_tes_substrate), float(g_substrate_bath)


def series_equivalent(g1: float, g2: float) -> float:
    g1 = float(g1)
    g2 = float(g2)
    if g1 <= 0.0 or g2 <= 0.0:
        raise ValueError("conductances must be positive")
    return float(1.0 / (1.0 / g1 + 1.0 / g2))


def _build_spec(
    *,
    config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
    c2_free: bool,
) -> dict:
    profile = config["shared_physical_profile"]

    ctes_cfg = profile["C_tes_material_multiplier"]
    csub_cfg = profile["C_substrate_J_per_K"]
    l_cfg = profile["L_H"]
    ratio_cfg = profile["G_tes_substrate_over_G_substrate_bath"]

    shared_bounds = {
        "C_tes": (
            material_C_tes * float(ctes_cfg["min"]),
            material_C_tes * float(ctes_cfg["max"]),
        ),
        "L": (
            float(l_cfg["min"]),
            float(l_cfg["max"]),
        ),
        "C_substrate": (
            float(csub_cfg["min"]),
            float(csub_cfg["max"]),
        ),
        "G_ratio": (
            float(ratio_cfg["min"]),
            float(ratio_cfg["max"]),
        ),
    }
    for name, (lower, upper) in shared_bounds.items():
        if not 0.0 < lower < upper:
            raise ValueError(f"invalid shared bound for {name}")

    day_names = tuple(config["day_specific_detector_parameters"])
    if day_names != DAY_PARAMETERS:
        raise ValueError(
            "this diagnostic requires day-specific alpha, beta, T_bath"
        )

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
    for name in ("C_tes", "L", "C_substrate", "G_ratio"):
        low, high = shared_bounds[name]
        names.append(f"shared_log10_{name}")
        bounds.append((np.log10(low), np.log10(high)))

    for prefix, problem in (
        ("reference", reference_problem),
        ("repeat", repeat_problem),
    ):
        for name in day_names:
            low, high = problem["detector_bounds"][name]
            names.append(f"{prefix}_{name}")
            bounds.append(
                (
                    competition.encode_value(name, low),
                    competition.encode_value(name, high),
                )
            )

    for prefix in ("reference", "repeat"):
        low, high = white_bounds[prefix]
        if not 0.0 < low < high:
            raise ValueError("invalid white scale bounds")
        names.append(f"{prefix}_log10_white_scale")
        bounds.append((np.log10(low), np.log10(high)))

    if c2_free:
        names.append("shared_c2")
        bounds.append((0.0, float(config["c2"]["upper_bound"])))

    return {
        "names": tuple(names),
        "bounds": tuple(bounds),
        "shared_bounds": shared_bounds,
        "white_bounds": white_bounds,
        "day_names": day_names,
        "c2_free": bool(c2_free),
    }


def _decode(
    vector,
    *,
    spec: dict,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
) -> dict:
    vector = np.asarray(vector, dtype=float)
    if vector.size != len(spec["names"]):
        raise ValueError("joint vector length mismatch")

    offset = 0
    shared_C_tes = float(10.0 ** vector[offset])
    offset += 1
    shared_L = float(10.0 ** vector[offset])
    offset += 1
    shared_C_substrate = float(10.0 ** vector[offset])
    offset += 1
    shared_G_ratio = float(10.0 ** vector[offset])
    offset += 1

    detectors = {}
    for prefix, problem in (
        ("reference", reference_problem),
        ("repeat", repeat_problem),
    ):
        detector = dict(problem["baseline_detector"])
        detector["C_tes"] = shared_C_tes
        detector["L"] = shared_L
        detector["thermal_extension"] = THERMAL_EXTENSION_TES_BATH_SERIES
        detector["C_substrate"] = shared_C_substrate

        g_equivalent = float(detector["G_tes-bath"])
        g_tes_substrate, g_substrate_bath = split_series_conductance(
            g_equivalent,
            shared_G_ratio,
        )
        detector["G_tes-substrate"] = g_tes_substrate
        detector["G_substrate-bath"] = g_substrate_bath

        for name in spec["day_names"]:
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

    c2 = 0.0
    if spec["c2_free"]:
        c2 = float(vector[offset])
        offset += 1

    if offset != vector.size:
        raise ValueError("joint vector decode mismatch")

    readout = dict(readout_reference)
    readout["c2"] = c2

    return {
        "shared_C_tes": shared_C_tes,
        "shared_L": shared_L,
        "shared_C_substrate": shared_C_substrate,
        "shared_G_ratio": shared_G_ratio,
        "detectors": detectors,
        "white_scales": white_scales,
        "readout": readout,
    }


def _encode_warm(
    *,
    spec: dict,
    reference_problem: dict,
    repeat_problem: dict,
    shared_C_tes: float,
    shared_L: float,
    shared_C_substrate: float,
    shared_G_ratio: float,
    reference_detector: dict,
    repeat_detector: dict,
    reference_white_scale: float,
    repeat_white_scale: float,
    c2: float = 0.0,
) -> np.ndarray:
    values = [
        np.log10(float(shared_C_tes)),
        np.log10(float(shared_L)),
        np.log10(float(shared_C_substrate)),
        np.log10(float(shared_G_ratio)),
    ]
    for detector in (reference_detector, repeat_detector):
        for name in spec["day_names"]:
            values.append(competition.encode_value(name, detector[name]))
    values.extend(
        [
            np.log10(float(reference_white_scale)),
            np.log10(float(repeat_white_scale)),
        ]
    )
    if spec["c2_free"]:
        values.append(float(c2))
    return np.asarray(values, dtype=float)


def _baseline_warm(
    spec: dict,
    reference_problem: dict,
    repeat_problem: dict,
) -> np.ndarray:
    bounds = spec["shared_bounds"]
    shared_C_tes = float(np.sqrt(np.prod(bounds["C_tes"])))
    shared_L = float(np.sqrt(np.prod(bounds["L"])))
    shared_C_substrate = float(np.sqrt(np.prod(bounds["C_substrate"])))
    shared_G_ratio = float(np.sqrt(np.prod(bounds["G_ratio"])))

    ref = dict(reference_problem["baseline_detector"])
    rep = dict(repeat_problem["baseline_detector"])
    for detector, problem in (
        (ref, reference_problem),
        (rep, repeat_problem),
    ):
        detector["C_tes"] = shared_C_tes
        detector["L"] = shared_L
        for name in spec["day_names"]:
            low, high = problem["detector_bounds"][name]
            detector[name] = float(min(max(detector[name], low), high))

    return _encode_warm(
        spec=spec,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        shared_C_tes=shared_C_tes,
        shared_L=shared_L,
        shared_C_substrate=shared_C_substrate,
        shared_G_ratio=shared_G_ratio,
        reference_detector=ref,
        repeat_detector=rep,
        reference_white_scale=1.0,
        repeat_white_scale=1.0,
        c2=50.0 if spec["c2_free"] else 0.0,
    )


def _fit_joint_branch(
    *,
    name: str,
    config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
    c2_free: bool,
    warm_from=None,
) -> dict:
    spec = _build_spec(
        config=config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        c2_free=c2_free,
    )
    bounds = spec["bounds"]
    lower = np.asarray([row[0] for row in bounds], dtype=float)
    upper = np.asarray([row[1] for row in bounds], dtype=float)
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
        decoded = _decode(
            vector,
            spec=spec,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
            readout_reference=readout_reference,
        )
        white_ref = (
            float(reference_problem["tracked_white"])
            * decoded["white_scales"]["reference"]
        )
        white_rep = (
            float(repeat_problem["tracked_white"])
            * decoded["white_scales"]["repeat"]
        )

        try:
            model_ref, point_ref = residual.model_for_candidate(
                decoded["detectors"]["reference"],
                decoded["readout"],
                reference_problem["frequency"],
                reference_problem["scale_hz"],
                white_ref,
            )
        except Exception:
            model_ref, point_ref = None, {}
        try:
            model_rep, point_rep = residual.model_for_candidate(
                decoded["detectors"]["repeat"],
                decoded["readout"],
                repeat_problem["frequency"],
                repeat_problem["scale_hz"],
                white_rep,
            )
        except Exception:
            model_rep, point_rep = None, {}

        if model_ref is None:
            invalid_reference += 1
        if model_rep is None:
            invalid_repeat += 1
        if model_ref is None or model_rep is None:
            return None, decoded, point_ref, point_rep

        return (
            {
                "reference": np.asarray(model_ref, dtype=float),
                "repeat": np.asarray(model_rep, dtype=float),
            },
            decoded,
            point_ref,
            point_rep,
        )

    def objective(vector):
        models, _, _, _ = evaluate(vector)
        if models is None:
            return penalty
        ref_score = opt.fit_score(
            models["reference"],
            reference_problem["target"],
            reference_problem["frequency"],
            reference_problem["fit_args"],
        )
        rep_score = opt.fit_score(
            models["repeat"],
            repeat_problem["target"],
            repeat_problem["frequency"],
            repeat_problem["fit_args"],
        )
        return float(0.5 * (ref_score + rep_score))

    def ls_residual(vector):
        models, _, _, _ = evaluate(vector)
        if models is None:
            return invalid_residual
        ref = opt.weighted_residual_vector(
            models["reference"],
            reference_problem["target"],
            reference_problem["frequency"],
            reference_problem["fit_args"],
        )
        rep = opt.weighted_residual_vector(
            models["repeat"],
            repeat_problem["target"],
            repeat_problem["frequency"],
            repeat_problem["fit_args"],
        )
        return np.concatenate((ref, rep)) / np.sqrt(2.0)

    optimizer = config["optimizer"]
    seed = (
        int(reference_problem["optimizer_cfg"]["seed"])
        + int(optimizer["seed_offset"])
        + (1000 if c2_free else 0)
    )
    de = differential_evolution(
        objective,
        bounds,
        seed=seed,
        maxiter=int(optimizer["DE_maxiter"]),
        popsize=int(optimizer["DE_popsize"]),
        tol=1.0e-8,
        polish=False,
        workers=1,
        updating="immediate",
    )
    candidates = [("de", np.asarray(de.x, dtype=float))]

    warm_vectors = [
        (
            "separated_substrate_midpoint",
            _baseline_warm(spec, reference_problem, repeat_problem),
        )
    ]

    if warm_from is not None:
        solution = warm_from["solution"]
        warm_zero = _encode_warm(
            spec=spec,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
            shared_C_tes=solution["shared_C_tes_J_per_K"],
            shared_L=solution["shared_L_H"],
            shared_C_substrate=solution["shared_C_substrate_J_per_K"],
            shared_G_ratio=solution[
                "shared_G_tes_substrate_over_G_substrate_bath"
            ],
            reference_detector=solution["reference_day"][
                "detector_candidate"
            ],
            repeat_detector=solution["repeat_day"]["detector_candidate"],
            reference_white_scale=solution["reference_day"]["white_scale"],
            repeat_white_scale=solution["repeat_day"]["white_scale"],
            c2=0.0,
        )
        warm_vectors.append(("c2_zero_joint_solution", warm_zero))
        if c2_free:
            warm_guess = warm_zero.copy()
            warm_guess[-1] = min(
                50.0,
                float(config["c2"]["upper_bound"]),
            )
            warm_vectors.append(("c2_zero_plus_c2_50", warm_guess))

    for label, vector in warm_vectors:
        if np.any(vector < lower) or np.any(vector > upper):
            continue
        candidates.append((f"{label}_exact", vector.copy()))

    initial = list(candidates)
    for label, vector in initial:
        score = objective(vector)
        if score >= penalty:
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
        raise RuntimeError(f"no stable solution for {name}")

    score, source, vector = min(scored, key=lambda row: row[0])
    models, decoded, point_ref, point_rep = evaluate(vector)
    if models is None:
        raise RuntimeError("selected joint solution became invalid")

    ref_metrics = shared._day_metrics(
        reference_problem,
        models["reference"],
    )
    rep_metrics = shared._day_metrics(
        repeat_problem,
        models["repeat"],
    )
    combined_rms = shared._combined_rms(ref_metrics, rep_metrics)

    sb = spec["shared_bounds"]
    shared_hits = {
        key: shared._encoded_boundary(
            decoded[
                {
                    "C_tes": "shared_C_tes",
                    "L": "shared_L",
                    "C_substrate": "shared_C_substrate",
                    "G_ratio": "shared_G_ratio",
                }[key]
            ],
            *sb[key],
            log=True,
        )
        for key in ("C_tes", "L", "C_substrate", "G_ratio")
    }
    white_hits = {
        prefix: shared._encoded_boundary(
            decoded["white_scales"][prefix],
            *spec["white_bounds"][prefix],
            log=True,
        )
        for prefix in ("reference", "repeat")
    }
    c2_hit = (
        shared._encoded_boundary(
            decoded["readout"]["c2"],
            0.0,
            float(config["c2"]["upper_bound"]),
            log=False,
        )
        if c2_free
        else None
    )

    def day_solution(prefix, problem, point, metrics):
        detector = decoded["detectors"][prefix]
        white_scale = decoded["white_scales"][prefix]
        return {
            "detector_candidate": {
                key: (
                    float(detector[key])
                    if key != "thermal_extension"
                    else str(detector[key])
                )
                for key in (
                    "alpha",
                    "beta",
                    "C_tes",
                    "L",
                    "T_bath",
                    "R",
                    "C_substrate",
                    "G_tes-substrate",
                    "G_substrate-bath",
                    "G_tes-bath",
                    "thermal_extension",
                )
                if key in detector
            },
            "day_specific_boundary_hits": (
                competition.parameter_boundary_hits(
                    detector,
                    spec["day_names"],
                    problem["detector_bounds"],
                )
            ),
            "white_scale": float(white_scale),
            "white_asd_A_rtHz": float(
                problem["tracked_white"] * white_scale
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
                "G_substrate_equivalent_W_per_K": (
                    float(point["G_substrate_equivalent_W_per_K"])
                    if point.get("G_substrate_equivalent_W_per_K")
                    is not None
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

    historical = float(
        config["shared_physical_profile"]["C_substrate_J_per_K"][
            "historical_PoST_total_reference_J_per_K"
        ]
    )
    local_total = (
        decoded["shared_C_tes"] + decoded["shared_C_substrate"]
    )

    solution = {
        "shared_C_tes_J_per_K": float(decoded["shared_C_tes"]),
        "shared_C_tes_over_material": float(
            decoded["shared_C_tes"] / material_C_tes
        ),
        "shared_C_substrate_J_per_K": float(
            decoded["shared_C_substrate"]
        ),
        "shared_local_C_total_J_per_K": float(local_total),
        "shared_local_C_total_over_historical_PoST": float(
            local_total / historical
        ),
        "historical_PoST_C_tes_reference_J_per_K": historical,
        "shared_L_H": float(decoded["shared_L"]),
        "shared_G_tes_substrate_over_G_substrate_bath": float(
            decoded["shared_G_ratio"]
        ),
        "shared_c2": float(decoded["readout"]["c2"]),
        "shared_physical_boundary_hits": shared_hits,
        "shared_c2_boundary_hit": c2_hit,
        "reference_day": day_solution(
            "reference",
            reference_problem,
            point_ref,
            ref_metrics,
        ),
        "repeat_day": day_solution(
            "repeat",
            repeat_problem,
            point_rep,
            rep_metrics,
        ),
    }

    return {
        "name": name,
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
        "_vector": vector,
    }


def _nested_summary(zero: dict, free: dict, screen: dict) -> dict:
    base = shared._nested_summary(zero, free, screen)
    material = bool(base["shared_c2_material_improvement"])

    c2_hit = free["solution"].get("shared_c2_boundary_hit")
    c2_upper = bool(c2_hit and c2_hit.get("at_upper"))

    profile_boundary_active = False
    boundary_rows = {}
    for branch_name, branch in (("c2_zero", zero), ("c2_free", free)):
        rows = branch["solution"]["shared_physical_boundary_hits"]
        boundary_rows[branch_name] = rows
        if any(
            bool(row.get("at_lower")) or bool(row.get("at_upper"))
            for row in rows.values()
        ):
            profile_boundary_active = True

    if not material:
        classification = (
            "separated_substrate_node_removes_material_c2_need"
        )
    elif c2_upper:
        classification = (
            "shared_c2_material_but_c2_upper_bound_limited_after_"
            "separated_substrate_profile"
        )
    elif profile_boundary_active:
        classification = (
            "shared_c2_material_but_separated_substrate_profile_"
            "boundary_active"
        )
    else:
        classification = (
            "shared_c2_material_after_separated_substrate_profile"
        )

    return {
        **base,
        "classification": classification,
        "shared_physical_profile_boundary_active": profile_boundary_active,
        "shared_physical_boundary_hits": boundary_rows,
    }


def run(config: dict, config_path: Path) -> dict:
    stack = compensation._load_stack(config, config_path)
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

    filter_cfg = config["magnicon_filter"]
    canonical = xxf1.second_order_bessel_canonical(
        float(filter_cfg["cutoff_Hz"]),
        str(filter_cfg["normalization"]),
    )
    readout_reference = {
        "pole_Hz": float(canonical["pole_Hz"]),
        "pole_Q": float(canonical["pole_Q"]),
        "c2": 0.0,
        "c4": float(filter_cfg["c4_fixed"]),
    }
    if readout_reference["c4"] != 0.0:
        raise ValueError("this diagnostic requires c4=0")

    material_C_tes = float(opt.C_TES_MATERIAL_J_PER_K)

    zero = _fit_joint_branch(
        name="separated_substrate_joint_c2_zero",
        config=config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
        c2_free=False,
    )
    free = _fit_joint_branch(
        name="separated_substrate_joint_c2_free",
        config=config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
        c2_free=True,
        warm_from=zero,
    )
    nested = _nested_summary(
        zero,
        free,
        config["materiality_screen"],
    )

    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "After separating a shared local substrate/membrane heat "
            "capacity from C_tes and placing it on a DC-preserving "
            "TES-to-bath series thermal path, is a shared residual c2 "
            "still materially required?"
        ),
        "fit_semantics": {
            "shared_across_days": [
                "L",
                "C_tes",
                "C_substrate",
                "G_tes-substrate/G_substrate-bath",
            ],
            "day_specific": [
                "alpha",
                "beta",
                "T_bath",
                "post_filter_white_asd",
            ],
            "thermal_topology": (
                "TES -- G_tes-substrate -- C_substrate -- "
                "G_substrate-bath -- bath, in parallel with the existing "
                "TES--Stycast--Pb branch"
            ),
            "DC_constraint": (
                "series_equivalent(G_tes-substrate, G_substrate-bath) "
                "equals existing G_tes-bath exactly"
            ),
            "nested_readout_test": (
                "shared c2=0 versus one shared c2 free"
            ),
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
            "objective": (
                "equal-weight mean of the two day-specific continuum "
                "shape scores over 1-200 kHz"
            ),
        },
        "physical_profile": {
            "C_tes_material_J_per_K": material_C_tes,
            **config["shared_physical_profile"],
        },
        "c2_zero": {
            key: value
            for key, value in zero.items()
            if not key.startswith("_")
        },
        "c2_free": {
            key: value
            for key, value in free.items()
            if not key.startswith("_")
        },
        "nested_test": nested,
        "interpretation": {
            "classification": nested["classification"],
            "shared_c2_material": bool(
                nested["shared_c2_material_improvement"]
            ),
            "shared_c2": float(free["solution"]["shared_c2"]),
            "c2_zero_shared_local_C_total_J_per_K": float(
                zero["solution"]["shared_local_C_total_J_per_K"]
            ),
            "c2_free_shared_local_C_total_J_per_K": float(
                free["solution"]["shared_local_C_total_J_per_K"]
            ),
            "c2_zero_shared_C_substrate_J_per_K": float(
                zero["solution"]["shared_C_substrate_J_per_K"]
            ),
            "c2_free_shared_C_substrate_J_per_K": float(
                free["solution"]["shared_C_substrate_J_per_K"]
            ),
            "c2_zero_shared_C_tes_J_per_K": float(
                zero["solution"]["shared_C_tes_J_per_K"]
            ),
            "c2_free_shared_C_tes_J_per_K": float(
                free["solution"]["shared_C_tes_J_per_K"]
            ),
            "c2_zero_shared_L_H": float(
                zero["solution"]["shared_L_H"]
            ),
            "c2_free_shared_L_H": float(
                free["solution"]["shared_L_H"]
            ),
            "c2_zero_shared_G_ratio": float(
                zero["solution"][
                    "shared_G_tes_substrate_over_G_substrate_bath"
                ]
            ),
            "c2_free_shared_G_ratio": float(
                free["solution"][
                    "shared_G_tes_substrate_over_G_substrate_bath"
                ]
            ),
            "guardrail": str(config["guardrail"]),
        },
        "inputs": {
            "config": str(config_path),
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
            "reference_free": free["_models"]["reference"].tolist(),
            "repeat_zero": zero["_models"]["repeat"].tolist(),
            "repeat_free": free["_models"]["repeat"].tolist(),
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
        figsize=(13.5, 8.5),
        sharex="col",
        gridspec_kw={"height_ratios": [2.0, 1.0]},
    )

    specs = [
        (
            0,
            "2024-12-06 reference",
            np.asarray(p["reference_target"], dtype=float),
            np.asarray(p["reference_raw"], dtype=float),
            np.asarray(p["reference_zero"], dtype=float),
            np.asarray(p["reference_free"], dtype=float),
            result["c2_zero"]["solution"]["reference_day"]["metrics"],
            result["c2_free"]["solution"]["reference_day"]["metrics"],
        ),
        (
            1,
            "2024-12-05 repeat",
            np.asarray(p["repeat_target"], dtype=float),
            np.asarray(p["repeat_raw"], dtype=float),
            np.asarray(p["repeat_zero"], dtype=float),
            np.asarray(p["repeat_free"], dtype=float),
            result["c2_zero"]["solution"]["repeat_day"]["metrics"],
            result["c2_free"]["solution"]["repeat_day"]["metrics"],
        ),
    ]

    for col, label, target, raw, zero, free, m0, m1 in specs:
        top = axes[0, col]
        bottom = axes[1, col]
        top.loglog(f, raw, linewidth=0.8, alpha=0.4, label="Raw ASD")
        top.loglog(f, target, linewidth=1.4, label="Continuum target")
        top.loglog(f, zero, linewidth=1.4, label="Separated node, c2=0")
        top.loglog(f, free, linewidth=1.4, label="Separated node, c2 free")
        top.set_title(
            f"{label}\nRMS {m0['continuum_rms_dB']:.3f} -> "
            f"{m1['continuum_rms_dB']:.3f} dB",
            fontsize=10,
        )
        top.set_ylabel("Normalized ASD")
        top.grid(True, which="both", alpha=0.2)
        top.legend(frameon=False, fontsize=8)

        bottom.semilogx(
            f,
            _residual_db(zero, target),
            linewidth=1.2,
            label="c2=0",
        )
        bottom.semilogx(
            f,
            _residual_db(free, target),
            linewidth=1.2,
            label="c2 free",
        )
        bottom.axhline(0.0, linewidth=1.0)
        bottom.set_xlabel("Frequency [Hz]")
        bottom.set_ylabel("Model / continuum [dB]")
        bottom.grid(True, which="both", alpha=0.2)
        bottom.legend(frameon=False, fontsize=8)

    z = result["c2_zero"]["solution"]
    q = result["c2_free"]["solution"]
    fig.suptitle(
        result["interpretation"]["classification"]
        + "\n"
        + (
            "c2=0: "
            f"Ctes={z['shared_C_tes_J_per_K']:.3e}, "
            f"Csub={z['shared_C_substrate_J_per_K']:.3e}, "
            f"L={z['shared_L_H']:.3e}; "
            "c2 free: "
            f"c2={q['shared_c2']:.3g}, "
            f"Ctes={q['shared_C_tes_J_per_K']:.3e}, "
            f"Csub={q['shared_C_substrate_J_per_K']:.3e}, "
            f"L={q['shared_L_H']:.3e}"
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
    return {
        key: value
        for key, value in result.items()
        if key != "_plot"
    }


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

    z = result["c2_zero"]["solution"]
    q = result["c2_free"]["solution"]
    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "classification": result["interpretation"][
                    "classification"
                ],
                "shared_c2_material": result["interpretation"][
                    "shared_c2_material"
                ],
                "c2_zero": {
                    "C_tes_J_per_K": z["shared_C_tes_J_per_K"],
                    "C_substrate_J_per_K": z[
                        "shared_C_substrate_J_per_K"
                    ],
                    "C_total_J_per_K": z[
                        "shared_local_C_total_J_per_K"
                    ],
                    "C_total_over_historical_PoST": z[
                        "shared_local_C_total_over_historical_PoST"
                    ],
                    "G_ratio": z[
                        "shared_G_tes_substrate_over_G_substrate_bath"
                    ],
                    "L_H": z["shared_L_H"],
                    "joint_rms_dB": result["c2_zero"][
                        "joint_continuum_rms_dB"
                    ],
                },
                "c2_free": {
                    "shared_c2": q["shared_c2"],
                    "C_tes_J_per_K": q["shared_C_tes_J_per_K"],
                    "C_substrate_J_per_K": q[
                        "shared_C_substrate_J_per_K"
                    ],
                    "C_total_J_per_K": q[
                        "shared_local_C_total_J_per_K"
                    ],
                    "C_total_over_historical_PoST": q[
                        "shared_local_C_total_over_historical_PoST"
                    ],
                    "G_ratio": q[
                        "shared_G_tes_substrate_over_G_substrate_bath"
                    ],
                    "L_H": q["shared_L_H"],
                    "joint_rms_dB": result["c2_free"][
                        "joint_continuum_rms_dB"
                    ],
                },
                "joint_rms_improvement_dB": result["nested_test"][
                    "joint_rms_improvement_dB"
                ],
                "joint_shape_score_ratio": result["nested_test"][
                    "joint_shape_score_ratio_c2_free_to_c2_zero"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
