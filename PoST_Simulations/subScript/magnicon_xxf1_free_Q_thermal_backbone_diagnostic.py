"""Second-order Magnicon 10 kHz free-Q diagnostic.

The XXF-1 manual documents a second-order Bessel output filter.  This test
therefore keeps the filter order at two and asks whether the residual can be
absorbed by the damping/shape parameter Q rather than by c2 or by deleting one
pole.

For the documented phase-normalized second-order Bessel,

    H(s) = 1 / ((s/w0)^2 + s/(Q*w0) + 1)

has Q = 1/sqrt(3) and w0/(2*pi) = 10 kHz.  The free-Q branch profiles Q while
keeping the natural denominator frequency f0 within 10 kHz +/-2.5%.

Both branches use c2=0 and the same relaxed thermal-backbone profile used by
the preceding filter-order diagnostic.
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
    CONFIG_DIR / "magnicon_xxf1_free_Q_thermal_backbone_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR / "magnicon_xxf1_free_Q_thermal_backbone_diagnostic.json"
)
DEFAULT_FIGURE = (
    WORK_DIR / "magnicon_xxf1_free_Q_thermal_backbone_diagnostic.png"
)
REFERENCE_HZ = 1000.0
BRANCH_BESSEL = "documented_bessel_Q"
BRANCH_FREE_Q = "free_Q"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_shared_L_Ctes_joint_c2_diagnostic as shared  # noqa: E402
from subScript import magnicon_xxf1_thermal_backbone_pure_zero_diagnostic as thermal  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402


def second_order_magnitude(frequency_hz, f0_hz: float, q: float):
    frequency = np.asarray(frequency_hz, dtype=float)
    f0 = float(f0_hz)
    q = float(q)
    if f0 <= 0.0 or q <= 0.0:
        raise ValueError("f0 and Q must be positive")
    x = frequency / f0
    denominator_squared = (1.0 - x**2) ** 2 + (x / q) ** 2
    return 1.0 / np.sqrt(denominator_squared)


def minus3db_frequency(f0_hz: float, q: float) -> float:
    """Return the positive frequency where |H|=1/sqrt(2)."""
    f0 = float(f0_hz)
    q = float(q)
    # Let y=(f/f0)^2.  Solve
    # (1-y)^2 + y/Q^2 = 2.
    b = -2.0 + 1.0 / (q * q)
    disc = b * b + 4.0
    y = (-b + np.sqrt(disc)) / 2.0
    if y <= 0.0:
        raise ValueError("no positive -3 dB solution")
    return float(f0 * np.sqrt(y))


def model_candidate(
    detector_candidate: dict,
    *,
    frequency,
    white_asd: float,
    magnicon_f0_hz: float,
    magnicon_q: float,
):
    point = opt.tes_operating_point(detector_candidate)
    if not point.get("valid") or not point.get("stable"):
        return None, point

    context = competition.fixed_white_context(
        detector_candidate,
        frequency,
        white_asd,
    )
    frequency = np.asarray(context["frequency_Hz"], dtype=float)
    alias_frequency = np.asarray(context["alias_frequency_Hz"], dtype=float)

    magnicon_main = second_order_magnitude(
        frequency,
        magnicon_f0_hz,
        magnicon_q,
    )
    magnicon_alias = second_order_magnitude(
        alias_frequency,
        magnicon_f0_hz,
        magnicon_q,
    )
    sim965_main = opt.hardware_filter_magnitude(
        frequency,
        cutoff_hz=opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=opt.TARGET_HARDWARE_BESSEL_ORDER,
        norm=opt.TARGET_HARDWARE_BESSEL_NORM,
    )
    sim965_alias = opt.hardware_filter_magnitude(
        alias_frequency,
        cutoff_hz=opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=opt.TARGET_HARDWARE_BESSEL_ORDER,
        norm=opt.TARGET_HARDWARE_BESSEL_NORM,
    )

    main = (
        np.asarray(context["main_intrinsic_asd"], dtype=float)
        * magnicon_main
        * sim965_main
    )
    alias = (
        np.asarray(context["alias_intrinsic_asd"], dtype=float)
        * magnicon_alias
        * sim965_alias
    )
    alias = np.where(context["same_bin"], 0.0, alias)
    white = float(context["post_filter_white_asd_A_rtHz"])
    absolute = np.sqrt(main**2 + alias**2 + white**2)
    model = opt.normalize_at(
        frequency,
        absolute,
        reference_hz=REFERENCE_HZ,
    )
    if np.any(~np.isfinite(model)) or np.any(model <= 0.0):
        return None, point
    return model, point


def build_spec(
    *,
    branch: str,
    config: dict,
    thermal_config: dict,
    separated_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    if branch not in {BRANCH_BESSEL, BRANCH_FREE_Q}:
        raise ValueError(f"unknown branch: {branch}")
    base_spec = thermal.build_spec(
        branch=thermal.BRANCH_ZERO,
        config=thermal_config,
        base_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    f0_cfg = config["magnicon_filter"]["natural_frequency_Hz"]
    f0_bounds = (float(f0_cfg["min"]), float(f0_cfg["max"]))
    names = list(base_spec["names"])
    bounds = list(base_spec["bounds"])
    names.append("shared_log10_magnicon_f0_Hz")
    bounds.append((np.log10(f0_bounds[0]), np.log10(f0_bounds[1])))

    q_bounds = None
    if branch == BRANCH_FREE_Q:
        q_cfg = config["magnicon_filter"]["free_Q"]
        q_bounds = (float(q_cfg["min"]), float(q_cfg["max"]))
        names.append("shared_log10_magnicon_Q")
        bounds.append((np.log10(q_bounds[0]), np.log10(q_bounds[1])))

    return {
        "branch": branch,
        "base_spec": base_spec,
        "names": tuple(names),
        "bounds": tuple(bounds),
        "f0_bounds": f0_bounds,
        "q_bounds": q_bounds,
    }


def decode_vector(
    vector,
    *,
    spec: dict,
    config: dict,
    thermal_config: dict,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    vector = np.asarray(vector, dtype=float)
    if vector.size != len(spec["names"]):
        raise ValueError("free-Q vector length mismatch")

    if spec["branch"] == BRANCH_FREE_Q:
        base_vector = vector[:-2]
        f0 = float(10.0 ** vector[-2])
        q = float(10.0 ** vector[-1])
    else:
        base_vector = vector[:-1]
        f0 = float(10.0 ** vector[-1])
        q = float(config["magnicon_filter"]["documented_bessel_Q"])

    dummy_readout = {
        "pole_Hz": 1.0,
        "pole_Q": 1.0,
        "c2": 0.0,
        "c4": 0.0,
    }
    decoded = thermal.decode_vector(
        base_vector,
        spec=spec["base_spec"],
        config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=dummy_readout,
    )
    decoded = dict(decoded)
    decoded["magnicon_f0_Hz"] = f0
    decoded["magnicon_Q"] = q
    return decoded


def midpoint_warm(
    *,
    spec: dict,
    config: dict,
    thermal_config: dict,
    reference_problem: dict,
    repeat_problem: dict,
) -> np.ndarray:
    base = thermal.warm_vector(
        spec=spec["base_spec"],
        config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        prior_solution=None,
    )
    f0 = float(config["magnicon_filter"]["natural_frequency_Hz"]["nominal"])
    values = list(base) + [np.log10(f0)]
    if spec["branch"] == BRANCH_FREE_Q:
        values.append(
            np.log10(float(config["magnicon_filter"]["documented_bessel_Q"]))
        )
    return np.asarray(values, dtype=float)


def fit_branch(
    *,
    branch: str,
    config: dict,
    thermal_config: dict,
    separated_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
    warm_from: dict | None = None,
) -> dict:
    spec = build_spec(
        branch=branch,
        config=config,
        thermal_config=thermal_config,
        separated_config=separated_config,
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
            thermal_config=thermal_config,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
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
                model, point = model_candidate(
                    decoded["detectors"][prefix],
                    frequency=problem["frequency"],
                    white_asd=white_asd,
                    magnicon_f0_hz=decoded["magnicon_f0_Hz"],
                    magnicon_q=decoded["magnicon_Q"],
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
        return float(
            np.mean(
                [
                    opt.fit_score(
                        models[prefix],
                        problem["target"],
                        problem["frequency"],
                        problem["fit_args"],
                    )
                    for prefix, problem in (
                        ("reference", reference_problem),
                        ("repeat", repeat_problem),
                    )
                ]
            )
        )

    def ls_residual(vector):
        models, _, _, _ = evaluate(vector)
        if models is None:
            return invalid_residual
        return np.concatenate(
            [
                opt.weighted_residual_vector(
                    models[prefix],
                    problem["target"],
                    problem["frequency"],
                    problem["fit_args"],
                )
                for prefix, problem in (
                    ("reference", reference_problem),
                    ("repeat", repeat_problem),
                )
            ]
        ) / np.sqrt(2.0)

    optimizer = config["optimizer"]
    seed = (
        int(reference_problem["optimizer_cfg"]["seed"])
        + int(optimizer["seed_offset"])
        + (1000 if branch == BRANCH_FREE_Q else 0)
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

    midpoint = midpoint_warm(
        spec=spec,
        config=config,
        thermal_config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    candidates.append(("physical_midpoint_warm", midpoint))

    if warm_from is not None and branch == BRANCH_FREE_Q:
        parent = np.asarray(warm_from["_vector"], dtype=float)
        q0 = np.log10(float(config["magnicon_filter"]["documented_bessel_Q"]))
        inherited = np.concatenate((parent, [q0]))
        if inherited.size == len(spec["names"]):
            candidates.append(("documented_Q_exact_warm", inherited))

    for label, vector in list(candidates):
        if np.any(vector < lower) or np.any(vector > upper):
            continue
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
        if np.any(vector < lower) or np.any(vector > upper):
            continue
        score = objective(vector)
        if score < penalty:
            scored.append((float(score), label, vector))
    if not scored:
        raise RuntimeError(f"no stable solution for {branch}")

    score, source, vector = min(scored, key=lambda row: row[0])
    models, decoded, point_ref, point_rep = evaluate(vector)
    if models is None:
        raise RuntimeError("selected free-Q solution became invalid")

    ref_metrics = shared._day_metrics(reference_problem, models["reference"])
    rep_metrics = shared._day_metrics(repeat_problem, models["repeat"])
    combined_rms = shared._combined_rms(ref_metrics, rep_metrics)

    base_spec = spec["base_spec"]
    shared_values = decoded["shared"]
    shared_hits = {
        name: shared._encoded_boundary(
            shared_values[name],
            *base_spec["shared_bounds"][name],
            log=True,
        )
        for name in base_spec["shared_names"]
    }
    f0_hit = shared._encoded_boundary(
        decoded["magnicon_f0_Hz"],
        *spec["f0_bounds"],
        log=True,
    )
    q_hit = None
    if branch == BRANCH_FREE_Q:
        q_hit = shared._encoded_boundary(
            decoded["magnicon_Q"],
            *spec["q_bounds"],
            log=True,
        )
    white_hits = {
        prefix: shared._encoded_boundary(
            decoded["white_scales"][prefix],
            *base_spec["white_bounds"][prefix],
            log=True,
        )
        for prefix in ("reference", "repeat")
    }

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
                    *base_spec["day_bounds"][prefix][name],
                    log=False,
                )
                for name in thermal.DAY_PARAMETERS
            },
            "white_scale": float(decoded["white_scales"][prefix]),
            "white_asd_A_rtHz": float(
                problem["tracked_white"] * decoded["white_scales"][prefix]
            ),
            "white_boundary_hit": white_hits[prefix],
            "operating_point": {
                "valid": bool(point.get("valid")),
                "stable": bool(point.get("stable")),
                "f_substrate_node_Hz": (
                    float(point["f_substrate_node_Hz"])
                    if point.get("f_substrate_node_Hz") is not None
                    else None
                ),
            },
            "metrics": metrics,
        }

    historical = float(
        separated_config["shared_physical_profile"]["C_substrate_J_per_K"][
            "historical_PoST_total_reference_J_per_K"
        ]
    )
    local_total = shared_values["C_tes"] + shared_values["C_substrate"]
    solution = {
        "magnicon_order": 2,
        "magnicon_f0_Hz": float(decoded["magnicon_f0_Hz"]),
        "magnicon_Q": float(decoded["magnicon_Q"]),
        "magnicon_minus3dB_Hz": minus3db_frequency(
            decoded["magnicon_f0_Hz"],
            decoded["magnicon_Q"],
        ),
        "magnicon_f0_boundary_hit": f0_hit,
        "magnicon_Q_boundary_hit": q_hit,
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
        "shared_physical_boundary_hits": shared_hits,
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
        "_vector": vector,
    }


def compare_branches(fixed: dict, free: dict, config: dict) -> dict:
    screen = config["materiality_screen"]
    score_ratio = float(
        free["joint_shape_score"] / fixed["joint_shape_score"]
    )
    rms_improvement = float(
        fixed["joint_continuum_rms_dB"] - free["joint_continuum_rms_dB"]
    )
    material = bool(
        score_ratio
        <= float(screen["max_freeQ_over_bessel_score_ratio"])
        and rms_improvement
        >= float(screen["min_freeQ_rms_improvement_dB"])
    )
    good = bool(
        free["joint_continuum_rms_dB"]
        <= float(screen["target_good_fit_rms_dB"])
    )
    q_hit = free["solution"]["magnicon_Q_boundary_hit"]
    q_censored = bool(
        q_hit and (q_hit.get("at_lower") or q_hit.get("at_upper"))
    )

    if material and good and not q_censored:
        classification = "free_Q_reaches_good_fit_with_interior_Q"
    elif material and q_censored:
        classification = "free_Q_material_but_Q_boundary_censored"
    elif material:
        classification = "free_Q_materially_improves_but_does_not_fully_close_residual"
    else:
        classification = "free_Q_does_not_materially_resolve_residual"

    return {
        "classification": classification,
        "free_Q_materially_better": material,
        "free_Q_reaches_target_good_fit": good,
        "free_Q_boundary_censored": q_censored,
        "freeQ_over_bessel_shape_score_ratio": score_ratio,
        "free_Q_rms_improvement_dB": rms_improvement,
        "screen": screen,
    }


def transfer_comparison(fixed: dict, free: dict) -> dict:
    frequencies = np.asarray(
        [2e3, 5e3, 10e3, 20e3, 40e3, 100e3],
        dtype=float,
    )
    fixed_mag = second_order_magnitude(
        frequencies,
        fixed["solution"]["magnicon_f0_Hz"],
        fixed["solution"]["magnicon_Q"],
    )
    free_mag = second_order_magnitude(
        frequencies,
        free["solution"]["magnicon_f0_Hz"],
        free["solution"]["magnicon_Q"],
    )
    fixed_ref = second_order_magnitude(
        np.asarray([REFERENCE_HZ]),
        fixed["solution"]["magnicon_f0_Hz"],
        fixed["solution"]["magnicon_Q"],
    )[0]
    free_ref = second_order_magnitude(
        np.asarray([REFERENCE_HZ]),
        free["solution"]["magnicon_f0_Hz"],
        free["solution"]["magnicon_Q"],
    )[0]
    ratio_db = 20.0 * np.log10(
        (free_mag / free_ref) / (fixed_mag / fixed_ref)
    )
    return {
        "reference_Hz": REFERENCE_HZ,
        "frequency_Hz": frequencies.tolist(),
        "free_Q_over_documented_bessel_dB": ratio_db.tolist(),
    }


def run(config: dict, config_path: Path) -> dict:
    thermal_path = thermal.resolve_config_path(
        config["base_thermal_backbone_config"],
        config_path,
    )
    thermal_config = json.loads(thermal_path.read_text(encoding="utf-8"))
    separated_path = thermal.resolve_config_path(
        thermal_config["base_separated_substrate_config"],
        thermal_path,
    )
    separated_config = json.loads(
        separated_path.read_text(encoding="utf-8")
    )
    stack = compensation._load_stack(separated_config, separated_path)
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

    material_C_tes = float(opt.C_TES_MATERIAL_J_PER_K)
    fixed = fit_branch(
        branch=BRANCH_BESSEL,
        config=config,
        thermal_config=thermal_config,
        separated_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    free = fit_branch(
        branch=BRANCH_FREE_Q,
        config=config,
        thermal_config=thermal_config,
        separated_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        warm_from=fixed,
    )
    comparison = compare_branches(fixed, free, config)

    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "With Magnicon order fixed at two and c2 fixed at zero, can a "
            "non-Bessel damping Q within a broad diagnostic range absorb the "
            "continuum residual under the same relaxed thermal backbone?"
        ),
        "manual_constraint": {
            "documented_order": int(
                config["magnicon_filter"]["documented_order"]
            ),
            "documented_bessel_Q": float(
                config["magnicon_filter"]["documented_bessel_Q"]
            ),
            "statement": str(config["magnicon_filter"]["manual_statement"]),
            "natural_frequency_profile_Hz": [
                float(config["magnicon_filter"]["natural_frequency_Hz"]["min"]),
                float(config["magnicon_filter"]["natural_frequency_Hz"]["max"]),
            ],
            "free_Q_profile": [
                float(config["magnicon_filter"]["free_Q"]["min"]),
                float(config["magnicon_filter"]["free_Q"]["max"]),
            ],
        },
        "fit_semantics": {
            "c2": 0.0,
            "filter_order": 2,
            "shared_thermal_parameters": [
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
            "external_analog_filter": {
                "model": "SRS SIM965",
                "front_panel_cutoff_Hz": float(
                    opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ
                ),
                "order": int(opt.TARGET_HARDWARE_BESSEL_ORDER),
                "scipy_norm": str(opt.TARGET_HARDWARE_BESSEL_NORM),
            },
        },
        "documented_bessel_Q": {
            key: value
            for key, value in fixed.items()
            if not key.startswith("_")
        },
        "free_Q": {
            key: value
            for key, value in free.items()
            if not key.startswith("_")
        },
        "comparison": comparison,
        "filter_transfer_comparison": transfer_comparison(fixed, free),
        "interpretation": {
            "classification": comparison["classification"],
            "guardrail": str(config["guardrail"]),
        },
        "inputs": {
            "config": str(config_path),
            "base_thermal_backbone_config": str(thermal_path),
            "base_separated_substrate_config": str(separated_path),
            "base_cross_day_config": str(stack["base_cross_path"]),
            "base_magnicon_config": str(stack["magnicon_config_path"]),
        },
        "_plot": {
            "frequency_Hz": reference_problem["frequency"].tolist(),
            "reference_target": reference_problem["target"].tolist(),
            "repeat_target": repeat_problem["target"].tolist(),
            "reference_raw": reference_problem["raw_target"].tolist(),
            "repeat_raw": repeat_problem["raw_target"].tolist(),
            "reference_fixed": fixed["_models"]["reference"].tolist(),
            "repeat_fixed": fixed["_models"]["repeat"].tolist(),
            "reference_free": free["_models"]["reference"].tolist(),
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
            np.asarray(p["reference_fixed"], dtype=float),
            np.asarray(p["reference_free"], dtype=float),
        ),
        (
            1,
            "2024-12-05 repeat",
            np.asarray(p["repeat_target"], dtype=float),
            np.asarray(p["repeat_raw"], dtype=float),
            np.asarray(p["repeat_fixed"], dtype=float),
            np.asarray(p["repeat_free"], dtype=float),
        ),
    ]
    for col, label, target, raw, fixed, free in rows:
        top = axes[0, col]
        bottom = axes[1, col]
        top.loglog(f, raw, linewidth=0.7, alpha=0.35, label="Raw ASD")
        top.loglog(f, target, linewidth=1.5, label="Continuum target")
        top.loglog(f, fixed, linewidth=1.3, label="2nd-order Bessel Q, c2=0")
        top.loglog(f, free, linewidth=1.4, label="2nd-order free Q, c2=0")
        top.set_title(label, fontsize=10)
        top.set_ylabel("Normalized ASD")
        top.grid(True, which="both", alpha=0.2)
        top.legend(frameon=False, fontsize=7)

        bottom.semilogx(
            f, _residual_db(fixed, target), linewidth=1.2, label="Bessel Q"
        )
        bottom.semilogx(
            f, _residual_db(free, target), linewidth=1.2, label="free Q"
        )
        bottom.axhline(0.0, linewidth=0.9)
        bottom.set_xlabel("Frequency [Hz]")
        bottom.set_ylabel("Model / continuum [dB]")
        bottom.grid(True, which="both", alpha=0.2)
        bottom.legend(frameon=False, fontsize=7)

    fixed = result["documented_bessel_Q"]
    free = result["free_Q"]
    fig.suptitle(
        result["interpretation"]["classification"]
        + "\n"
        + (
            f"Bessel: f0={fixed['solution']['magnicon_f0_Hz']:.1f} Hz, "
            f"Q={fixed['solution']['magnicon_Q']:.4g}, "
            f"RMS={fixed['joint_continuum_rms_dB']:.3f} dB; "
            f"free: f0={free['solution']['magnicon_f0_Hz']:.1f} Hz, "
            f"Q={free['solution']['magnicon_Q']:.4g}, "
            f"RMS={free['joint_continuum_rms_dB']:.3f} dB"
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
                "documented_bessel_Q": {
                    "joint_rms_dB": result["documented_bessel_Q"][
                        "joint_continuum_rms_dB"
                    ],
                    "f0_Hz": result["documented_bessel_Q"]["solution"][
                        "magnicon_f0_Hz"
                    ],
                    "Q": result["documented_bessel_Q"]["solution"]["magnicon_Q"],
                },
                "free_Q": {
                    "joint_rms_dB": result["free_Q"]["joint_continuum_rms_dB"],
                    "f0_Hz": result["free_Q"]["solution"]["magnicon_f0_Hz"],
                    "Q": result["free_Q"]["solution"]["magnicon_Q"],
                    "minus3dB_Hz": result["free_Q"]["solution"][
                        "magnicon_minus3dB_Hz"
                    ],
                },
                "comparison": result["comparison"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
