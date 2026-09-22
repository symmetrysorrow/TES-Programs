"""Counterfactual Magnicon 10 kHz filter-order test.

The XXF-1 manual explicitly documents the connector-box anti-alias filters as
second-order Bessel, 10 kHz +/-2.5%.  This diagnostic does not dispute that
documentation.  It asks a narrower residual-shape question:

    If the measured chain attenuated by one pole less than the present
    second-order model, could the c2/zero-like residual disappear?

Both branches keep c2=0 and use exactly the same relaxed thermal-backbone
profile:
  * shared L, C_tes, C_substrate, substrate conductance split;
  * broad shared effective G_tes-bath scale;
  * shared Pb thickness 0.4--0.6 mm, with C_abs and G_abs-abs recomputed from
    the current Elmer Pb material table;
  * day-specific alpha, beta, T_bath, and white ASD scale.

The only branch distinction is Magnicon Bessel order 1 versus 2.  Each branch
profiles one shared cutoff over the documented 9.75--10.25 kHz tolerance.
The external SRS SIM965 model remains unchanged.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.optimize import differential_evolution, least_squares

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
WORK_DIR = ROOT / ".noise_optimization_work_rsh_sweep"

DEFAULT_CONFIG = (
    CONFIG_DIR
    / "magnicon_xxf1_filter_order_thermal_backbone_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR
    / "magnicon_xxf1_filter_order_thermal_backbone_diagnostic.json"
)
DEFAULT_FIGURE = (
    WORK_DIR
    / "magnicon_xxf1_filter_order_thermal_backbone_diagnostic.png"
)
REFERENCE_HZ = 1000.0

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


def analog_bessel_magnitude(
    frequency_hz,
    *,
    order: int,
    cutoff_hz: float,
    norm: str,
):
    """Return the analog low-pass Bessel magnitude for one Magnicon branch."""
    frequency = np.asarray(frequency_hz, dtype=float)
    order = int(order)
    cutoff = float(cutoff_hz)
    if order < 1:
        raise ValueError("Bessel order must be >= 1")
    if cutoff <= 0.0:
        raise ValueError("cutoff_hz must be positive")
    if norm not in {"mag", "phase", "delay"}:
        raise ValueError("unsupported Bessel normalization")
    if np.any(frequency < 0.0):
        raise ValueError("frequency must be non-negative")

    b, a = signal.bessel(
        order,
        2.0 * np.pi * cutoff,
        btype="low",
        analog=True,
        output="ba",
        norm=norm,
    )
    _, response = signal.freqs(
        np.asarray(b, dtype=float),
        np.asarray(a, dtype=float),
        worN=2.0 * np.pi * frequency,
    )
    return np.abs(response)


def normalized_magnicon_magnitude(
    frequency_hz,
    *,
    order: int,
    cutoff_hz: float,
    norm: str,
    reference_hz: float = REFERENCE_HZ,
):
    frequency = np.asarray(frequency_hz, dtype=float)
    response = analog_bessel_magnitude(
        frequency,
        order=order,
        cutoff_hz=cutoff_hz,
        norm=norm,
    )
    reference = float(
        analog_bessel_magnitude(
            np.asarray([reference_hz], dtype=float),
            order=order,
            cutoff_hz=cutoff_hz,
            norm=norm,
        )[0]
    )
    return response / reference


def model_candidate(
    detector_candidate: dict,
    *,
    frequency,
    white_asd: float,
    magnicon_order: int,
    magnicon_cutoff_hz: float,
    magnicon_norm: str,
):
    """Expected pre-analysis normalized ASD for one filter-order branch."""
    point = opt.tes_operating_point(detector_candidate)
    if not point.get("valid") or not point.get("stable"):
        return None, point

    context = competition.fixed_white_context(
        detector_candidate,
        frequency,
        white_asd,
    )
    frequency = np.asarray(context["frequency_Hz"], dtype=float)
    alias_frequency = np.asarray(
        context["alias_frequency_Hz"],
        dtype=float,
    )

    hardware_main = opt.hardware_filter_magnitude(
        frequency,
        cutoff_hz=opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=opt.TARGET_HARDWARE_BESSEL_ORDER,
        norm=opt.TARGET_HARDWARE_BESSEL_NORM,
    )
    hardware_alias = opt.hardware_filter_magnitude(
        alias_frequency,
        cutoff_hz=opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ,
        order=opt.TARGET_HARDWARE_BESSEL_ORDER,
        norm=opt.TARGET_HARDWARE_BESSEL_NORM,
    )
    magnicon_main = analog_bessel_magnitude(
        frequency,
        order=magnicon_order,
        cutoff_hz=magnicon_cutoff_hz,
        norm=magnicon_norm,
    )
    magnicon_alias = analog_bessel_magnitude(
        alias_frequency,
        order=magnicon_order,
        cutoff_hz=magnicon_cutoff_hz,
        norm=magnicon_norm,
    )

    main = (
        np.asarray(context["main_intrinsic_asd"], dtype=float)
        * magnicon_main
        * hardware_main
    )
    alias = (
        np.asarray(context["alias_intrinsic_asd"], dtype=float)
        * magnicon_alias
        * hardware_alias
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
    config: dict,
    thermal_config: dict,
    separated_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    base_spec = thermal.build_spec(
        branch=thermal.BRANCH_ZERO,
        config=thermal_config,
        base_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    cutoff_cfg = config["magnicon_filter"]["cutoff_Hz"]
    cutoff_bounds = (
        float(cutoff_cfg["min"]),
        float(cutoff_cfg["max"]),
    )
    if not 0.0 < cutoff_bounds[0] < cutoff_bounds[1]:
        raise ValueError("invalid Magnicon cutoff bounds")

    return {
        "base_spec": base_spec,
        "names": tuple(base_spec["names"])
        + ("shared_log10_magnicon_cutoff_Hz",),
        "bounds": tuple(base_spec["bounds"])
        + (
            (
                np.log10(cutoff_bounds[0]),
                np.log10(cutoff_bounds[1]),
            ),
        ),
        "cutoff_bounds": cutoff_bounds,
    }


def decode_vector(
    vector,
    *,
    spec: dict,
    thermal_config: dict,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    vector = np.asarray(vector, dtype=float)
    if vector.size != len(spec["names"]):
        raise ValueError("filter-order vector length mismatch")

    dummy_readout = {
        "pole_Hz": 1.0,
        "pole_Q": 1.0,
        "c2": 0.0,
        "c4": 0.0,
    }
    decoded = thermal.decode_vector(
        vector[:-1],
        spec=spec["base_spec"],
        config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=dummy_readout,
    )
    decoded = dict(decoded)
    decoded["magnicon_cutoff_Hz"] = float(10.0 ** vector[-1])
    return decoded


def warm_vector(
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
    cutoff = float(config["magnicon_filter"]["cutoff_Hz"]["nominal"])
    return np.concatenate((base, [np.log10(cutoff)]))


def fit_order(
    *,
    order: int,
    config: dict,
    thermal_config: dict,
    separated_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    spec = build_spec(
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
    norm = str(config["magnicon_filter"]["normalization"])

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
                    magnicon_order=order,
                    magnicon_cutoff_hz=decoded["magnicon_cutoff_Hz"],
                    magnicon_norm=norm,
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
        pieces = []
        for prefix, problem in (
            ("reference", reference_problem),
            ("repeat", repeat_problem),
        ):
            pieces.append(
                opt.weighted_residual_vector(
                    models[prefix],
                    problem["target"],
                    problem["frequency"],
                    problem["fit_args"],
                )
            )
        return np.concatenate(pieces) / np.sqrt(2.0)

    optimizer = config["optimizer"]
    seed = (
        int(reference_problem["optimizer_cfg"]["seed"])
        + int(optimizer["seed_offset"])
        + 1000 * int(order)
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

    warm = warm_vector(
        spec=spec,
        config=config,
        thermal_config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
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
        raise RuntimeError(f"no stable solution for Magnicon order {order}")

    score, source, vector = min(scored, key=lambda row: row[0])
    models, decoded, point_ref, point_rep = evaluate(vector)
    if models is None:
        raise RuntimeError("selected filter-order solution became invalid")

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
    cutoff_hit = shared._encoded_boundary(
        decoded["magnicon_cutoff_Hz"],
        *spec["cutoff_bounds"],
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

    historical = float(
        separated_config["shared_physical_profile"]["C_substrate_J_per_K"][
            "historical_PoST_total_reference_J_per_K"
        ]
    )
    local_total = shared_values["C_tes"] + shared_values["C_substrate"]
    solution = {
        "magnicon_order": int(order),
        "magnicon_cutoff_Hz": float(decoded["magnicon_cutoff_Hz"]),
        "magnicon_normalization": norm,
        "magnicon_cutoff_boundary_hit": cutoff_hit,
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
        "name": f"magnicon_order_{int(order)}_c2_zero",
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


def transfer_ratio_diagnostics(
    first: dict,
    second: dict,
    *,
    norm: str,
) -> dict:
    frequencies = np.asarray(
        [2_000.0, 5_000.0, 10_000.0, 20_000.0, 40_000.0, 100_000.0],
        dtype=float,
    )
    first_mag = normalized_magnicon_magnitude(
        frequencies,
        order=1,
        cutoff_hz=first["solution"]["magnicon_cutoff_Hz"],
        norm=norm,
    )
    second_mag = normalized_magnicon_magnitude(
        frequencies,
        order=2,
        cutoff_hz=second["solution"]["magnicon_cutoff_Hz"],
        norm=norm,
    )
    ratio_db = 20.0 * np.log10(first_mag / second_mag)
    return {
        "reference_Hz": REFERENCE_HZ,
        "frequency_Hz": frequencies.tolist(),
        "first_order_over_second_order_dB": ratio_db.tolist(),
    }


def compare_orders(first: dict, second: dict, config: dict) -> dict:
    screen = config["materiality_screen"]
    score_ratio = float(
        first["joint_shape_score"] / second["joint_shape_score"]
    )
    rms_improvement = float(
        second["joint_continuum_rms_dB"]
        - first["joint_continuum_rms_dB"]
    )
    first_materially_better = bool(
        score_ratio <= float(screen["max_first_over_second_score_ratio"])
        and rms_improvement
        >= float(screen["min_first_order_rms_improvement_dB"])
    )
    first_good = bool(
        first["joint_continuum_rms_dB"]
        <= float(screen["target_good_fit_rms_dB"])
    )

    if first_materially_better and first_good:
        classification = (
            "counterfactual_first_order_reaches_good_fit_and_beats_documented_second_order"
        )
    elif first_materially_better:
        classification = (
            "counterfactual_first_order_materially_improves_but_does_not_fully_close_residual"
        )
    else:
        classification = (
            "counterfactual_first_order_does_not_materially_resolve_residual"
        )

    return {
        "classification": classification,
        "first_order_materially_better": first_materially_better,
        "first_order_reaches_target_good_fit": first_good,
        "first_over_second_shape_score_ratio": score_ratio,
        "first_order_rms_improvement_dB": rms_improvement,
        "screen": screen,
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
    orders = tuple(int(value) for value in config["orders"])
    if orders != (1, 2):
        raise ValueError("this diagnostic expects orders [1, 2]")

    first = fit_order(
        order=1,
        config=config,
        thermal_config=thermal_config,
        separated_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    second = fit_order(
        order=2,
        config=config,
        thermal_config=thermal_config,
        separated_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    comparison = compare_orders(first, second, config)
    norm = str(config["magnicon_filter"]["normalization"])
    transfer = transfer_ratio_diagnostics(first, second, norm=norm)

    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "With c2 fixed to zero and the same relaxed thermal-backbone "
            "profile, does replacing the documented second-order Magnicon "
            "10 kHz Bessel response by a counterfactual first-order response "
            "remove the continuum mismatch?"
        ),
        "manual_constraint": {
            "documented_order": int(
                config["magnicon_filter"]["documented_order"]
            ),
            "statement": str(
                config["magnicon_filter"]["manual_statement"]
            ),
            "cutoff_profile_Hz": [
                float(config["magnicon_filter"]["cutoff_Hz"]["min"]),
                float(config["magnicon_filter"]["cutoff_Hz"]["max"]),
            ],
            "normalization_used_for_isolated_order_test": norm,
        },
        "fit_semantics": {
            "c2": 0.0,
            "same_number_of_free_parameters_per_order": bool(
                first["n_free_parameters"] == second["n_free_parameters"]
            ),
            "shared_thermal_parameters": [
                "C_tes",
                "L",
                "C_substrate",
                "G_tes-substrate/G_substrate-bath",
                "G_tes-bath scale to inherited snapshot",
                "Pb absorber thickness",
            ],
            "shared_filter_parameter": "Magnicon cutoff within 10 kHz +/-2.5%",
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
        "counterfactual_first_order": {
            key: value
            for key, value in first.items()
            if not key.startswith("_")
        },
        "documented_second_order": {
            key: value
            for key, value in second.items()
            if not key.startswith("_")
        },
        "comparison": comparison,
        "filter_transfer_comparison": transfer,
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
            "reference_order1": first["_models"]["reference"].tolist(),
            "repeat_order1": first["_models"]["repeat"].tolist(),
            "reference_order2": second["_models"]["reference"].tolist(),
            "repeat_order2": second["_models"]["repeat"].tolist(),
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
            np.asarray(p["reference_order1"], dtype=float),
            np.asarray(p["reference_order2"], dtype=float),
        ),
        (
            1,
            "2024-12-05 repeat",
            np.asarray(p["repeat_target"], dtype=float),
            np.asarray(p["repeat_raw"], dtype=float),
            np.asarray(p["repeat_order1"], dtype=float),
            np.asarray(p["repeat_order2"], dtype=float),
        ),
    ]

    for col, label, target, raw, order1, order2 in rows:
        top = axes[0, col]
        bottom = axes[1, col]
        top.loglog(f, raw, linewidth=0.7, alpha=0.35, label="Raw ASD")
        top.loglog(f, target, linewidth=1.5, label="Continuum target")
        top.loglog(
            f,
            order1,
            linewidth=1.4,
            label="counterfactual Magnicon order 1, c2=0",
        )
        top.loglog(
            f,
            order2,
            linewidth=1.4,
            label="documented Magnicon order 2, c2=0",
        )
        top.set_title(label, fontsize=10)
        top.set_ylabel("Normalized ASD")
        top.grid(True, which="both", alpha=0.2)
        top.legend(frameon=False, fontsize=7)

        bottom.semilogx(
            f,
            _residual_db(order1, target),
            linewidth=1.2,
            label="order 1",
        )
        bottom.semilogx(
            f,
            _residual_db(order2, target),
            linewidth=1.2,
            label="order 2",
        )
        bottom.axhline(0.0, linewidth=0.9)
        bottom.set_xlabel("Frequency [Hz]")
        bottom.set_ylabel("Model / continuum [dB]")
        bottom.grid(True, which="both", alpha=0.2)
        bottom.legend(frameon=False, fontsize=7)

    first = result["counterfactual_first_order"]
    second = result["documented_second_order"]
    fig.suptitle(
        result["interpretation"]["classification"]
        + "\n"
        + (
            f"order1: fc={first['solution']['magnicon_cutoff_Hz']:.1f} Hz, "
            f"RMS={first['joint_continuum_rms_dB']:.3f} dB; "
            f"order2: fc={second['solution']['magnicon_cutoff_Hz']:.1f} Hz, "
            f"RMS={second['joint_continuum_rms_dB']:.3f} dB"
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
                "counterfactual_first_order": {
                    "joint_rms_dB": result["counterfactual_first_order"][
                        "joint_continuum_rms_dB"
                    ],
                    "joint_shape_score": result["counterfactual_first_order"][
                        "joint_shape_score"
                    ],
                    "cutoff_Hz": result["counterfactual_first_order"][
                        "solution"
                    ]["magnicon_cutoff_Hz"],
                },
                "documented_second_order": {
                    "joint_rms_dB": result["documented_second_order"][
                        "joint_continuum_rms_dB"
                    ],
                    "joint_shape_score": result["documented_second_order"][
                        "joint_shape_score"
                    ],
                    "cutoff_Hz": result["documented_second_order"][
                        "solution"
                    ]["magnicon_cutoff_Hz"],
                },
                "comparison": result["comparison"],
                "filter_transfer_comparison": result[
                    "filter_transfer_comparison"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
