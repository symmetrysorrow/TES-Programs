"""Replace phenomenological c2 with concrete source-class candidates.

This diagnostic keeps the separated TES/substrate thermal topology and the
same adjacent-day joint constraints, then compares four branches:

  1. shared c2 = 0 baseline;
  2. shared c2 free reference;
  3. c2 fixed to zero + one finite real lead/lag readout pair;
  4. c2 fixed to zero + one shared additive electronics-noise ASD injected
     downstream of the Magnicon 10 kHz LPF but upstream of the known 100 kHz
     hardware Bessel and first-alias fold.

The two candidate families are deliberately narrow.  A good fit is evidence
only that the source class can reproduce the required shape; it is not a
component identification.
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
    / "magnicon_xxf1_c2_physical_replacement_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR
    / "magnicon_xxf1_c2_physical_replacement_diagnostic.json"
)
DEFAULT_FIGURE = (
    WORK_DIR
    / "magnicon_xxf1_c2_physical_replacement_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_separated_substrate_joint_c2_diagnostic as separated  # noqa: E402
from subScript import magnicon_xxf1_shared_L_Ctes_joint_c2_diagnostic as shared  # noqa: E402
from subScript import readout_detector_state_competition_diagnostic as competition  # noqa: E402


REFERENCE_HZ = 1_000.0
FAMILY_LEAD_LAG = "real_lead_lag"
FAMILY_BYPASS_NOISE = "bypass_electronics_noise"
FAMILIES = (FAMILY_LEAD_LAG, FAMILY_BYPASS_NOISE)


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def magnicon_bessel_magnitude(
    frequency_hz,
    pole_hz: float,
    pole_q: float,
):
    frequency = np.asarray(frequency_hz, dtype=float)
    pole_hz = float(pole_hz)
    pole_q = float(pole_q)
    if pole_hz <= 0.0 or pole_q <= 0.0:
        raise ValueError("Magnicon pole and Q must be positive")
    x = frequency / pole_hz
    denominator_squared = (1.0 - x**2) ** 2 + (x / pole_q) ** 2
    return 1.0 / np.sqrt(denominator_squared)


def lead_lag_magnitude(
    frequency_hz,
    zero_hz: float,
    pole_hz: float,
):
    frequency = np.asarray(frequency_hz, dtype=float)
    zero_hz = float(zero_hz)
    pole_hz = float(pole_hz)
    if zero_hz <= 0.0 or pole_hz <= 0.0:
        raise ValueError("lead/lag zero and pole must be positive")
    if pole_hz <= zero_hz:
        raise ValueError("lead candidate requires pole_Hz > zero_Hz")
    numerator = 1.0 + (frequency / zero_hz) ** 2
    denominator = 1.0 + (frequency / pole_hz) ** 2
    return np.sqrt(numerator / denominator)


def equivalent_c2_zero_hz(c2: float, scale_hz: float) -> float | None:
    c2 = float(c2)
    scale_hz = float(scale_hz)
    if c2 <= 0.0:
        return None
    return float(scale_hz / np.sqrt(c2))


def _base_spec(
    *,
    base_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    return separated._build_spec(
        config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        c2_free=False,
    )


def build_candidate_spec(
    *,
    family: str,
    config: dict,
    base_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    spec = copy.deepcopy(
        _base_spec(
            base_config=base_config,
            material_C_tes=material_C_tes,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
        )
    )
    names = list(spec["names"])
    bounds = list(spec["bounds"])
    base_size = len(names)

    family_cfg = config["candidate_families"][family]
    if family == FAMILY_LEAD_LAG:
        zero_cfg = family_cfg["zero_Hz"]
        pole_cfg = family_cfg["pole_Hz"]
        names.extend(
            (
                "shared_log10_lead_zero_Hz",
                "shared_log10_lead_pole_Hz",
            )
        )
        bounds.extend(
            (
                (
                    np.log10(float(zero_cfg["min"])),
                    np.log10(float(zero_cfg["max"])),
                ),
                (
                    np.log10(float(pole_cfg["min"])),
                    np.log10(float(pole_cfg["max"])),
                ),
            )
        )
    elif family == FAMILY_BYPASS_NOISE:
        scale_cfg = family_cfg[
            "shared_asd_scale_to_geometric_mean_tracked_white"
        ]
        names.append("shared_log10_bypass_noise_scale")
        bounds.append(
            (
                np.log10(float(scale_cfg["min"])),
                np.log10(float(scale_cfg["max"])),
            )
        )
    else:
        raise ValueError(f"unknown candidate family: {family}")

    spec["names"] = tuple(names)
    spec["bounds"] = tuple(bounds)
    spec["base_size"] = int(base_size)
    spec["family"] = family
    return spec


def decode_candidate(
    vector,
    *,
    spec: dict,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
) -> dict:
    vector = np.asarray(vector, dtype=float)
    if vector.size != len(spec["names"]):
        raise ValueError("candidate vector length mismatch")

    base_vector = vector[: spec["base_size"]]
    decoded = separated._decode(
        base_vector,
        spec={
            **spec,
            "names": tuple(spec["names"][: spec["base_size"]]),
            "bounds": tuple(spec["bounds"][: spec["base_size"]]),
            "c2_free": False,
        },
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
    )
    tail = vector[spec["base_size"] :]

    if spec["family"] == FAMILY_LEAD_LAG:
        zero_hz = float(10.0 ** tail[0])
        pole_hz = float(10.0 ** tail[1])
        decoded["family_parameters"] = {
            "zero_Hz": zero_hz,
            "pole_Hz": pole_hz,
        }
        decoded["candidate_valid"] = bool(pole_hz > zero_hz)
    elif spec["family"] == FAMILY_BYPASS_NOISE:
        decoded["family_parameters"] = {
            "asd_scale_to_geomean_tracked_white": float(10.0 ** tail[0]),
        }
        decoded["candidate_valid"] = True
    else:
        raise ValueError(f"unknown family: {spec['family']}")

    return decoded


def _base_warm_vector(
    *,
    spec: dict,
    reference_problem: dict,
    repeat_problem: dict,
    solution: dict | None,
) -> np.ndarray:
    base_names = tuple(spec["names"][: spec["base_size"]])
    base_spec = {
        **spec,
        "names": base_names,
        "bounds": tuple(spec["bounds"][: spec["base_size"]]),
        "c2_free": False,
    }
    if solution is None:
        return separated._baseline_warm(
            base_spec,
            reference_problem,
            repeat_problem,
        )

    return separated._encode_warm(
        spec=base_spec,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        shared_C_tes=solution["shared_C_tes_J_per_K"],
        shared_L=solution["shared_L_H"],
        shared_C_substrate=solution["shared_C_substrate_J_per_K"],
        shared_G_ratio=solution[
            "shared_G_tes_substrate_over_G_substrate_bath"
        ],
        reference_detector=solution["reference_day"]["detector_candidate"],
        repeat_detector=solution["repeat_day"]["detector_candidate"],
        reference_white_scale=solution["reference_day"]["white_scale"],
        repeat_white_scale=solution["repeat_day"]["white_scale"],
        c2=0.0,
    )


def warm_candidate_vector(
    *,
    spec: dict,
    reference_problem: dict,
    repeat_problem: dict,
    zero_solution: dict | None,
    c2_reference_solution: dict | None,
) -> np.ndarray:
    base = _base_warm_vector(
        spec=spec,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        solution=zero_solution,
    )

    if spec["family"] == FAMILY_LEAD_LAG:
        zero_hz = 2_000.0
        if c2_reference_solution is not None:
            candidate = equivalent_c2_zero_hz(
                c2_reference_solution["shared_c2"],
                reference_problem["scale_hz"],
            )
            if candidate is not None:
                zero_hz = candidate
        zlo, zhi = spec["bounds"][spec["base_size"]]
        plo, phi = spec["bounds"][spec["base_size"] + 1]
        zero_hz = float(
            np.clip(zero_hz, 10.0**zlo, 10.0**zhi)
        )
        pole_hz = max(1_000_000.0, 20.0 * zero_hz)
        pole_hz = float(
            np.clip(pole_hz, 10.0**plo, 10.0**phi)
        )
        if pole_hz <= zero_hz:
            pole_hz = min(10.0**phi, 1.1 * zero_hz)
        tail = np.log10([zero_hz, pole_hz])
    else:
        tail = np.asarray([0.0], dtype=float)  # scale = 1

    return np.concatenate((base, np.asarray(tail, dtype=float)))


def _hardware_magnitudes(context):
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
    return frequency, alias_frequency, hardware_main, hardware_alias


def model_candidate(
    detector_candidate: dict,
    *,
    family: str,
    family_parameters: dict,
    frequency,
    white_asd: float,
    magnicon_pole_hz: float,
    magnicon_pole_q: float,
    bypass_noise_asd: float | None = None,
):
    point = opt.tes_operating_point(detector_candidate)
    if not point.get("valid") or not point.get("stable"):
        return None, point

    context = competition.fixed_white_context(
        detector_candidate,
        frequency,
        white_asd,
    )
    (
        frequency,
        alias_frequency,
        hardware_main,
        hardware_alias,
    ) = _hardware_magnitudes(context)

    magnicon_main = magnicon_bessel_magnitude(
        frequency,
        magnicon_pole_hz,
        magnicon_pole_q,
    )
    magnicon_alias = magnicon_bessel_magnitude(
        alias_frequency,
        magnicon_pole_hz,
        magnicon_pole_q,
    )

    if family == FAMILY_LEAD_LAG:
        extra_main = lead_lag_magnitude(
            frequency,
            family_parameters["zero_Hz"],
            family_parameters["pole_Hz"],
        )
        extra_alias = lead_lag_magnitude(
            alias_frequency,
            family_parameters["zero_Hz"],
            family_parameters["pole_Hz"],
        )
    elif family == FAMILY_BYPASS_NOISE:
        extra_main = np.ones_like(frequency)
        extra_alias = np.ones_like(alias_frequency)
    else:
        raise ValueError(f"unknown family: {family}")

    main = (
        context["main_intrinsic_asd"]
        * hardware_main
        * magnicon_main
        * extra_main
    )
    alias = (
        context["alias_intrinsic_asd"]
        * hardware_alias
        * magnicon_alias
        * extra_alias
    )
    alias = np.where(context["same_bin"], 0.0, alias)

    total_psd = main**2 + alias**2
    white = float(context["post_filter_white_asd_A_rtHz"])
    total_psd = total_psd + white**2

    if family == FAMILY_BYPASS_NOISE:
        if bypass_noise_asd is None or bypass_noise_asd <= 0.0:
            raise ValueError("bypass noise ASD must be positive")
        bypass_main = float(bypass_noise_asd) * hardware_main
        bypass_alias = float(bypass_noise_asd) * hardware_alias
        bypass_alias = np.where(
            context["same_bin"],
            0.0,
            bypass_alias,
        )
        total_psd = (
            total_psd
            + bypass_main**2
            + bypass_alias**2
        )

    absolute = np.sqrt(total_psd)
    model = opt.normalize_at(
        frequency,
        absolute,
        reference_hz=REFERENCE_HZ,
    )
    if np.any(~np.isfinite(model)) or np.any(model <= 0.0):
        return None, point
    return model, point


def fit_candidate_branch(
    *,
    family: str,
    config: dict,
    base_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
    zero_branch: dict,
    c2_reference_branch: dict,
) -> dict:
    spec = build_candidate_spec(
        family=family,
        config=config,
        base_config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    lower = np.asarray(
        [row[0] for row in spec["bounds"]],
        dtype=float,
    )
    upper = np.asarray(
        [row[1] for row in spec["bounds"]],
        dtype=float,
    )
    penalty = float(
        reference_problem["optimizer_cfg"]["instability_penalty"]
    )

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

    tracked_geomean = float(
        np.sqrt(
            float(reference_problem["tracked_white"])
            * float(repeat_problem["tracked_white"])
        )
    )

    evaluation_count = 0
    invalid_reference = 0
    invalid_repeat = 0

    def evaluate(vector):
        nonlocal evaluation_count, invalid_reference, invalid_repeat
        evaluation_count += 1
        decoded = decode_candidate(
            vector,
            spec=spec,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
            readout_reference=readout_reference,
        )
        if not decoded["candidate_valid"]:
            return None, decoded, {}, {}

        bypass_asd = None
        if family == FAMILY_BYPASS_NOISE:
            bypass_asd = (
                tracked_geomean
                * decoded["family_parameters"][
                    "asd_scale_to_geomean_tracked_white"
                ]
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
                    family=family,
                    family_parameters=decoded["family_parameters"],
                    frequency=problem["frequency"],
                    white_asd=white_asd,
                    magnicon_pole_hz=readout_reference["pole_Hz"],
                    magnicon_pole_q=readout_reference["pole_Q"],
                    bypass_noise_asd=bypass_asd,
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

        return (
            models,
            decoded,
            points["reference"],
            points["repeat"],
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
    family_index = FAMILIES.index(family)
    seed = (
        int(reference_problem["optimizer_cfg"]["seed"])
        + int(optimizer["seed_offset"])
        + 1000 * family_index
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

    warm = warm_candidate_vector(
        spec=spec,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        zero_solution=zero_branch["solution"],
        c2_reference_solution=c2_reference_branch["solution"],
    )
    if np.all(warm >= lower) and np.all(warm <= upper):
        candidates.append(("c2_zero_based_warm", warm.copy()))

    initial = list(candidates)
    for label, vector in initial:
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
        raise RuntimeError(f"no stable solution for {family}")

    score, source, vector = min(scored, key=lambda row: row[0])
    models, decoded, point_ref, point_rep = evaluate(vector)
    if models is None:
        raise RuntimeError("selected candidate solution became invalid")

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

    family_hits = {}
    offset = spec["base_size"]
    if family == FAMILY_LEAD_LAG:
        family_hits = {
            "zero_Hz": shared._encoded_boundary(
                decoded["family_parameters"]["zero_Hz"],
                10.0 ** spec["bounds"][offset][0],
                10.0 ** spec["bounds"][offset][1],
                log=True,
            ),
            "pole_Hz": shared._encoded_boundary(
                decoded["family_parameters"]["pole_Hz"],
                10.0 ** spec["bounds"][offset + 1][0],
                10.0 ** spec["bounds"][offset + 1][1],
                log=True,
            ),
        }
    else:
        family_hits = {
            "asd_scale_to_geomean_tracked_white": (
                shared._encoded_boundary(
                    decoded["family_parameters"][
                        "asd_scale_to_geomean_tracked_white"
                    ],
                    10.0 ** spec["bounds"][offset][0],
                    10.0 ** spec["bounds"][offset][1],
                    log=True,
                )
            )
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
            "white_scale": float(decoded["white_scales"][prefix]),
            "white_asd_A_rtHz": float(
                problem["tracked_white"]
                * decoded["white_scales"][prefix]
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

    solution = {
        "shared_C_tes_J_per_K": float(decoded["shared_C_tes"]),
        "shared_C_tes_over_material": float(
            decoded["shared_C_tes"] / material_C_tes
        ),
        "shared_C_substrate_J_per_K": float(
            decoded["shared_C_substrate"]
        ),
        "shared_L_H": float(decoded["shared_L"]),
        "shared_G_tes_substrate_over_G_substrate_bath": float(
            decoded["shared_G_ratio"]
        ),
        "family": family,
        "family_parameters": dict(decoded["family_parameters"]),
        "family_boundary_hits": family_hits,
        "shared_physical_boundary_hits": shared_hits,
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
    if family == FAMILY_BYPASS_NOISE:
        solution["shared_bypass_noise_asd_A_rtHz"] = float(
            tracked_geomean
            * decoded["family_parameters"][
                "asd_scale_to_geomean_tracked_white"
            ]
        )

    return {
        "name": family,
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


def replacement_summary(
    candidate: dict,
    c2_zero: dict,
    c2_reference: dict,
    screen: dict,
) -> dict:
    score_ratio_to_zero = float(
        candidate["joint_shape_score"]
        / c2_zero["joint_shape_score"]
    )
    rms_improvement_vs_zero = float(
        c2_zero["joint_continuum_rms_dB"]
        - candidate["joint_continuum_rms_dB"]
    )
    score_ratio_to_c2 = float(
        candidate["joint_shape_score"]
        / c2_reference["joint_shape_score"]
    )
    rms_delta_to_c2 = float(
        candidate["joint_continuum_rms_dB"]
        - c2_reference["joint_continuum_rms_dB"]
    )

    material_vs_zero = bool(
        score_ratio_to_zero
        <= float(
            screen[
                "max_score_ratio_to_c2_zero_for_material_improvement"
            ]
        )
        and rms_improvement_vs_zero
        >= float(screen["min_rms_improvement_vs_c2_zero_dB"])
    )
    near_c2 = bool(
        score_ratio_to_c2
        <= float(screen["max_score_ratio_to_c2_reference"])
        and rms_delta_to_c2
        <= float(screen["max_rms_delta_from_c2_reference_dB"])
    )

    if near_c2 and material_vs_zero:
        classification = "candidate_replaces_c2_within_screen"
    elif material_vs_zero:
        classification = (
            "candidate_materially_improves_c2_zero_but_not_c2_reference"
        )
    else:
        classification = "candidate_does_not_materially_replace_c2"

    return {
        "classification": classification,
        "material_improvement_vs_c2_zero": material_vs_zero,
        "within_c2_reference_screen": near_c2,
        "score_ratio_to_c2_zero": score_ratio_to_zero,
        "rms_improvement_vs_c2_zero_dB": rms_improvement_vs_zero,
        "score_ratio_to_c2_reference": score_ratio_to_c2,
        "rms_delta_from_c2_reference_dB": rms_delta_to_c2,
        "screen": screen,
    }


def run(config: dict, config_path: Path) -> dict:
    base_path = resolve_config_path(
        config["base_separated_substrate_config"],
        config_path,
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

    c2_zero = separated._fit_joint_branch(
        name="replacement_reference_c2_zero",
        config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
        c2_free=False,
    )
    c2_reference = separated._fit_joint_branch(
        name="replacement_reference_c2_free",
        config=base_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
        c2_free=True,
        warm_from=c2_zero,
    )

    candidates = {}
    comparisons = {}
    for family in FAMILIES:
        row = fit_candidate_branch(
            family=family,
            config=config,
            base_config=base_config,
            material_C_tes=material_C_tes,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
            readout_reference=readout_reference,
            zero_branch=c2_zero,
            c2_reference_branch=c2_reference,
        )
        candidates[family] = row
        comparisons[family] = replacement_summary(
            row,
            c2_zero,
            c2_reference,
            config["replacement_screen"],
        )

    lead = candidates[FAMILY_LEAD_LAG]
    bypass = candidates[FAMILY_BYPASS_NOISE]
    lead_solution = lead["solution"]
    lead_pole = float(
        lead_solution["family_parameters"]["pole_Hz"]
    )
    lead_pole_upper = float(
        config["candidate_families"][FAMILY_LEAD_LAG]["pole_Hz"]["max"]
    )

    if any(
        row["within_c2_reference_screen"]
        for row in comparisons.values()
    ):
        classification = "at_least_one_physical_candidate_replaces_c2"
    elif any(
        row["material_improvement_vs_c2_zero"]
        for row in comparisons.values()
    ):
        classification = (
            "candidate_physics_improves_c2_zero_but_c2_shape_still_unmatched"
        )
    else:
        classification = "tested_physical_candidates_do_not_replace_c2"

    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "Can either one finite real lead/lag readout section or one "
            "downstream additive electronics-noise source replace the "
            "phenomenological shared c2 after separating the substrate node?"
        ),
        "candidate_semantics": {
            FAMILY_LEAD_LAG: config["candidate_families"][
                FAMILY_LEAD_LAG
            ]["semantics"],
            FAMILY_BYPASS_NOISE: config["candidate_families"][
                FAMILY_BYPASS_NOISE
            ]["semantics"],
        },
        "c2_zero_reference": {
            key: value
            for key, value in c2_zero.items()
            if not key.startswith("_")
        },
        "c2_free_reference": {
            key: value
            for key, value in c2_reference.items()
            if not key.startswith("_")
        },
        "candidates": {
            name: {
                key: value
                for key, value in row.items()
                if not key.startswith("_")
            }
            for name, row in candidates.items()
        },
        "replacement_comparisons": comparisons,
        "interpretation": {
            "classification": classification,
            "c2_reference_value": float(
                c2_reference["solution"]["shared_c2"]
            ),
            "c2_reference_equivalent_zero_Hz": (
                equivalent_c2_zero_hz(
                    c2_reference["solution"]["shared_c2"],
                    reference_problem["scale_hz"],
                )
            ),
            "lead_lag_zero_Hz": float(
                lead_solution["family_parameters"]["zero_Hz"]
            ),
            "lead_lag_pole_Hz": lead_pole,
            "lead_lag_pole_over_zero": float(
                lead_pole
                / lead_solution["family_parameters"]["zero_Hz"]
            ),
            "lead_lag_pole_within_1pct_of_upper_bound": bool(
                lead_pole >= 0.99 * lead_pole_upper
            ),
            "bypass_noise_asd_A_rtHz": float(
                bypass["solution"]["shared_bypass_noise_asd_A_rtHz"]
            ),
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
            "reference_c2_zero": c2_zero["_models"]["reference"].tolist(),
            "repeat_c2_zero": c2_zero["_models"]["repeat"].tolist(),
            "reference_c2_free": c2_reference["_models"]["reference"].tolist(),
            "repeat_c2_free": c2_reference["_models"]["repeat"].tolist(),
            "reference_lead_lag": lead["_models"]["reference"].tolist(),
            "repeat_lead_lag": lead["_models"]["repeat"].tolist(),
            "reference_bypass_noise": bypass["_models"]["reference"].tolist(),
            "repeat_bypass_noise": bypass["_models"]["repeat"].tolist(),
        },
    }


def _residual_db(model, target):
    return 20.0 * np.log10(
        np.asarray(model, dtype=float)
        / np.asarray(target, dtype=float)
    )


def make_plot(result: dict, output: Path, show=False) -> None:
    import matplotlib.pyplot as plt

    p = result["_plot"]
    f = np.asarray(p["frequency_Hz"], dtype=float)
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(14.0, 8.8),
        sharex="col",
        gridspec_kw={"height_ratios": [2.0, 1.0]},
    )
    specs = [
        (
            0,
            "2024-12-06 reference",
            np.asarray(p["reference_target"], dtype=float),
            np.asarray(p["reference_raw"], dtype=float),
            np.asarray(p["reference_c2_zero"], dtype=float),
            np.asarray(p["reference_c2_free"], dtype=float),
            np.asarray(p["reference_lead_lag"], dtype=float),
            np.asarray(p["reference_bypass_noise"], dtype=float),
        ),
        (
            1,
            "2024-12-05 repeat",
            np.asarray(p["repeat_target"], dtype=float),
            np.asarray(p["repeat_raw"], dtype=float),
            np.asarray(p["repeat_c2_zero"], dtype=float),
            np.asarray(p["repeat_c2_free"], dtype=float),
            np.asarray(p["repeat_lead_lag"], dtype=float),
            np.asarray(p["repeat_bypass_noise"], dtype=float),
        ),
    ]

    for col, label, target, raw, zero, c2, lead, bypass in specs:
        top = axes[0, col]
        bottom = axes[1, col]
        top.loglog(f, raw, linewidth=0.7, alpha=0.35, label="Raw ASD")
        top.loglog(f, target, linewidth=1.5, label="Continuum target")
        top.loglog(f, zero, linewidth=1.1, label="c2=0")
        top.loglog(f, c2, linewidth=1.4, label="c2 free reference")
        top.loglog(f, lead, linewidth=1.3, label="real lead-lag")
        top.loglog(
            f,
            bypass,
            linewidth=1.3,
            label="bypass electronics noise",
        )
        top.set_title(label, fontsize=10)
        top.set_ylabel("Normalized ASD")
        top.grid(True, which="both", alpha=0.2)
        top.legend(frameon=False, fontsize=7)

        for model, name in (
            (zero, "c2=0"),
            (c2, "c2 free"),
            (lead, "lead-lag"),
            (bypass, "bypass noise"),
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

    lead = result["candidates"][FAMILY_LEAD_LAG]
    bypass = result["candidates"][FAMILY_BYPASS_NOISE]
    c2 = result["c2_free_reference"]
    fig.suptitle(
        result["interpretation"]["classification"]
        + "\n"
        + (
            f"c2={c2['solution']['shared_c2']:.3g}, "
            f"lead zero/pole="
            f"{lead['solution']['family_parameters']['zero_Hz']:.3g}/"
            f"{lead['solution']['family_parameters']['pole_Hz']:.3g} Hz, "
            f"bypass ASD="
            f"{bypass['solution']['shared_bypass_noise_asd_A_rtHz']:.3e} "
            "A/rtHz"
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
        json.dumps(cleaned_result(result), indent=2, allow_nan=False)
        + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "figure": str(args.figure),
                "classification": result["interpretation"][
                    "classification"
                ],
                "c2_reference": {
                    "c2": result["interpretation"]["c2_reference_value"],
                    "joint_rms_dB": result["c2_free_reference"][
                        "joint_continuum_rms_dB"
                    ],
                },
                "real_lead_lag": {
                    "parameters": result["candidates"][
                        FAMILY_LEAD_LAG
                    ]["solution"]["family_parameters"],
                    "joint_rms_dB": result["candidates"][
                        FAMILY_LEAD_LAG
                    ]["joint_continuum_rms_dB"],
                    "comparison": result["replacement_comparisons"][
                        FAMILY_LEAD_LAG
                    ],
                },
                "bypass_electronics_noise": {
                    "shared_asd_A_rtHz": result["interpretation"][
                        "bypass_noise_asd_A_rtHz"
                    ],
                    "joint_rms_dB": result["candidates"][
                        FAMILY_BYPASS_NOISE
                    ]["joint_continuum_rms_dB"],
                    "comparison": result["replacement_comparisons"][
                        FAMILY_BYPASS_NOISE
                    ],
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
