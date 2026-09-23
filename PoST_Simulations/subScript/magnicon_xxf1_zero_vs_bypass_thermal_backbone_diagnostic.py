"""Compare zero-like transfer against Magnicon-bypass electronics noise.

All branches use the same relaxed thermal-backbone profile introduced by
magnicon_xxf1_thermal_backbone_pure_zero_diagnostic:

  * shared C_tes, L, C_substrate, conductance split;
  * broad shared G_tes-bath scale;
  * shared Pb absorber thickness with C/G recomputed from the Elmer material;
  * day-specific alpha, beta, T_bath, and post-filter white ASD.

Four branches are reported:
  1. c2=0;
  2. c2 free reference;
  3. one shared pure zero (the c2 numerator reparameterized by f_z);
  4. c2=0 plus one shared additive electronics-noise ASD injected after the
     Magnicon 10 kHz LPF and before the fixed 100 kHz SIM965/sampling stage.

The documented Magnicon response remains second-order Bessel at nominal 10 kHz
for all four branches.  This diagnostic does not identify a component; it asks
which source class has the required spectral leverage once the thermal backbone
is allowed the same freedoms in every branch.
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
    / "magnicon_xxf1_zero_vs_bypass_thermal_backbone_diagnostic_config.json"
)
DEFAULT_OUTPUT = (
    WORK_DIR
    / "magnicon_xxf1_zero_vs_bypass_thermal_backbone_diagnostic.json"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

import Opt_noise as opt  # noqa: E402
from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_c2_physical_replacement_diagnostic as physical  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import magnicon_xxf1_lpf_continuum_diagnostic as xxf1  # noqa: E402
from subScript import magnicon_xxf1_shared_L_Ctes_joint_c2_diagnostic as shared  # noqa: E402
from subScript import magnicon_xxf1_thermal_backbone_pure_zero_diagnostic as thermal  # noqa: E402

FAMILY_BYPASS = "bypass_electronics_noise"


def build_bypass_spec(
    *,
    config: dict,
    thermal_config: dict,
    separated_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
) -> dict:
    base = thermal.build_spec(
        branch=thermal.BRANCH_ZERO,
        config=thermal_config,
        base_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
    )
    scale_cfg = config[FAMILY_BYPASS][
        "shared_asd_scale_to_geometric_mean_tracked_white"
    ]
    lower = float(scale_cfg["min"])
    upper = float(scale_cfg["max"])
    if not 0.0 < lower < upper:
        raise ValueError("invalid bypass-noise scale bounds")
    return {
        "base_spec": base,
        "names": tuple(base["names"]) + ("shared_log10_bypass_noise_scale",),
        "bounds": tuple(base["bounds"])
        + ((np.log10(lower), np.log10(upper)),),
        "scale_bounds": (lower, upper),
    }


def decode_bypass(
    vector,
    *,
    spec: dict,
    thermal_config: dict,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
) -> dict:
    vector = np.asarray(vector, dtype=float)
    if vector.size != len(spec["names"]):
        raise ValueError("bypass vector length mismatch")
    decoded = thermal.decode_vector(
        vector[:-1],
        spec=spec["base_spec"],
        config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
    )
    decoded = dict(decoded)
    decoded["bypass_noise_scale"] = float(10.0 ** vector[-1])
    return decoded


def warm_bypass(
    *,
    spec: dict,
    thermal_config: dict,
    reference_problem: dict,
    repeat_problem: dict,
    zero_solution: dict | None,
) -> np.ndarray:
    base = thermal.warm_vector(
        spec=spec["base_spec"],
        config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        prior_solution=zero_solution,
    )
    lo, hi = spec["scale_bounds"]
    scale = float(np.clip(1.0, lo, hi))
    return np.concatenate((base, [np.log10(scale)]))


def fit_bypass(
    *,
    config: dict,
    thermal_config: dict,
    separated_config: dict,
    material_C_tes: float,
    reference_problem: dict,
    repeat_problem: dict,
    readout_reference: dict,
    zero_branch: dict,
) -> dict:
    spec = build_bypass_spec(
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
        decoded = decode_bypass(
            vector,
            spec=spec,
            thermal_config=thermal_config,
            reference_problem=reference_problem,
            repeat_problem=repeat_problem,
            readout_reference=readout_reference,
        )
        bypass_asd = tracked_geomean * decoded["bypass_noise_scale"]
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
                model, point = physical.model_candidate(
                    decoded["detectors"][prefix],
                    family=physical.FAMILY_BYPASS_NOISE,
                    family_parameters={
                        "asd_scale_to_geomean_tracked_white": decoded[
                            "bypass_noise_scale"
                        ]
                    },
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

    warm = warm_bypass(
        spec=spec,
        thermal_config=thermal_config,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        zero_solution=zero_branch["solution"],
    )
    if np.all(warm >= lower) and np.all(warm <= upper):
        candidates.append(("c2_zero_based_warm", warm))

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
        raise RuntimeError("no stable solution for bypass electronics noise")

    score, source, vector = min(scored, key=lambda row: row[0])
    models, decoded, point_ref, point_rep = evaluate(vector)
    if models is None:
        raise RuntimeError("selected bypass solution became invalid")

    ref_metrics = shared._day_metrics(reference_problem, models["reference"])
    rep_metrics = shared._day_metrics(repeat_problem, models["repeat"])
    combined_rms = shared._combined_rms(ref_metrics, rep_metrics)

    base_spec = spec["base_spec"]
    shared_hits = {
        name: shared._encoded_boundary(
            decoded["shared"][name],
            *base_spec["shared_bounds"][name],
            log=True,
        )
        for name in base_spec["shared_names"]
    }
    white_hits = {
        prefix: shared._encoded_boundary(
            decoded["white_scales"][prefix],
            *base_spec["white_bounds"][prefix],
            log=True,
        )
        for prefix in ("reference", "repeat")
    }
    bypass_hit = shared._encoded_boundary(
        decoded["bypass_noise_scale"],
        *spec["scale_bounds"],
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

    shared_values = decoded["shared"]
    historical = float(
        separated_config["shared_physical_profile"]["C_substrate_J_per_K"][
            "historical_PoST_total_reference_J_per_K"
        ]
    )
    local_total = shared_values["C_tes"] + shared_values["C_substrate"]
    bypass_asd = tracked_geomean * decoded["bypass_noise_scale"]

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
        "shared_bypass_noise_scale_to_geomean_tracked_white": float(
            decoded["bypass_noise_scale"]
        ),
        "shared_bypass_noise_asd_A_rtHz": float(bypass_asd),
        "shared_physical_boundary_hits": shared_hits,
        "bypass_noise_boundary_hit": bypass_hit,
        "reference_day": day_solution(
            "reference", reference_problem, point_ref, ref_metrics
        ),
        "repeat_day": day_solution(
            "repeat", repeat_problem, point_rep, rep_metrics
        ),
    }
    return {
        "name": FAMILY_BYPASS,
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


def replacement_metrics(
    candidate: dict,
    *,
    zero: dict,
    c2: dict,
    screen: dict,
) -> dict:
    score_ratio_c2 = float(
        candidate["joint_shape_score"] / c2["joint_shape_score"]
    )
    rms_delta_c2 = float(
        candidate["joint_continuum_rms_dB"] - c2["joint_continuum_rms_dB"]
    )
    score_ratio_zero = float(
        candidate["joint_shape_score"] / zero["joint_shape_score"]
    )
    rms_gain_zero = float(
        zero["joint_continuum_rms_dB"] - candidate["joint_continuum_rms_dB"]
    )
    within_c2 = bool(
        score_ratio_c2 <= float(screen["max_score_ratio_to_c2_reference"])
        and rms_delta_c2 <= float(screen["max_rms_delta_from_c2_reference_dB"])
    )
    material_vs_zero = bool(
        score_ratio_zero
        <= float(screen[
            "max_score_ratio_to_c2_zero_for_material_improvement"
        ])
        and rms_gain_zero
        >= float(screen["min_rms_improvement_vs_c2_zero_dB"])
    )
    return {
        "within_c2_reference_screen": within_c2,
        "material_improvement_vs_c2_zero": material_vs_zero,
        "score_ratio_to_c2_reference": score_ratio_c2,
        "rms_delta_from_c2_reference_dB": rms_delta_c2,
        "score_ratio_to_c2_zero": score_ratio_zero,
        "rms_improvement_vs_c2_zero_dB": rms_gain_zero,
    }


def compare_pure_zero_and_bypass(
    pure: dict,
    bypass: dict,
    *,
    screen: dict,
) -> dict:
    pure_over_bypass = float(
        pure["joint_shape_score"] / bypass["joint_shape_score"]
    )
    bypass_over_pure = float(1.0 / pure_over_bypass)
    pure_rms_advantage = float(
        bypass["joint_continuum_rms_dB"] - pure["joint_continuum_rms_dB"]
    )
    max_ratio = float(screen["direct_pair_max_score_ratio"])
    min_advantage = float(screen["direct_pair_min_rms_advantage_dB"])

    pure_material = bool(
        pure_over_bypass <= max_ratio
        and pure_rms_advantage >= min_advantage
    )
    bypass_material = bool(
        bypass_over_pure <= max_ratio
        and -pure_rms_advantage >= min_advantage
    )
    if pure_material:
        classification = "pure_zero_materially_better_than_bypass_noise"
    elif bypass_material:
        classification = "bypass_noise_materially_better_than_pure_zero"
    else:
        classification = "pure_zero_and_bypass_noise_not_discriminated_by_screen"

    return {
        "classification": classification,
        "pure_zero_materially_better": pure_material,
        "bypass_noise_materially_better": bypass_material,
        "pure_over_bypass_shape_score_ratio": pure_over_bypass,
        "pure_zero_rms_advantage_dB": pure_rms_advantage,
        "screen": {
            "max_score_ratio": max_ratio,
            "min_rms_advantage_dB": min_advantage,
        },
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

    filter_cfg = separated_config["magnicon_filter"]
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

    common = {
        "config": thermal_config,
        "base_config": separated_config,
        "material_C_tes": material_C_tes,
        "reference_problem": reference_problem,
        "repeat_problem": repeat_problem,
        "readout_reference": readout_reference,
    }
    zero = thermal.fit_branch(branch=thermal.BRANCH_ZERO, **common)
    c2 = thermal.fit_branch(
        branch=thermal.BRANCH_C2,
        warm_from=zero,
        **common,
    )
    pure = thermal.fit_branch(
        branch=thermal.BRANCH_PURE_ZERO,
        warm_from=c2,
        **common,
    )
    bypass = fit_bypass(
        config=config,
        thermal_config=thermal_config,
        separated_config=separated_config,
        material_C_tes=material_C_tes,
        reference_problem=reference_problem,
        repeat_problem=repeat_problem,
        readout_reference=readout_reference,
        zero_branch=zero,
    )

    screen = config["replacement_screen"]
    pure_replacement = replacement_metrics(
        pure, zero=zero, c2=c2, screen=screen
    )
    bypass_replacement = replacement_metrics(
        bypass, zero=zero, c2=c2, screen=screen
    )
    direct = compare_pure_zero_and_bypass(
        pure, bypass, screen=screen
    )

    return {
        "diagnostic_only": True,
        "production_default_topology_unchanged": True,
        "tested_question": (
            "Under the same relaxed thermal backbone, is the residual better "
            "represented by a multiplicative zero-like response or by an "
            "additive electronics-noise source that bypasses the Magnicon "
            "10 kHz filter?"
        ),
        "fit_semantics": {
            "magnicon": {
                "order": 2,
                "cutoff_Hz": float(filter_cfg["cutoff_Hz"]),
                "normalization": str(filter_cfg["normalization"]),
                "pole_Q": float(canonical["pole_Q"]),
            },
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
            "bypass_source_location": (
                "downstream of Magnicon 10 kHz LPF; upstream of fixed "
                "100 kHz SIM965 and sampling/first-alias fold"
            ),
            "external_analog_filter": {
                "model": "SRS SIM965",
                "front_panel_cutoff_Hz": float(
                    opt.TARGET_HARDWARE_BESSEL_CUTOFF_HZ
                ),
                "order": int(opt.TARGET_HARDWARE_BESSEL_ORDER),
                "scipy_norm": str(opt.TARGET_HARDWARE_BESSEL_NORM),
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
        "bypass_electronics_noise": {
            key: value
            for key, value in bypass.items()
            if not key.startswith("_")
        },
        "replacement_comparisons": {
            "pure_zero": pure_replacement,
            "bypass_electronics_noise": bypass_replacement,
        },
        "direct_source_class_comparison": direct,
        "interpretation": {
            "classification": direct["classification"],
            "pure_zero_Hz": float(
                pure["solution"]["shared_pure_zero_Hz"]
            ),
            "pure_zero_equivalent_c2": float(
                pure["solution"]["shared_c2"]
            ),
            "bypass_noise_asd_A_rtHz": float(
                bypass["solution"]["shared_bypass_noise_asd_A_rtHz"]
            ),
            "bypass_noise_scale_to_geomean_tracked_white": float(
                bypass["solution"][
                    "shared_bypass_noise_scale_to_geomean_tracked_white"
                ]
            ),
            "lab_if_pure_zero_better": (
                "Prioritize swept injected-signal/complex transfer measurement "
                "through the Magnicon/FLL/connector-box path."
            ),
            "lab_if_bypass_noise_better": (
                "Prioritize off-TES/reference-input electronics-noise spectra "
                "and measurements at nodes before/after the Magnicon output "
                "filter."
            ),
            "guardrail": str(config["guardrail"]),
        },
        "inputs": {
            "config": str(config_path),
            "base_thermal_backbone_config": str(thermal_path),
            "base_separated_substrate_config": str(separated_path),
            "base_cross_day_config": str(stack["base_cross_path"]),
            "base_magnicon_config": str(stack["magnicon_config_path"]),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    result = run(config, args.config)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )

    print(
        json.dumps(
            {
                "output": str(args.output),
                "classification": result[
                    "direct_source_class_comparison"
                ]["classification"],
                "joint_rms_dB": {
                    "c2_zero": result["c2_zero"][
                        "joint_continuum_rms_dB"
                    ],
                    "c2_reference": result["c2_reference"][
                        "joint_continuum_rms_dB"
                    ],
                    "pure_zero": result["pure_zero"][
                        "joint_continuum_rms_dB"
                    ],
                    "bypass_electronics_noise": result[
                        "bypass_electronics_noise"
                    ]["joint_continuum_rms_dB"],
                },
                "pure_zero_Hz": result["interpretation"][
                    "pure_zero_Hz"
                ],
                "bypass_noise_asd_A_rtHz": result["interpretation"][
                    "bypass_noise_asd_A_rtHz"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
