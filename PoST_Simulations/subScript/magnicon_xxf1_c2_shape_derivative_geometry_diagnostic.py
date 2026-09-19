"""Local shape-derivative geometry for c2 versus C_tes and L.

This diagnostic follows the c2 detector-compensation experiment.  Instead of
asking which optimizer branch wins after large parameter moves, it measures the
local tangent directions of the normalized ASD model at each day's fitted
Magnicon+c2 solution.

For each day, Magnicon pole/Q and the profiled white floor are held fixed.  The
diagnostic estimates

    d ln ASD / d c2,
    d ln ASD / d ln C_tes,
    d ln ASD / d ln L,

using bound-aware finite differences.  It then compares those vectors with
signed/absolute cosine similarity and a least-squares one-direction projection.

A second comparison uses the finite model change produced by forcing c2 to zero
with the detector held fixed.  This reveals whether the local derivative
picture remains informative for the large c2 removal used by the preceding
compensation diagnostic.

Derivative alignment is an identifiability statement only.  It is not a claim
that C_tes, L, or c2 is the physical origin of the residual.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
REPOSITORY_ROOT = ROOT.parent
CONFIG_DIR = ROOT / "config"
DEFAULT_CONFIG = (
    CONFIG_DIR / "magnicon_xxf1_c2_shape_derivative_geometry_config.json"
)
DEFAULT_OUTPUT = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_c2_shape_derivative_geometry_diagnostic.json"
)
DEFAULT_FIGURE = (
    ROOT
    / ".noise_optimization_work_rsh_sweep"
    / "magnicon_xxf1_c2_shape_derivative_geometry_diagnostic.png"
)

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from subScript import magnicon_xxf1_c2_detector_compensation_diagnostic as compensation  # noqa: E402
from subScript import magnicon_xxf1_lpf_c2_cross_day_repeatability_diagnostic as crossday  # noqa: E402
from subScript import readout_residual_dof_competition_diagnostic as residual  # noqa: E402


def resolve_config_path(value, config_path: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (config_path.parent / path).resolve()


def cosine_similarity(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.shape != b.shape:
        raise ValueError("cosine vectors must have matching shape")
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))
    if na <= 0.0 or nb <= 0.0:
        return None
    return float(np.dot(a, b) / (na * nb))


def projection_summary(target, basis):
    target = np.asarray(target, dtype=float)
    basis = np.asarray(basis, dtype=float)
    if target.shape != basis.shape:
        raise ValueError("projection vectors must have matching shape")
    target_norm = float(np.linalg.norm(target))
    basis_power = float(np.dot(basis, basis))
    if target_norm <= 0.0 or basis_power <= 0.0:
        return {
            "scale": None,
            "fractional_residual_norm": None,
            "explained_norm_fraction": None,
        }
    scale = float(np.dot(target, basis) / basis_power)
    residual_vector = target - scale * basis
    residual_fraction = float(
        np.linalg.norm(residual_vector) / target_norm
    )
    return {
        "scale": scale,
        "fractional_residual_norm": residual_fraction,
        "explained_norm_fraction": float(
            max(0.0, 1.0 - residual_fraction**2)
        ),
    }


def vector_comparison(target, basis):
    signed = cosine_similarity(target, basis)
    projection = projection_summary(target, basis)
    return {
        "cosine_similarity": signed,
        "absolute_cosine_similarity": (
            None if signed is None else float(abs(signed))
        ),
        "projection": projection,
    }


def band_comparisons(frequency, target, basis, bands):
    frequency = np.asarray(frequency, dtype=float)
    target = np.asarray(target, dtype=float)
    basis = np.asarray(basis, dtype=float)
    rows = {}
    for band in bands:
        low = float(band["min"])
        high = float(band["max"])
        mask = (frequency >= low) & (frequency <= high)
        if not np.any(mask):
            continue
        rows[str(band["name"])] = vector_comparison(
            target[mask],
            basis[mask],
        )
    return rows


def _model(
    problem: dict,
    detector: dict,
    readout: dict,
    white_asd: float,
):
    model, point = residual.model_for_candidate(
        detector,
        readout,
        problem["frequency"],
        problem["scale_hz"],
        white_asd,
    )
    if model is None:
        raise RuntimeError(
            "finite-difference model is invalid/unstable: "
            + json.dumps(point, default=float)
        )
    return np.asarray(model, dtype=float)


def _log_shape(model):
    model = np.asarray(model, dtype=float)
    if np.any(model <= 0.0) or np.any(~np.isfinite(model)):
        raise ValueError("model must be finite and positive")
    return np.log(model)


def _within(value, bounds):
    low, high = (float(x) for x in bounds)
    return bool(float(value) >= low and float(value) <= high)


def log_parameter_derivative(
    *,
    problem: dict,
    detector: dict,
    readout: dict,
    white_asd: float,
    parameter: str,
    log_step: float,
):
    """Estimate d ln(model) / d ln(parameter), bound-aware."""
    step = float(log_step)
    if step <= 0.0:
        raise ValueError("log_step must be positive")
    base_value = float(detector[parameter])
    if base_value <= 0.0:
        raise ValueError(f"{parameter} must be positive for log derivative")
    bounds = problem["detector_bounds"][parameter]

    base_model = _model(problem, detector, readout, white_asd)
    plus_value = float(base_value * np.exp(step))
    minus_value = float(base_value * np.exp(-step))

    plus_model = None
    minus_model = None
    if _within(plus_value, bounds):
        plus_detector = dict(detector)
        plus_detector[parameter] = plus_value
        try:
            plus_model = _model(
                problem,
                plus_detector,
                readout,
                white_asd,
            )
        except RuntimeError:
            plus_model = None
    if _within(minus_value, bounds):
        minus_detector = dict(detector)
        minus_detector[parameter] = minus_value
        try:
            minus_model = _model(
                problem,
                minus_detector,
                readout,
                white_asd,
            )
        except RuntimeError:
            minus_model = None

    if plus_model is not None and minus_model is not None:
        derivative = (
            _log_shape(plus_model) - _log_shape(minus_model)
        ) / (2.0 * step)
        scheme = "central_log"
    elif plus_model is not None:
        derivative = (
            _log_shape(plus_model) - _log_shape(base_model)
        ) / step
        scheme = "forward_log"
    elif minus_model is not None:
        derivative = (
            _log_shape(base_model) - _log_shape(minus_model)
        ) / step
        scheme = "backward_log"
    else:
        raise RuntimeError(
            f"no stable in-bounds finite difference for {parameter}"
        )

    return {
        "parameter": parameter,
        "derivative_definition": f"d ln ASD / d ln {parameter}",
        "base_value": base_value,
        "step": step,
        "scheme": scheme,
        "lower_bound": float(bounds[0]),
        "upper_bound": float(bounds[1]),
        "_vector": np.asarray(derivative, dtype=float),
    }


def c2_derivative(
    *,
    problem: dict,
    detector: dict,
    readout: dict,
    white_asd: float,
    relative_step: float,
    absolute_min_step: float,
    c2_upper_bound: float,
):
    """Estimate d ln(model) / d c2, bound-aware."""
    base_value = float(readout["c2"])
    if base_value < 0.0:
        raise ValueError("c2 must be non-negative")
    step = max(
        abs(base_value) * float(relative_step),
        float(absolute_min_step),
    )
    if step <= 0.0:
        raise ValueError("c2 finite-difference step must be positive")

    base_model = _model(problem, detector, readout, white_asd)
    plus_value = base_value + step
    minus_value = base_value - step

    plus_model = None
    minus_model = None
    if plus_value <= float(c2_upper_bound):
        plus_readout = dict(readout)
        plus_readout["c2"] = plus_value
        try:
            plus_model = _model(
                problem,
                detector,
                plus_readout,
                white_asd,
            )
        except RuntimeError:
            plus_model = None
    if minus_value >= 0.0:
        minus_readout = dict(readout)
        minus_readout["c2"] = minus_value
        try:
            minus_model = _model(
                problem,
                detector,
                minus_readout,
                white_asd,
            )
        except RuntimeError:
            minus_model = None

    if plus_model is not None and minus_model is not None:
        derivative = (
            _log_shape(plus_model) - _log_shape(minus_model)
        ) / (2.0 * step)
        scheme = "central_linear"
    elif plus_model is not None:
        derivative = (
            _log_shape(plus_model) - _log_shape(base_model)
        ) / step
        scheme = "forward_linear"
    elif minus_model is not None:
        derivative = (
            _log_shape(base_model) - _log_shape(minus_model)
        ) / step
        scheme = "backward_linear"
    else:
        raise RuntimeError("no stable in-bounds finite difference for c2")

    return {
        "parameter": "c2",
        "derivative_definition": "d ln ASD / d c2",
        "base_value": base_value,
        "step": step,
        "scheme": scheme,
        "lower_bound": 0.0,
        "upper_bound": float(c2_upper_bound),
        "_vector": np.asarray(derivative, dtype=float),
    }


def _finite_c2_removal_compensation(
    *,
    problem: dict,
    detector: dict,
    readout: dict,
    white_asd: float,
):
    base_model = _model(problem, detector, readout, white_asd)
    zero_readout = dict(readout)
    zero_readout["c2"] = 0.0
    zero_readout["c4"] = 0.0
    zero_model = _model(
        problem,
        detector,
        zero_readout,
        white_asd,
    )
    removal_change = _log_shape(zero_model) - _log_shape(base_model)
    needed_compensation = -removal_change
    return {
        "definition": (
            "-[ln ASD(c2=0, detector fixed) - ln ASD(c2=best, detector fixed)]"
        ),
        "_vector": needed_compensation,
        "rms_ln_shape": float(
            np.sqrt(np.mean(needed_compensation**2))
        ),
        "max_abs_ln_shape": float(
            np.max(np.abs(needed_compensation))
        ),
    }


def _clean_derivative(row: dict):
    return {
        key: value
        for key, value in row.items()
        if key != "_vector"
    }


def _top_parameter(comparisons: dict):
    valid = [
        (name, row["absolute_cosine_similarity"])
        for name, row in comparisons.items()
        if row["absolute_cosine_similarity"] is not None
    ]
    if not valid:
        return None
    return max(valid, key=lambda item: item[1])[0]


def _classify_day(finite_comparisons: dict, screen: dict):
    top = _top_parameter(finite_comparisons)
    top_value = (
        None
        if top is None
        else finite_comparisons[top]["absolute_cosine_similarity"]
    )
    if top_value is None:
        classification = "c2_compensation_geometry_unresolved"
    elif top_value >= float(screen["strong_abs_cosine"]):
        classification = "strong_single_local_shape_alignment"
    elif top_value >= float(screen["moderate_abs_cosine"]):
        classification = "moderate_single_local_shape_alignment"
    else:
        classification = "weak_single_local_shape_alignment"
    return {
        "classification": classification,
        "top_parameter": top,
        "top_absolute_cosine_similarity": top_value,
        "similarity_screen": screen,
    }


def _run_day(
    *,
    problem: dict,
    c2_config: dict,
    magnicon_config: dict,
    derivative_parameters,
    fd_config: dict,
    bands,
    screen: dict,
    seed_offset: int,
):
    local = crossday._fit_local_case(
        problem=problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        seed_offset=seed_offset,
    )
    free = local["branch"]["profiled_c2_free"]
    detector = dict(free["_detector_full"])
    readout = dict(free["_readout_full"])
    white = float(free["profiled_white_asd_A_rtHz"])

    c2_row = c2_derivative(
        problem=problem,
        detector=detector,
        readout=readout,
        white_asd=white,
        relative_step=float(fd_config["c2_relative_step"]),
        absolute_min_step=float(fd_config["c2_absolute_min_step"]),
        c2_upper_bound=float(c2_config["c2_upper_bound"]),
    )
    detector_rows = {}
    for name in derivative_parameters:
        detector_rows[name] = log_parameter_derivative(
            problem=problem,
            detector=detector,
            readout=readout,
            white_asd=white,
            parameter=name,
            log_step=float(fd_config["detector_log_step"]),
        )

    finite = _finite_c2_removal_compensation(
        problem=problem,
        detector=detector,
        readout=readout,
        white_asd=white,
    )

    local_comparisons = {}
    finite_comparisons = {}
    for name, row in detector_rows.items():
        local_comparisons[name] = {
            **vector_comparison(c2_row["_vector"], row["_vector"]),
            "bands": band_comparisons(
                problem["frequency"],
                c2_row["_vector"],
                row["_vector"],
                bands,
            ),
        }
        finite_comparisons[name] = {
            **vector_comparison(finite["_vector"], row["_vector"]),
            "bands": band_comparisons(
                problem["frequency"],
                finite["_vector"],
                row["_vector"],
                bands,
            ),
        }

    interpretation = _classify_day(
        finite_comparisons,
        screen,
    )

    return {
        "case": {
            "label": problem["repeat_case"]["label"],
            "role": problem["repeat_case"].get("role"),
            "accepted_records": int(problem["accepted_records"]),
            "comparison_source": problem["comparison_source"],
        },
        "base_fit": {
            "readout": {
                key: float(value)
                for key, value in readout.items()
            },
            "detector": {
                name: float(detector[name])
                for name in derivative_parameters
            },
            "white_asd_A_rtHz": white,
            "continuum_rms_dB": float(
                free["continuum_metrics_full"][
                    "residual_metrics"
                ]["rms_residual_dB"]
            ),
            "detector_boundary_hits": {
                name: free.get("detector_boundary_hits", {}).get(name)
                for name in derivative_parameters
            },
            "readout_boundary_hits": free.get(
                "readout_boundary_hits", {}
            ),
        },
        "derivatives": {
            "c2": _clean_derivative(c2_row),
            "detector": {
                name: _clean_derivative(row)
                for name, row in detector_rows.items()
            },
        },
        "local_derivative_similarity_to_c2": local_comparisons,
        "finite_c2_removal_compensation": {
            key: value
            for key, value in finite.items()
            if key != "_vector"
        },
        "finite_c2_removal_similarity": finite_comparisons,
        "interpretation": interpretation,
        "_vectors": {
            "frequency_Hz": problem["frequency"],
            "c2": c2_row["_vector"],
            "finite_compensation": finite["_vector"],
            **{
                name: row["_vector"]
                for name, row in detector_rows.items()
            },
        },
    }


def _cross_day_vector_similarity(reference: dict, repeat: dict, parameters):
    rows = {
        "c2": vector_comparison(
            reference["_vectors"]["c2"],
            repeat["_vectors"]["c2"],
        ),
        "finite_c2_removal_compensation": vector_comparison(
            reference["_vectors"]["finite_compensation"],
            repeat["_vectors"]["finite_compensation"],
        ),
    }
    for name in parameters:
        rows[name] = vector_comparison(
            reference["_vectors"][name],
            repeat["_vectors"][name],
        )
    return rows


def run(config: dict, config_path: Path):
    stack = compensation._load_stack(config, config_path)
    base_cross = stack["base_cross"]
    c2_config = dict(stack["c2_config"])
    c2_config["primary_normalization"] = str(
        base_cross["normalization"]
    )
    magnicon_config = stack["magnicon_config"]

    derivative_parameters = tuple(
        str(name) for name in config["derivative_parameters"]
    )
    allowed = {"C_tes", "L"}
    unknown = sorted(set(derivative_parameters) - allowed)
    if unknown:
        raise ValueError(
            "this diagnostic currently supports only C_tes and L: "
            + ", ".join(unknown)
        )

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

    reference = _run_day(
        problem=reference_problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        derivative_parameters=derivative_parameters,
        fd_config=config["finite_difference"],
        bands=config["bands_Hz"],
        screen=config["similarity_screen"],
        seed_offset=15100,
    )
    repeat = _run_day(
        problem=repeat_problem,
        c2_config=c2_config,
        magnicon_config=magnicon_config,
        derivative_parameters=derivative_parameters,
        fd_config=config["finite_difference"],
        bands=config["bands_Hz"],
        screen=config["similarity_screen"],
        seed_offset=15700,
    )

    ref_top = reference["interpretation"]["top_parameter"]
    rep_top = repeat["interpretation"]["top_parameter"]
    cross_day = _cross_day_vector_similarity(
        reference,
        repeat,
        derivative_parameters,
    )

    if ref_top is not None and rep_top is not None and ref_top != rep_top:
        classification = "c2_compensation_geometry_is_day_dependent"
    elif (
        reference["interpretation"]["classification"]
        == "strong_single_local_shape_alignment"
        and repeat["interpretation"]["classification"]
        == "strong_single_local_shape_alignment"
    ):
        classification = "same_detector_direction_tracks_c2_on_both_days"
    else:
        classification = "c2_compensation_geometry_not_strongly_resolved"

    result = {
        "diagnostic_only": True,
        "production_noise_model_unchanged": True,
        "tested_question": (
            "Do local normalized-ASD shape derivatives explain why C_tes "
            "compensates removed c2 on one day while L compensates it on the "
            "other?"
        ),
        "derivative_semantics": {
            "c2": "d ln(normalized ASD) / d c2",
            "C_tes": "d ln(normalized ASD) / d ln C_tes",
            "L": "d ln(normalized ASD) / d ln L",
            "magnicon_pole_Q_frozen": True,
            "white_floor_frozen": True,
            "local_to_each_day_best_fit": True,
        },
        "reference_day": {
            key: value
            for key, value in reference.items()
            if key != "_vectors"
        },
        "repeat_day": {
            key: value
            for key, value in repeat.items()
            if key != "_vectors"
        },
        "cross_day_derivative_similarity": cross_day,
        "interpretation": {
            "classification": classification,
            "reference_top_parameter": ref_top,
            "repeat_top_parameter": rep_top,
            "same_top_parameter": bool(
                ref_top is not None and ref_top == rep_top
            ),
            "guardrail": config["guardrail"],
        },
        "inputs": {
            "config": str(config_path),
            "base_cross_day_config": str(stack["base_cross_path"]),
            "base_c2_config": str(stack["c2_config_path"]),
            "base_magnicon_config": str(
                stack["magnicon_config_path"]
            ),
        },
        "_plot": {
            "reference": {
                key: (
                    value.tolist()
                    if isinstance(value, np.ndarray)
                    else value
                )
                for key, value in reference["_vectors"].items()
            },
            "repeat": {
                key: (
                    value.tolist()
                    if isinstance(value, np.ndarray)
                    else value
                )
                for key, value in repeat["_vectors"].items()
            },
        },
    }
    return result


def make_plot(result: dict, output: Path, show=False):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 8.0), sharex=True)

    for ax, day_key, title in [
        (axes[0], "reference", "12/06 reference"),
        (axes[1], "repeat", "12/05 repeat"),
    ]:
        p = result["_plot"][day_key]
        frequency = np.asarray(p["frequency_Hz"], dtype=float)
        c2 = np.asarray(p["c2"], dtype=float)
        finite = np.asarray(p["finite_compensation"], dtype=float)

        def norm_shape(values):
            values = np.asarray(values, dtype=float)
            norm = np.linalg.norm(values)
            return values / norm if norm > 0.0 else values

        ax.semilogx(
            frequency,
            norm_shape(c2),
            label="local c2 derivative",
        )
        ax.semilogx(
            frequency,
            norm_shape(finite),
            label="finite c2-removal compensation",
        )
        for name in ("C_tes", "L"):
            if name in p:
                ax.semilogx(
                    frequency,
                    norm_shape(np.asarray(p[name], dtype=float)),
                    label=name,
                )
        ax.axhline(0.0, linewidth=1.0)
        ax.set_ylabel("Unit-norm shape direction")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.2)
        ax.legend(frameon=False, fontsize=8)

    axes[1].set_xlabel("Frequency [Hz]")
    fig.suptitle(result["interpretation"]["classification"], fontsize=10)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)


def cleaned_result(result: dict):
    return {
        key: value
        for key, value in result.items()
        if key != "_plot"
    }


def main():
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
                "reference_top_parameter": result["interpretation"][
                    "reference_top_parameter"
                ],
                "repeat_top_parameter": result["interpretation"][
                    "repeat_top_parameter"
                ],
                "cross_day_derivative_similarity": result[
                    "cross_day_derivative_similarity"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
