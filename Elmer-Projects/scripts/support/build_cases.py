"""Generate self-contained Elmer case SIFs from the `cases` section of
elmer_project.json (redesign plan, Phase 1).

Every runnable SIF lives in generated/cases/ and inlines all values — no
INCLUDEs, no MATC $-variables, no hand-edited copies. Values given as strings
in the case definitions are dimensioned expressions evaluated against the
project parameters (e.g. "20.02[ms]", "1332[keV]"); plain numbers are SI.

Templates:
- steady:    steady-state circuit case (heat_source: circuit_local uses the
             legacy per-sweep TESHeatSource, circuit_implicit uses the
             transient UDF's steady mode). Optionally writes a .result.
- transient: uniform- or staged-dt transient with the implicit circuit UDF.
- pulse:     transient restarted from another case's .result plus a
             rectangular-window Gaussian pulse in the absorber. The pulse
             center defaults to the absorber centroid of the case's mesh and
             the discrete norm is computed from the mesh at build time.
"""
from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from scripts.support.mesh_names import MeshNames, parse_mesh_names
from scripts.support.mesh_quantities import absorber_centroid, gaussian_discrete_norm
from scripts.support.vendored.dimensioned_expression import (
    dimension_name_of,
    evaluate_dimensioned_expression,
)

# Upstream unit tables use keV as the energy base unit.
KEV_TO_JOULE = 1.602176634e-16

STEFAN_BOLTZMANN = "5.670374419e-8"

# Material index (matches materials_block's fixed Material 1..7 ordering),
# keyed by a body's base name (its mesh.names name with any dual-TES `_L`/
# `_R` stack suffix stripped).
MATERIAL_BY_BASE_NAME = {
    "abs": 1,
    "TES": 2,
    "Stycast": 3,
    "Membrane_SiNx": 7,
    "SiO2_1": 4,
    "Si_1": 5,
    "SiNx": 6,
    "Si_2": 5,
    "SiO2_2": 4,
    "Membrane_Si1": 7,
}

# Slave -> master mortar pairs, by (base name, face) on each side. The abs
# side of the third pair is never suffixed: in a dual-TES mesh both stacks'
# Stycast slaves mortar onto the same absorber boundary.
_MORTAR_PAIRS = [
    ("TES", "zmin", "Membrane_SiNx", "zmax", True),
    ("Stycast", "zmin", "TES", "zmax", True),
    ("Stycast", "zmax", "abs", "zmin", False),
]

_FACE_LABEL = {"zmin": "bottom", "zmax": "top"}

# Display name used in mortar BC labels; only "abs" differs from its body
# name (the absorber's material, and the legacy SIF's label for it, is "Pb").
_MORTAR_LABEL_BASE = {"abs": "Pb"}

HEAT_SOURCES = {
    "circuit_local": ("tes_heat_source_t0", "TESHeatSource", "TES shunt power from local T"),
    "circuit_implicit": (
        "tes_transient_heat_source_t0",
        "TESTransientHeatSource",
        "TES R(T,I) plus branch inductance",
    ),
    "circuit_parallel": (
        "tes_parallel_circuit",
        "TESParallelHeatSource",
        "MPI-safe lumped TES circuit power",
    ),
    "circuit_inner": (
        "tes_parallel_circuit",
        "TESParallelHeatSource",
        "MPI-safe circuit updated inside HeatSolve nonlinear iterations",
    ),
}

PULSE_PROCEDURE = ("tes_transient_heat_source_t0", "AbsorberWindowPulseHeatSource")


def fmt(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.16g}"
    return str(value)


def fmt_real(value: Any) -> str:
    """Format a float so Elmer's SIF parser types it as Real even when the
    value is integral (a bare "1" becomes an integer entry, which silently
    changes e.g. which relaxation code path the solver takes). Matches the
    old chain, where these values arrived as MATC variables (always Real)."""
    text = fmt(float(value))
    if re.fullmatch(r"-?\d+", text):
        text += ".0"
    return text


def _round_noise(value: float) -> float:
    # Round away float64 arithmetic noise (10*1e-6 != 1e-5) while keeping far
    # more precision than the simulation needs.
    return float(f"{value:.15g}")


def eval_si(value: Any, params: dict[str, float]) -> float:
    """Evaluate a case-definition value: numbers are SI already, strings are
    dimensioned expressions over the project parameters."""
    if isinstance(value, bool):
        raise ValueError(f"Unexpected boolean where a number was expected: {value!r}")
    if isinstance(value, (int, float)):
        return float(value)
    quantity = evaluate_dimensioned_expression(str(value), variables=params)
    if dimension_name_of(quantity) == "energy":
        return _round_noise(quantity.value * KEV_TO_JOULE)
    return _round_noise(float(quantity.value))


def matc_guarded_membrane_expr(model: dict) -> str:
    membrane_expr = str(model["materials"]["Membrane"]["k"]["expression"])
    params = model["parameters"]
    # Inline parameter values (the generated SIFs are self-contained, so no
    # $-variables exist to resolve identifiers), longest names first so that
    # e.g. 'membrane_dx' is not clobbered while replacing a shorter name.
    names = sorted(
        (n for n in set(re.findall(r"[A-Za-z_]\w*", membrane_expr)) if n in params),
        key=len,
        reverse=True,
    )
    inlined = membrane_expr
    for name in names:
        inlined = re.sub(rf"\b{re.escape(name)}\b", f"({fmt(params[name])})", inlined)
    # Guard the non-integer power law against negative or zero iterates in
    # the nonlinear solve while preserving the original k(T) for T > 0.
    t_guard = "(0.5*((tx+1.0e-12)+abs(tx-1.0e-12)))"
    guarded = re.sub(r"\bT\b", t_guard, inlined)
    return guarded.replace("**", "^")


def materials_block(model: dict) -> list[str]:
    materials = model["materials"]
    lines: list[str] = []
    order = [
        (1, "Pb", "Pb"),
        (2, "TES", "TES"),
        (3, "Stycast", "Stycast"),
        (4, "SiO2", "SiO2"),
        (5, "Si", "Si"),
        (6, "SiNx", "SiNx"),
    ]
    for idx, sif_name, key in order:
        mat = materials[key]
        lines += [
            f"Material {idx}",
            f'  Name = "{sif_name}"',
            f"  Density = {fmt(mat['rho']['nominal'])}",
            f"  Heat Capacity = {fmt(mat['cp']['nominal'])}",
            f"  Heat Conductivity = {fmt(mat['k']['nominal'])}",
            "End",
            "",
        ]
    membrane = materials["Membrane"]
    lines += [
        "Material 7",
        '  Name = "Membrane"',
        f"  Density = {fmt(membrane['rho']['nominal'])}",
        f"  Heat Capacity = {fmt(membrane['cp']['nominal'])}",
        "  Heat Conductivity = Variable Temperature",
        f'    Real MATC "{matc_guarded_membrane_expr(model)}"',
        "End",
        "",
    ]
    return lines


def tes_constants_block(
    params: dict[str, float], series_file: str | None, state_file: str | None = None
) -> list[str]:
    lines = [
        "Constants",
        f"  Stefan Boltzmann = Real {STEFAN_BOLTZMANN}",
    ]
    if series_file:
        lines.append(f'  TES Series File = String "{series_file}"')
    if state_file:
        lines.append(f'  TES State File = String "{state_file}"')
    lines += [
        f"  TES Bias Current = Real {fmt(params['I_bias'])}",
        f"  TES Shunt Resistance = Real {fmt(params['R_sh'])}",
        f"  TES R0 = Real {fmt(params['R_0'])}",
        f"  TES Rmin = Real {fmt(params['R_min'])}",
        f"  TES Alpha = Real {fmt(params['alpha'])}",
        f"  TES Beta = Real {fmt(params['beta'])}",
        f"  TES I0 = Real {fmt(params['I_0'])}",
        f"  TES Tc = Real {fmt(params['T_c'])}",
        f"  TES T0 = Real {fmt(params['T_0'])}",
        f"  TES Volume = Real {fmt(params['TES_volume'])}",
        f"  TES Inductance = Real {fmt(params['L_tes'])}",
    ]
    return lines


def _side_series_file(series_file: str, side: str) -> str:
    """Insert a dual-TES side tag ('L'/'R') into a series_file name,
    following the existing '..._series.csv' naming convention (e.g.
    'tes_dual_steady_series.csv' -> 'tes_dual_steady_L_series.csv')."""
    suffix = "_series.csv"
    if series_file.endswith(suffix):
        return f"{series_file[: -len(suffix)]}_{side}{suffix}"
    stem, dot, ext = series_file.rpartition(".")
    return f"{stem}_{side}.{ext}" if dot else f"{series_file}_{side}"


def tes_side_constants_block(
    side: str, params: dict[str, float], series_file: str, state_file: str | None = None
) -> list[str]:
    """Constants entries (no `Constants`/Stefan-Boltzmann header) for one
    dual-TES circuit instance, prefixed 'TES <side> ' to match
    TESTransientHeatSource<side> in tes_transient_heat_source.f90. Takes
    *params* and *series_file* as explicit arguments (rather than closing
    over a single shared dict) so a future case spec can override one
    side's circuit parameters or series file independently of the other.

    *state_file*, if given, emits 'TES <side> State File' -- the persisted
    steady-state circuit state (docs/dual_tes_plan.md). A steady dual case
    writes it every nonlinear iteration (last write = converged state); a
    restarted transient/pulse dual case reads it once at initialization
    instead of falling back to the T0 analytic estimate, which sits ~1.3 mK
    off the true steady operating point and otherwise causes a restart
    transient. Omitted entirely for single-pixel cases (no 'side'), which
    keeps their SIFs byte-identical."""
    prefix = f"TES {side} "
    lines = [
        f"  {prefix}Bias Current = Real {fmt(params['I_bias'])}",
        f"  {prefix}Shunt Resistance = Real {fmt(params['R_sh'])}",
        f"  {prefix}Inductance = Real {fmt(params['L_tes'])}",
        f"  {prefix}R0 = Real {fmt(params['R_0'])}",
        f"  {prefix}Rmin = Real {fmt(params['R_min'])}",
        f"  {prefix}Alpha = Real {fmt(params['alpha'])}",
        f"  {prefix}Beta = Real {fmt(params['beta'])}",
        f"  {prefix}I0 = Real {fmt(params['I_0'])}",
        f"  {prefix}T0 = Real {fmt(params['T_0'])}",
        f"  {prefix}Volume = Real {fmt(params['TES_volume'])}",
        f'  {prefix}Series File = String "{series_file}"',
    ]
    if state_file:
        lines.append(f'  {prefix}State File = String "{state_file}"')
    return lines


def dual_tes_constants_block(
    sides: list[str],
    params: dict[str, float],
    series_file_base: str,
    state_files: dict[str, str] | None = None,
) -> list[str]:
    """Full `Constants` block for a dual-TES case: Stefan-Boltzmann plus a
    'TES <side> ...' set per side, each side's series file derived from
    *series_file_base* via _side_series_file. *state_files*, if given, maps
    side -> state file path (see tes_side_constants_block)."""
    lines = ["Constants", f"  Stefan Boltzmann = Real {STEFAN_BOLTZMANN}"]
    for side in sides:
        state_file = state_files.get(side) if state_files else None
        lines += tes_side_constants_block(
            side, params, _side_series_file(series_file_base, side), state_file
        )
    return lines


def _side_state_file(mesh_dir_name: str, case_name: str, side: str) -> str:
    """Path (relative to the repo root -- ElmerSolver's cwd, see run.py) for
    one side's persisted circuit state file. Lives in the mesh directory
    alongside the .result restart interface (never collected into
    results/<case>/ by run.py) so it survives across runs the same way."""
    return f"{mesh_dir_name}/{case_name}_{side}.state"


def solver1_block(
    solver: dict[str, Any],
    *,
    solver_index: int = 1,
    equation_name: str = "Heat Equation",
    inner_circuit: bool = False,
    calculate_loads: bool = False,
    lumped_mass: bool = False,
    transient_restart: bool = False,
    comment: str | None = None,
) -> list[str]:
    lines = [
        f"Solver {solver_index}",
        f"  Equation = {equation_name}",
        '  Procedure = "HeatSolve" "HeatSolver"',
        "  Variable = Temperature",
        "  Variable DOFs = 1",
    ]
    if calculate_loads:
        lines.append("  Calculate Loads = True")
    if comment:
        lines.append(f"! {comment}")
    if inner_circuit:
        # Implemented in the custom HeatSolve module.  Unlike an external
        # slave solver, this hook executes within HeatSolve's nonlinear loop
        # and is collective-safe under MPI.
        lines += [
            '  "TES Inner Circuit Update" = Logical True',
            '  "TES Body ID" = Integer 2',
        ]
    lines += [
        f"  Nonlinear System Max Iterations = {solver['nonlinear_max_iterations']}",
        f"  Nonlinear System Convergence Tolerance = {fmt_real(solver['nonlinear_convergence_tolerance'])}",
        f"  Nonlinear System Relaxation Factor = {fmt_real(solver['nonlinear_relaxation_factor'])}",
    ]
    if lumped_mass:
        lines += [
            "! Lumped mass + BDF1 keeps the discrete maximum principle when dt is",
            "! below the element diffusion limit (h^2*rho*c/6k): prevents unphysical",
            "! undershoot around the sharp pulse deposition.",
            "  Lumped Mass Matrix = True",
        ]
    linear_system = solver.get("linear_system", "direct")
    if linear_system == "iterative":
        lines += [
            "  Linear System Solver = Iterative",
            "  Linear System Iterative Method = BiCGStabl",
            "  Linear System Preconditioning = ILU2",
            "  Linear System Max Iterations = 2000",
            "  Linear System Convergence Tolerance = 1.0e-10",
        ]
    elif linear_system == "iterative_hypre_boomeramg":
        # HYPRE is an optional Elmer build dependency.  This configuration is
        # based on Elmer's upstream BoomerAMG regression cases and is intended
        # for the nonsymmetric heat-equation matrix of the TES model.
        lines += [
            "  Linear System Use Hypre = True",
            "  Linear System Solver = Iterative",
            "  Linear System Iterative Method = BiCGStab",
            "  Linear System Preconditioning = BoomerAMG",
            "  Linear System Max Iterations = 1000",
            "  Linear System Convergence Tolerance = 1.0e-10",
            "  Linear System Abort Not Converged = True",
            "  Linear System Residual Output = 20",
            "  BoomerAMG Relax Type = 3",
            "  BoomerAMG Coarsen Type = 0",
            "  BoomerAMG Num Sweeps = 1",
            "  BoomerAMG Max Levels = 25",
            "  BoomerAMG Interpolation Type = 0",
            "  BoomerAMG Smooth Type = 6",
            "  BoomerAMG Cycle Type = 1",
            "  BoomerAMG Num Functions = 1",
            "  BoomerAMG Strong Threshold = 0.25",
        ]
    elif linear_system == "mumps":
        lines += [
            "  Linear System Solver = Direct",
            "  Linear System Direct Method = MUMPS",
        ]
    else:
        lines += [
            "  Linear System Solver = Direct",
            "  Linear System Direct Method = Umfpack",
        ]
    lines += [
        "  Apply Mortar BCs = True",
        f"  Steady State Convergence Tolerance = {fmt_real(solver['steady_state_convergence_tolerance'])}",
    ]
    if transient_restart:
        lines.append("  Transient Restart = Logical True")
    lines.append("End")
    return lines


def vtu_solver_block(case_name: str, vtu: Any, *, solver_index: int = 2) -> list[str]:
    """Build the optional result-output solver.

    ``circuit_inner`` reserves Solver 1 for its nonlinear pre-solver and
    Solver 2 for HeatSolve, so its output solver must use a distinct index.
    """
    lines = [f"Solver {solver_index}"]
    if vtu == "after_simulation":
        lines.append("  Exec Solver = After Simulation")
    elif vtu == "after_timestep":
        lines.append("  Exec Solver = After Timestep")
    elif isinstance(vtu, dict) and "exec_intervals" in vtu:
        intervals = vtu["exec_intervals"]
        lines.append("  Exec Solver = After Timestep")
        lines.append("! per timestep-stage save frequency (same convention as Output Intervals)")
        if len(intervals) == 1:
            lines.append(f"  Exec Intervals = {intervals[0]}")
        else:
            lines.append(
                f"  Exec Intervals({len(intervals)}) = " + " ".join(str(v) for v in intervals)
            )
    else:
        raise ValueError(f"Unsupported vtu spec: {vtu!r}")
    lines += [
        "  Equation = Result Output",
        '  Procedure = "ResultOutputSolve" "ResultOutputSolver"',
        f'  Output File Name = "{case_name}"',
        "  Vtu Format = True",
        "End",
    ]
    return lines


# Body Force name for a dual-TES side's circuit instance (TESTransientHeatSourceL/R
# in tes_transient_heat_source.f90). Natural per-side extension of the
# single-instance HEAT_SOURCES["circuit_implicit"] label.
_SIDE_LABEL = {
    "L": "TES L shunt power from local T (implicit circuit, R(T,I) + inductance)",
    "R": "TES R shunt power from local T (implicit circuit, R(T,I) + inductance)",
}


def body_force_blocks(heat_source: str, tes_body_names: list[str], with_pulse: bool) -> list[str]:
    """Body Force 1..N for the TES bodies (in the same order Body Force
    numbers are assigned to them in bodies_and_bcs -- ascending target id),
    then, if `with_pulse`, one more Body Force for the absorber pulse.

    A single-pixel mesh has exactly one TES body ("TES", unsuffixed) and
    this reduces to the original single Body Force 1 block. A dual-TES mesh
    has two ("TES_L", "TES_R"), each wired to its own UDF instance
    (TESTransientHeatSourceL/R) via the shared circuit_implicit dll.
    """
    dll, unprefixed_proc, unprefixed_label = HEAT_SOURCES[heat_source]
    lines: list[str] = []
    for i, name in enumerate(tes_body_names, start=1):
        side = name[len(_base_body_name(name)):].lstrip("_")
        if side:
            if heat_source != "circuit_implicit":
                raise ValueError(
                    f"{name}: heat_source '{heat_source}' has no per-side TES procedure"
                )
            proc = f"{unprefixed_proc}{side}"
            label = _SIDE_LABEL[side]
        else:
            proc = unprefixed_proc
            label = unprefixed_label
        lines += [
            f"Body Force {i}",
            f'  Name = "{label}"',
            "  Volumetric Heat Source = Variable Temperature",
            f'    Real Procedure "{dll}" "{proc}"',
            "End",
            "",
        ]
    if with_pulse:
        pulse_dll, pulse_proc = PULSE_PROCEDURE
        lines += [
            f"Body Force {len(tes_body_names) + 1}",
            '  Name = "Rectangular window pulse at the absorber center"',
            "  Volumetric Heat Source = Variable Temperature",
            f'    Real Procedure "{pulse_dll}" "{pulse_proc}"',
            "End",
            "",
        ]
    return lines


def _base_body_name(name: str) -> str:
    """Strip a dual-TES `_L`/`_R` stack suffix, if present, from a body
    name (`TES_L` -> `TES`; `abs` -> `abs`, it is never suffixed)."""
    for suffix in ("_L", "_R"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def resolve_bodies(mesh_names: MeshNames) -> list[tuple[int, str, int]]:
    """(target id, body name, material index) for every body in
    *mesh_names*, ordered by target id ascending (mesh.names' own order,
    which is also the order the legacy hand-written SIFs list bodies in)."""
    bodies = []
    for name, target in mesh_names.bodies.items():
        base = _base_body_name(name)
        material = MATERIAL_BY_BASE_NAME.get(base)
        if material is None:
            raise ValueError(f"No material role for body '{name}' (base name '{base}')")
        bodies.append((target, name, material))
    bodies.sort(key=lambda item: item[0])
    return bodies


def resolve_tes_body_names(mesh_names: MeshNames) -> list[str]:
    """TES-role body names ordered by ascending target id -- the order Body
    Force 1..N are assigned to them in bodies_and_bcs (and that
    body_force_blocks must match). ["TES"] for a single-pixel mesh,
    ["TES_L", "TES_R"] for dual-TES."""
    return [name for _, name, _ in resolve_bodies(mesh_names) if _base_body_name(name) == "TES"]


def _present_suffixes(mesh_names: MeshNames, base: str) -> list[str]:
    """Which of the stack suffixes "", "_L", "_R" actually have a
    `<base><suffix>` body in this mesh."""
    return [sfx for sfx in ("", "_L", "_R") if f"{base}{sfx}" in mesh_names.bodies]


def resolve_bath_boundaries(mesh_names: MeshNames) -> list[int]:
    """Target IDs of every `SiO2_2<sfx>__zmin` boundary (the bath contact),
    ascending. A single-pixel mesh has one; a dual-TES mesh has one per
    stack."""
    targets = [
        mesh_names.boundaries[f"SiO2_2{sfx}__zmin"]
        for sfx in _present_suffixes(mesh_names, "SiO2_2")
    ]
    return sorted(targets)


def resolve_mortar_pairs(mesh_names: MeshNames) -> list[tuple[int, str, int, str]]:
    """(slave target, slave label, master target, master label) for every
    mortar pair, expanded over whichever of the "", "_L", "_R" stack
    suffixes are actually present in the mesh -- looked up by boundary
    *name*, not by the mesh's ID numbering convention or its spurious
    `_free` boundary groups (see docs/dual_tes_plan.md Phase A notes).

    Ordered by pair (TES/Membrane_SiNx, Stycast/TES, abs/Stycast) then by
    suffix within a pair, matching the legacy single-pixel SIF's BC order.
    """
    suffixes = _present_suffixes(mesh_names, "TES")
    pairs = []
    for slave_base, slave_face, master_base, master_face, master_suffixed in _MORTAR_PAIRS:
        for sfx in suffixes:
            master_sfx = sfx if master_suffixed else ""
            slave_name = f"{slave_base}{sfx}__{slave_face}"
            master_name = f"{master_base}{master_sfx}__{master_face}"
            slave_label = (
                f"{_MORTAR_LABEL_BASE.get(slave_base, slave_base)}{sfx} "
                f"{_FACE_LABEL[slave_face]} mortar"
            )
            master_label = (
                f"{_MORTAR_LABEL_BASE.get(master_base, master_base)}{master_sfx} "
                f"{_FACE_LABEL[master_face]} mortar"
            )
            pairs.append(
                (
                    mesh_names.boundaries[slave_name],
                    slave_label,
                    mesh_names.boundaries[master_name],
                    master_label,
                )
            )
    return pairs


def _target_boundaries_line(targets: list[int]) -> str:
    return f"  Target Boundaries({len(targets)}) = " + " ".join(str(t) for t in targets)


def bodies_and_bcs(mesh_names: MeshNames, with_pulse: bool) -> list[str]:
    bodies = resolve_bodies(mesh_names)
    # TES-role bodies get Body Force 1..N (ascending target id, i.e. one per
    # stack); the pulse (if any) is appended after them on the abs body.
    tes_body_names = resolve_tes_body_names(mesh_names)
    body_force_by_name = {name: i for i, name in enumerate(tes_body_names, start=1)}
    pulse_body_force = len(tes_body_names) + 1

    lines: list[str] = []
    for i, (target, name, material) in enumerate(bodies, start=1):
        lines += [
            f"Body {i}",
            f"  Target Bodies(1) = {target}",
            f'  Name = "{name}"',
            "  Equation = 1",
            f"  Material = {material}",
        ]
        if name in body_force_by_name:
            lines.append(f"  Body Force = {body_force_by_name[name]}")
        if _base_body_name(name) == "abs" and with_pulse:
            lines.append(f"  Body Force = {pulse_body_force}")
        lines += ["  Initial Condition = 1", "End", ""]

    bc_blocks: list[list[str]] = [
        [
            "Boundary Condition 1",
            _target_boundaries_line(resolve_bath_boundaries(mesh_names)),
            '  Name = "bath on SiO2_2__zmin"',
            "  Temperature = __T_BATH__",
            "End",
        ]
    ]
    bc_number = 2
    master_bc_of: dict[int, int] = {}
    for slave_target, slave_label, master_target, master_label in resolve_mortar_pairs(mesh_names):
        slave_bc = bc_number
        bc_number += 1
        is_new_master = master_target not in master_bc_of
        if is_new_master:
            master_bc_of[master_target] = bc_number
            bc_number += 1
        master_bc = master_bc_of[master_target]

        bc_blocks.append(
            [
                f"Boundary Condition {slave_bc}",
                _target_boundaries_line([slave_target]),
                f'  Name = "{slave_label}"',
                f"  Mortar BC = {master_bc}",
                "  Galerkin Projector = True",
                "  Plane Projector = True",
                "End",
            ]
        )
        if is_new_master:
            bc_blocks.append(
                [
                    f"Boundary Condition {master_bc}",
                    _target_boundaries_line([master_target]),
                    f'  Name = "{master_label}"',
                    "End",
                ]
            )

    for i, block in enumerate(bc_blocks):
        if i > 0:
            lines.append("")
        lines += block
    return lines


def _initial_temperature(spec: dict, params: dict[str, float]) -> float:
    value = spec.get("initial_temperature", "T_0")
    if isinstance(value, str) and value in params:
        return float(params[value])
    return eval_si(value, params)


def _timestep_lines(spec: dict, params: dict[str, float]) -> list[str]:
    stages = spec["timesteps"]
    sizes = [eval_si(size, params) for size, _ in stages]
    counts = [int(count) for _, count in stages]
    n = len(stages)
    lines = [
        f"  Timestep Sizes({n}) = " + " ".join(fmt(s) for s in sizes),
        f"  Timestep Intervals({n}) = " + " ".join(str(c) for c in counts),
    ]
    outputs = spec.get("output_intervals")
    if outputs is not None:
        if len(outputs) != n:
            raise ValueError(
                f"output_intervals has {len(outputs)} stages but timesteps has {n}"
            )
        if n == 1:
            lines.append(f"  Output Intervals = {outputs[0]}")
        else:
            lines.append(f"  Output Intervals({n}) = " + " ".join(str(v) for v in outputs))
    return lines


def _resolve_pulse_center(
    center_spec: Any, params: dict[str, float], mesh_dir: Path, mesh_name: str
) -> tuple[tuple[float, float, float], str]:
    """Resolve a pulse `center` spec to (x, y, z) plus a SIF comment noting
    provenance. Accepts:
    - "auto": the absorber centroid of the mesh (all three components).
    - {"x": ..., "y": ..., "z": ...}: per-component, where each value is
      either "auto" (that component of the absorber centroid) or a
      dimensioned expression string evaluated against the project
      parameters (dual-TES offset-pulse cases use this to fix x while
      leaving y/z at their auto centroid value).
    """
    axes = ("x", "y", "z")

    def centroid() -> tuple[float, float, float]:
        # Round to 12 significant digits: the centroid is geometry, not noise.
        return tuple(float(f"{c:.12g}") for c in absorber_centroid(mesh_dir))

    if center_spec == "auto":
        return centroid(), f"auto (absorber centroid of {mesh_name})"

    if not isinstance(center_spec, dict):
        raise ValueError(f"Unsupported pulse center spec: {center_spec!r}")

    cached_centroid: tuple[float, float, float] | None = None
    comps: list[float] = []
    notes: list[str] = []
    for i, axis in enumerate(axes):
        value = center_spec.get(axis, "auto")
        if value == "auto":
            if cached_centroid is None:
                cached_centroid = centroid()
            comps.append(cached_centroid[i])
            notes.append(f"{axis}=auto (absorber centroid)")
        else:
            comps.append(eval_si(value, params))
            notes.append(f"{axis}={value}")
    return (comps[0], comps[1], comps[2]), "; ".join(notes)


def _pulse_constants(
    case_name: str, spec: dict, params: dict[str, float], root: Path, mesh_dir_name: str
) -> list[str]:
    pulse = spec["pulse"]
    mesh_dir = root / mesh_dir_name
    sigma = eval_si(pulse["sigma"], params)
    center, center_note = _resolve_pulse_center(
        pulse.get("center", "auto"), params, mesh_dir, spec["mesh"]
    )
    norm = gaussian_discrete_norm(mesh_dir, center, sigma)
    energy = eval_si(pulse["energy"], params)
    return [
        f"  Pulse Energy = Real {fmt(energy)}",
        f"  Pulse Start Time = Real {fmt(eval_si(pulse['start'], params))}",
        f"  Pulse Duration = Real {fmt(eval_si(pulse['duration'], params))}",
        f"  Pulse Sigma = Real {fmt(sigma)}",
        f"! Pulse center: {center_note}",
        f"  Pulse Center X = Real {fmt(center[0])}",
        f"  Pulse Center Y = Real {fmt(center[1])}",
        f"  Pulse Center Z = Real {fmt(center[2])}",
        f"! FE integral of the nodal Gaussian over the absorber of {spec['mesh']}",
        "! (recomputed automatically at build time)",
        f"  Pulse Discrete Norm = Real {fmt(norm)}",
    ]


def build_case(case_name: str, spec: dict, model: dict, root: Path) -> str:
    params = model["parameters"]
    template = spec["template"]
    meshes = model.get("meshes", {})
    if spec["mesh"] not in meshes:
        raise ValueError(f"{case_name}: mesh '{spec['mesh']}' is not in the meshes registry")
    mesh_dir_name = meshes[spec["mesh"]]["dir"]
    if not (root / mesh_dir_name / "mesh.header").exists():
        raise ValueError(
            f"{case_name}: mesh directory '{mesh_dir_name}' has no mesh.header "
            "(build it with: python build_mesh.py " + spec["mesh"] + ")"
        )
    heat_source = spec.get(
        "heat_source", "circuit_implicit" if template in ("transient", "pulse") else "circuit_local"
    )
    with_pulse = template == "pulse"
    series_file = spec.get("series_file")
    if heat_source == "circuit_implicit" and not series_file:
        raise ValueError(f"{case_name}: circuit_implicit cases need a series_file")

    mesh_names = parse_mesh_names(root / mesh_dir_name / "mesh.names")
    tes_body_names = resolve_tes_body_names(mesh_names)
    # "" for a single-pixel mesh's unsuffixed "TES" body, "L"/"R" for a
    # dual-TES mesh's "TES_L"/"TES_R" bodies (see body_force_blocks).
    tes_sides = [name[len(_base_body_name(name)):].lstrip("_") for name in tes_body_names]
    is_dual_tes = tes_sides != [""]

    lines: list[str] = [
        f"! Auto-generated by scripts/support/build_cases.py from elmer_project.json",
        f"! (case '{case_name}'). Do not edit; edit the case definition instead.",
        "",
        "Header",
        "  CHECK KEYWORDS Warn",
        f'  Mesh DB "." "{mesh_dir_name}"',
        "End",
        "",
        "Simulation",
        "  Max Output Level = 5",
        "  Coordinate System = Cartesian 3D",
    ]

    if template == "steady":
        lines += [
            "  Simulation Type = Steady State",
            f"  Steady State Max Iterations = {spec.get('steady_state_max_iterations', 1)}",
            f"  Output Intervals = {spec.get('output_intervals', 1)}",
        ]
    elif template in ("transient", "pulse"):
        lines += [
            "  Simulation Type = Transient",
            "  Timestepping Method = BDF",
            "  BDF Order = 1",
        ]
        lines += _timestep_lines(spec, params)
        lines.append(
            f"  Steady State Max Iterations = {spec.get('steady_state_max_iterations', 1)}"
        )
    else:
        raise ValueError(f"{case_name}: unknown template {template!r}")

    lines.append(f"  Solver Input File = generated/cases/{case_name}.sif")

    restart_from = spec.get("restart_from")
    if restart_from:
        dep = model["cases"].get(restart_from)
        if not dep or not dep.get("output_result"):
            raise ValueError(
                f"{case_name}: restart_from '{restart_from}' must exist and set output_result"
            )
        lines += [
            f"! Run {restart_from} first to produce the restart field.",
            f"  Restart File = {restart_from}.result",
            f"  Restart Position = {spec.get('restart_position', 0)}",
        ]
        if "restart_time" in spec:
            lines.append(f"  Restart Time = Real {fmt(float(spec['restart_time']))}")
    if spec.get("output_result"):
        lines.append(f"  Output File = {case_name}.result")
    if spec.get("post_file"):
        lines.append(f"  Post File = {case_name}.ep")
    lines += ["End", ""]

    lines += [
        "Initial Condition 1",
        f"  Temperature = {fmt(_initial_temperature(spec, params))}",
        "End",
        "",
    ]

    if is_dual_tes:
        # State File Constants are dual-only (docs/dual_tes_plan.md): a dual
        # steady case writes its converged circuit state; a dual transient/
        # pulse case that restarts from another case reads that case's
        # written state. Anything else (a dual case with neither trait) gets
        # none, which keeps the UDF on its T0-fallback init path.
        state_files: dict[str, str] | None = None
        if restart_from:
            state_files = {
                side: _side_state_file(mesh_dir_name, restart_from, side) for side in tes_sides
            }
        elif template == "steady":
            state_files = {
                side: _side_state_file(mesh_dir_name, case_name, side) for side in tes_sides
            }
        lines += dual_tes_constants_block(tes_sides, params, series_file, state_files)
    else:
        lines += tes_constants_block(params, series_file, spec.get("state_file"))
    if with_pulse:
        lines += _pulse_constants(case_name, spec, params, root, mesh_dir_name)
    lines += ["End", ""]

    parallel_circuit_iterations = int(spec.get("parallel_circuit_iterations", 1))
    if parallel_circuit_iterations < 1:
        raise ValueError(f"{case_name}: parallel_circuit_iterations must be >= 1")
    if heat_source == "circuit_parallel":
        if is_dual_tes:
            raise ValueError(f"{case_name}: circuit_parallel prototype supports one TES only")
        circuit_relaxation = float(spec.get("parallel_circuit_relaxation", 0.04))
        for coupling_iter in range(parallel_circuit_iterations):
            circuit_index = 2 * coupling_iter + 1
            heat_index = circuit_index + 1
            lines += [
                f"Solver {circuit_index}",
                f'  Equation = "TES Parallel Circuit {coupling_iter + 1}"',
                '  Procedure = "tes_parallel_circuit" "TESParallelCircuitSolver"',
                '  "TES Body ID" = Integer 2',
                f'  "TES Circuit Relaxation" = Real {fmt_real(circuit_relaxation)}',
                f'  "TES Write Series" = Logical {"True" if coupling_iter == parallel_circuit_iterations - 1 else "False"}',
                "  Exec Solver = Always",
                "End",
                "",
            ]
            lines += solver1_block(
                spec["solver"],
                solver_index=heat_index,
                equation_name=f"Heat Equation Coupling {coupling_iter + 1}",
                calculate_loads=bool(spec.get("calculate_loads")),
                lumped_mass=bool(spec.get("lumped_mass")),
                transient_restart=bool(spec.get("transient_restart")),
                comment=spec.get("solver_comment"),
            )
            lines.append("")
    elif heat_source == "circuit_inner":
        if is_dual_tes:
            raise ValueError(f"{case_name}: circuit_inner supports one TES only")
        lines += solver1_block(
            spec["solver"],
            inner_circuit=True,
            calculate_loads=bool(spec.get("calculate_loads")),
            lumped_mass=bool(spec.get("lumped_mass")),
            transient_restart=bool(spec.get("transient_restart")),
            comment=spec.get("solver_comment"),
        )
        lines.append("")
    else:
        lines += solver1_block(
            spec["solver"],
            calculate_loads=bool(spec.get("calculate_loads")),
            lumped_mass=bool(spec.get("lumped_mass")),
            transient_restart=bool(spec.get("transient_restart")),
            comment=spec.get("solver_comment"),
        )
        lines.append("")

    vtu_default = "after_simulation" if template == "steady" else "after_timestep"
    vtu_spec = spec.get("vtu", vtu_default)
    if vtu_spec:
        lines += vtu_solver_block(case_name, vtu_spec)
        lines.append("")
    if heat_source == "circuit_parallel":
        active = " ".join(str(i) for i in range(1, 2 * parallel_circuit_iterations + 1))
        lines += ["Equation 1", '  Name = "Heat"', f"  Active Solvers({2 * parallel_circuit_iterations}) = {active}", "End", ""]
    else:
        lines += ["Equation 1", '  Name = "Heat"', "  Active Solvers(1) = 1", "End", ""]

    lines += materials_block(model)
    lines += body_force_blocks(heat_source, tes_body_names, with_pulse)
    body_lines = bodies_and_bcs(mesh_names, with_pulse)
    body_lines = [
        line.replace("__T_BATH__", fmt(params["T_bath"])) for line in body_lines
    ]
    lines += body_lines
    lines.append("")
    return "\n".join(lines)


def build_all_cases(model: dict, out_dir: Path, root: Path) -> list[str]:
    cases = model.get("cases", {})
    if not cases:
        raise ValueError("elmer_project.json has no 'cases' section")
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for case_name, spec in cases.items():
        text = build_case(case_name, spec, model, root)
        (out_dir / f"{case_name}.sif").write_text(text, encoding="utf-8", newline="\n")
        written.append(case_name)
    return written
