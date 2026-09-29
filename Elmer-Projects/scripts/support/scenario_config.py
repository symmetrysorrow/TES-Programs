"""Compile a user project TOML and its reusable model into solver projects.

The user-facing format deliberately has two levels: projects/*.toml owns run
conditions, while projects/models/*.toml owns geometry, materials, layout and
mesh. Generated JSON remains an internal adapter for the existing SIF/run
pipeline.
"""

from __future__ import annotations

import copy
import hashlib
import json
import re
import tomllib
from pathlib import Path

from scripts.support.reconcile_project import reconcile_project
from scripts.support.vendored.dimensioned_expression import evaluate_dimensioned_expression

ROOT = Path(__file__).resolve().parents[2]
BASE_PROJECT = ROOT / "artifacts" / "phase24_conformal_hybrid" / "project.json"
PROJECT_DIR = ROOT / "projects"
MODEL_DIR = PROJECT_DIR / "models"
ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
CIRCUIT_KEYS = {"I_bias", "R_sh", "L_tes", "R_0", "R_min", "alpha", "beta", "I_0", "T_c", "T_0", "TES_volume"}
PROJECT_CIRCUIT_KEYS = {"I_bias", "R_sh", "L_tes", "R_0", "alpha", "beta"}
MODEL_FORBIDDEN_PARAMETERS = PROJECT_CIRCUIT_KEYS | {"I_0"}


def read_toml(path: Path) -> dict:
    if not path.is_file():
        raise ValueError(f"TOML file not found: {path}")
    return tomllib.loads(path.read_text(encoding="utf-8"))


def require_id(value: object, label: str) -> str:
    if not isinstance(value, str) or not ID_RE.fullmatch(value):
        raise ValueError(f"{label} must be an identifier (letters, digits, underscore): {value!r}")
    return value


def reject_unknown(data: dict, allowed: set[str], label: str) -> None:
    extra = set(data) - allowed
    if extra:
        raise ValueError(f"{label} has unknown keys: {sorted(extra)}")


def quantity(value: object, context: dict[str, float], label: str, *, positive: bool = False) -> float:
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        raise ValueError(f"{label} must be a dimensioned expression")
    try:
        result = float(evaluate_dimensioned_expression(str(value), variables=context).value)
    except Exception as exc:
        raise ValueError(f"{label}: cannot evaluate {value!r}: {exc}") from exc
    if positive and result <= 0:
        raise ValueError(f"{label} must be positive")
    return result


def digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()[:12]


def resolve_model_file(project_file: Path, value: object) -> Path:
    """Resolve the explicit model path written in a project."""
    if not isinstance(value, str) or not value.strip():
        raise ValueError("model must be a relative TOML path such as models/single_pixel.toml")
    ref = Path(value)
    if ref.is_absolute() or ref.suffix.lower() != ".toml":
        raise ValueError("model must be a relative .toml path")
    model_file = (project_file.parent / ref).resolve()
    if not model_file.is_file():
        raise ValueError(f"model TOML not found: {model_file}")
    return model_file


def load_scenario(project_file: Path) -> dict:
    """Resolve one user project into stable TES positions and final run values."""
    project = read_toml(project_file)
    reject_unknown(
        project,
        {"schema_version", "name", "model", "circuits", "pulse", "solver", "timesteps", "window", "save_field"},
        "project",
    )
    if project.get("schema_version") != 3:
        raise ValueError("new project TOML needs schema_version = 3")

    name = require_id(project.get("name", project_file.stem), "project name")
    model_file = resolve_model_file(project_file, project.get("model"))
    config = read_toml(model_file)
    reject_unknown(config, {"model_version", "parameters", "materials", "components", "layout", "mesh"}, "model")
    if config.get("model_version") != 2:
        raise ValueError(f"{model_file} needs model_version = 2")

    layout, mesh = config.get("layout", {}), config.get("mesh", {})
    defaults = {"parameters": config.get("parameters", {}), "materials": config.get("materials", {})}
    components = config.get("components", {})
    if not all(isinstance(v, dict) for v in (layout, mesh, defaults["parameters"], defaults["materials"], components)):
        raise ValueError("model sections must be tables")

    duplicated = MODEL_FORBIDDEN_PARAMETERS & set(defaults["parameters"])
    if duplicated:
        raise ValueError(
            "model contains project-owned circuit keys; move them to projects/*.toml: "
            f"{sorted(duplicated)}"
        )

    reject_unknown(layout, {"absorber_length", "absorber_width", "stycast_diameter", "dimensions", "groups"}, "layout")
    reject_unknown(mesh, {"source", "global_size", "tes_local_size", "absorber_local_size", "absorber_local_radius", "stycast_layers", "si_1_layers", "sio2_1_layers", "sinx_layers"}, "mesh")

    circuits = project.get("circuits")
    if not isinstance(circuits, dict) or not circuits:
        raise ValueError("project needs a [circuits.<TES ID>] table for every TES")
    if not isinstance(project.get("solver", {}), dict):
        raise ValueError("solver must be a table")

    base = json.loads(BASE_PROJECT.read_text(encoding="utf-8"))
    original_geometry_parameters = reconcile_project(base)["parameters"]
    expressions = {key: str(value) for key, value in defaults["parameters"].items()}
    base["parameter_expressions"] = dict(expressions)
    for material, props in defaults["materials"].items():
        if material not in base["materials"]:
            raise ValueError(f"unknown model material: {material}")
        for prop, expression in props.items():
            base["materials"][material][prop]["expression"] = str(expression)

    dimensions = layout.get("dimensions", {})
    if not isinstance(dimensions, dict) or any(not k.endswith(("_dx", "_dy", "_dz")) for k in dimensions):
        raise ValueError("layout [dimensions] may contain only *_dx, *_dy, *_dz expressions")
    expressions.update({k: str(v) for k, v in dimensions.items()})
    context = reconcile_project({"parameter_expressions": expressions})["parameters"]

    groups = layout.get("groups", [])
    if not isinstance(groups, list) or not groups:
        raise ValueError("layout needs at least one [[groups]] entry")

    tes = []
    group_ids = set()
    for group in groups:
        reject_unknown(group, {"id", "component", "origin"}, "group")
        group_id = require_id(group.get("id"), "group id")
        if group_id in group_ids:
            raise ValueError(f"duplicate group id: {group_id}")
        group_ids.add(group_id)
        component_name = require_id(group.get("component"), f"group {group_id} component")
        if component_name not in components:
            raise ValueError(f"unknown component: {component_name}")
        component = components[component_name]
        reject_unknown(component, {"tes"}, f"component {component_name}")
        if not isinstance(component.get("tes"), list) or not component["tes"]:
            raise ValueError(f"component {component_name} needs at least one [[tes]]")
        origin = group.get("origin", ["0[mm]", "0[mm]"])
        if not isinstance(origin, list) or len(origin) != 2:
            raise ValueError(f"group {group_id} origin needs [x, y]")
        gx, gy = (quantity(v, context, f"group {group_id} origin {axis}") for v, axis in zip(origin, "xy"))
        for item in component.get("tes", []):
            reject_unknown(item, {"id", "position"}, f"component {component_name} TES")
            local_id = require_id(item.get("id"), f"component {component_name} TES id")
            position = item.get("position")
            if not isinstance(position, list) or len(position) != 2:
                raise ValueError(f"component {component_name}.{local_id} position needs [x, y]")
            x = gx + quantity(position[0], context, f"{group_id}.{local_id} x")
            y = gy + quantity(position[1], context, f"{group_id}.{local_id} y")
            tes.append({"id": f"{group_id}.{local_id}", "x": x, "y": y})

    if not 1 <= len(tes) <= 8:
        raise ValueError("the current numbered circuit UDF supports 1 to 8 TES")
    if len({(t["x"], t["y"]) for t in tes}) != len(tes):
        raise ValueError("two TES stacks have the same position")

    expected_ids = {t["id"] for t in tes}
    actual_ids = set(circuits)
    unknown_ids = actual_ids - expected_ids
    missing_ids = expected_ids - actual_ids
    if unknown_ids:
        raise ValueError(
            f"unknown circuit TES IDs: {sorted(unknown_ids)}; available IDs: {sorted(expected_ids)}"
        )
    if missing_ids:
        raise ValueError(f"missing circuit TES IDs: {sorted(missing_ids)}")

    for tes_item in tes:
        tes_id = tes_item["id"]
        values = circuits[tes_id]
        if not isinstance(values, dict):
            raise ValueError(f"circuits.{tes_id} must be a table")
        reject_unknown(values, PROJECT_CIRCUIT_KEYS, f"circuits.{tes_id}")
        missing = PROJECT_CIRCUIT_KEYS - set(values)
        if missing:
            raise ValueError(f"circuits.{tes_id} is missing final values: {sorted(missing)}")
        for key, value in values.items():
            quantity(value, context, f"circuits.{tes_id}.{key}")
        tes_item["circuit"] = dict(values)
        tes_item["sources"] = {key: f"project {project_file.name}" for key in values}

    pulse = project.get("pulse", {})
    if not isinstance(pulse, dict):
        raise ValueError("pulse must be a table")
    if pulse:
        reject_unknown(pulse, {"energy", "start", "duration", "sigma", "center", "radius"}, "pulse")
        missing = {"energy", "start", "duration", "sigma", "center"} - set(pulse)
        if missing:
            raise ValueError(f"pulse is missing: {sorted(missing)}")

    tes.sort(key=lambda t: (t["x"], t["y"], t["id"]))
    chip_dx, chip_dy = context["Si_dx"], context["Si_dy"]
    for i, left in enumerate(tes):
        for right in tes[i + 1:]:
            if abs(left["x"] - right["x"]) < chip_dx and abs(left["y"] - right["y"]) < chip_dy:
                raise ValueError(f"TES chip stacks overlap: {left['id']} and {right['id']}")

    geometry = {
        "tes": [{"id": t["id"], "x": t["x"], "y": t["y"]} for t in tes],
        "absorber_length": quantity(layout.get("absorber_length"), context, "absorber_length", positive=True),
        "absorber_width": quantity(layout.get("absorber_width", "abs_dy"), context, "absorber_width", positive=True),
    }
    contact_diameter = quantity(layout.get("stycast_diameter", "498[um]"), context, "stycast_diameter", positive=True)
    if geometry["absorber_length"] < tes[-1]["x"] - tes[0]["x"] + contact_diameter:
        raise ValueError("absorber_length does not cover every TES contact")
    ys = [item["y"] for item in tes]
    if geometry["absorber_width"] < max(ys) - min(ys) + contact_diameter:
        raise ValueError("absorber_width does not cover every TES contact")
    geometry["stycast_diameter"] = contact_diameter

    source_mesh = mesh.get("source")
    if source_mesh is not None:
        if source_mesh != "mesh_hybrid_fullconf_h8" or source_mesh not in base["meshes"]:
            raise ValueError(f"unsupported prebuilt mesh: {source_mesh!r}")
        dimension_keys = {k for k in context if k.endswith(("_dx", "_dy", "_dz"))}
        if any(abs(context[k] - original_geometry_parameters[k]) > 1e-12 for k in dimension_keys):
            raise ValueError("prebuilt single-pixel mesh requires its original chip dimensions")
        if (
            len(tes) != 1
            or abs(tes[0]["x"]) > 1e-12
            or abs(tes[0]["y"]) > 1e-12
            or abs(geometry["absorber_length"] - 1e-3) > 1e-12
            or abs(geometry["absorber_width"] - 1e-3) > 1e-12
            or abs(contact_diameter - 498e-6) > 1e-12
        ):
            raise ValueError("prebuilt single-pixel mesh requires one centered TES and a 1 mm absorber")

    mesh_values = {
        "global_size": quantity(mesh.get("global_size", "50[um]"), context, "global_size", positive=True),
        "tes_local_size": quantity(mesh.get("tes_local_size", "8[um]"), context, "tes_local_size", positive=True),
        "absorber_local_size": quantity(mesh.get("absorber_local_size", "35[um]"), context, "absorber_local_size", positive=True),
        "absorber_local_radius": quantity(mesh.get("absorber_local_radius", "50[um]"), context, "absorber_local_radius", positive=True),
        "stycast_layers": int(mesh.get("stycast_layers", 32)),
        "si_1_layers": int(mesh.get("si_1_layers", 4)),
        "sio2_1_layers": int(mesh.get("sio2_1_layers", 2)),
        "sinx_layers": int(mesh.get("sinx_layers", 2)),
    }
    if mesh_values["tes_local_size"] > mesh_values["global_size"]:
        raise ValueError("tes_local_size must not exceed global_size")
    if mesh_values["absorber_local_size"] > mesh_values["global_size"]:
        raise ValueError("absorber_local_size must not exceed global_size")
    if any(mesh_values[k] < 1 for k in ("stycast_layers", "si_1_layers", "sio2_1_layers", "sinx_layers")):
        raise ValueError("mesh layer counts must be positive")

    physics = {"parameters": expressions, "tes_circuits": [t["circuit"] for t in tes], "pulse": project.get("pulse", {})}
    geometry_params = {k: v for k, v in context.items() if k.endswith(("_dx", "_dy", "_dz"))}
    geometry_identity = {k: v for k, v in geometry.items() if k != "stycast_diameter"}
    if contact_diameter != 498.0e-6:
        geometry_identity["stycast_diameter"] = contact_diameter
    mesh_hash = digest({"geometry": geometry_identity, "mesh": mesh_values if source_mesh is None else {"source": source_mesh}, "dimensions": geometry_params})
    steady_hash = digest({"geometry_mesh": mesh_hash, "parameters": expressions, "tes_circuits": physics["tes_circuits"], "materials": base["materials"], "solver": project.get("solver", {})})
    pulse_hash = digest({"steady": steady_hash, "pulse": physics["pulse"], "timesteps": project.get("timesteps"), "window": project.get("window"), "save_field": bool(project.get("save_field", False))})
    return {"name": name, "project_file": project_file, "model_file": model_file, "project": project, "base": base, "context": context, "tes": tes, "geometry": geometry, "mesh": mesh_values, "source_mesh": source_mesh, "physics": physics, "mesh_hash": mesh_hash, "steady_hash": steady_hash, "pulse_hash": pulse_hash}

def compile_project(scenario: dict) -> tuple[dict, str, str, str]:
    name = scenario["name"]
    mesh_name = scenario["source_mesh"] or f"mesh_toml_{scenario['mesh_hash']}"
    steady_name = f"case_toml_{scenario['steady_hash']}_steady"
    pulse_name = f"case_{name}_{scenario['pulse_hash']}_pulse"
    project = copy.deepcopy(scenario["base"])
    project["parameter_expressions"] = scenario["physics"]["parameters"]
    project["geometries"]["toml_layout"] = {"name": "toml_layout", "kind": "layout", "children": []}
    if scenario["source_mesh"]:
        project["meshes"] = {mesh_name: copy.deepcopy(scenario["base"]["meshes"][mesh_name])}
    else:
        project["meshes"] = {mesh_name: {"geometry": "toml_layout", "dir": mesh_name, "notes": f"TOML mesh {scenario['mesh_hash']}", "recipe": {"generator": "generate_hybrid_prism_geometry.py"}}}
    steady = copy.deepcopy(scenario["base"]["cases"]["case_conformal_hybrid_steady_base"])
    circuits = []
    for circuit in scenario["physics"]["tes_circuits"]:
        values = dict(circuit)
        p = dict(scenario["context"])
        p.update({key: quantity(value, p, key) for key, value in values.items()})
        values["I_0"] = p["I_bias"] * p["R_sh"] / (p["R_0"] + p["R_sh"])
        circuits.append(values)
    steady.update(mesh=mesh_name, series_file=f"{steady_name}_series.csv", iteration_series_file=f"{steady_name}_iterations.csv", state_file=f"work/meshes/{mesh_name}/{steady_name}.state", output_file_path=f"../work/meshes/{mesh_name}/{steady_name}.result", tes_circuits=circuits)
    steady["solver"].update(scenario["project"].get("solver", {}))
    project["cases"] = {steady_name: steady}
    if scenario["physics"]["pulse"]:
        pulse = copy.deepcopy(scenario["base"]["cases"]["case_conformal_hybrid_steady_base"])
        if scenario["project"].get("window") == "75ms_bdf2h15":
            from scripts.support.run_phase24_gate4_5_nomortar import trim_schedule
            timesteps = trim_schedule([], "75ms_bdf2h15")
        elif scenario["project"].get("window"):
            raise ValueError(f"unknown time window: {scenario['case']['window']}")
        else:
            timesteps = scenario["project"].get("timesteps", [["1[ms]", 20], ["10[us]", 2], ["1[ns]", 1], ["10[ns]", 10], ["100[ns]", 9], ["1[us]", 9], ["5[us]", 20]])
        pulse.update(template="pulse", mesh=mesh_name, restart_file_base=steady_name, restart_time=0.0, preexisting_restart=True, heat_source="circuit_inner", series_file=f"{pulse_name}_series.csv", iteration_series_file=f"{pulse_name}_iterations.csv", state_file=f"work/meshes/{mesh_name}/{steady_name}.state", output_file_path=f"../work/meshes/{mesh_name}/{pulse_name}.result", tes_circuits=circuits, pulse=scenario["physics"]["pulse"], lumped_mass=True, vtu=False, timesteps=timesteps, output_intervals=[1] * len(timesteps))
        pulse["output_result"] = bool(scenario["project"].get("save_field", False))
        if not pulse["output_result"]:
            pulse.pop("output_file_path", None)
        pulse["solver"].update(scenario["project"].get("solver", {}))
        if scenario["project"].get("window"):
            pulse["bdf_order"] = 2
            pulse["phase24_bdf2_predictor"] = True
        project["cases"][pulse_name] = pulse
    return project, mesh_name, steady_name, pulse_name


def project_output(scenario: dict) -> Path:
    return ROOT / "generated" / "scenarios" / scenario["name"] / f"{scenario['pulse_hash']}.json"


def write_project(scenario: dict, *, steady_only: bool = False) -> tuple[Path, str, str, str]:
    project, mesh, steady, pulse = compile_project(scenario)
    if steady_only:
        project["cases"] = {steady: project["cases"][steady]}
        path = ROOT / "generated" / "scenarios" / f"steady_{scenario['steady_hash']}.json"
    else:
        project["cases"] = {pulse: project["cases"][pulse]} if pulse in project["cases"] else project["cases"]
        path = project_output(scenario)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(project, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return path, mesh, steady, pulse
