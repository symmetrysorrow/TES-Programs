"""Create an isolated project for the mixed tetra/prism single-pixel mesh.

The established all-tetrahedral meshes and cases are left untouched.  The new
steady and pulse cases use `mesh_hybrid_abs_tet_layers_prism_pilot` only.
"""
from __future__ import annotations

import copy
import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "elmer_project_comsol_timegrid.json"
OUT = ROOT / "elmer_project_hybrid_prism.json"

MESH = "mesh_hybrid_abs_tet_layers_prism_conformal"
MESH_MPI = "mesh_hybrid_abs_tet_layers_prism_conformal_repart_x"
MESH_MPI_NOBC = "mesh_hybrid_abs_tet_layers_prism_conformal_repart_x_nobcoptim"
STEADY = "case_tes_steady_hybrid_prism"
PULSE = "case_tes_pulse_hybrid_prism_fast_compare"
STEADY_MPI = "case_tes_steady_hybrid_prism_mpi4"
PULSE_MPI = "case_tes_pulse_hybrid_prism_mpi4_regression"
STEADY_INNER_1 = "case_tes_steady_hybrid_prism_inner_1rank"
PULSE_INNER_1 = "case_tes_pulse_hybrid_prism_inner_1rank_regression"
STEADY_MPI_TIGHT = "case_tes_steady_hybrid_prism_mpi4_tight"
STEADY_INNER_1_TIGHT = "case_tes_steady_hybrid_prism_inner_1rank_tight"
PULSE_MPI_TIGHT = "case_tes_pulse_hybrid_prism_mpi4_tight"
PULSE_INNER_1_TIGHT = "case_tes_pulse_hybrid_prism_inner_1rank_tight"
STEADY_MPI_NOBC_TIGHT = "case_tes_steady_hybrid_prism_mpi4_nobcoptim_tight"
PULSE_SERIAL_MAP_1 = "case_tes_pulse_hybrid_prism_from_20ms_serial_tight_1rank"
PULSE_SERIAL_MAP_4 = "case_tes_pulse_hybrid_prism_from_20ms_serial_tight_mpi4"
PULSE_SERIAL_MAP_1_MIN3 = PULSE_SERIAL_MAP_1 + "_min3"
PULSE_SERIAL_MAP_4_MIN3 = PULSE_SERIAL_MAP_4 + "_min3"
PULSE_SERIAL_MAP_1_MIN5 = PULSE_SERIAL_MAP_1 + "_min5"
PULSE_SERIAL_MAP_4_MIN5 = PULSE_SERIAL_MAP_4 + "_min5"
PULSE_SERIAL_MAP_1_COMMIT = PULSE_SERIAL_MAP_1 + "_step_commit"
PULSE_SERIAL_MAP_4_COMMIT = PULSE_SERIAL_MAP_4 + "_step_commit"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT, help="project JSON to create (default: elmer_project_hybrid_prism.json)")
    args = parser.parse_args()
    project = json.loads(SOURCE.read_text(encoding="utf-8"))
    project["meshes"][MESH] = {
        "geometry": "single_pixel",
        "dir": MESH,
        "recipe": {
            "generator": "generate_hybrid_prism_geometry.py",
            "mesh_overrides": {},
            "elmergrid_args": [
                "14", "2", "gmsh/project_hybrid_prism.msh", "-merge", "1e-10", "-out", MESH,
            ],
        },
        "notes": (
            "Pilot hybrid mesh: absorber uses tetrahedra; all support-stack "
            "volumes are triangular footprints extruded as prisms. "
            "The all-tetra mesh_refined_3x remains the regression reference."
        ),
    }
    # ElmerGrid creates this directory by partitioning MESH.  It is a new
    # target so the serial hybrid mesh and all established partitions remain
    # untouched.
    project["meshes"][MESH_MPI] = {
        "geometry": "single_pixel",
        "dir": MESH_MPI,
        "recipe": {
            "partition_from": MESH,
            "elmergrid_args": ["2", "2", MESH, "-partition", "4", "1", "1", "-out", MESH_MPI],
        },
        "notes": "Four-rank x-slab partition of the conformal hybrid mesh.",
    }
    project["meshes"][MESH_MPI_NOBC] = {"geometry":"single_pixel","dir":MESH_MPI_NOBC,"recipe":{"partition_from":MESH,"elmergrid_args":["2","2",MESH,"-partition","4","1","1","-partnobcoptim","-out",MESH_MPI_NOBC]},"notes":"4-rank hybrid partition without boundary ownership optimization."}

    steady = copy.deepcopy(project["cases"]["case_tes_steady_3x_comsol_grid"])
    steady.update({
        "mesh": MESH,
        "series_file": "tes_steady_hybrid_prism_series.csv",
        "state_file": f"{MESH}/{STEADY}.state",
        "output_result": True,
        "vtu": False,
    })
    project["cases"][STEADY] = steady

    pulse = copy.deepcopy(project["cases"]["case_tes_pulse_20ms_3x_fast_compare"])
    pulse.update({
        "mesh": MESH,
        "restart_from": STEADY,
        "restart_time": 0.0,
        "state_file": f"{MESH}/{STEADY}.state",
        "series_file": "tes_pulse_hybrid_prism_fast_compare_series.csv",
        "vtu": False,
    })
    project["cases"][PULSE] = pulse

    steady_mpi = copy.deepcopy(project["cases"]["case_tes_steady_3x_mumps_inner_circuit"])
    steady_mpi.update({
        "mesh": MESH_MPI,
        "heat_source": "circuit_inner",
        "series_file": "tes_steady_hybrid_prism_mpi4_series.csv",
        "iteration_series_file": "tes_steady_hybrid_prism_mpi4_iterations.csv",
        "state_file": f"{MESH_MPI}/{STEADY_MPI}.state",
        "output_result": True,
        "vtu": False,
    })
    project["cases"][STEADY_MPI] = steady_mpi

    pulse_mpi = copy.deepcopy(project["cases"]["case_tes_mpi_legacy_regression"])
    pulse_mpi.update({
        "mesh": MESH_MPI,
        "restart_from": STEADY_MPI,
        "restart_time": 0.0,
        "heat_source": "circuit_inner",
        "series_file": "tes_hybrid_prism_mpi4_regression_series.csv",
        "iteration_series_file": "tes_hybrid_prism_mpi4_regression_iterations.csv",
        "state_file": f"{MESH_MPI}/{STEADY_MPI}.state",
        "vtu": False,
    })
    project["cases"][PULSE_MPI] = pulse_mpi

    steady_inner_1 = copy.deepcopy(steady_mpi)
    steady_inner_1.update({"mesh": MESH, "series_file": "tes_steady_hybrid_prism_inner_1rank_series.csv", "iteration_series_file": "tes_steady_hybrid_prism_inner_1rank_iterations.csv", "state_file": f"{MESH}/{STEADY_INNER_1}.state"})
    project["cases"][STEADY_INNER_1] = steady_inner_1
    pulse_inner_1 = copy.deepcopy(pulse_mpi)
    pulse_inner_1.update({"mesh": MESH, "restart_from": STEADY_INNER_1, "state_file": f"{MESH}/{STEADY_INNER_1}.state", "series_file": "tes_hybrid_prism_inner_1rank_regression_series.csv", "iteration_series_file": "tes_hybrid_prism_inner_1rank_regression_iterations.csv"})
    project["cases"][PULSE_INNER_1] = pulse_inner_1
    for name, base, mesh in ((STEADY_MPI_TIGHT, steady_mpi, MESH_MPI), (STEADY_INNER_1_TIGHT, steady_inner_1, MESH)):
        item = copy.deepcopy(base); item.update({"mesh":mesh,"state_file":f"{mesh}/{name}.state","series_file":f"tes_{name}_series.csv","iteration_series_file":f"tes_{name}_iterations.csv"}); item["solver"]["nonlinear_convergence_tolerance"]=1e-8; project["cases"][name]=item
    for name, base, restart, mesh in ((PULSE_MPI_TIGHT,pulse_mpi,STEADY_MPI_TIGHT,MESH_MPI),(PULSE_INNER_1_TIGHT,pulse_inner_1,STEADY_INNER_1_TIGHT,MESH)):
        item=copy.deepcopy(base); item.update({"mesh":mesh,"restart_from":restart,"state_file":f"{mesh}/{restart}.state","series_file":f"tes_{name}_series.csv","iteration_series_file":f"tes_{name}_iterations.csv"}); project["cases"][name]=item
    item=copy.deepcopy(project["cases"][STEADY_MPI_TIGHT]); item.update({"mesh":MESH_MPI_NOBC,"state_file":f"{MESH_MPI_NOBC}/{STEADY_MPI_NOBC_TIGHT}.state","series_file":f"tes_{STEADY_MPI_NOBC_TIGHT}_series.csv","iteration_series_file":f"tes_{STEADY_MPI_NOBC_TIGHT}_iterations.csv"}); project["cases"][STEADY_MPI_NOBC_TIGHT]=item
    for name, base, mesh, restart in ((PULSE_SERIAL_MAP_1,pulse_inner_1,MESH,STEADY_INNER_1_TIGHT),(PULSE_SERIAL_MAP_4,pulse_mpi,MESH_MPI,"case_tes_steady_hybrid_prism_serial_tight_mapped_v2")):
        item=copy.deepcopy(base); item.update({"mesh":mesh,"restart_from":None,"preexisting_restart":True,"restart_time":0.02,"state_file":f"{mesh}/{restart}.state","series_file":f"tes_{name}_series.csv","iteration_series_file":f"tes_{name}_iterations.csv","restart_file_base":restart})
        item["timesteps"]=item["timesteps"][1:]; item["output_intervals"]=item["output_intervals"][1:]; project["cases"][name]=item
    for name, base in ((PULSE_SERIAL_MAP_1_MIN3, PULSE_SERIAL_MAP_1), (PULSE_SERIAL_MAP_4_MIN3, PULSE_SERIAL_MAP_4)):
        item=copy.deepcopy(project["cases"][base]); item.update({"series_file":f"tes_{name}_series.csv","iteration_series_file":f"tes_{name}_iterations.csv"}); item["solver"]["nonlinear_min_iterations"]=3; project["cases"][name]=item
    for name, base in ((PULSE_SERIAL_MAP_1_MIN5, PULSE_SERIAL_MAP_1), (PULSE_SERIAL_MAP_4_MIN5, PULSE_SERIAL_MAP_4)):
        item=copy.deepcopy(project["cases"][base]); item.update({"series_file":f"tes_{name}_series.csv","iteration_series_file":f"tes_{name}_iterations.csv"}); item["solver"]["nonlinear_min_iterations"]=5; project["cases"][name]=item
    for name, base in ((PULSE_SERIAL_MAP_1_COMMIT,PULSE_SERIAL_MAP_1),(PULSE_SERIAL_MAP_4_COMMIT,PULSE_SERIAL_MAP_4)):
        item=copy.deepcopy(project["cases"][base]); item.update({"series_file":f"tes_{name}_series.csv","iteration_series_file":f"tes_{name}_iterations.csv","inner_circuit_step_commit":True}); project["cases"][name]=item

    args.output.write_text(json.dumps(project, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}")
    print(f"New meshes: {MESH}, {MESH_MPI}; cases: {STEADY}, {PULSE}, {STEADY_MPI}, {PULSE_MPI}")


if __name__ == "__main__":
    main()
