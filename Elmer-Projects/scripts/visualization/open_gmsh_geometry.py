from pathlib import Path

import gmsh


geometry = Path(__file__).resolve().parents[2] / "gmsh" / "project_shifted.brep"

gmsh.initialize()
try:
    gmsh.open(str(geometry))
    gmsh.fltk.run()
finally:
    gmsh.finalize()
