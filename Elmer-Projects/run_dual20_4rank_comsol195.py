from pathlib import Path
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "case_dual20_pos30_4rank_comsol195"
OUT.mkdir(parents=True, exist_ok=True)
mesh_local = OUT / "partitioning.4"
if not mesh_local.exists():
    shutil.copytree(ROOT / "mesh_dual_20mm_repart_x" / "partitioning.4", mesh_local)
partitioned_root = ROOT / "mesh_dual_20mm_repart_x"
for mesh_name in ("mesh.header", "mesh.nodes", "mesh.elements", "mesh.boundary", "mesh.names"):
    target = partitioned_root / mesh_name
    if not target.exists():
        shutil.copy2(ROOT / "mesh_dual_20mm" / mesh_name, target)

src = ROOT / "generated" / "cases" / "case_dual20_pos30.sif"
sif = OUT / "case_dual20_pos30_4rank_comsol195.sif"
text = src.read_text(encoding="utf-8")
text = text.replace('"mesh_dual_20mm"', '"mesh_dual_20mm_repart_x"', 1)
text = text.replace("mesh_dual_20mm/", "../mesh_dual_20mm/")
text = text.replace("Restart File = case_dual20_steady.result",
                    "Restart File = ../mesh_dual_20mm/case_dual20_steady.result")
text = text.replace("Timestep Intervals(8) = 20 1 2 1 10 9 2000 700",
                    "Timestep Intervals(9) = 20 1 2 1 10 9 2000 1649 1")
text = text.replace("Timestep Sizes(8) =", "Timestep Sizes(9) =")
text = re.sub(r"(Timestep Sizes\(9\)\s*=.*?\s)0\.0001\s*$",
              r"\g<1>0.0001 7.8999999e-05", text, flags=re.MULTILINE)
text = text.replace("Output Intervals(8) = 10 1 1 1 5 3 10 10",
                    "Output Intervals(9) = 10 1 1 1 5 3 10 10 1")
text = text.replace("Linear System Solver = Direct\n  Linear System Direct Method = Umfpack",
                    "Linear System Solver = Direct\n  Linear System Direct Method = Umfpack")
text = text.replace("generated/cases/case_dual20_pos30.sif",
                    str(sif).replace("\\", "/"))
text = text.replace("tes_dual20_pos30_L_series.csv", "tes_dual20_pos30_4rank_comsol195_L_series.csv")
text = text.replace("tes_dual20_pos30_R_series.csv", "tes_dual20_pos30_4rank_comsol195_R_series.csv")
text = text.replace('"case_dual20_pos30"', '"case_dual20_pos30_4rank_comsol195"')
sif.write_text(text, encoding="utf-8")

solver = Path(r"C:\Program Files\Elmer 26.1-Release\bin\ElmerSolver.exe")
mpiexec = Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe")
log = OUT / "solver.log"
with log.open("w", encoding="utf-8") as handle:
    proc = subprocess.run([str(mpiexec), "-n", "4", str(solver), str(sif)],
                          cwd=ROOT, stdout=handle, stderr=subprocess.STDOUT)

for name in ("tes_dual20_pos30_4rank_comsol195_L_series.csv",
             "tes_dual20_pos30_4rank_comsol195_R_series.csv"):
    p = ROOT / name
    if p.exists():
        shutil.move(str(p), OUT / name)

print(f"exit_code={proc.returncode}")
print(f"output={OUT}")
