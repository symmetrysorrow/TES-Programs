from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results" / "case_tes_single_pixel_4rank_comsol200"
OUT.mkdir(parents=True, exist_ok=True)
src = ROOT / "generated" / "cases" / "case_tes_mpi_comsol_grid.sif"
sif = ROOT / "generated" / "cases" / "case_tes_single_pixel_4rank_comsol200.sif"
text = src.read_text(encoding="utf-8")
sizes = re.search(r"Timestep Sizes\(227\) = (.*)", text).group(1).split()
sizes += ["0.005"] * 35 + ["0.00438"]
intervals = ["1"] * len(sizes)
text = re.sub(r"Timestep Sizes\(227\) = .*", "Timestep Sizes(%d) = %s" % (len(sizes), " ".join(sizes)), text)
text = re.sub(r"Timestep Intervals\(227\) = .*", "Timestep Intervals(%d) = %s" % (len(sizes), " ".join(intervals)), text)
text = re.sub(r"Output Intervals\(227\) = .*", "Output Intervals(%d) = %s" % (len(sizes), " ".join(intervals)), text)
text = text.replace("case_tes_mpi_comsol_grid", "case_tes_single_pixel_4rank_comsol200")
text = text.replace("tes_mpi_comsol_grid_series.csv", "tes_single_pixel_4rank_comsol200_series.csv")
text = text.replace("tes_mpi_comsol_grid_iterations.csv", "tes_single_pixel_4rank_comsol200_iterations.csv")
text = text.replace("Linear System Direct Method = MUMPS", "Linear System Direct Method = Umfpack")
text = text.replace("generated/cases/case_tes_single_pixel_4rank_comsol200.sif",
                    "generated/cases/case_tes_single_pixel_4rank_comsol200.sif")
sif.write_text(text, encoding="utf-8")
log = OUT / "solver.log"
solver = Path(r"C:\Program Files\Elmer 26.1-Release\bin\ElmerSolver.exe")
mpiexec = Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe")
with log.open("w", encoding="utf-8") as f:
    p = subprocess.run([str(mpiexec), "-n", "4", str(solver), str(sif)], cwd=ROOT,
                       stdout=f, stderr=subprocess.STDOUT)
for name in ("tes_single_pixel_4rank_comsol200_series.csv",
             "tes_single_pixel_4rank_comsol200_iterations.csv"):
    path = ROOT / name
    if path.exists(): shutil.move(str(path), OUT / name)
print(f"exit_code={p.returncode}")
