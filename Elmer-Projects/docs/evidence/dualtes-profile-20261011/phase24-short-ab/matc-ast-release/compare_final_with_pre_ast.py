from pathlib import Path
import json, math

HERE = Path(__file__).resolve().parent
RUNTIME = Path(r"C:\Users\Nao\AppData\Local\Temp\phase24-cache-runtime")

def read_result(path):
    lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    result = {}
    i = 0
    while i < len(lines):
        if not lines[i].startswith("Time:"):
            i += 1
            continue
        step = int(lines[i].split()[1]); i += 1
        variable = lines[i].strip(); i += 1
        count = int(lines[i].split()[2]); i += 1
        ids = [int(lines[i+j].split()[1]) for j in range(count)]; i += count
        vals = [float(lines[i+j].strip()) for j in range(count)]; i += count
        result[step] = (variable, dict(zip(ids, vals)))
    return result

old = json.loads((HERE / "runtime_baseline_pre_ast.json").read_text())
new = json.loads((HERE / "runtime_runs.json").read_text())
rows = []
for case, mesh, n_ranks in (("dualtes", "mesh_toml_95333801e1a3", 6),
                           ("singlepixel", "mesh_hybrid_fullconf_h8", 6)):
    for mode in ("off", "on"):
        prior = next(r for r in old["records"] if r["case"] == case and r["mode"] == mode)
        current = next(r for r in new["records"] if r["case"] == case and r["mode"] == mode)
        per = []
        for rank in range(n_ranks):
            base = RUNTIME / mesh
            a = read_result(base / f"recon_{case}_{mode}.result.{rank}")
            b = read_result(base / f"final_{case}_{mode}.result.{rank}")
            if a.keys() != b.keys(): raise SystemExit(f"timestep mismatch: {case} {mode} rank {rank}")
            for step in a:
                av, bv = a[step][1], b[step][1]
                if av.keys() != bv.keys(): raise SystemExit(f"node mismatch: {case} {mode} rank {rank} step {step}")
                diffs = [abs(av[k] - bv[k]) for k in av]
                per.append({"rank": rank, "timestep": step, "nodes": len(av),
                            "max_abs_K": max(diffs, default=0.0),
                            "max_rel": max((d/max(abs(av[k]), abs(bv[k]), 1e-300)
                                            for k, d in zip(av, diffs)), default=0.0),
                            "changed_values": sum(d != 0.0 for d in diffs)})
        rows.append({"case": case, "mode": mode,
                     "pre_ast_wall_seconds": prior["wall_seconds"],
                     "final_wall_seconds": current["wall_seconds"],
                     "pre_ast_solver_sha256": prior["solver_sha256"],
                     "final_solver_sha256": current["solver_sha256"],
                     "result_norms_identical": prior["result_norms"] == current["result_norms"],
                     "hypre_records_identical": prior["hypre_pcg"] == current["hypre_pcg"],
                     "tes_csv_hashes_identical": [x["sha256"] for x in prior["outputs"] if x["file"].endswith(".csv")] ==
                                                 [x["sha256"] for x in current["outputs"] if x["file"].endswith(".csv")],
                     "all_nodes_all_four_steps": per,
                     "max_abs_K": max(p["max_abs_K"] for p in per),
                     "max_rel": max(p["max_rel"] for p in per),
                     "changed_values": sum(p["changed_values"] for p in per)})
out = {"baseline": "saved four-step runtime on reconstructed 01-only Release before AST API",
       "candidate": "same Release configuration after 02 MATC AST purity patch",
       "note": "solver/HeatSolve SHA are compared; the candidate MATC DLL hash is separately recorded in manifest",
       "comparisons": rows}
(HERE / "pre_ast_runtime_comparison.json").write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
print(json.dumps([{k:r[k] for k in ("case","mode","pre_ast_wall_seconds","final_wall_seconds","result_norms_identical","hypre_records_identical","tes_csv_hashes_identical","max_abs_K","max_rel","changed_values")} for r in rows], indent=2))
