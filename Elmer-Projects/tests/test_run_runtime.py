from __future__ import annotations

import os
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import run


class RuntimeSifTests(unittest.TestCase):
    def test_pinned_udf_creates_copy_and_preserves_source(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "case.sif"
            source_text = 'Real Procedure "tes_transient_heat_source_t0" "TESTransientHeatSource"\n'
            source.write_text(source_text)
            udf = root / "matching.dll"
            udf.write_bytes(b"dll")
            runtime_sif, applied = run.write_runtime_sif(source, root / "results", str(udf))
            self.assertTrue(applied)
            self.assertEqual(source.read_text(), source_text)
            runtime_text = runtime_sif.read_text()
            self.assertIn(f'"{udf.resolve().as_posix()}"', runtime_text)
            self.assertNotIn('"tes_transient_heat_source_t0"', runtime_text)
            self.assertTrue(udf.resolve().suffix == ".dll")

    def test_pinned_udf_requires_expected_token(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source, udf = root / "case.sif", root / "matching.dll"
            source.write_text('Real Procedure "other" "TESTransientHeatSource"\n')
            udf.write_bytes(b"dll")
            with self.assertRaisesRegex(ValueError, "token"):
                run.write_runtime_sif(source, root / "results", str(udf))

    def test_udf_is_not_applied_to_non_udf_dependency(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source, udf = root / "steady.sif", root / "matching.dll"
            source.write_text('Procedure "HeatSolve" "HeatSolver"\n')
            udf.write_bytes(b"dll")
            runtime_sif, applied = run.write_runtime_sif(source, root / "results", str(udf))
            self.assertFalse(applied)
            self.assertEqual(runtime_sif, source)

    def test_runtime_environment_orders_paths(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            solver = root / "install" / "bin" / "ElmerSolver.exe"
            runtime = root / "runtime"
            solver.parent.mkdir(parents=True)
            runtime.mkdir()
            solver.write_bytes(b"solver")
            old_path = os.environ.get("PATH", "")
            env, found_solver, prefix = run.runtime_environment(str(solver), str(runtime))
            self.assertEqual(found_solver, solver.resolve())
            self.assertEqual(prefix, (root / "install").resolve())
            self.assertEqual(env["PATH"].split(os.pathsep)[:3], [
                str(solver.parent.resolve()),
                str((prefix / "share" / "elmersolver" / "lib").resolve()),
                str(runtime.resolve()),
            ])
            self.assertTrue(env["PATH"].endswith(old_path))

    def test_mpi_launcher_uses_pinned_path(self) -> None:
        with mock.patch.object(run.shutil, "which", return_value="C:/runtime/mpiexec.exe") as which:
            self.assertEqual(run.mpi_launcher({"PATH": "C:/runtime"}), str(Path("C:/runtime/mpiexec.exe").resolve()))
            which.assert_called_once_with("mpiexec", path="C:/runtime")

    def test_preexisting_restart_requires_state_file(self) -> None:
        model = {
            "meshes": {"m": {"dir": "mesh"}},
            "cases": {
                "pulse": {
                    "mesh": "m",
                    "preexisting_restart": True,
                    "restart_file_base": "mapped",
                }
            },
        }
        with self.assertRaisesRegex(ValueError, "state_file"):
            run.preexisting_restart_paths(model, "pulse", 4)

    def test_preexisting_mpi_restart_validates_every_rank_and_state(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            mesh = root / "mesh"
            mesh.mkdir()
            model = {
                "meshes": {"m": {"dir": "mesh"}},
                "cases": {
                    "pulse": {
                        "mesh": "m",
                        "preexisting_restart": True,
                        "restart_file_base": "mapped",
                        "state_file": "mesh/mapped.state",
                    }
                },
            }
            for rank in range(4):
                (mesh / f"mapped.result.{rank}").write_bytes(str(rank).encode())
            (mesh / "mapped.state").write_bytes(b"state")
            with mock.patch.object(run, "ROOT", root):
                paths = run.validate_preexisting_restart(model, "pulse", 4)
                self.assertEqual(
                    [path.name for path in paths],
                    [
                        "mapped.result.0",
                        "mapped.result.1",
                        "mapped.result.2",
                        "mapped.result.3",
                        "mapped.state",
                    ],
                )
                (mesh / "mapped.result.2").unlink()
                with self.assertRaisesRegex(FileNotFoundError, r"mapped\.result\.2"):
                    run.validate_preexisting_restart(model, "pulse", 4)


if __name__ == "__main__":
    unittest.main()
