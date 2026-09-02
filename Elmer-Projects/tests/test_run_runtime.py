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

    def test_amgx_override_replaces_direct_solver_in_runtime_copy(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "case.sif"
            source.write_text(
                "Solver 1\n"
                "  Linear System Solver = Direct\n"
                "  Linear System Direct Method = MUMPS\n"
                "End\n"
            )
            config = root / "amgx.json"
            config.write_text("{}")
            runtime_sif, applied = run.write_runtime_sif(
                source, root / "results", None, str(config)
            )
            self.assertTrue(applied)
            runtime_text = runtime_sif.read_text()
            self.assertIn('Linear System Solver = "AMGX"', runtime_text)
            self.assertIn('Linear System Iterative Method = "FGMRES"', runtime_text)
            self.assertIn('Linear System Preconditioning = "AMG"', runtime_text)
            self.assertIn(
                'Linear System Scaling Method = "row equilibration"', runtime_text
            )
            self.assertIn(f'AMGX Config = String "{config.resolve().as_posix()}"', runtime_text)
            self.assertNotIn("Linear System Direct Method", runtime_text)
            self.assertIn("Linear System Solver = Direct", source.read_text())

    def test_amgx_eliminates_mortar_multiplier_rows(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "case.sif"
            source.write_text(
                "Solver 1\n"
                "  Linear System Solver = Direct\n"
                "  Apply Mortar BCs = True\n"
                "End\n"
            )
            config = root / "amgx.json"
            config.write_text("{}")
            runtime_sif, _ = run.write_runtime_sif(
                source, root / "results", None, str(config)
            )
            self.assertIn(
                "Eliminate Linear Constraints = Logical True",
                runtime_sif.read_text(),
            )

    def test_amgx_slave_constraint_mode(self) -> None:
        text = (
            "Solver 1\n"
            "  Linear System Solver = Direct\n"
            "  Apply Mortar BCs = True\n"
            "End\n"
        )
        configured = run.configure_amgx_sif(
            text, Path("amgx.json"), constraint_mode="slave"
        )
        self.assertIn("Eliminate Slave = Logical True", configured)

    def test_amgx_no_scaling_mode_omits_row_equilibration(self) -> None:
        text = (
            "Solver 1\n"
            "  Linear System Solver = Direct\n"
            "  Apply Mortar BCs = True\n"
            "End\n"
        )
        configured = run.configure_amgx_sif(
            text, Path("amgx.json"), constraint_mode="no-scaling"
        )
        self.assertNotIn("Linear System Scaling Method", configured)
        self.assertIn("Linear System Scaling = Logical False", configured)

    def test_amgx_dual_lagrange_constraint_mode(self) -> None:
        text = (
            "Solver 1\n"
            "  Linear System Solver = Direct\n"
            "  Apply Mortar BCs = True\n"
            "End\n"
            "Boundary Condition 1\n"
            "  Mortar BC = 2\n"
            "End\n"
            "Boundary Condition 2\n"
            "  Name = not_mortar\n"
            "End\n"
        )
        configured = run.configure_amgx_sif(
            text, Path("amgx.json"), constraint_mode="dual-lagrange"
        )
        self.assertIn("Biorthogonal Dual Lagrange Coefficients = Logical True", configured)
        self.assertIn("Biorthogonal Dual Slave = Logical False", configured)
        self.assertIn("Biorthogonal Dual Master = Logical False", configured)
        self.assertEqual(configured.count("Use Biorthogonal Basis"), 1)

    def test_amgx_penalty_constraint_mode(self) -> None:
        text = (
            "Solver 1\n"
            "  Linear System Solver = Direct\n"
            "  Apply Mortar BCs = True\n"
            "End\n"
        )
        configured = run.configure_amgx_sif(
            text, Path("amgx.json"), constraint_mode="penalty",
            constraint_penalty=2500.0,
        )
        self.assertNotIn("Eliminate Linear Constraints = Logical True", configured)
        self.assertIn("Penalty Linear Constraints = Logical True", configured)
        self.assertIn("Linear Constraint Penalty = Real 2500", configured)

    def test_amgx_schur_constraint_mode(self) -> None:
        text = (
            "Solver 1\n"
            "  Linear System Solver = Direct\n"
            "  Apply Mortar BCs = True\n"
            "End\n"
        )
        configured = run.configure_amgx_sif(
            text, Path("amgx.json"), constraint_mode="schur"
        )
        self.assertIn('Linear System Solver = "AMGX Schur"', configured)
        self.assertIn('Linear System Iterative Method = "GMRES"', configured)
        self.assertNotIn("Eliminate Linear Constraints", configured)
        self.assertNotIn("Penalty Linear Constraints", configured)
        self.assertIn("AMGX Schur Augmentation = Real 10000", configured)
        self.assertIn("Linear System Min Iterations = Integer 1", configured)
        self.assertIn("Linear System Convergence Tolerance = Real 1.0e-6", configured)
        self.assertIn("Linear System GMRES Restart = Integer 100", configured)
        self.assertIn("Linear System Scaling = Logical True", configured)
        self.assertIn("AMGX Allow Not Converged = Logical True", configured)

    def test_amgx_stabilized_constraint_mode(self) -> None:
        text = (
            "Solver 1\n"
            "  Linear System Solver = Direct\n"
            "  Apply Mortar BCs = True\n"
            "End\n"
        )
        configured = run.configure_amgx_sif(
            text, Path("amgx.json"), constraint_mode="stabilized",
            constraint_penalty=4.0,
        )
        self.assertIn('Linear System Solver = "AMGX Stabilized"', configured)
        self.assertIn('Linear System Iterative Method = "GCR"', configured)
        self.assertNotIn("Eliminate Linear Constraints", configured)
        self.assertIn("AMGX Constraint Stabilization = Real 0.25", configured)
        self.assertIn("Linear System Scaling = Logical True", configured)

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
            mesh = root / "work" / "meshes" / "mesh"
            mesh.mkdir(parents=True)
            model = {
                "meshes": {"m": {"dir": "mesh"}},
                "cases": {
                    "pulse": {
                        "mesh": "m",
                        "preexisting_restart": True,
                        "restart_file_base": "mapped",
                        "state_file": "work/meshes/mesh/mapped.state",
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
