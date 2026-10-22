import importlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
module = importlib.import_module("tes_cpp.dump2event")


class PbPhysicsApiTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        root = Path(tmp.name)
        self.dump = root / "dump.dat"
        self.dump.write_text("fixture")
        self.config = root / "config.json"
        self.config.write_text("{}")
        self.table = root / "pb.csv"
        self.table.write_text("fixture")
        self.output = root / "event.h5"

    def test_default_api_is_pb_egs5_hdf5_v2(self):
        with patch.object(module, "_executable", return_value="dump2event"), \
             patch.object(module.subprocess, "run") as run:
            result = module.dump2event(self.dump, self.output)
        command = run.call_args.args[0]
        self.assertEqual(command[:3], ["dump2event", str(self.dump), str(self.output)])
        self.assertNotIn("--physics-config", command)
        self.assertNotIn("--stopping-table", command)
        self.assertNotIn("--input-energy", command)
        self.assertNotIn("--normalize-deposits", command)
        self.assertEqual(result.output_path, self.output)
        self.assertTrue(run.call_args.kwargs["check"])

    def test_bundled_table_can_be_selected_without_explicit_argument(self):
        with patch.object(module, "_executable", return_value="dump2event"), \
             patch.object(module.subprocess, "run") as run:
            module.dump2event(self.dump, self.output, physics_config=self.config)
        command = run.call_args.args[0]
        self.assertIn("--physics-config", command)
        self.assertNotIn("--stopping-table", command)
        self.assertTrue(run.call_args.kwargs["check"])

    def test_pb_physics_explicit_alias_is_accepted(self):
        with patch.object(module, "_executable", return_value="dump2event"), \
             patch.object(module.subprocess, "run") as run:
            module.dump2event(self.dump, self.output, physics="pb-egs5",
                              physics_config=self.config, stopping_table=self.table)
        self.assertTrue(run.called)

    def test_two_argument_api_is_supported(self):
        with patch.object(module, "_executable", return_value="dump2event"), \
             patch.object(module.subprocess, "run") as run:
            module.dump2event(self.dump, self.output)
        self.assertEqual(run.call_args.args[0][1:3], [str(self.dump), str(self.output)])

    def test_rejects_legacy_input_energy_and_filtering(self):
        for obsolete in ({"input_energy": 2.0}, {"save_all": True},
                         {"normalize_deposits": True}, {"full_energy_only": True},
                         {"reaction_audit": "audit.csv"}):
            with self.subTest(obsolete=obsolete):
                with self.assertRaisesRegex(ValueError, "old dump2event options"):
                    module.dump2event(self.dump, self.output,
                                      stopping_table=self.table,
                                      **obsolete)

    def test_rejects_unsupported_physics(self):
        with self.assertRaisesRegex(ValueError, "only PHITS"):
            module.dump2event(self.dump, self.output, physics="legacy",
                              physics_config=self.config, stopping_table=self.table)

    def test_rejects_missing_stopping_table(self):
        with self.assertRaises(FileNotFoundError):
            module.dump2event(self.dump, self.output,
                              physics_config=self.config,
                              stopping_table=self.table.with_name("missing.csv"))

    def test_rejects_legacy_output_extensions(self):
        with self.assertRaisesRegex(ValueError, "HDF5 v2"):
            module.dump2event(self.dump, self.output.with_suffix(".json"),
                              physics_config=self.config, stopping_table=self.table)

    def test_existing_output_is_never_replaced(self):
        self.output.write_text("must stay unchanged")
        with self.assertRaises(FileExistsError):
            module.dump2event(self.dump, self.output,
                              physics_config=self.config, stopping_table=self.table)
        self.assertEqual(self.output.read_text(), "must stay unchanged")


if __name__ == "__main__":
    unittest.main()
