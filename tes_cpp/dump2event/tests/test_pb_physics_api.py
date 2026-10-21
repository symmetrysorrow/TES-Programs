import sys
import importlib
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
module = importlib.import_module("tes_cpp.dump2event")


class PbPhysicsApiTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        self.dump = root / "dump.dat"
        self.dump.write_text("fixture")
        self.config = root / "config.json"
        self.config.write_text("{}")
        self.table = root / "pb.csv"
        self.table.write_text("fixture")
        self.output = root / "event.h5"

    def test_default_api_keeps_input_energy_path(self):
        with patch.object(module, "_executable", return_value="dump2event"), \
             patch.object(module.subprocess, "run") as run:
            module.dump2event(self.dump, self.output, input_energy=1.17)
        args = run.call_args.args[0]
        self.assertIn("--input-energy", args)
        self.assertNotIn("--physics", args)

    def test_pb_mode_is_explicit_hdf5_v2_opt_in(self):
        with patch.object(module, "_executable", return_value="dump2event"), \
             patch.object(module.subprocess, "run") as run:
            module.dump2event(self.dump, self.output, physics="pb-egs5",
                              physics_config=self.config, stopping_table=self.table)
        args = run.call_args.args[0]
        self.assertIn("--physics", args)
        self.assertIn("pb-egs5", args)
        self.assertIn("--hdf5-version", args)
        self.assertIn("2", args)

    def test_pb_mode_requires_explicit_config_and_table(self):
        with self.assertRaisesRegex(ValueError, "requires physics_config"):
            module.dump2event(self.dump, self.output, physics="pb-egs5")

    def test_pb_mode_rejects_legacy_normalization(self):
        with self.assertRaisesRegex(ValueError, "not supported"):
            module.dump2event(self.dump, self.output, physics="pb-egs5",
                              physics_config=self.config, stopping_table=self.table,
                              normalize_deposits=True)

    def test_pb_mode_rejects_missing_file_before_launch(self):
        with self.assertRaises(FileNotFoundError):
            module.dump2event(self.dump, self.output, physics="pb-egs5",
                              physics_config=self.config, stopping_table=self.table.with_name("missing.csv"))


if __name__ == "__main__":
    unittest.main()
