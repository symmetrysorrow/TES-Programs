import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import build_physical_hdf5_v2 as build


VALID = """[Parameters]
icntl = 0
bitrseed = 0101101011001010010111000011100101101001110000111011001010010111
maxcas = 10000
maxbch = 1
negs = 1
e-mode = 2
emin(12) = 0.001
emin(13) = 0.001
emin(14) = 0.001
[Source]
totfact = 1.0
1.173 1
1.332 1
[Material]
mat[2] Pb 1
[Cell]
1000 2 -9.78 -100
[t-deposit]
mesh = reg
reg = 1000
part = all
deposit = 0
unit = 3
output = deposit
e-type = 2
ne = 150
emin = 0.0
emax = 1.5
"""


class PhysicalInputGuardTests(unittest.TestCase):
    def write(self, text):
        f = tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False)
        f.write(text)
        f.close()
        self.addCleanup(Path(f.name).unlink, missing_ok=True)
        return Path(f.name)

    def test_accepts_explicit_validated_pb_input_slice(self):
        build.validate_source_input(self.write(VALID))

    def test_rejects_input_without_explicit_positron_cutoff(self):
        with self.assertRaisesRegex(ValueError, r"emin\(13\)"):
            build.validate_source_input(self.write(VALID.replace("emin(13) = 0.001\n", "")))

    def test_rejects_wrong_region(self):
        with self.assertRaisesRegex(ValueError, "reg must be '1000'"):
            build.validate_source_input(self.write(VALID.replace("reg = 1000", "reg = 1")))

    def test_rejects_wrong_pb_region_material(self):
        with self.assertRaisesRegex(ValueError, "region 1000"):
            build.validate_source_input(self.write(VALID.replace("1000 2 -9.78", "1000 4 -1.0")))

    def test_rejects_non_pb_material(self):
        with self.assertRaisesRegex(ValueError, "pure Pb"):
            build.validate_source_input(self.write(VALID.replace("mat[2] Pb 1", "mat[2] Cu 1")))

    def test_example_json_contract_resolves_and_validates_source_input(self):
        config_path = Path(__file__).resolve().parents[1] / "validation" / "pb_physics_cli_integration" / "pb_config.example.json"
        data, source = build.load_physics_config(config_path)
        self.assertEqual(data["phits_version"], "3.341")
        self.assertTrue(source.is_file())

    def test_rejects_non_unit_weight_contract(self):
        import json
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "trace.inp"
            source.write_text(VALID, encoding="utf-8")
            config = {"phits_version": "3.341", "phits_input": "trace.inp", "region_id": 1000,
                      "region_material": "Pb", "dumpall_material_id": 2, "egs5_medium": 2,
                      "electron_cutoff_MeV": 0.001, "positron_cutoff_MeV": 0.001,
                      "weight_mode": "0.5", "energy_unit": "MeV",
                      "deposit_mode": "t-deposit-output-deposit"}
            path = root / "config.json"
            path.write_text(json.dumps(config), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "weight_mode"):
                build.load_physics_config(path)


if __name__ == "__main__":
    unittest.main()
