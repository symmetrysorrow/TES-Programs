import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "validation"))
import build_physical_hdf5_v2 as build


VALID = """[Parameters]
icntl = 0
maxcas = 10000
maxbch = 1
negs = 1
e-mode = 2
emin(12) = 0.001
emin(13) = 0.001
emin(14) = 0.001
[Source]
1.173 1
1.332 1
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


if __name__ == "__main__":
    unittest.main()
