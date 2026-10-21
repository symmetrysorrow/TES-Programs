import importlib.util
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "validation" / "pb_physics_cli_integration" / "scripts" / "audit_addition_count_closure.py"
spec = importlib.util.spec_from_file_location("audit_addition_count_closure", SCRIPT)
audit_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit_module)


class AdditionClosureTests(unittest.TestCase):
    def test_exclusive_category_sum_is_checked(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "ledger.csv"
            path.write_text("process,status,E_deposit_MeV\nexact_zero,reconstructed,0\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "addition closure failed"):
                audit_module.audit(path)


if __name__ == "__main__":
    unittest.main()
