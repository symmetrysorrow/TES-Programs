"""Synthetic tests: these are validator tests, not PHITS physics validation."""

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / 'validation' / 'compare_physical_deposits.py'
FIELDS = ('addition_id', 'history_id', 'region', 'x', 'y', 'z', 'energy')


def row(key='1', energy='0.1', history='1', region='1000', x='0'):
    return dict(zip(FIELDS, (key, history, region, x, '0', '0', energy)))


class TallyLedgerValidatorTests(unittest.TestCase):
    def run_case(self, ref, cand, *, add_args=()):
        with tempfile.TemporaryDirectory() as directory:
            ref_path, cand_path = Path(directory) / 'reference.csv', Path(directory) / 'candidate.csv'
            for path, rows in ((ref_path, ref), (cand_path, cand)):
                with path.open('w', encoding='utf-8', newline='') as f:
                    writer = csv.DictWriter(f, fieldnames=FIELDS)
                    writer.writeheader()
                    writer.writerows(rows)
            cmd = [sys.executable, str(SCRIPT), str(ref_path), str(cand_path),
                   '--key-column', 'addition_id', '--reference-energy-column', 'energy',
                   '--candidate-energy-column', 'energy', '--history-column', 'history_id',
                   '--region-column', 'region', '--position-columns', 'x', 'y', 'z',
                   *add_args]
            return subprocess.run(cmd, capture_output=True, text=True)

    def test_exact_match_and_zeros(self):
        records = [row(), row('2', '0', history='2')]
        got = self.run_case(records, records)
        self.assertEqual(got.returncode, 0, got.stderr)
        self.assertIn('"reference_zero_count": 1', got.stdout)

    def test_near_machine_precision(self):
        got = self.run_case([row(energy='0.1')], [row(energy='0.10000000000000006')])
        self.assertEqual(got.returncode, 0, got.stderr)

    def test_energy_mismatch_cannot_be_masked_by_total(self):
        got = self.run_case([row('1', '0.2'), row('2', '0.3')],
                            [row('1', '0.1'), row('2', '0.4')])
        self.assertEqual(got.returncode, 1)
        self.assertIn('"mismatched_count": 2', got.stdout)

    def test_missing_contribution_fails(self):
        got = self.run_case([row(), row('2')], [row()])
        self.assertEqual(got.returncode, 1)
        self.assertIn('"missing_count": 1', got.stdout)

    def test_extra_contribution_fails(self):
        got = self.run_case([row()], [row(), row('2')])
        self.assertEqual(got.returncode, 1)
        self.assertIn('"unexpected_count": 1', got.stdout)

    def test_duplicate_key_is_invalid(self):
        got = self.run_case([row(), row()], [row()])
        self.assertEqual(got.returncode, 2)
        self.assertIn('duplicate', got.stderr)

    def test_negative_is_invalid(self):
        got = self.run_case([row()], [row(energy='-0.01')])
        self.assertEqual(got.returncode, 2)
        self.assertIn('negative', got.stderr)

    def test_nan_is_invalid(self):
        got = self.run_case([row()], [row(energy='NaN')])
        self.assertEqual(got.returncode, 2)
        self.assertIn('non-finite', got.stderr)

    def test_region_mismatch_fails(self):
        got = self.run_case([row()], [row(region='2000')])
        self.assertEqual(got.returncode, 1)

    def test_history_mismatch_fails(self):
        got = self.run_case([row()], [row(history='2')])
        self.assertEqual(got.returncode, 1)

    def test_position_mismatch_fails(self):
        got = self.run_case([row()], [row(x='1.0')])
        self.assertEqual(got.returncode, 1)

    def test_empty_ledger_is_invalid(self):
        got = self.run_case([], [])
        self.assertEqual(got.returncode, 2)

    def test_zero_tolerance(self):
        got = self.run_case([row()], [row(energy='0.10000000000000006')],
                            add_args=['--energy-abs-tol', '0'])
        self.assertEqual(got.returncode, 1)


if __name__ == '__main__':
    unittest.main()
