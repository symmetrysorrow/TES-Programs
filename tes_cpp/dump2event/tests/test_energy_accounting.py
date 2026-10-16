"""Synthetic PHITS dumpall regression tests for reaction energy bookkeeping.

Set TES_CPP_DUMP2EVENT_EXECUTABLE to a built binary or put dump2event on PATH.
These tests check parser/accounting invariants, NOT physics validation against
real PHITS simulations or an independent [t-deposit] tally.
"""

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


def record(ncol, particle_id, kind, e_before, e_at, secondaries=()):
    # Matches the line-count state machine in read_dump(): electron records
    # omit a photon-only data line.
    lines = [str(ncol), f"{particle_id} 0 {kind}"]
    lines.extend(["0"] * (8 if kind in (12, 13) else 7))
    lines.extend([f"{e_before} 0", "0", f"{e_at} 0 0", "0 0", "0", "0"])
    if ncol in (13, 14):
        lines.append(str(len(secondaries)))
        if secondaries:
            lines.append("0 0 13")  # MATHZ, MATHN, JCOLL
            for secondary_kind, status, energy in secondaries:
                lines.extend([
                    "0", f"0 0 0 {secondary_kind} {status}",
                    "0", "0", "0", f"0 {energy}", "0", "0", "0",
                ])
    return lines


class DumpEnergyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.exe = os.environ.get("TES_CPP_DUMP2EVENT_EXECUTABLE") or shutil.which("dump2event")
        if not cls.exe:
            raise unittest.SkipTest("build dump2event and set TES_CPP_DUMP2EVENT_EXECUTABLE")

    def convert(self, records, *, normalize=False, source_energy=1.0, input_energy=1.0):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            dump = root / "dumpall.dat"
            output = root / "events.json"
            summary = root / "summary.csv"
            # PHITS opens each history with a complete NCOL=4 particle
            # record. Only then can the next NCOL header be parsed.
            initial = ["0", "4", "101", "1 0 14", *["0"] * 7,
                       f"{source_energy} 0", "0", f"{source_energy} 0 0", "0 0", "0", "0"]
            dump.write_text("\n".join([*initial, *records, "2"]) + "\n")
            command = [self.exe, str(dump), str(output), "--input-energy", str(input_energy),
                       "--save-all", "--history-summary", str(summary)]
            if normalize:
                command.append("--normalize-deposits")
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            doc = json.loads(output.read_text())
            with summary.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 1)
            return doc["101"], rows[0]

    def test_transported_secondary_is_not_deposited_twice(self):
        records = [
            *record(13, 1, 14, 1.0, 0.0, [(12, 0, 0.7)]),
            *record(11, 2, 12, 0.7, 0.7),
        ]
        event, diag = self.convert(records)
        self.assertAlmostEqual(event["1"]["E_deposit"][0], 0.3)
        self.assertAlmostEqual(event["2"]["E_deposit"][0], 0.7)
        self.assertAlmostEqual(float(diag["raw_total_deposit"]), 1.0)
        self.assertEqual(diag["energy_consistent"], "1")

    def test_zero_secondary_reaction_does_not_consume_next_ncol(self):
        event, diag = self.convert(record(13, 1, 14, 0.6, 0.0))
        self.assertAlmostEqual(event["1"]["E_deposit"][0], 0.6)
        self.assertEqual(diag["energy_consistent"], "0")

    def test_no_implicit_scaling_and_opt_in_legacy_mode(self):
        records = [
            *record(13, 1, 14, 0.6, 0.0, [(12, 0, 0.4)]),
            *record(11, 2, 12, 0.4, 0.4),
        ]
        event, diag = self.convert(records)
        self.assertAlmostEqual(sum(sum(p["E_deposit"]) for p in event.values()), 0.6)
        self.assertAlmostEqual(float(diag["raw_total_deposit"]), 0.6)
        self.assertEqual(diag["energy_consistent"], "0")
        normalized, n_diag = self.convert(records, normalize=True)
        self.assertAlmostEqual(sum(sum(p["E_deposit"]) for p in normalized.values()), 1.0)
        self.assertAlmostEqual(float(n_diag["raw_total_deposit"]), 0.6)
        self.assertEqual(n_diag["energy_consistent"], "0")

    def test_negative_parent_deposit_is_reported_not_hidden(self):
        records = [
            *record(13, 1, 14, 0.4, 0.0, [(12, 0, 0.6)]),
            *record(11, 2, 12, 0.6, 0.6),
        ]
        event, diag = self.convert(records)
        self.assertAlmostEqual(event["1"]["E_deposit"][0], -0.2)
        self.assertEqual(diag["negative_deposits"], "1")

    def test_source_line_energy_overrides_fixed_configured_energy_for_audit(self):
        # A 60Co-like spectrum with one sampled gamma per source history.
        # The CLI energy is 1.332 MeV, but the NCOL=4 photon is 1.173 MeV.
        event, diag = self.convert(
            record(13, 1, 14, 1.173, 0.0),
            source_energy=1.173,
            input_energy=1.332,
        )
        self.assertAlmostEqual(sum(sum(p["E_deposit"]) for p in event.values()), 1.173)
        self.assertAlmostEqual(float(diag["source_energy"]), 1.173)
        self.assertEqual(diag["source_energy_available"], "1")
        self.assertAlmostEqual(float(diag["energy_balance_reference"]), 1.173)
        self.assertEqual(diag["energy_consistent"], "1")
        # This would have failed under the old global input-energy comparison.
        self.assertAlmostEqual(1.332 - float(diag["raw_total_deposit"]), 0.159)

    def test_parent_type_preserved_for_non_em_secondary(self):
        event, _ = self.convert(record(13, 1, 14, 0.8, 0.0, [(1, 0, 0.3)]))
        self.assertAlmostEqual(event["1"]["E_deposit"][0], 0.5)


if __name__ == "__main__":
    unittest.main()
