from __future__ import annotations

import unittest
import tempfile
from pathlib import Path

from scripts.prep.run_singlepixel_prod_v2_original_timegrid import (
    steady_state_from_iterations,
    truncate_timesteps,
)


class SmokeScheduleTests(unittest.TestCase):
    def test_truncate_preserves_grouped_stages(self) -> None:
        schedule = [["1[us]", 2], ["2[us]", 4], ["3[us]", 1]]
        self.assertEqual(
            truncate_timesteps(schedule, 5),
            [["1[us]", 2], ["2[us]", 3]],
        )

    def test_truncate_rejects_invalid_or_oversize_request(self) -> None:
        with self.assertRaisesRegex(ValueError, ">= 1"):
            truncate_timesteps([["1[us]", 1]], 0)
        with self.assertRaisesRegex(ValueError, "only 1"):
            truncate_timesteps([["1[us]", 1]], 2)

    def test_known_steady_state_uses_final_converged_values(self) -> None:
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "iterations.csv"
            path.write_text(
                "time_s,time_step,nonlinear_iter,tes_temperature_K,"
                "previous_current_A,raw_current_A,tes_resistance_ohm,"
                "raw_power_W,residual_W,omega,omega_cap,relaxed_power_W\n"
                "1,1,1,0.16,1e-4,2e-4,0.01,3e-10,0,0.5,0.5,4e-10\n"
                "1,1,2,0.17,5e-4,6e-4,0.02,7e-10,0,0.5,0.5,8e-10\n",
                encoding="utf-8",
            )
            self.assertEqual(
                steady_state_from_iterations(path),
                (0.17, 6e-4, 0.02, 8e-10, 5e-4),
            )


if __name__ == "__main__":
    unittest.main()
