"""Manual smoke test for a locally built posi2pulse executable."""

from pathlib import Path

from tes_cpp import posi2pulse


root = Path(__file__).parent
pulses = posi2pulse(root / "input_smoke.json", [1, 2])
assert [pulse.position for pulse in pulses] == [1, 2]
assert all(len(pulse.time) == len(pulse.ch0) == len(pulse.ch1) == 8 for pulse in pulses)
