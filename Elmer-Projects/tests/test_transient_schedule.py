import pytest

from support.run_phase24_gate4_5_nomortar import trim_schedule
from support.transient_schedule import project_schedule


@pytest.mark.parametrize(
    "window",
    [
        "75ms_bdf2h15",
        "1ms_bdf2h15",
        "1ms_bdf2fixed2p5",
        "1ms_bdf2fixed1p25",
        "1ms_bdf2fixed2p5_5",
        "singlepixel_hybrid_1p48ms",
    ],
)
def test_project_windows_match_archived_validated_schedule(window):
    assert project_schedule(window) == trim_schedule([], window)
