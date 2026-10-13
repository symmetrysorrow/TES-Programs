import pytest

from support.build_cases import bdf2_predictor_lines


@pytest.mark.parametrize(
    "key", ("bdf2_predictor", "phase24_bdf2_predictor"),
)
def test_predictor_produces_the_same_stable_sif_keyword(key):
    assert bdf2_predictor_lines({"bdf_order": 2, key: True}) == [
        "  Phase24 BDF2 Predictor = Logical True"
    ]


def test_predictor_disabled_and_consistent_aliases():
    assert bdf2_predictor_lines({"bdf_order": 2}) == []
    assert bdf2_predictor_lines({"bdf_order": 2, "bdf2_predictor": False}) == []
    assert bdf2_predictor_lines({
        "bdf_order": 2, "bdf2_predictor": True, "phase24_bdf2_predictor": True,
    }) == ["  Phase24 BDF2 Predictor = Logical True"]


def test_predictor_rejects_conflicting_aliases_and_bdf1():
    with pytest.raises(ValueError, match="conflicting"):
        bdf2_predictor_lines({
            "bdf_order": 2, "bdf2_predictor": False, "phase24_bdf2_predictor": True,
        })
    for key in ("bdf2_predictor", "phase24_bdf2_predictor"):
        with pytest.raises(ValueError, match="bdf_order"):
            bdf2_predictor_lines({"bdf_order": 1, key: True})
