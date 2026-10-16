"""Regression coverage for two-channel event-key alignment."""

import numpy as np
import pandas as pd
import pytest

from Analyze_Experimental_Data.tes_analysis import analysis_utils as general


def test_two_channel_values_are_joined_by_key_after_row_shuffle():
    ch0 = pd.DataFrame({"key": [30, 10, 20], "Peak": [3.0, 1.0, 2.0]})
    ch1 = pd.DataFrame({"key": [20, 30, 10], "Peak": [200.0, 300.0, 100.0]})

    keys, peak0, peak1 = general._align_keyed_columns(ch0, ch1, "Peak", "Peak")

    np.testing.assert_array_equal(keys, [10, 20, 30])
    np.testing.assert_array_equal(peak0, [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(peak1, [100.0, 200.0, 300.0])


def test_keyisin_returns_rows_in_the_same_key_order_after_shuffle():
    ch0 = pd.DataFrame({"key": [30, 10, 20], "Peak": [3.0, 1.0, 2.0]})
    ch1 = pd.DataFrame({"key": [20, 30, 10], "Peak": [200.0, 300.0, 100.0]})

    left, right = general.KeyIsin(ch0, ch1)

    np.testing.assert_array_equal(left["key"], right["key"])
    np.testing.assert_array_equal(left["key"], [10, 20, 30])
    np.testing.assert_array_equal(left["Peak"], [1.0, 2.0, 3.0])
    np.testing.assert_array_equal(right["Peak"], [100.0, 200.0, 300.0])


@pytest.mark.parametrize("bad_frame", [
    pd.DataFrame({"key": [1, 1], "Peak": [0.1, 0.2]}),
    pd.DataFrame({"key": [1, np.nan], "Peak": [0.1, 0.2]}),
])
def test_keyed_pairing_rejects_duplicate_or_missing_event_keys(bad_frame):
    other = pd.DataFrame({"key": [1, 2], "Peak": [1.0, 2.0]})

    with pytest.raises(ValueError):
        general._align_keyed_columns(bad_frame, other, "Peak", "Peak")

