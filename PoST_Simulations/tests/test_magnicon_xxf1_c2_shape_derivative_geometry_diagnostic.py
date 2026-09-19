from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from subScript import magnicon_xxf1_c2_shape_derivative_geometry_diagnostic as diag


def test_cosine_similarity_parallel_and_antiparallel():
    a = np.array([1.0, 2.0, 3.0])
    assert diag.cosine_similarity(a, 2.0 * a) == pytest.approx(1.0)
    assert diag.cosine_similarity(a, -3.0 * a) == pytest.approx(-1.0)


def test_cosine_similarity_orthogonal():
    assert diag.cosine_similarity(
        np.array([1.0, 0.0]),
        np.array([0.0, 2.0]),
    ) == pytest.approx(0.0)


def test_projection_summary_exact_direction():
    target = np.array([2.0, 4.0, 6.0])
    basis = np.array([1.0, 2.0, 3.0])
    row = diag.projection_summary(target, basis)
    assert row["scale"] == pytest.approx(2.0)
    assert row["fractional_residual_norm"] == pytest.approx(0.0)
    assert row["explained_norm_fraction"] == pytest.approx(1.0)


def test_projection_summary_partial_direction():
    target = np.array([1.0, 1.0])
    basis = np.array([1.0, 0.0])
    row = diag.projection_summary(target, basis)
    assert row["scale"] == pytest.approx(1.0)
    assert row["fractional_residual_norm"] == pytest.approx(
        1.0 / np.sqrt(2.0)
    )
    assert row["explained_norm_fraction"] == pytest.approx(0.5)


def test_band_comparisons():
    frequency = np.array([1_000.0, 2_000.0, 10_000.0, 20_000.0])
    target = np.array([1.0, 2.0, 1.0, -1.0])
    basis = np.array([2.0, 4.0, -1.0, 1.0])
    rows = diag.band_comparisons(
        frequency,
        target,
        basis,
        [
            {"name": "low", "min": 1_000.0, "max": 2_000.0},
            {"name": "high", "min": 10_000.0, "max": 20_000.0},
        ],
    )
    assert rows["low"]["cosine_similarity"] == pytest.approx(1.0)
    assert rows["high"]["cosine_similarity"] == pytest.approx(-1.0)
    assert rows["high"]["absolute_cosine_similarity"] == pytest.approx(1.0)


def test_top_parameter_uses_absolute_cosine():
    rows = {
        "C_tes": {"absolute_cosine_similarity": 0.95},
        "L": {"absolute_cosine_similarity": 0.80},
    }
    assert diag._top_parameter(rows) == "C_tes"


def test_classify_day_strong():
    rows = {
        "C_tes": {"absolute_cosine_similarity": 0.95},
        "L": {"absolute_cosine_similarity": 0.60},
    }
    result = diag._classify_day(
        rows,
        {"strong_abs_cosine": 0.90, "moderate_abs_cosine": 0.70},
    )
    assert result["top_parameter"] == "C_tes"
    assert result["classification"] == "strong_single_local_shape_alignment"


def test_classify_day_moderate():
    rows = {
        "C_tes": {"absolute_cosine_similarity": 0.75},
        "L": {"absolute_cosine_similarity": 0.65},
    }
    result = diag._classify_day(
        rows,
        {"strong_abs_cosine": 0.90, "moderate_abs_cosine": 0.70},
    )
    assert result["classification"] == "moderate_single_local_shape_alignment"


def test_classify_day_weak():
    rows = {
        "C_tes": {"absolute_cosine_similarity": 0.50},
        "L": {"absolute_cosine_similarity": 0.40},
    }
    result = diag._classify_day(
        rows,
        {"strong_abs_cosine": 0.90, "moderate_abs_cosine": 0.70},
    )
    assert result["classification"] == "weak_single_local_shape_alignment"
