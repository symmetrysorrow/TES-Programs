from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

sp = pytest.importorskip("scipy.sparse")

from scripts.support.analyze_phase24_native_branch_operator_audit import (  # noqa: E402
    LocalOperator, annulus_ntd_shape_factor,
)
from scripts.support.phase24_element_dissipation import _tet, _wedge  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "artifacts" / "phase24_native_branch_operator_audit"
GHIST = 1.8782944616314638e-8
GBEST = 1.959713948e-8


def load(name: str):
    return json.loads((AUDIT / name).read_text(encoding="utf-8"))


def rows(name: str) -> list[dict[str, str]]:
    with (AUDIT / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


# ---------------------------------------------------------------- unit level


@pytest.mark.parametrize("builder,points", [
    (_tet, np.array([[0, 0, 0], [2, 0, 0], [0, 1, 0], [0.3, 0.2, 0.5]], dtype=float)),
    (_wedge, np.array([[0, 0, 0], [2, 0, 0], [0, 1, 0], [0, 0, 0.5], [2, 0, 0.5], [0, 1, 0.5]], dtype=float)),
])
def test_element_conductance_is_exact_for_linear_fields(builder, points) -> None:
    k = builder(points)
    assert np.allclose(k.sum(axis=1), 0.0, atol=1e-12)
    grad = np.array([0.7, -1.3, 2.0])
    t = points @ grad
    volume = abs(np.linalg.det(np.c_[np.ones(4), points[:4]])) / 6.0 if len(points) == 4 else 0.5 * 2 * 1 * 0.5
    assert t @ k @ t == pytest.approx(volume * grad @ grad, rel=1e-12)


def test_local_operator_condenses_a_series_chain_exactly() -> None:
    # 5-node chain with unit link conductances; trace node 0, ground node 4
    n = 5
    rows_, cols, vals = [], [], []
    for i in range(n - 1):
        for a, b, v in ((i, i, 1.0), (i + 1, i + 1, 1.0), (i, i + 1, -1.0), (i + 1, i, -1.0)):
            rows_.append(a); cols.append(b); vals.append(v)
    K = sp.csr_matrix((vals, (rows_, cols)), shape=(n, n))
    op = LocalOperator(K, np.array([0]), np.array([4]), np.arange(n))
    assert op.apply(np.array([1.0])).item() == pytest.approx(0.25)


def test_local_operator_drops_uninformative_hanging_constraint() -> None:
    n = 3
    K = sp.csr_matrix(np.array([[1.0, -1.0, 0.0], [-1.0, 2.0, -1.0], [0.0, -1.0, 1.0]]))
    # node 0 trace, node 2 ground; a constraint touching only the trace
    C = sp.csr_matrix(np.array([[1.0, 0.0, 0.0]]))
    op = LocalOperator(K, np.array([0]), np.array([2]), np.arange(n), C)
    assert op.dropped_constraint_rows == 1
    assert op.apply(np.array([1.0])).item() == pytest.approx(0.5)


def test_annulus_reference_converges() -> None:
    assert annulus_ntd_shape_factor(8.0) == pytest.approx(11.770, abs=0.01)


# ---------------------------------------------------------------- artifacts


def test_native_capture_reproduces_reported_conductance() -> None:
    acc = load("remaining_gap_accounting.json")["best_phase24_minus_historical"]
    assert acc["G_reference"] == pytest.approx(GHIST, rel=1e-9)
    assert acc["G_case"] == pytest.approx(GBEST, rel=5e-5)
    assert acc["ratio"] == pytest.approx(1.043348, abs=5e-5)


def test_element_reassembly_matches_native_matrix() -> None:
    ver = load("provenance.json")["element_reassembly_verification"]
    assert ver["historical"]["max_entry_diff_relative_to_row_diagonal"] < 1e-10
    for key in ("best_phase24", "refined_parent"):
        assert ver[key]["max_entry_diff_relative_to_row_diagonal"] < 2e-6


def test_native_branch_flux_closes() -> None:
    by = {}
    for row in rows("native_branch_flux.csv"):
        by.setdefault((row["case"], row["P_fraction"]), {})[row["branch_id"]] = row
    for (case, _), branches in by.items():
        P = float(branches["B0"]["P_W"])
        assert float(branches["B0"]["Q_over_P"]) == pytest.approx(1.0, abs=5e-5)
        assert abs(float(branches["B0s"]["Q_W"])) < 1e-6 * P
        assert float(branches["M"]["Q_over_P"]) == pytest.approx(1.0, abs=5e-5)
        exits = sum(float(branches[b]["Q_W"]) for b in ("A", "B", "C"))
        assert exits == pytest.approx(float(branches["B0"]["Q_W"]), rel=1e-6)
        for b in ("A", "B", "C"):
            assert abs(float(branches[b]["closure_W"])) < 1e-6 * P


def test_gap_is_closed_by_direct_evidence() -> None:
    acc = load("remaining_gap_accounting.json")["best_phase24_minus_historical"]
    assert acc["explained_fraction"] > 0.99
    cats = acc["categories_fraction_of_gap_eff"]
    assert abs(cats["merge_frame_to_bath_network"]) < 0.01
    assert abs(cats["TES_to_Membrane_operator (TES body + coupling)"]) < 0.05
    assert cats["branch_C_SiO2_1_window_sheet"] == max(v for k, v in cats.items())


def test_operator_split_localizes_gap_to_tes_edge_loading() -> None:
    dec = load("operator_gap_decomposition.json")["best_phase24_over_historical"]
    assert dec["bulk_stack_operator_ratio(uniform flux)"] == pytest.approx(1.0082, abs=5e-4)
    assert dec["log_share_edge"] > 0.75
    assert dec["G_ratio"] == pytest.approx(dec["bulk_stack_operator_ratio(uniform flux)"] * dec["TES_edge_loading_ratio"], rel=1e-12)


def test_historical_mortar_slave_trace_is_flagged_ill_posed() -> None:
    tes = load("tes_membrane_condensed_operator.json")
    assert tes["historical"]["mortar_slave_Dirichlet_wellposedness"]["conclusion"].endswith("ill-posed")
