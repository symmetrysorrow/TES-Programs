import json
from pathlib import Path
import sys

import h5py
import numpy as np

SIMULATION_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SIMULATION_ROOT))
import PoST_Simulation as simulation  # noqa: E402


def _write_input(tmp_path, **updates):
    source = Path(simulation.__file__).with_name("input.json")
    para = json.loads(source.read_text(encoding="utf-8"))
    para.update({"rate": 100000.0, "samples": 256.0, "position": [1]})
    para.update(updates)
    (tmp_path / "input.json").write_text(json.dumps(para), encoding="utf-8")


def test_hanging_noise_has_seven_states_and_two_itfn_sources(tmp_path, monkeypatch):
    _write_input(
        tmp_path,
        tes_internal_model="hanging",
        C_tes_hanging=2.0e-13,
        **{"G_tes-hanging": 2.0e-7},
    )
    monkeypatch.setattr(simulation, "output", str(tmp_path))
    simulation.MakeNoise()

    with h5py.File(tmp_path / "noise.h5", "r") as h5:
        assert h5.attrs["noise_model"] == "seven_state_hanging_tes"
        assert h5.attrs["state_count"] == 7
        assert h5["transfer_ch0"].shape[1] == 10
        assert h5["transfer_ch1"].shape == h5["transfer_ch0"].shape
        assert json.loads(h5.attrs["source_names"])[-2:] == [
            "phonon_tes1_hanging",
            "phonon_tes2_hanging",
        ]
        for name in ("phonon_tes1_hanging", "phonon_tes2_hanging"):
            assert name in h5["components_ch0"]
            assert np.any(h5["components_ch0"][name][:] > 0)


def test_hanging_source_matrix_preserves_left_right_symmetry(tmp_path, monkeypatch):
    """A state-index mistake in matrix_N breaks this symmetry immediately."""
    _write_input(
        tmp_path,
        tes_internal_model="hanging",
        C_tes_hanging=2.0e-13,
        **{"G_tes-hanging": 2.0e-7},
    )
    monkeypatch.setattr(simulation, "output", str(tmp_path))
    simulation.MakeNoise()

    pairs = [
        ("johnson_tes1", "johnson_tes2"),
        ("johnson_load1", "johnson_load2"),
        ("phonon_tes1_bath", "phonon_tes2_bath"),
        ("phonon_tes1_absorber_effective", "phonon_tes2_absorber_effective"),
        ("phonon_tes1_hanging", "phonon_tes2_hanging"),
    ]
    with h5py.File(tmp_path / "noise.h5", "r") as h5:
        for left, right in pairs:
            np.testing.assert_allclose(
                h5["components_ch0"][left][:],
                h5["components_ch1"][right][:],
                rtol=2.0e-12,
                atol=1.0e-24,
            )


def test_hanging_zero_conductance_is_exact_five_state_limit(tmp_path, monkeypatch):
    none_dir = tmp_path / "none"
    zero_dir = tmp_path / "zero"
    none_dir.mkdir()
    zero_dir.mkdir()
    _write_input(none_dir)
    _write_input(
        zero_dir,
        tes_internal_model="hanging",
        tes_johnson_model="constant_M",
        C_tes_hanging=2.0e-13,
        **{"G_tes-hanging": 0.0},
    )
    monkeypatch.setattr(simulation, "output", str(none_dir))
    simulation.MakeNoise()
    monkeypatch.setattr(simulation, "output", str(zero_dir))
    simulation.MakeNoise()

    with h5py.File(none_dir / "noise.h5", "r") as lhs, h5py.File(
        zero_dir / "noise.h5", "r"
    ) as rhs:
        assert rhs.attrs["noise_model"] == lhs.attrs["noise_model"]
        np.testing.assert_allclose(rhs["total"][:], lhs["total"][:])
        np.testing.assert_allclose(rhs["transfer_ch0"][:], lhs["transfer_ch0"][:])


def test_lorentzian_resistance_fluctuation_is_independent_source(tmp_path, monkeypatch):
    _write_input(
        tmp_path,
        tes_resistance_fluctuation_model="lorentzian",
        resistance_fluctuation_M0=0.5,
        resistance_fluctuation_tau_s=1.0e-5,
    )
    monkeypatch.setattr(simulation, "output", str(tmp_path))
    simulation.MakeNoise()
    with h5py.File(tmp_path / "noise.h5", "r") as h5:
        assert h5.attrs["tes_resistance_fluctuation_model"] == "lorentzian"
        assert h5["transfer_ch0"].shape[1] == 10
        assert h5["components_ch0"]["resistance_fluctuation_tes1"][0] > 0
        assert h5["components_ch0"]["resistance_fluctuation_tes1"][-1] < h5[
            "components_ch0"
        ]["resistance_fluctuation_tes1"][0]


def test_hanging_and_lorentzian_models_compose(tmp_path, monkeypatch):
    _write_input(
        tmp_path,
        tes_internal_model="hanging",
        C_tes_hanging=2.0e-13,
        **{
            "G_tes-hanging": 2.0e-7,
            "tes_resistance_fluctuation_model": "lorentzian",
            "resistance_fluctuation_M0": 0.5,
            "resistance_fluctuation_tau_s": 1.0e-5,
        },
    )
    monkeypatch.setattr(simulation, "output", str(tmp_path))
    simulation.MakeNoise()
    with h5py.File(tmp_path / "noise.h5", "r") as h5:
        assert h5["transfer_ch0"].shape[1] == 12
        assert h5.attrs["noise_model"] == "seven_state_hanging_tes"
        assert "resistance_fluctuation_tes2" in h5["components_ch1"]


def test_resistance_fluctuation_none_is_legacy_output(tmp_path, monkeypatch):
    legacy = tmp_path / "legacy"
    explicit = tmp_path / "explicit"
    legacy.mkdir()
    explicit.mkdir()
    _write_input(legacy)
    _write_input(explicit, tes_resistance_fluctuation_model="none")
    monkeypatch.setattr(simulation, "output", str(legacy))
    simulation.MakeNoise()
    monkeypatch.setattr(simulation, "output", str(explicit))
    simulation.MakeNoise()
    with h5py.File(legacy / "noise.h5", "r") as lhs, h5py.File(
        explicit / "noise.h5", "r"
    ) as rhs:
        assert lhs["transfer_ch0"].shape[1] == rhs["transfer_ch0"].shape[1] == 8
        assert json.loads(lhs.attrs["source_names"]) == json.loads(
            rhs.attrs["source_names"]
        )
        np.testing.assert_array_equal(lhs["total"][:], rhs["total"][:])
        np.testing.assert_array_equal(
            lhs["transfer_ch0"][:], rhs["transfer_ch0"][:]
        )
