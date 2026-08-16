import argparse
import json

from scripts.prep import run_singlepixel_resolution_autotune as autotune
from scripts.prep import run_singlepixel_resolution_pilot as pilot


def test_spatial_selection_uses_first_accepted_refinement() -> None:
    manifest = {"comparisons": [{"accepted": False}, {"accepted": True}]}
    assert autotune.selected_spatial_level(manifest, [1.0, 2.0, 4.0]) == (4.0, "accepted")


def test_spatial_selection_marks_unresolved_finest_level() -> None:
    manifest = {"comparisons": [{"accepted": False}, {"accepted": False}]}
    assert autotune.selected_spatial_level(manifest, [1.0, 2.0, 4.0]) == (4.0, "unresolved_at_finest_candidate")


def test_time_selection_keeps_finest_when_no_coarser_step_passes() -> None:
    manifest = {"comparisons": [{"accepted": False}]}
    assert autotune.selected_time_step(manifest, [0.625, 1.25]) == (0.625, "finest_candidate_required")


def test_context_tag_is_short_stable_and_distinguishes_settings() -> None:
    selected = {"sinx_layers": 1.0, "stycast_layers": 32.0, "tes_layers": 1.0}
    tag = autotune.context_tag(selected)
    assert tag == autotune.context_tag(dict(reversed(list(selected.items()))))
    assert tag.startswith("auto_")
    assert len(tag) <= 20
    assert tag != autotune.context_tag({"sinx_layers": 2.0, "stycast_layers": 32.0, "tes_layers": 1.0})


def test_autotune_manifest_path_uses_pilot_single_source(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(autotune, "ROOT", tmp_path)
    monkeypatch.setattr(pilot, "OUT_DIR", tmp_path / "pilots")
    axis = "stack_local_size_um"
    tag = "auto_base"
    _, expected = pilot.pilot_output_paths(axis, tag)
    expected.parent.mkdir(parents=True)

    def fake_run(*args, **kwargs) -> None:
        expected.write_text(json.dumps({"comparisons": [{"accepted": True}]}), encoding="utf-8")

    monkeypatch.setattr(autotune.subprocess, "run", fake_run)
    args = argparse.Namespace(
        end_us=105.0,
        early_step_us=0.625,
        linear_system="direct",
        tail_start_us=None,
        elmer_solver=None,
        runtime_bin=None,
        execute=True,
    )
    config = {"spatial_axes": {axis: [25.0, 16.6667]}}
    result = autotune.run_axis(args, config, axis, {})
    assert result["manifest"] == str(expected.relative_to(tmp_path))
    assert expected.name.endswith(f"_{pilot.PILOT_IMPLEMENTATION_SUFFIX}_pilot_manifest.json")
