import json
from pathlib import Path

from scripts.support.run_phase24_short_difference_campaign import build_project


def test_campaign_project_has_three_one_microsecond_variants(tmp_path: Path, monkeypatch) -> None:
    # Keep the test independent of the repository's generated campaign file.
    import scripts.support.run_phase24_short_difference_campaign as campaign

    monkeypatch.setattr(campaign, "CAMPAIGN_DIR", tmp_path)
    monkeypatch.setattr(campaign, "PROJECT_PATH", tmp_path / "campaign.json")
    project_path, names = build_project()
    project = json.loads(project_path.read_text(encoding="utf-8"))

    assert set(names) == {"bdf2_reuse_on", "bdf2_reuse_off", "bdf1_reuse_off"}
    for name in names.values():
        case = project["cases"][name]
        assert len(case["timesteps"]) == 5
        assert len(case["output_intervals"]) == 5
        assert case["timesteps"][-1][0] == "100[ns]"
