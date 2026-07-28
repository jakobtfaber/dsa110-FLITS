"""Check that the batch runner resolves the existing hand-tuned scint configs."""

from pathlib import Path

from flits.batch.batch_runner import discover_scint_configs

def test_discovers_hand_tuned_scint_configs(tmp_path):
    for telescope in ("chime", "dsa"):
        for nickname in ("casey", "johndoeII"):
            (tmp_path / f"{nickname}_{telescope}.yaml").write_text("{}\n")
    (tmp_path / "casey_chime_hi.yaml").write_text("{}\n")
    (tmp_path / "casey_dsa_temp.yaml").write_text("{}\n")

    found = discover_scint_configs(tmp_path, ["chime", "dsa"])
    dsa_bursts = {b for b, tels in found.items() if "dsa" in tels}
    assert dsa_bursts == {"casey", "johndoeii"}
    assert found["casey"].keys() >= {"chime", "dsa"}
    p = found["casey"]["dsa"]
    assert p.exists() and p.name == "casey_dsa.yaml"
    # variant configs (_hi, _temp) must not be swallowed as bursts
    assert "casey_chime_hi" not in found
    assert all("_temp" not in b and not b.endswith("_hi") for b in found)
    # keys are lowercased to match BurstInfo.from_filename burst names
    assert "johndoeii" in dsa_bursts
    assert "johndoeII" not in found
    assert found["johndoeii"]["dsa"].name == "johndoeII_dsa.yaml"


def test_missing_telescope_dir_is_empty_not_error(tmp_path):
    found = discover_scint_configs(tmp_path / "nope", ["chime", "dsa"])
    assert found == {}
