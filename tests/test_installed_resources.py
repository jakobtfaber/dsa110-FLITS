from pathlib import Path

from flits.resources import path
from galaxies.foreground import sightline_budget


def test_installed_generic_resources_are_available():
    for name in (
        "matplotlibrc",
        "scattering_sampler.yaml",
        "scattering_telescopes.yaml",
        "scintillation_chime.yaml",
        "scintillation_dsa.yaml",
    ):
        resource = path(name)
        assert isinstance(resource, Path)
        assert resource.is_file()
        assert resource.stat().st_size > 0


def test_budget_forwards_explicit_census_inputs(monkeypatch, tmp_path):
    seen = {}

    def fake_budget(*args, **kwargs):
        seen.update(kwargs)
        return {"name": args[0]}

    monkeypatch.setattr(sightline_budget, "build_sightline_budget", fake_budget)
    registry = tmp_path / "intervening_census_registry.csv"
    frame = sightline_budget.build_all_budgets(
        targets=[("synthetic", "00h00m00s", "+00d00m00s", 0.1)],
        registry_path=registry,
        census_data_dir=tmp_path,
    )
    assert list(frame["name"]) == ["synthetic"]
    assert seen["registry_path"] == registry
    assert seen["census_data_dir"] == tmp_path
