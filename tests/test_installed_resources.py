from pathlib import Path

from flits.resources import path


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
