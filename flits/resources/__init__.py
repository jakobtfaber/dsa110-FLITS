"""Installed FLITS resource lookup."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path


def path(name: str) -> Path:
    """Return an installed generic FLITS resource."""
    resource = files(__package__).joinpath(name)
    if not resource.is_file():
        raise FileNotFoundError(f"unknown FLITS resource: {name}")
    return Path(str(resource))
