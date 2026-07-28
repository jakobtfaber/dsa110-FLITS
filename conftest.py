"""Keep repository-root helper modules importable during broad test collection."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


_ROOT = Path(__file__).parent
_SCRIPTS_INIT = _ROOT / "scripts" / "__init__.py"
_loaded_scripts = sys.modules.get("scripts")
if _loaded_scripts is None or Path(getattr(_loaded_scripts, "__file__", "")).resolve() != _SCRIPTS_INIT:
    spec = importlib.util.spec_from_file_location(
        "scripts", _SCRIPTS_INIT, submodule_search_locations=[str(_SCRIPTS_INIT.parent)]
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules["scripts"] = module
    spec.loader.exec_module(module)
