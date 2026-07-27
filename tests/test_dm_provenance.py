import csv
import pathlib

ROOT = pathlib.Path(__file__).parents[1]
TAB = ROOT / "crossmatching/dm_provenance.csv"

REQUIRED = {"nickname",
            "dm_dsa", "dm_dsa_err", "dm_dsa_method", "dm_dsa_source",
            "dm_chime", "dm_chime_err", "dm_chime_method",
            "dm_chime_source", "delta_dm", "delta_dm_sigma"}


def test_dm_provenance_covers_all_twelve():
    rows = list(csv.DictReader(TAB.open()))
    assert len(rows) == 12
    assert REQUIRED <= set(rows[0])
    assert all(r["dm_dsa_method"] and r["dm_chime_method"]
               for r in rows)
    assert all(r["dm_dsa_source"] and r["dm_chime_source"]
               for r in rows)


def _load_builder():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "build_dm_provenance", ROOT / "scripts/build_dm_provenance.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_builder_raises_without_dsa_catalog():
    import pytest
    mod = _load_builder()
    with pytest.raises(RuntimeError, match="no longer carries the independent DSA"):
        mod.build_rows()


def test_builder_sources_dsa_from_catalog(tmp_path, monkeypatch):
    import json
    mod = _load_builder()
    cat = tmp_path / "catalog.csv"
    cat.write_text(
        "nick,dsa_dm,dsa_sigma\n"
        "zach,100.5,0.02\n"
        "johndoeII,200.25,0.03\n"
    )
    inputs = tmp_path / "chime_side_inputs.json"
    inputs.write_text(json.dumps([
        {"name": "zach", "dm_chime": 100.1, "dm_chime_err": 0.01,
         "method": "arrival regression"},
        {"name": "johndoeii", "dm_chime": None, "dm_chime_err": None,
         "dm_status": "unconstrained"},
    ]))
    monkeypatch.setattr(mod, "CHIME_INPUTS", inputs)
    rows = mod.build_rows(cat)
    assert [r["dm_dsa"] for r in rows] == ["100.5000", "200.2500"]
    assert [r["dm_dsa_err"] for r in rows] == ["0.0200", "0.0300"]
