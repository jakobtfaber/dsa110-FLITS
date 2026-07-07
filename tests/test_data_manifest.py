import csv
import pathlib

import pytest

MANIFEST = pathlib.Path(__file__).parents[1] / "data-manifest.csv"


def _rows():
    with MANIFEST.open() as fh:
        return list(csv.DictReader(fh))


def test_manifest_parses_and_is_nonempty():
    rows = _rows()
    assert len(rows) == 24


# xfail removed by P2.2 once the h17-resident rows are hashed.
@pytest.mark.xfail(reason="h17-resident rows pending P2.2", strict=False)
def test_manifest_has_no_pending_checksums():
    pending = [r["filename"] for r in _rows()
               if r["sha256"].strip().upper().startswith("PENDING")]
    assert pending == []
