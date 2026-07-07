import os
import pathlib

import numpy as np
import pytest

DATA = pathlib.Path(os.environ.get("FLITS_DATA", "~/Data/Faber2026")).expanduser()
CUBES = sorted(DATA.rglob("*_32000b_cntr_bpc.npy"))

pytestmark = pytest.mark.skipif(not CUBES, reason="data root absent")


def _snr_profile(cube):
    arr = np.load(cube, mmap_mode="r")
    prof = np.nanmean(np.asarray(arr, float), axis=0)
    base = np.nanmedian(prof)
    noise = 1.4826 * np.nanmedian(np.abs(prof - base))
    return (prof - base) / max(noise, 1e-12)


@pytest.mark.parametrize("cube", CUBES, ids=lambda p: p.stem)
def test_no_edge_significance_and_centered(cube):
    snr = _snr_profile(cube)
    n = snr.size
    hot = snr > 5.0
    assert hot.any(), "no burst detected above 5 sigma"
    edge = n // 50
    assert not hot[:edge].any() and not hot[-edge:].any()
    idx = np.flatnonzero(hot)
    centroid = float((snr[idx] * idx).sum() / snr[idx].sum())
    assert 0.4 * n < centroid < 0.6 * n
