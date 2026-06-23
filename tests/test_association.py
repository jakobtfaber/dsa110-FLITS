"""Tests for crossmatching/association.py — CHIME-DSA association significance (pillars 1-4)."""

import math

import numpy as np
import pytest

from crossmatching.association import (
    OMEGA_WIN_BASELINE_DEG2,
    chance_mu,
    chance_probability,
    dm_agreement,
    expected_chance_associations,
    f_dm,
)

BASE = dict(rate_per_day=1000.0, omega_win_deg2=OMEGA_WIN_BASELINE_DEG2, dt_s=1.0, ddm=5.0)


# --- Pillar 1: chance-coincidence probability ---------------------------------
def test_chance_mu_regression_dm500():
    # pinned from the validated experiment (.agents/experiment-chance-coincidence-falsealarm.md)
    assert chance_mu(500.0, **BASE) == pytest.approx(5.023345e-09, rel=1e-4)


def test_chance_mu_scales_linearly_in_small_window():
    base = chance_mu(500.0, **BASE)
    assert chance_mu(500.0, **{**BASE, "dt_s": 2.0}) == pytest.approx(2 * base, rel=1e-9)
    assert chance_mu(500.0, **{**BASE, "ddm": 10.0}) == pytest.approx(2 * base, rel=1e-9)


def test_chance_matches_monte_carlo_in_measurable_regime():
    # analytic must equal a direct background MC where the MC has enough hits (mu ~ 0.046)
    infl = dict(rate_per_day=1000.0, omega_win_deg2=200.0, dt_s=3600.0, ddm=50.0)
    p_an = chance_probability(500.0, **infl)
    rng = np.random.default_rng(7)
    lam = chance_mu(500.0, **infl) / f_dm(500.0, 50.0)  # mean events in pos+time box
    n = 2_000_000
    counts = rng.poisson(lam, size=n)
    total = int(counts.sum())
    dms = np.exp(rng.normal(math.log(500.0), 0.7, size=total))
    hit = np.zeros(n, bool)
    hit[np.repeat(np.arange(n), counts)[np.abs(dms - 500.0) <= 50.0]] = True
    p_mc = hit.mean()
    assert p_mc == pytest.approx(p_an, rel=0.05)


def test_expected_chance_associations_sums_mu():
    dms = [262.4, 500.0, 960.1]
    assert expected_chance_associations(dms, **BASE) == pytest.approx(
        sum(chance_mu(d, **BASE) for d in dms), rel=1e-12
    )


# --- Pillar 2: independent DM agreement ---------------------------------------
def test_dm_agreement_consistent():
    r = dm_agreement(dm_chime=500.0, dm_chime_err=2.0, dm_dsa=502.0, dm_dsa_err=1.0)
    assert r["delta"] == pytest.approx(2.0)
    assert r["sigma"] == pytest.approx(math.sqrt(5.0))
    assert r["n_sigma"] == pytest.approx(2.0 / math.sqrt(5.0))
    assert r["consistent"] is True


def test_dm_agreement_inconsistent_beyond_3sigma():
    r = dm_agreement(dm_chime=500.0, dm_chime_err=1.0, dm_dsa=510.0, dm_dsa_err=1.0)
    assert r["consistent"] is False


def test_dm_agreement_missing_chime_dm_returns_null_reason():
    r = dm_agreement(dm_chime=None, dm_chime_err=None, dm_dsa=502.0, dm_dsa_err=1.0)
    assert r["consistent"] is None and "no CHIME DM" in r["reason"]
