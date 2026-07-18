from __future__ import annotations

import numpy as np
import pytest

from scintillation.scint_analysis.rigorous_campaign import (
    REQUIRED_QUALIFICATION_GATES,
    bootstrap_acf_fit,
    combine_uncertainties,
    fit_acf_contract,
    generalized_lorentzian_sensitivity,
    moving_block_resample,
    qualify_gates,
    summarize_components,
)


def _lorentz(lags: np.ndarray, gamma: float = 0.8, m: float = 0.7) -> np.ndarray:
    return m**2 / (1.0 + (lags / gamma) ** 2) + 0.01


def test_components_are_ordered_and_modulation_semantics_are_explicit() -> None:
    report = summarize_components(
        [
            {"dnu_mhz": 5.0, "dnu_err": 0.2, "m": 0.8, "m_err": 0.03},
            {"dnu_mhz": 0.5, "dnu_err": 0.05, "m": 0.6, "m_err": 0.02},
        ],
        channel_width_mhz=0.05,
        fit_range_mhz=20.0,
        model_stable=True,
    )

    assert [row["role"] for row in report["components"]] == ["narrow", "broad"]
    assert report["bandwidth"]["gamma_mhz"] == pytest.approx(0.5)
    assert report["m_narrow"]["value"] == pytest.approx(0.6)
    assert report["m_broad"]["value"] == pytest.approx(0.8)
    assert report["m_total"]["value"] == pytest.approx(1.0)
    assert report["m_total"]["eligible"] is True


def test_unresolved_broad_component_blocks_broad_and_total_but_not_narrow() -> None:
    report = summarize_components(
        [
            {"dnu_mhz": 0.5, "dnu_err": 0.05, "m": 0.6, "m_err": 0.02},
            {"dnu_mhz": 19.0, "dnu_err": 1.0, "m": 0.8, "m_err": 0.03},
        ],
        channel_width_mhz=0.05,
        fit_range_mhz=20.0,
        model_stable=True,
    )

    assert report["bandwidth"]["eligible"] is True
    assert report["m_narrow"]["eligible"] is True
    assert report["m_broad"]["eligible"] is False
    assert report["m_total"]["eligible"] is False
    assert "width_near_fit_limit" in report["m_broad"]["reasons"]


def test_missing_qualification_gate_fails_closed() -> None:
    gates = {name: True for name in REQUIRED_QUALIFICATION_GATES}
    gates.pop("matched_injection")
    result = qualify_gates(gates)
    assert result["qualified"] is False
    assert result["failed"] == ["matched_injection:missing"]


def test_moving_block_resample_is_deterministic_and_length_preserving() -> None:
    values = np.arange(17.0)
    first = moving_block_resample(values, block_length=4, rng=np.random.default_rng(19))
    second = moving_block_resample(values, block_length=4, rng=np.random.default_rng(19))
    np.testing.assert_array_equal(first, second)
    assert first.shape == values.shape


def test_common_fit_recovers_known_lorentzian() -> None:
    lags = np.arange(-10.0, 10.0001, 0.05)
    acf = _lorentz(lags)
    err = np.full_like(acf, 0.01)
    fit = fit_acf_contract(
        lags,
        acf,
        err,
        channel_width_mhz=0.05,
        fit_range_mhz=8.0,
        first_positive_lag=1,
        max_components=2,
    )
    assert fit["fit_ok"] is True
    assert fit["n_preferred"] == 1
    assert fit["components"]["bandwidth"]["gamma_mhz"] == pytest.approx(0.8, rel=0.02)
    assert fit["components"]["m_narrow"]["value"] == pytest.approx(0.7, rel=0.02)


def test_block_bootstrap_reselects_models_deterministically() -> None:
    rng = np.random.default_rng(7)
    lags = np.arange(-8.0, 8.0001, 0.05)
    acf = _lorentz(lags) + rng.normal(0.0, 0.004, lags.size)
    err = np.full_like(acf, 0.01)
    kwargs = dict(
        channel_width_mhz=0.05,
        fit_range_mhz=6.0,
        first_positive_lag=1,
        max_components=2,
        n_bootstrap=24,
        block_length=6,
        seed=1234,
    )
    first = bootstrap_acf_fit(lags, acf, err, **kwargs)
    second = bootstrap_acf_fit(lags, acf, err, **kwargs)
    assert first == second
    assert first["n_success"] >= 20
    assert first["gamma_mhz"]["q16"] < first["gamma_mhz"]["q84"]
    assert sum(first["model_counts"].values()) == first["n_success"]


def test_uncertainties_are_combined_without_hiding_terms() -> None:
    result = combine_uncertainties(
        covariance_sigma=0.1,
        bootstrap_q16=0.7,
        bootstrap_q84=1.1,
        systematic_half_range=0.3,
    )
    assert result["bootstrap_sigma"] == pytest.approx(0.2)
    assert result["total_sigma"] == pytest.approx(np.sqrt(0.1**2 + 0.2**2 + 0.3**2))


def test_alternative_shape_keeps_central_broad_component_in_two_scale_acf() -> None:
    lags = np.arange(-20.0, 20.0001, 0.05)
    acf = (
        0.6**2 / (1.0 + (lags / 0.5) ** 2)
        + 0.4**2 / (1.0 + (lags / 5.0) ** 2)
        + 0.01
    )
    error = np.full_like(acf, 0.005)
    central = fit_acf_contract(
        lags,
        acf,
        error,
        channel_width_mhz=0.05,
        fit_range_mhz=15.0,
        max_components=2,
    )
    assert central["n_preferred"] == 2

    alternative = generalized_lorentzian_sensitivity(central)

    assert alternative["pass"] is True
    assert alternative["gamma_mhz"] == pytest.approx(0.5, rel=0.05)
