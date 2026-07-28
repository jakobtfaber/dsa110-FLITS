"""Known-truth checks for factor-aware ACF fitting contracts."""

from __future__ import annotations

import numpy as np
import pytest

from scint_analysis.acf_fitting import (
    CHIME_COARSE_CHANNEL_WIDTH_MHZ,
    SUPPORTED_CHIME_UPCHANNEL_FACTORS,
    build_subband_plan,
    fit_lorentzian_components,
    lorentzian_component,
    multiplicative_lorentzian_acf,
    noise_corrected_modulation_index,
    spectrum_difference_noise_rms,
    summed_lorentzian_acf,
    total_modulation_index,
    total_modulation_index_uncertainty,
    validate_chime_factor_grid,
)


def _uniform_grid(n_channels: int, channel_width_mhz: float) -> np.ndarray:
    return 400.0 + np.arange(n_channels, dtype=float) * channel_width_mhz


@pytest.mark.parametrize("factor", SUPPORTED_CHIME_UPCHANNEL_FACTORS)
def test_factor_grid_and_single_lorentzian_recovery(factor):
    channel_width = CHIME_COARSE_CHANNEL_WIDTH_MHZ / factor
    validate_chime_factor_grid(factor, channel_width)
    positive_lags = np.arange(1, int(1.5 / channel_width) + 1) * channel_width
    lags = np.r_[-positive_lags[::-1], positive_lags]
    true_gamma = 0.2
    true_m = 0.7
    acf = lorentzian_component(lags, true_gamma, true_m) + 0.01
    errors = np.full_like(acf, 0.002)

    result = fit_lorentzian_components(
        lags,
        acf,
        max_components=1,
        acf_err=errors,
        channel_width_mhz=channel_width,
        upchannel_factor=factor,
    )
    component = result["fits"][0]["components"][0]

    assert result["upchannel_factor"] == factor
    assert result["fit_domain"] == "positive lags only; zero lag excluded"
    assert result["selection_scope"].startswith("diagnostic_only")
    assert component["dnu_mhz"] == pytest.approx(true_gamma, rel=2.0e-4)
    assert component["m"] == pytest.approx(true_m, rel=2.0e-4)


def test_invalid_factor_grid_is_rejected():
    with pytest.raises(ValueError, match="unsupported CHIME upchannel factor"):
        validate_chime_factor_grid(8, CHIME_COARSE_CHANNEL_WIDTH_MHZ / 8)
    with pytest.raises(ValueError, match="does not match"):
        validate_chime_factor_grid(32, 0.02)


def test_equal_channel_subbands_cover_grid_once():
    frequencies = _uniform_grid(103, 0.02)
    plan = build_subband_plan(frequencies, 4, mode="equal_channels")
    slices = [(item.start, item.stop) for item in plan.subbands]
    counts = [item.channel_count for item in plan.subbands]

    assert slices[0][0] == 0
    assert slices[-1][1] == frequencies.size
    assert all(left[1] == right[0] for left, right in zip(slices[:-1], slices[1:], strict=True))
    assert max(counts) - min(counts) <= 1
    assert plan.weight_definition == "channel count"


def test_equal_snr_subbands_use_snr_squared_not_total_signal():
    frequencies = _uniform_grid(100, 0.02)
    signal = np.r_[np.full(50, 4.0), np.ones(50)]
    noise = np.ones(100)
    plan = build_subband_plan(
        frequencies,
        2,
        mode="equal_snr",
        signal=signal,
        noise_rms=noise,
    )

    first, second = plan.subbands
    # Bright channels carry 16 times the S/N² weight, so half the information
    # is reached before the midpoint even though the intervals remain contiguous.
    assert first.stop < 50
    channel_weights = np.square(signal / noise)
    largest_channel_weight = float(channel_weights.max())
    assert abs(first.allocation_weight - second.allocation_weight) <= largest_channel_weight
    assert plan.weight_definition == "sum of independent-channel S/N squared"


def test_equal_snr_requires_explicit_noise():
    frequencies = _uniform_grid(32, 0.02)
    with pytest.raises(ValueError, match="requires per-channel signal and noise_rms"):
        build_subband_plan(frequencies, 4, mode="equal_snr", signal=np.ones(32))


def test_spectrum_difference_noise_uses_exact_mean_weights():
    off = np.tile(np.array([-1.0, 1.0, -1.0, 1.0]), (2, 1))
    on = np.ones((2, 4))
    on[1, -1] = np.nan
    weights = np.array([1.0, 2.0, 1.0, 0.0])
    measured = spectrum_difference_noise_rms(off, on, time_weights=weights)
    off_std = np.std(off[0], ddof=1)
    expected_first = off_std * np.sqrt((1 + 4 + 1) / 4**2 + 1 / 4)
    expected_second = expected_first  # the masked sample already had zero weight

    assert measured[0] == pytest.approx(expected_first)
    assert measured[1] == pytest.approx(expected_second)


def test_fixed_subbands_must_cover_every_channel():
    frequencies = _uniform_grid(32, 0.02)
    with pytest.raises(ValueError, match="stop at"):
        build_subband_plan(
            frequencies,
            2,
            mode="fixed",
            fixed_slices=[(0, 8), (8, 24)],
        )


def test_two_lorentzian_fit_recovers_separated_scales():
    positive_lags = np.arange(1, 401, dtype=float) * 0.005
    lags = np.r_[-positive_lags[::-1], positive_lags]
    true_gammas = np.array([0.06, 0.7])
    true_mods = np.array([0.5, 0.35])
    rng = np.random.default_rng(20260728)
    acf = summed_lorentzian_acf(lags, true_gammas, true_mods, baseline=0.01)
    acf += rng.normal(0.0, 5.0e-4, lags.size)
    result = fit_lorentzian_components(
        lags,
        acf,
        max_components=2,
        acf_err=np.full_like(acf, 5.0e-4),
    )

    assert result["n_preferred"] == 2
    fitted = sorted(result["fits"][1]["components"], key=lambda item: item["dnu_mhz"])
    assert [item["dnu_mhz"] for item in fitted] == pytest.approx(true_gammas, rel=0.08)
    assert [item["m"] for item in fitted] == pytest.approx(true_mods, rel=0.08)
    assert result["fits"][1]["modulation_parameterization"] == "phenomenological_sum"
    assert result["fits"][1]["m_total_is_physical_screen_model"] is False


def test_single_lorentzian_does_not_gain_unjustified_component():
    positive_lags = np.arange(1, 301, dtype=float) * 0.01
    lags = np.r_[-positive_lags[::-1], positive_lags]
    rng = np.random.default_rng(4)
    acf = lorentzian_component(lags, 0.3, 0.8) + 0.02
    acf += rng.normal(0.0, 0.002, lags.size)
    result = fit_lorentzian_components(
        lags,
        acf,
        max_components=2,
        acf_err=np.full_like(acf, 0.002),
    )

    assert result["n_preferred"] == 1


def test_fit_span_railed_width_is_a_limit_not_measurement():
    positive_lags = np.arange(1, 201, dtype=float) * 0.005
    lags = np.r_[-positive_lags[::-1], positive_lags]
    acf = lorentzian_component(lags, 10.0, 0.8) + 0.01
    result = fit_lorentzian_components(
        lags,
        acf,
        max_components=1,
        acf_err=np.full_like(acf, 0.001),
    )
    component = result["fits"][0]["components"][0]

    assert component["width_status"] == "lower_limit_from_fit_span"
    assert component["measurement_admissible"] is False


def test_multiplicative_two_screen_model_keeps_product_term():
    lags = np.array([0.0, 0.2, 1.0])
    gammas = [0.1, 1.0]
    mods = [1.0, 1.0]
    summed = summed_lorentzian_acf(lags, gammas, mods)
    physical = multiplicative_lorentzian_acf(lags, gammas, mods)

    assert summed[0] == pytest.approx(2.0)
    assert physical[0] == pytest.approx(3.0)
    assert total_modulation_index(
        mods, parameterization="phenomenological_sum"
    ) == pytest.approx(np.sqrt(2.0))
    assert total_modulation_index(
        mods, parameterization="multiplicative_screens"
    ) == pytest.approx(np.sqrt(3.0))
    assert np.all(physical >= summed)


def test_total_modulation_uncertainty_uses_full_covariance():
    mods = [0.6, 0.8]
    covariance = np.array([[0.01, -0.002], [-0.002, 0.04]])
    gradient = np.array(mods) / np.hypot(*mods)
    expected = np.sqrt(gradient @ covariance @ gradient)

    measured = total_modulation_index_uncertainty(
        mods,
        covariance,
        parameterization="phenomenological_sum",
    )

    assert measured == pytest.approx(expected)


def test_noise_corrected_modulation_index_matches_variance_definition():
    result = noise_corrected_modulation_index(
        mean_signal=10.0,
        observed_variance=13.0,
        noise_variance=4.0,
        mean_signal_err=0.1,
        observed_variance_err=0.5,
        noise_variance_err=0.3,
    )
    nondetection = noise_corrected_modulation_index(
        mean_signal=10.0,
        observed_variance=3.0,
        noise_variance=4.0,
    )

    assert result["modulation_index"] == pytest.approx(0.3)
    expected_variance_err = np.hypot(0.5, 0.3)
    expected_m_err = np.hypot(
        expected_variance_err / (2.0 * 10.0 * 3.0),
        3.0 * 0.1 / 10.0**2,
    )
    assert result["modulation_index_err"] == pytest.approx(expected_m_err)
    assert result["intrinsic_variance"] == pytest.approx(9.0)
    assert result["detected"] is True
    assert nondetection["modulation_index"] is None
    assert nondetection["detected"] is False
    assert nondetection["intrinsic_variance"] == pytest.approx(-1.0)
    assert nondetection["nonnegative_intrinsic_variance"] == pytest.approx(0.0)


def test_covariance_evidence_contract_is_factor_tagged_and_physical(monkeypatch):
    from scint_analysis import acf_evidence

    def fake_nested(_loglike, _prior, ndim, *_args):
        median = [0.2, 0.49, 0.0] if ndim == 3 else [0.1, 0.36, 5.0, 0.49, 0.0]
        return {
            "logz": 1.0 if ndim == 3 else 2.0,
            "logz_err": 0.1,
            "samples": np.tile(median, (100, 1)),
            "weights": np.full(100, 0.01),
            "median": median,
            "n_like_calls": 100,
        }

    monkeypatch.setattr(acf_evidence, "_run_nested", fake_nested)
    factor = 32
    channel_width = CHIME_COARSE_CHANNEL_WIDTH_MHZ / factor
    lags = np.arange(1, 21, dtype=float) * channel_width
    acf = lorentzian_component(lags, 0.2, 0.7)
    result = acf_evidence.compare_acf_evidence(
        lags,
        acf,
        np.eye(lags.size) * 0.01,
        channel_width,
        band_width_mhz=6.0,
        upchannel_factor=factor,
    )

    assert result["fit_contract"]["upchannel_factor"] == factor
    assert result["m1"]["params_median"]["m"] == pytest.approx(0.7)
    expected_total = total_modulation_index(
        [0.6, 0.7], parameterization="multiplicative_screens"
    )
    assert result["m2"]["params_median"]["m_total"] == pytest.approx(expected_total)
    assert result["m2"]["params_median"]["m_total_err_lower"] == pytest.approx(0.0)
    assert result["m2"]["params_median"]["m_total_err_upper"] == pytest.approx(0.0)
    assert (
        result["m2"]["modulation_parameterization"]
        == "multiplicative_screens_with_cross_terms"
    )


def test_covariance_evidence_rejects_factor_grid_mismatch(monkeypatch):
    from scint_analysis import acf_evidence

    monkeypatch.setattr(
        acf_evidence,
        "_run_nested",
        lambda *_args, **_kwargs: pytest.fail("sampler should not run"),
    )
    lags = np.arange(1, 21, dtype=float) * 0.02
    with pytest.raises(ValueError, match="does not match"):
        acf_evidence.compare_acf_evidence(
            lags,
            np.zeros_like(lags),
            np.eye(lags.size),
            channel_width_mhz=0.02,
            band_width_mhz=6.0,
            upchannel_factor=32,
        )
