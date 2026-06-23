"""Tests for analysis/flux_cal.py — radiometer flux calibration (S/N -> Jy)."""

import numpy as np

from analysis.flux_cal import calibrated_band_integral_jy_ms_hz, radiometer_sigma_jy


def test_sigma_jy_analytic():
    # SEFD=2000 Jy, n_pol=2, dnu=1e6 Hz, dt=1e-3 s, G=1 -> 2000/sqrt(2*1e6*1e-3) = sqrt(2000)
    s = radiometer_sigma_jy(2000.0, n_pol=2, dnu_hz=1e6, dt_s=1e-3, g=1.0)
    assert abs(s - 2000.0 / np.sqrt(2000.0)) < 1e-9
    # beam attenuation G=0.5 doubles the noise
    assert abs(radiometer_sigma_jy(2000.0, 2, 1e6, 1e-3, 0.5) - 2.0 * s) < 1e-9


def test_band_integral_flat_oracle():
    # flat per-channel S/N integral A, flat sigma_S=s0, band [nu1,nu2]:
    # integral = A*s0*dt_ms*(nu2-nu1)
    nf = 64
    freq_hz = np.linspace(1.311e9, 1.499e9, nf)
    sn_integrated = np.full(nf, 3.0)  # per-channel sum_onpulse(S/N)
    sigma_jy = np.full(nf, 5.0)  # per-channel sigma_S [Jy]
    dt_ms = 0.131072
    i_band = calibrated_band_integral_jy_ms_hz(sn_integrated, sigma_jy, freq_hz, dt_ms)
    oracle = 3.0 * 5.0 * dt_ms * (freq_hz[-1] - freq_hz[0])
    assert abs(i_band - oracle) / oracle < 1e-9


def test_sn_spectrum_synthetic(tmp_path):
    # synthetic (n_freq, n_time) S/N data: unit noise + a centred Gaussian pulse on every channel
    from analysis.flux_cal import sn_spectrum_from_npy

    rng = np.random.default_rng(0)
    nf, nt = 64, 512
    t = np.arange(nt)
    prof = np.exp(-0.5 * ((t - nt // 2) / 4.0) ** 2)
    data = rng.standard_normal((nf, nt)) + 8.0 * prof[None, :]  # S/N peak ~8
    p = tmp_path / "synth.npy"
    np.save(p, data)
    freq_hz, sn_int, dt_ms, dnu_hz = sn_spectrum_from_npy(
        p, telescope="dsa", f_factor=1, t_factor=1
    )
    assert sn_int.shape == (nf,)
    assert np.all(np.isfinite(sn_int))
    assert np.median(sn_int) > 5.0  # integrated on-pulse S/N well above zero
    assert freq_hz[0] < freq_hz[-1] and dt_ms > 0 and dnu_hz > 0


def test_dsa_sigma_jy_constant_at_boresight():
    from analysis.flux_cal import dsa_sigma_jy

    freq_hz = np.linspace(1.311e9, 1.499e9, 32)
    dnu_hz, dt_s = 30517.578, 1.31072e-4
    sig = dsa_sigma_jy(
        freq_hz,
        dnu_hz,
        sefd_jy=4000.0,
        dt_s=dt_s,
        theta_deg=0.0,
        phi_deg=0.0,
        beam_gain_fn=lambda th, ph, f: 1.0,
    )
    expect = 4000.0 / np.sqrt(2 * dnu_hz * dt_s)  # G=1, n_pol=2 -> constant analytic value
    assert np.allclose(sig, expect)


def test_dsa_beam_offset_is_dec_difference():
    from analysis.flux_cal import dsa_beam_offset

    # transit geometry: same-RA boresight -> separation is exactly the dec difference
    theta, phi = dsa_beam_offset(dec_src=70.31, dec_pointing=72.0)
    assert abs(theta - 1.69) < 1e-6 and phi == 0.0


def test_burst_epoch_position_from_yaml():
    from analysis.flux_cal import burst_epoch_position

    mjd, ra, dec = burst_epoch_position("casey")
    assert abs(mjd - 60369.371) < 1e-3
    assert abs(ra - 169.983542) < 1e-6 and abs(dec - 70.676222) < 1e-6
