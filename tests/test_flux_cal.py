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
