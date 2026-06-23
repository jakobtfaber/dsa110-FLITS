#!/usr/bin/env python
"""Radiometer flux calibration: per-channel S/N -> Jy for FLITS dynamic spectra.

FLITS data is per-channel z-scored S/N (scattering/scat_analysis/pipeline/io.py:131-146:
_bandpass_correct divides each channel by its own off-pulse mean/std and the loader keeps
"units as S/N"). Physical flux density is therefore one radiometer multiply away:

    S_nu(nu,t) [Jy] = (S/N)(nu,t) * sigma_S(nu),
    sigma_S(nu) = SEFD(nu) / (sqrt(n_pol * dnu * dt) * G(theta,phi,nu)),

with SEFD = 2 k_B T_sys / A_eff and G the beam gain (boresight=1). Dividing by the
per-channel off-pulse noise also cancels the bandpass gain, so the S/N spectrum (and hence
this calibration) is free of the band-edge rolloff that rails the fitted gamma_D.

The band fluence integral (Jy*ms*Hz) returned here feeds
analysis/calculate_burst_energies.band_energy_erg with flux_scale=1 -- the per-channel scale
is already folded in. See docs/rse/specs/plan-radiometer-flux-cal.md.
"""

from __future__ import annotations

import numpy as np


def radiometer_sigma_jy(sefd_jy, n_pol, dnu_hz, dt_s, g):
    """Per-sample radiometer noise [Jy]: SEFD / (sqrt(n_pol*dnu*dt) * G).

    sefd_jy: system-equivalent flux density [Jy]; n_pol: summed polarisations;
    dnu_hz: channel bandwidth [Hz]; dt_s: sample time [s]; g: beam gain (boresight=1).
    """
    return sefd_jy / (np.sqrt(n_pol * dnu_hz * dt_s) * g)


def calibrated_band_integral_jy_ms_hz(sn_integrated, sigma_jy, freq_hz, dt_ms):
    """int_band [ sigma_S(nu) * dt_ms * sum_onpulse(S/N)(nu) ] dnu   [Jy*ms*Hz].

    sn_integrated: per-channel sum over the on-pulse window of (data/noise_std) [dimensionless];
    sigma_jy: per-channel sigma_S [Jy]; freq_hz: ascending channel centres [Hz]; dt_ms: sample [ms].
    """
    chan_fluence_jy_ms = sigma_jy * dt_ms * sn_integrated  # [Jy*ms] per channel
    return float(np.trapezoid(chan_fluence_jy_ms, freq_hz))  # [Jy*ms*Hz]


def _check() -> None:
    # 1. radiometer noise vs analytic: 2000/sqrt(2*1e6*1e-3) == sqrt(2000); G=0.5 doubles it
    s = radiometer_sigma_jy(2000.0, 2, 1e6, 1e-3, 1.0)
    assert abs(s - np.sqrt(2000.0)) < 1e-9, s
    assert abs(radiometer_sigma_jy(2000.0, 2, 1e6, 1e-3, 0.5) - 2.0 * s) < 1e-9
    # 2. flat-band integral oracle: trapz(const) == const*(nu2-nu1)
    freq_hz = np.linspace(1.311e9, 1.499e9, 64)
    i_band = calibrated_band_integral_jy_ms_hz(
        np.full(64, 3.0), np.full(64, 5.0), freq_hz, 0.131072
    )
    oracle = 3.0 * 5.0 * 0.131072 * (freq_hz[-1] - freq_hz[0])
    assert abs(i_band - oracle) / oracle < 1e-9, (i_band, oracle)
    print("self-check OK: radiometer noise + flat-band integral match analytic oracles")


if __name__ == "__main__":
    _check()
