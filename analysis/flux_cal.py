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


def sn_spectrum_from_npy(inpath, telescope, f_factor=1, t_factor=1, onpulse_thresh=3.0):
    """Per-channel on-pulse S/N spectrum from a burst .npy via the FLITS loader.

    Returns (freq_hz ascending, sn_integrated [per channel, sum over on-pulse of data/noise_std],
    dt_ms, dnu_hz). BurstDataset z-scores each channel (io.py:131-146), so model.data is already
    S/N; dividing by the full-window per-channel noise makes the unit explicit and robust to the
    downsample rescaling.
    """
    import sys
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo / "scattering"))  # io.py uses package-relative imports (..burstfit)
    from scat_analysis.config_utils import load_telescope_block
    from scat_analysis.pipeline.io import BurstDataset

    tel = load_telescope_block(str(repo / "scattering" / "configs" / "telescopes.yaml"), telescope)
    ds = BurstDataset(
        inpath,
        inpath,
        telescope=tel,
        f_factor=f_factor,
        t_factor=t_factor,
        onpulse_crop=True,
        onpulse_thresh=onpulse_thresh,
    )
    m = ds.model
    noise = np.clip(m.noise_std, 1e-9, None)  # per-channel full-window noise (n_freq,)
    sn = m.data / noise[:, None]  # (n_freq, n_time_onpulse) signal-to-noise
    sn_integrated = np.nansum(sn, axis=1)  # sum over the cropped on-pulse window
    return m.freq * 1e9, sn_integrated, ds.dt_ms, ds.df_MHz * 1e6


def dsa_sigma_jy(freq_hz, dnu_hz, sefd_jy, dt_s, theta_deg, phi_deg, beam_gain_fn):
    """Per-channel radiometer noise sigma_S(nu) [Jy] for DSA (n_pol=2) at a beam offset.

    beam_gain_fn(theta_deg, phi_deg, freq_ghz) -> normalized power-beam gain (boresight=1);
    dnu_hz is the channel bandwidth (downsampled ds.df_MHz*1e6), dt_s the sample time [s].
    """
    g = np.array([beam_gain_fn(theta_deg, phi_deg, f / 1e9) for f in freq_hz])
    return radiometer_sigma_jy(sefd_jy, 2, dnu_hz, dt_s, g)


def dsa_beam_offset(dec_src, dec_pointing):
    """(theta_deg, phi_deg) of the source from the DSA primary-beam boresight.

    DSA-110 is a transit array: at the burst's meridian transit (HA~0) the boresight sits at the
    source RA, so the angular separation is exactly the declination difference; phi=0 (meridian).
    HA refinement (small E-W term) is out of scope -- see the plan.
    """
    return abs(dec_src - dec_pointing), 0.0


def burst_epoch_position(nick):
    """(mjd, ra_deg, dec_deg) for a burst nickname from configs/bursts.yaml."""
    from pathlib import Path

    import yaml

    b = yaml.safe_load(
        (Path(__file__).resolve().parents[1] / "configs" / "bursts.yaml").read_text()
    )
    e = b["bursts"][nick]
    return float(e["mjd"]), float(e["ra_deg"]), float(e["dec_deg"])


def _csv_lookup(filename, key, key_col="burst", val_col=None):
    import csv
    from pathlib import Path

    p = Path(__file__).resolve().parents[1] / "analysis" / "burst_energies" / filename
    if not p.exists():
        raise FileNotFoundError(f"{p} missing -- run the Phase 4 acquisition")
    for row in csv.DictReader(p.open()):
        if row[key_col] == key:
            return row[val_col]
    raise KeyError(f"{key} not in {p}")


def load_dsa_sefd(nick):
    """Measured DSA SEFD [Jy] nearest the burst epoch, from dsa_sefd.csv (Phase 4)."""
    return float(_csv_lookup("dsa_sefd.csv", nick, val_col="sefd_jy"))


def dsa_pointing_dec(nick):
    """DSA primary-beam pointing declination [deg] for a burst, from dsa_pointing.csv (Phase 4)."""
    return float(_csv_lookup("dsa_pointing.csv", nick, val_col="pointing_dec_deg"))


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
