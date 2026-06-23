"""Phase 5: DSA-band flux calibration wired into the energetics gate.

The "fluxcal" sentinel routes a band through the data-driven radiometer integral
(analysis.flux_cal.dsa_band_fluence_jy_ms_hz) instead of native-units x scalar, while the
both-bands-or-nothing gate stays closed until CHIME is calibrated too (Phase 6).
"""

import analysis.calculate_burst_energies as E


def test_dsa_calibrated_but_gate_closed(monkeypatch):
    monkeypatch.setattr(E, "dsa_band_fluence_jy_ms_hz", lambda nick: 1.234e6)
    rows = E.compute(scales={"C": None, "D": "fluxcal"})  # D via flux_cal, C still uncalibrated
    assert rows
    assert all("E_iso_erg" not in r for r in rows)  # gate closed -> no summed energy
    assert any(abs(r.get("I_DSA_jy_ms_hz", 0) - 1.234e6) < 1 for r in rows)  # DSA Jy on the row


def test_both_bands_fluxcal_opens_gate(monkeypatch):
    # when CHIME is also flux-calibrated, the gate opens and E_iso is summed from both Jy integrals
    monkeypatch.setattr(E, "dsa_band_fluence_jy_ms_hz", lambda nick: 2.0e6)
    rows = E.compute(scales={"C": 5.0, "D": "fluxcal"})  # C legacy scalar, D data-driven
    assert rows and all("E_iso_erg" in r for r in rows)
    r = rows[0]
    # (1+z) k-correction relates the two totals exactly; both per-band Jy energies are present
    assert abs(r["E_iso_erg"] - r["E_iso_erg_no_kcorr"] / (1.0 + r["z"])) < 1e-9 * r["E_iso_erg"]
    assert "E_iso_CHIME_erg" in r and "E_iso_DSA_erg" in r and r["I_DSA_jy_ms_hz"] == 2.0e6


def test_dsa_burst_config_resolves():
    # data-independent: the batch config maps nick -> the canonical .npy name + the fit's binning
    from analysis.flux_cal import _dsa_burst_config

    npy, f_factor, t_factor = _dsa_burst_config("chromatica")
    assert npy.name == "chromatica_dsa_I_272_368_2500b_cntr_bpc.npy"
    assert f_factor == 384 and t_factor == 2
