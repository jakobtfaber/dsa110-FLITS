"""Gate logic for the committed joint fits: Level-1 bounds, rail flag,
Level-2 (reused classify_fit_quality), Level-3 alpha-physics, worst-of-band."""

from gate_joint_committed import gate_one


def _fit(alpha, tau=0.3, bounds=(1.0, 6.0), err=None):
    a = {"median": alpha}
    if err is not None:
        a["err_minus"] = a["err_plus"] = err
    return {
        "alpha": a,
        "tau_1ghz": {"median": tau},
        "alpha_bounds": list(bounds),
    }


def test_alpha_subkolmogorov_is_marginal_not_fail():
    # ADR-0004: 1.0 <= alpha < 2.0 is sub-Kolmogorov -> L3 MARGINAL (inspect), not L1 FAIL.
    # alpha=1.4 with tight err sits 0.4 (>> 3*0.05 sigma) above the 1.0 floor -> not railed.
    v = gate_one("x", _fit(1.4, err=0.05), {"chi2_chime": 1.0, "chi2_dsa": 1.0})
    assert v["final"] == "MARGINAL"
    assert "sub-kolmogorov" in v["reason"].lower()
    assert v["rail"] is False


def test_alpha_below_hard_floor_fails_level1():
    # ADR-0004: alpha < 1.0 is achromatic (tau prop nu^-a meaningless) -> hard L1 FAIL.
    v = gate_one("x", _fit(0.9, err=0.05), {"chi2_chime": 1.0, "chi2_dsa": 1.0})
    assert v["final"] == "FAIL" and "alpha" in v["reason"].lower()


def test_alpha_at_floor_is_admitted_not_l1_fail():
    # alpha == 1.0 is admitted (L1 PASS), but sits on the prior bound -> rail-MARGINAL, never FAIL.
    v = gate_one("x", _fit(1.0, err=0.05), {"chi2_chime": 1.0, "chi2_dsa": 1.0})
    assert v["l1"] == "PASS" and v["final"] == "MARGINAL"


def test_wide_posterior_within_3sigma_of_bound_is_rail_marginal():
    # ADR-0004: oran case -- median 1.44 with a wide err (0.33) is within 3 sigma of the 1.0
    # floor -> rail-MARGINAL regardless of value (unconstrained, not a measurement).
    v = gate_one("x", _fit(1.44, err=0.33), {"chi2_chime": 1.0, "chi2_dsa": 1.0})
    assert v["rail"] is True and v["final"] == "MARGINAL"


def test_alpha_at_ceiling_fails_level1():
    v = gate_one("x", _fit(6.0), {"chi2_chime": 1.0, "chi2_dsa": 1.0})
    assert v["final"] == "FAIL"


def test_all_evaluable_levels_pass_is_capped_at_marginal():
    # Physical alpha, Kolmogorov window, good chi2, no rail: every EVALUABLE level
    # passes, but tau x dnu (L3) is not evaluable here, so the gate caps at MARGINAL
    # rather than certifying PASS on an incomplete contract check.
    v = gate_one("x", _fit(4.0), {"chi2_chime": 1.0, "chi2_dsa": 1.2})
    assert v["final"] == "MARGINAL"
    assert "tau x dnu" in v["reason"].lower()


def test_catastrophic_chi2_fails():
    v = gate_one("x", _fit(4.0), {"chi2_chime": 1.0, "chi2_dsa": 12.0})
    assert v["final"] == "FAIL"


def test_elevated_chi2_is_marginal():
    v = gate_one("x", _fit(4.0), {"chi2_chime": 1.0, "chi2_dsa": 2.0})
    assert v["final"] == "MARGINAL"


def test_alpha_off_kolmogorov_is_marginal():
    # physical (1.5<a<6) and good chi2, but a=2.7 deviates from 3.5-4.5 -> MARGINAL (L3)
    v = gate_one("x", _fit(2.7), {"chi2_chime": 1.0, "chi2_dsa": 1.0})
    assert v["final"] == "MARGINAL"


def test_rail_flagged_at_ceiling():
    # freya/chromatica/hamilton case: alpha pinned within EDGE of the prior ceiling.
    v = gate_one("x", _fit(5.99), {"chi2_chime": 1.1, "chi2_dsa": 1.1})
    assert v["rail"] is True


def test_no_rail_when_interior():
    v = gate_one("x", _fit(4.0), {"chi2_chime": 1.0, "chi2_dsa": 1.0})
    assert v["rail"] is False


def test_missing_ppc_is_marginal_unknown_l2():
    v = gate_one("x", _fit(4.0), None)
    assert v["final"] == "MARGINAL" and "chi2" in v["reason"].lower()


def test_incomplete_ppc_is_marginal_not_crash():
    # PPC present but missing a chi2_* key -> unknown chi2 (MARGINAL), not a format crash.
    v = gate_one("x", _fit(4.0), {"chi2_chime": 1.0})  # chi2_dsa absent
    assert v["final"] == "MARGINAL" and "incomplete" in v["reason"].lower()
