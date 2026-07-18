from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from fit_cross_band import (
    Point,
    build_result,
    intrinsic_scatter_power_law,
    load_dsa_points,
    load_rigorous_points,
    weighted_power_law,
)


def _point(nu: float, gamma: float, error: float = 0.01) -> Point:
    return Point("test", nu, gamma, error, True, None, "synthetic", 0)


def test_weighted_power_law_recovers_known_parameters() -> None:
    alpha = 4.0
    gamma_ref = 0.3
    points = [_point(nu, gamma_ref * (nu / 1000.0) ** alpha) for nu in (600, 700, 1300, 1450)]
    result = weighted_power_law(points)
    assert result["alpha"] == pytest.approx(alpha, abs=1e-10)
    assert result["gamma_ref_mhz"] == pytest.approx(gamma_ref, abs=1e-10)
    assert result["chi_square"] == pytest.approx(0.0, abs=1e-18)


def test_intrinsic_scatter_fit_recovers_slope_with_scatter() -> None:
    rng = np.random.default_rng(42)
    frequencies = np.geomspace(550.0, 1550.0, 40)
    log_gamma = (
        np.log(0.35) + 3.7 * np.log(frequencies / 1000.0) + rng.normal(0.0, 0.15, frequencies.size)
    )
    points = [
        _point(float(nu), float(np.exp(gamma)), 0.01)
        for nu, gamma in zip(frequencies, log_gamma, strict=True)
    ]
    result = intrinsic_scatter_power_law(points)
    assert result["alpha"] == pytest.approx(3.7, abs=0.15)
    assert result["intrinsic_log_scatter"] == pytest.approx(0.15, abs=0.05)


def test_dsa_gate_rejects_failed_off_pulse_null(tmp_path: Path) -> None:
    payload = {
        "subbands": [
            {
                "center_freq_mhz": 1400.0,
                "selected_components": [{"dnu_mhz": 1.0, "dnu_err": 0.1, "quality_flags": []}],
                "off_pulse_null": {"null_pass": False},
                "low_lag_stability": {"stable": True},
            }
        ]
    }
    path = tmp_path / "dsa.json"
    path.write_text(json.dumps(payload))
    point = load_dsa_points(path)[0]
    assert not point.accepted
    assert point.exclusion_reason == "off-pulse null did not pass"


def test_nonpositive_input_fails_closed() -> None:
    with pytest.raises(ValueError, match="at least two accepted points"):
        weighted_power_law([Point("test", 1000.0, 1.0, 0.1, False, "bad", "test", 0)])


def test_near_zero_intrinsic_scatter_has_no_spurious_infinite_error() -> None:
    points = [_point(nu, 0.3 * (nu / 1000.0) ** 4.0) for nu in (600, 700, 900, 1300, 1450)]
    result = intrinsic_scatter_power_law(points)
    assert result["intrinsic_log_scatter"] < 1e-4
    assert result["intrinsic_log_scatter_uncertainty_identifiable"] is False
    assert result["intrinsic_log_scatter_err_plus"] is None


def _rigorous_payload(*, accepted: bool = False) -> dict:
    failed = [] if accepted else ["matched_injection:failed"]
    return {
        "schema": "flits.rigorous-scintillation-campaign/v1",
        "bands": {
            band: {
                "subbands": [
                    {
                        "index": 0,
                        "center_frequency_mhz": frequency,
                        "accepted_for_cross_band": accepted,
                        "qualification": {"qualified": accepted, "failed": failed},
                        "central_fit": {
                            "fit_ok": True,
                            "components": {
                                "bandwidth": {
                                    "gamma_mhz": gamma,
                                    "total_sigma_mhz": 0.1 * gamma,
                                    "admitted": accepted,
                                }
                            },
                        },
                    }
                ]
            }
            for band, frequency, gamma in (
                ("CHIME/FRB", 700.0, 0.1),
                ("DSA-110", 1400.0, 1.0),
            )
        },
    }


def test_rigorous_loader_rejects_legacy_dsa_schema(tmp_path: Path) -> None:
    path = tmp_path / "legacy.json"
    path.write_text(json.dumps({"subbands": []}))
    with pytest.raises(ValueError, match="rigorous-scintillation-campaign"):
        load_rigorous_points(path)


def test_rigorous_loader_preserves_failed_gate_reasons(tmp_path: Path) -> None:
    path = tmp_path / "rigorous.json"
    path.write_text(json.dumps(_rigorous_payload()))
    points = load_rigorous_points(path)
    assert len(points) == 2
    assert all(not point.accepted for point in points)
    assert all(point.exclusion_reason == "matched_injection:failed" for point in points)


def test_build_result_refuses_joint_fit_without_two_qualified_dsa_points(
    tmp_path: Path,
) -> None:
    path = tmp_path / "rigorous.json"
    path.write_text(json.dumps(_rigorous_payload()))
    result = build_result(path)
    assert result["joint_fit"]["available"] is False
    assert result["joint_fit"]["reason"] == "requires at least two qualified points per band"
    assert result["formal_no_extra_scatter"] is None
    assert result["intrinsic_scatter"] is None
