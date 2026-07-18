"""Common fit, uncertainty, and qualification contract for scintillation ACFs.

Telescope adapters are responsible for calibration, masking, window selection, and
subband preparation.  Everything after an ACF is formed passes through this module so
CHIME/FRB and DSA-110 measurements have identical component semantics and gates.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np
from scipy.optimize import curve_fit

from .analysis import harmonic_lag_mask
from .revalidation import compare_lorentzian_components

REQUIRED_QUALIFICATION_GATES = (
    "normalization",
    "fit_quality",
    "off_pulse_null",
    "low_lag_stability",
    "bootstrap_stability",
    "variant_stability",
    "alternative_shape",
    "matched_injection",
)


def _finite_positive(value: Any) -> bool:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return False
    return bool(np.isfinite(number) and number > 0)


def _component_reasons(
    component: Mapping[str, Any], *, channel_width_mhz: float, fit_range_mhz: float
) -> tuple[list[str], list[str]]:
    width_reasons: list[str] = []
    modulation_reasons: list[str] = []
    gamma = float(component.get("dnu_mhz", np.nan))
    gamma_err = float(component.get("dnu_err", np.nan))
    modulation = float(component.get("m", np.nan))
    modulation_err = float(component.get("m_err", np.nan))

    if not _finite_positive(gamma):
        width_reasons.append("invalid_width")
    else:
        if gamma <= 2.0 * channel_width_mhz:
            width_reasons.append("unresolved_width")
        if gamma >= 0.8 * fit_range_mhz:
            width_reasons.append("width_near_fit_limit")
        if not np.isfinite(gamma_err) or gamma_err <= 0:
            width_reasons.append("missing_width_error")
        elif gamma_err / gamma >= 1.0:
            width_reasons.append("fractional_width_error_ge_1")

    if not _finite_positive(modulation):
        modulation_reasons.append("invalid_modulation")
    else:
        if modulation > 1.2:
            modulation_reasons.append("modulation_gt_1p2")
        if not np.isfinite(modulation_err) or modulation_err <= 0:
            modulation_reasons.append("missing_modulation_error")
        elif modulation_err / modulation >= 1.0:
            modulation_reasons.append("fractional_modulation_error_ge_1")
    return width_reasons, modulation_reasons


def _quadrature_value(components: Sequence[Mapping[str, Any]]) -> tuple[float, float]:
    values = np.asarray([float(row["m"]) for row in components], dtype=float)
    errors = np.asarray([float(row["m_err"]) for row in components], dtype=float)
    total = float(np.sqrt(np.sum(np.square(values))))
    if total <= 0 or not np.all(np.isfinite(errors)):
        return total, float("nan")
    error = float(np.sqrt(np.sum(np.square(values * errors))) / total)
    return total, error


def summarize_components(
    components: Sequence[Mapping[str, Any]],
    *,
    channel_width_mhz: float,
    fit_range_mhz: float,
    model_stable: bool,
    minimum_scale_separation: float = 4.0,
) -> dict[str, Any]:
    """Order fitted components and define bandwidth/modulation eligibility.

    Width and modulation eligibility are separate.  An unresolved broad component
    therefore blocks ``m_broad`` and ``m_total`` without silently deleting a usable
    narrow width or narrow modulation measurement.
    """
    ordered = [dict(row) for row in sorted(components, key=lambda row: float(row["dnu_mhz"]))]
    if not ordered:
        return {
            "components": [],
            "bandwidth": {"eligible": False, "reasons": ["no_components"]},
            "m_narrow": {"eligible": False, "reasons": ["no_components"]},
            "m_broad": {"eligible": False, "reasons": ["no_broad_component"]},
            "m_total": {"eligible": False, "reasons": ["no_components"]},
        }

    separation_ok = True
    for left, right in zip(ordered, ordered[1:], strict=False):
        separation = float(right["dnu_mhz"]) / float(left["dnu_mhz"])
        separation_ok &= bool(np.isfinite(separation) and separation >= minimum_scale_separation)

    reports = []
    for index, component in enumerate(ordered):
        width_reasons, modulation_reasons = _component_reasons(
            component,
            channel_width_mhz=channel_width_mhz,
            fit_range_mhz=fit_range_mhz,
        )
        if len(ordered) > 1 and not separation_ok:
            width_reasons.append("component_scale_separation_lt_4")
        if not model_stable:
            width_reasons.append("model_unstable")
        modulation_all = list(dict.fromkeys(width_reasons + modulation_reasons))
        reports.append(
            {
                **component,
                "role": "narrow" if index == 0 else "broad",
                "width_eligible": not width_reasons,
                "width_reasons": width_reasons,
                "modulation_eligible": not modulation_all,
                "modulation_reasons": modulation_all,
            }
        )

    narrow = reports[0]
    bandwidth = {
        "gamma_mhz": float(narrow["dnu_mhz"]),
        "covariance_sigma_mhz": float(narrow.get("dnu_err", np.nan)),
        "eligible": bool(narrow["width_eligible"]),
        "reasons": list(narrow["width_reasons"]),
    }
    m_narrow = {
        "value": float(narrow["m"]),
        "covariance_sigma": float(narrow.get("m_err", np.nan)),
        "eligible": bool(narrow["modulation_eligible"]),
        "reasons": list(narrow["modulation_reasons"]),
    }

    if len(reports) > 1:
        broad_value, broad_error = _quadrature_value(reports[1:])
        broad_reasons = list(
            dict.fromkeys(
                reason for row in reports[1:] for reason in row["modulation_reasons"]
            )
        )
        m_broad = {
            "value": broad_value,
            "covariance_sigma": broad_error,
            "eligible": not broad_reasons,
            "reasons": broad_reasons,
        }
    else:
        m_broad = {
            "value": None,
            "covariance_sigma": None,
            "eligible": False,
            "reasons": ["no_broad_component"],
        }

    total_value, total_error = _quadrature_value(reports)
    total_reasons = list(
        dict.fromkeys(reason for row in reports for reason in row["modulation_reasons"])
    )
    m_total = {
        "value": total_value,
        "covariance_sigma": total_error,
        "eligible": not total_reasons,
        "reasons": total_reasons,
    }
    return {
        "components": reports,
        "bandwidth": bandwidth,
        "m_narrow": m_narrow,
        "m_broad": m_broad,
        "m_total": m_total,
        "minimum_scale_separation": minimum_scale_separation,
        "scale_separation_pass": separation_ok,
    }


def qualify_gates(gates: Mapping[str, bool | None]) -> dict[str, Any]:
    """Require an explicit ``True`` for every scientific qualification gate."""
    failed = []
    normalized = {}
    for name in REQUIRED_QUALIFICATION_GATES:
        value = gates.get(name)
        normalized[name] = value
        if name not in gates:
            failed.append(f"{name}:missing")
        elif value is not True:
            failed.append(f"{name}:{'inconclusive' if value is None else 'failed'}")
    return {"qualified": not failed, "failed": failed, "gates": normalized}


def moving_block_resample(
    values: Sequence[float], *, block_length: int, rng: np.random.Generator
) -> np.ndarray:
    """Circular moving-block resample with the same length as ``values``."""
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size == 0:
        raise ValueError("values must be a non-empty one-dimensional array")
    block = int(block_length)
    if block < 1 or block > array.size:
        raise ValueError("block_length must be between 1 and len(values)")
    n_blocks = int(np.ceil(array.size / block))
    starts = rng.integers(0, array.size, size=n_blocks)
    offsets = np.arange(block)
    sampled = [array[(start + offsets) % array.size] for start in starts]
    return np.concatenate(sampled)[: array.size]


def _selected_fit(verdict: Mapping[str, Any]) -> Mapping[str, Any] | None:
    preferred = int(verdict.get("n_preferred", 1))
    return next(
        (
            row
            for row in verdict.get("fits", [])
            if int(row.get("n", -1)) == preferred and row.get("success")
        ),
        None,
    )


def _model_curve(lags: np.ndarray, fit: Mapping[str, Any]) -> np.ndarray:
    model = np.full_like(lags, float(fit.get("constant", 0.0)), dtype=float)
    for component in fit.get("components", []):
        gamma = float(component["dnu_mhz"])
        modulation = float(component["m"])
        model += modulation**2 / (1.0 + np.square(lags / gamma))
    return model


def fit_acf_contract(
    lags: Sequence[float],
    acf: Sequence[float],
    acf_err: Sequence[float] | None,
    *,
    channel_width_mhz: float,
    fit_range_mhz: float,
    first_positive_lag: int = 1,
    max_components: int = 3,
    harmonic_spacing_mhz: float | None = None,
    harmonic_halfwidth_mhz: float | None = None,
) -> dict[str, Any]:
    """Apply the common weighted fit and component-report contract to one ACF."""
    x = np.asarray(lags, dtype=float)
    y = np.asarray(acf, dtype=float)
    error = None if acf_err is None else np.asarray(acf_err, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y) & (x > 0) & (x <= fit_range_mhz)
    minimum_lag = (max(int(first_positive_lag), 1) - 0.5) * channel_width_mhz
    keep &= x >= minimum_lag
    if error is not None:
        keep &= np.isfinite(error) & (error > 0)
    if harmonic_spacing_mhz is not None and harmonic_halfwidth_mhz is not None:
        keep &= harmonic_lag_mask(x, harmonic_spacing_mhz, harmonic_halfwidth_mhz)
    x_fit = x[keep]
    y_fit = y[keep]
    error_fit = None if error is None else error[keep]
    if x_fit.size < 8:
        return {
            "fit_ok": False,
            "reason": "fewer_than_8_fit_lags",
            "n_fit_points": int(x_fit.size),
        }

    verdict = compare_lorentzian_components(
        x_fit,
        y_fit,
        max_components=max_components,
        acf_err=error_fit,
    )
    selected = _selected_fit(verdict)
    if selected is None:
        return {
            "fit_ok": False,
            "reason": "preferred_model_fit_failed",
            "n_fit_points": int(x_fit.size),
        }
    components = summarize_components(
        selected.get("components", []),
        channel_width_mhz=channel_width_mhz,
        fit_range_mhz=fit_range_mhz,
        model_stable=True,
    )
    model = _model_curve(x_fit, selected)
    return {
        "fit_ok": True,
        "n_preferred": int(verdict["n_preferred"]),
        "criterion": verdict.get("criterion"),
        "delta_bic": {str(k): float(v) for k, v in verdict.get("delta_bic", {}).items()},
        "f_test_p": {str(k): float(v) for k, v in verdict.get("f_test", {}).items()},
        "selected_bic": float(selected["bic"]),
        "selected_redchi": float(selected["redchi"]),
        "n_fit_points": int(x_fit.size),
        "fit_lags_mhz": x_fit.tolist(),
        "fit_acf": y_fit.tolist(),
        "fit_err": None if error_fit is None else error_fit.tolist(),
        "model_acf": model.tolist(),
        "residual_acf": (y_fit - model).tolist(),
        "components": components,
    }


def _quantiles(values: Sequence[float]) -> dict[str, float] | None:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return None
    q16, q50, q84 = np.quantile(array, [0.16, 0.5, 0.84])
    return {"q16": float(q16), "q50": float(q50), "q84": float(q84)}


def bootstrap_acf_fit(
    lags: Sequence[float],
    acf: Sequence[float],
    acf_err: Sequence[float] | None,
    *,
    channel_width_mhz: float,
    fit_range_mhz: float,
    first_positive_lag: int = 1,
    max_components: int = 3,
    n_bootstrap: int = 200,
    block_length: int = 8,
    seed: int = 0,
    harmonic_spacing_mhz: float | None = None,
    harmonic_halfwidth_mhz: float | None = None,
) -> dict[str, Any]:
    """Moving-block residual bootstrap with model reselection on every draw."""
    fit_kwargs = {
        "channel_width_mhz": channel_width_mhz,
        "fit_range_mhz": fit_range_mhz,
        "first_positive_lag": first_positive_lag,
        "max_components": max_components,
        "harmonic_spacing_mhz": harmonic_spacing_mhz,
        "harmonic_halfwidth_mhz": harmonic_halfwidth_mhz,
    }
    central = fit_acf_contract(lags, acf, acf_err, **fit_kwargs)
    if not central.get("fit_ok"):
        return {
            "status": "fail",
            "reason": central.get("reason", "central_fit_failed"),
            "seed": int(seed),
            "n_bootstrap": int(n_bootstrap),
        }

    x = np.asarray(central["fit_lags_mhz"], dtype=float)
    model = np.asarray(central["model_acf"], dtype=float)
    residual = np.asarray(central["residual_acf"], dtype=float)
    error = None if central["fit_err"] is None else np.asarray(central["fit_err"], dtype=float)
    rng = np.random.default_rng(seed)
    model_counts: Counter[int] = Counter()
    gammas: list[float] = []
    narrow_modulations: list[float] = []
    broad_modulations: list[float] = []
    total_modulations: list[float] = []

    for _ in range(int(n_bootstrap)):
        synthetic = model + moving_block_resample(residual, block_length=block_length, rng=rng)
        draw = fit_acf_contract(x, synthetic, error, **fit_kwargs)
        if not draw.get("fit_ok"):
            continue
        model_counts[int(draw["n_preferred"])] += 1
        component_report = draw["components"]
        gammas.append(float(component_report["bandwidth"]["gamma_mhz"]))
        narrow_modulations.append(float(component_report["m_narrow"]["value"]))
        if component_report["m_broad"]["value"] is not None:
            broad_modulations.append(float(component_report["m_broad"]["value"]))
        total_modulations.append(float(component_report["m_total"]["value"]))

    n_success = int(sum(model_counts.values()))
    preferred = int(central["n_preferred"])
    preferred_fraction = model_counts.get(preferred, 0) / n_success if n_success else 0.0
    gamma_quantiles = _quantiles(gammas)
    if gamma_quantiles is None:
        relative_half_interval = float("inf")
    else:
        relative_half_interval = (
            0.5 * (gamma_quantiles["q84"] - gamma_quantiles["q16"])
            / gamma_quantiles["q50"]
        )
    success_fraction = n_success / int(n_bootstrap) if n_bootstrap else 0.0
    gamma_stable = bool(success_fraction >= 0.8 and relative_half_interval <= 0.5)
    model_stable = bool(success_fraction >= 0.8 and preferred_fraction >= 0.7)
    return {
        "status": "pass" if gamma_stable else "fail",
        "stable": gamma_stable,
        "gamma_stable": gamma_stable,
        "model_stable": model_stable,
        "seed": int(seed),
        "n_bootstrap": int(n_bootstrap),
        "block_length": int(block_length),
        "n_success": n_success,
        "success_fraction": float(success_fraction),
        "central_n_preferred": preferred,
        "preferred_model_fraction": float(preferred_fraction),
        "model_counts": {str(k): int(v) for k, v in sorted(model_counts.items())},
        "gamma_relative_half_interval": float(relative_half_interval),
        "gamma_mhz": gamma_quantiles,
        "m_narrow": _quantiles(narrow_modulations),
        "m_broad": _quantiles(broad_modulations),
        "m_total": _quantiles(total_modulations),
    }


def combine_uncertainties(
    *,
    covariance_sigma: float,
    bootstrap_q16: float,
    bootstrap_q84: float,
    systematic_half_range: float,
) -> dict[str, float]:
    """Combine named independent uncertainty terms in quadrature."""
    bootstrap_sigma = 0.5 * (float(bootstrap_q84) - float(bootstrap_q16))
    terms = {
        "covariance_sigma": float(covariance_sigma),
        "bootstrap_sigma": float(bootstrap_sigma),
        "systematic_half_range": float(systematic_half_range),
    }
    if any(not np.isfinite(value) or value < 0 for value in terms.values()):
        raise ValueError(f"uncertainty terms must be finite and non-negative: {terms}")
    return {**terms, "total_sigma": float(np.sqrt(sum(value**2 for value in terms.values())))}


def generalized_lorentzian_sensitivity(
    fit: Mapping[str, Any], *, maximum_fractional_shift: float = 0.35
) -> dict[str, Any]:
    """Fit a generalized narrow component while retaining central broad terms.

    Fitting a single generalized curve to a selected two-scale ACF conflates the
    broad envelope with the narrow shape.  Broad Lorentzians are therefore held at
    their central values while the narrow width, amplitude, shape, and shared
    constant are varied.
    """
    if not fit.get("fit_ok"):
        return {"pass": False, "reason": "central_fit_failed"}
    x = np.asarray(fit["fit_lags_mhz"], dtype=float)
    y = np.asarray(fit["fit_acf"], dtype=float)
    error_values = fit.get("fit_err")
    sigma = None if error_values is None else np.asarray(error_values, dtype=float)
    central_gamma = float(fit["components"]["bandwidth"]["gamma_mhz"])
    central_m = float(fit["components"]["m_narrow"]["value"])
    broad_components = fit["components"]["components"][1:]

    def generalized(x_values, gamma, modulation, shape, constant):
        model = modulation**2 / (1.0 + np.power(x_values / gamma, shape)) + constant
        for component in broad_components:
            broad_gamma = float(component["dnu_mhz"])
            broad_m = float(component["m"])
            model += broad_m**2 / (1.0 + np.square(x_values / broad_gamma))
        return model

    spacing = float(np.median(np.diff(x)))
    try:
        params, covariance = curve_fit(
            generalized,
            x,
            y,
            p0=[central_gamma, central_m, 2.0, float(np.median(y[-max(3, x.size // 5) :]))],
            sigma=sigma,
            absolute_sigma=sigma is not None,
            bounds=([spacing / 10.0, 0.0, 0.5, -10.0], [2.0 * x.max(), 3.0, 8.0, 10.0]),
            maxfev=40000,
        )
    except Exception as exc:
        return {"pass": False, "reason": f"alternative_fit_failed:{type(exc).__name__}"}
    errors = np.sqrt(np.diag(covariance))
    shift = abs(float(params[0]) / central_gamma - 1.0)
    passed = bool(
        np.all(np.isfinite(params))
        and np.all(np.isfinite(errors))
        and shift <= maximum_fractional_shift
        and 0.55 < float(params[2]) < 7.95
    )
    return {
        "pass": passed,
        "gamma_mhz": float(params[0]),
        "gamma_err_mhz": float(errors[0]),
        "m": float(params[1]),
        "shape": float(params[2]),
        "shape_err": float(errors[2]),
        "fractional_gamma_shift": float(shift),
        "maximum_fractional_shift": float(maximum_fractional_shift),
        "broad_components_fixed": len(broad_components),
    }
