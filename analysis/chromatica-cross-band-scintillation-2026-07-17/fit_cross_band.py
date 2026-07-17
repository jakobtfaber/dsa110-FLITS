#!/usr/bin/env python3
"""Fit Chromatica's CHIME + DSA-110 scintillation bandwidth law."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2


@dataclass(frozen=True)
class Point:
    band: str
    frequency_mhz: float
    gamma_mhz: float
    gamma_err_mhz: float
    accepted: bool
    exclusion_reason: str | None
    uncertainty_policy: str
    subband_index: int


def _positive(value: Any, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and positive, got {value!r}")
    return result


def load_chime_points(path: Path) -> list[Point]:
    payload = json.loads(path.read_text())
    points = []
    for index, subband in enumerate(payload.get("subbands", [])):
        nu = _positive(subband.get("center_mhz"), f"CHIME subband {index} frequency")
        gamma = _positive(subband.get("gamma"), f"CHIME subband {index} gamma")
        fit_err = _positive(subband.get("gamma_err"), f"CHIME subband {index} fit error")
        scintle_err = _positive(
            subband.get("gamma_scintle_err"), f"CHIME subband {index} scintle error"
        )
        window_err = _positive(
            subband.get("gamma_win_sys"), f"CHIME subband {index} window systematic"
        )
        accepted = bool(subband.get("ok")) and bool(subband.get("resolved"))
        reason = None if accepted else "CHIME campaign point is not both ok and resolved"
        points.append(
            Point(
                band="CHIME/FRB",
                frequency_mhz=nu,
                gamma_mhz=gamma,
                gamma_err_mhz=math.sqrt(fit_err**2 + scintle_err**2 + window_err**2),
                accepted=accepted,
                exclusion_reason=reason,
                uncertainty_policy="quadrature(fit, finite_scintle, window_systematic)",
                subband_index=index,
            )
        )
    if not points:
        raise ValueError("CHIME input contains no subbands")
    return points


def load_dsa_points(path: Path) -> list[Point]:
    payload = json.loads(path.read_text())
    if "results" in payload:
        results = payload["results"]
        if len(results) != 1:
            raise ValueError("DSA aggregate input must contain exactly one burst")
        payload = results[0]
    points = []
    for index, subband in enumerate(payload.get("subbands", [])):
        components = sorted(
            subband.get("selected_components", []),
            key=lambda component: float(component.get("dnu_mhz", np.inf)),
        )
        if not components:
            raise ValueError(f"DSA subband {index} has no selected component")
        gamma_1 = components[0]
        nu = _positive(subband.get("center_freq_mhz"), f"DSA subband {index} frequency")
        gamma = _positive(gamma_1.get("dnu_mhz"), f"DSA subband {index} gamma_1")
        error = _positive(gamma_1.get("dnu_err"), f"DSA subband {index} gamma_1 error")

        reasons = []
        if gamma_1.get("quality_flags"):
            reasons.append("component flags: " + ", ".join(gamma_1["quality_flags"]))
        if subband.get("off_pulse_null", {}).get("null_pass") is not True:
            reasons.append("off-pulse null did not pass")
        if subband.get("low_lag_stability", {}).get("stable") is not True:
            reasons.append("low-lag stability did not pass")
        points.append(
            Point(
                band="DSA-110",
                frequency_mhz=nu,
                gamma_mhz=gamma,
                gamma_err_mhz=error,
                accepted=not reasons,
                exclusion_reason="; ".join(reasons) if reasons else None,
                uncertainty_policy="Lorentzian covariance error",
                subband_index=index,
            )
        )
    if not points:
        raise ValueError("DSA input contains no subbands")
    return points


def weighted_power_law(points: list[Point], nu_ref_mhz: float = 1000.0) -> dict[str, Any]:
    accepted = [point for point in points if point.accepted]
    if len(accepted) < 2:
        raise ValueError("at least two accepted points are required")
    nu = np.asarray([point.frequency_mhz for point in accepted])
    gamma = np.asarray([point.gamma_mhz for point in accepted])
    gamma_err = np.asarray([point.gamma_err_mhz for point in accepted])
    x = np.log(nu / nu_ref_mhz)
    y = np.log(gamma)
    sigma_y = gamma_err / gamma
    design = np.column_stack((np.ones_like(x), x))
    precision = 1.0 / np.square(sigma_y)
    normal = design.T @ (precision[:, None] * design)
    covariance = np.linalg.inv(normal)
    coefficients = covariance @ (design.T @ (precision * y))
    residual = (y - design @ coefficients) / sigma_y
    chi_square = float(residual @ residual)
    dof = len(accepted) - 2
    return {
        "model": "gamma_ref * (nu / nu_ref)^alpha",
        "nu_ref_mhz": nu_ref_mhz,
        "gamma_ref_mhz": float(np.exp(coefficients[0])),
        "gamma_ref_err_mhz": float(np.exp(coefficients[0]) * np.sqrt(covariance[0, 0])),
        "alpha": float(coefficients[1]),
        "alpha_err": float(np.sqrt(covariance[1, 1])),
        "covariance_log_gamma_ref_alpha": covariance.tolist(),
        "chi_square": chi_square,
        "dof": dof,
        "reduced_chi_square": chi_square / dof if dof > 0 else None,
        "goodness_of_fit_p": float(chi2.sf(chi_square, dof)) if dof > 0 else None,
        "n_points": len(accepted),
    }


def _intrinsic_nll(
    parameters: np.ndarray, x: np.ndarray, y: np.ndarray, sigma: np.ndarray
) -> float:
    log_gamma_ref, alpha, log_scatter = parameters
    variance = np.square(sigma) + np.exp(2.0 * log_scatter)
    residual = y - log_gamma_ref - alpha * x
    return float(0.5 * np.sum(np.square(residual) / variance + np.log(2.0 * np.pi * variance)))


def _numerical_hessian(function, location: np.ndarray) -> np.ndarray:
    location = np.asarray(location, dtype=float)
    steps = 1e-4 * np.maximum(1.0, np.abs(location))
    n = location.size
    hessian = np.empty((n, n), dtype=float)
    f0 = function(location)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = steps[i]
        hessian[i, i] = (function(location + ei) - 2.0 * f0 + function(location - ei)) / steps[
            i
        ] ** 2
        for j in range(i):
            ej = np.zeros(n)
            ej[j] = steps[j]
            value = (
                function(location + ei + ej)
                - function(location + ei - ej)
                - function(location - ei + ej)
                + function(location - ei - ej)
            ) / (4.0 * steps[i] * steps[j])
            hessian[i, j] = hessian[j, i] = value
    return hessian


def intrinsic_scatter_power_law(points: list[Point], nu_ref_mhz: float = 1000.0) -> dict[str, Any]:
    accepted = [point for point in points if point.accepted]
    if len(accepted) < 3:
        raise ValueError("at least three accepted points are required for intrinsic scatter")
    formal = weighted_power_law(accepted, nu_ref_mhz)
    nu = np.asarray([point.frequency_mhz for point in accepted])
    gamma = np.asarray([point.gamma_mhz for point in accepted])
    gamma_err = np.asarray([point.gamma_err_mhz for point in accepted])
    x = np.log(nu / nu_ref_mhz)
    y = np.log(gamma)
    sigma = gamma_err / gamma
    initial = np.array([np.log(formal["gamma_ref_mhz"]), formal["alpha"], np.log(0.2)], dtype=float)

    def objective(parameters: np.ndarray) -> float:
        return _intrinsic_nll(parameters, x, y, sigma)

    fit = minimize(
        objective,
        initial,
        method="L-BFGS-B",
        bounds=[(None, None), (None, None), (np.log(1e-6), np.log(10.0))],
    )
    if not fit.success:
        raise RuntimeError(f"intrinsic-scatter optimization failed: {fit.message}")
    hessian = _numerical_hessian(objective, fit.x)
    covariance = np.linalg.inv(hessian)
    errors = np.sqrt(np.diag(covariance))
    scatter = float(np.exp(fit.x[2]))
    scatter_uncertainty_identifiable = bool(
        scatter >= 1e-4 and np.isfinite(errors[2]) and fit.x[2] + errors[2] < 700.0
    )
    if scatter_uncertainty_identifiable:
        scatter_err_minus = float(scatter - np.exp(fit.x[2] - errors[2]))
        scatter_err_plus = float(np.exp(fit.x[2] + errors[2]) - scatter)
    else:
        scatter_err_minus = None
        scatter_err_plus = None
    return {
        "model": "gamma_ref * (nu / nu_ref)^alpha with Gaussian intrinsic log scatter",
        "nu_ref_mhz": nu_ref_mhz,
        "gamma_ref_mhz": float(np.exp(fit.x[0])),
        "gamma_ref_err_mhz": float(np.exp(fit.x[0]) * errors[0]),
        "alpha": float(fit.x[1]),
        "alpha_err": float(errors[1]),
        "intrinsic_log_scatter": scatter,
        "intrinsic_log_scatter_err_minus": scatter_err_minus,
        "intrinsic_log_scatter_err_plus": scatter_err_plus,
        "intrinsic_log_scatter_uncertainty_identifiable": scatter_uncertainty_identifiable,
        "covariance_log_gamma_ref_alpha_log_scatter": covariance.tolist(),
        "negative_log_likelihood": float(fit.fun),
        "n_points": len(accepted),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _dsa_raw_input_provenance(dsa_path: Path) -> dict[str, Any]:
    payload = json.loads(dsa_path.read_text())
    if "results" in payload:
        results = payload["results"]
        payload = results[0] if len(results) == 1 else {}
    raw_value = payload.get("input_data_path")
    if not raw_value:
        return {"dsa_raw_data_path": None, "dsa_raw_data_sha256": None}
    raw_path = Path(raw_value).expanduser()
    if not raw_path.is_file():
        raise FileNotFoundError(f"DSA raw input recorded by fit does not exist: {raw_path}")
    return {
        "dsa_raw_data_path": str(raw_path.resolve()),
        "dsa_raw_data_sha256": _sha256(raw_path),
    }


def build_result(chime_path: Path, dsa_path: Path) -> dict[str, Any]:
    points = load_chime_points(chime_path) + load_dsa_points(dsa_path)
    accepted = [point for point in points if point.accepted]
    formal = weighted_power_law(accepted)
    intrinsic = intrinsic_scatter_power_law(accepted)
    chime_only = weighted_power_law([point for point in points if point.band == "CHIME/FRB"])
    sensitivities = {}
    for point in accepted:
        if point.band != "DSA-110":
            continue
        retained = [candidate for candidate in accepted if candidate is not point]
        key = f"leave_out_dsa_subband_{point.subband_index}_{point.frequency_mhz:.3f}_mhz"
        sensitivities[key] = intrinsic_scatter_power_law(retained)

    residuals = []
    for point in accepted:
        x = math.log(point.frequency_mhz / formal["nu_ref_mhz"])
        observed = math.log(point.gamma_mhz)
        formal_expected = math.log(formal["gamma_ref_mhz"]) + formal["alpha"] * x
        measurement_sigma = point.gamma_err_mhz / point.gamma_mhz
        intrinsic_expected = math.log(intrinsic["gamma_ref_mhz"]) + intrinsic["alpha"] * x
        intrinsic_sigma = math.sqrt(measurement_sigma**2 + intrinsic["intrinsic_log_scatter"] ** 2)
        residuals.append(
            {
                "band": point.band,
                "subband_index": point.subband_index,
                "frequency_mhz": point.frequency_mhz,
                "formal_standardized_log_residual": (observed - formal_expected)
                / measurement_sigma,
                "intrinsic_scatter_standardized_log_residual": (observed - intrinsic_expected)
                / intrinsic_sigma,
            }
        )

    exact_law_rejected = bool(
        formal["goodness_of_fit_p"] is not None and formal["goodness_of_fit_p"] < 0.01
    )
    return {
        "analysis": "Chromatica CHIME/FRB + DSA-110 scintillation bandwidth power law",
        "selection_policy": {
            "chime": "campaign ok and resolved; narrow component",
            "chime_uncertainty": "quadrature of fit, finite-scintle, and window systematic",
            "dsa": "narrowest selected component with no component flags, off-pulse null pass, and low-lag stability pass",
            "dsa_uncertainty": "Lorentzian covariance error; no DSA window-campaign systematic available",
        },
        "inputs": {
            "chime_json": str(chime_path.resolve()),
            "chime_sha256": _sha256(chime_path),
            "dsa_json": str(dsa_path.resolve()),
            "dsa_sha256": _sha256(dsa_path),
            "fitter_sha256": _sha256(Path(__file__)),
            **_dsa_raw_input_provenance(dsa_path),
        },
        "points": [asdict(point) for point in points],
        "formal_no_extra_scatter": formal,
        "intrinsic_scatter": intrinsic,
        "chime_only_with_total_uncertainty": chime_only,
        "residuals": residuals,
        "sensitivity": sensitivities,
        "interpretation": {
            "primary_model": "intrinsic_scatter",
            "exact_single_power_law_rejected_at_p_lt_0p01": exact_law_rejected,
            "reason": (
                "The no-extra-scatter fit fails its chi-square goodness-of-fit test; "
                "the intrinsic-scatter model is the primary cross-band characterization."
            ),
        },
    }


def _write_csv(points: list[dict[str, Any]], path: Path) -> None:
    columns = list(points[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, lineterminator="\n")
        writer.writeheader()
        writer.writerows(points)


def _plot(result: dict[str, Any], output_stem: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from flits.plotting import use_flits_style

    use_flits_style()
    plt.rcParams.update({"pdf.fonttype": 42, "savefig.dpi": 300})
    fig, ax = plt.subplots(figsize=(6.2, 4.5), constrained_layout=True)
    colors = {"CHIME/FRB": "#4c78a8", "DSA-110": "#f58518"}
    markers = {"CHIME/FRB": "o", "DSA-110": "s"}
    for band in ("CHIME/FRB", "DSA-110"):
        accepted = [
            point for point in result["points"] if point["band"] == band and point["accepted"]
        ]
        excluded = [
            point for point in result["points"] if point["band"] == band and not point["accepted"]
        ]
        if accepted:
            ax.errorbar(
                [point["frequency_mhz"] for point in accepted],
                [point["gamma_mhz"] for point in accepted],
                yerr=[point["gamma_err_mhz"] for point in accepted],
                fmt=markers[band],
                color=colors[band],
                capsize=2.5,
                label=band,
                zorder=3,
            )
        if excluded:
            ax.scatter(
                [point["frequency_mhz"] for point in excluded],
                [point["gamma_mhz"] for point in excluded],
                marker="x",
                s=58,
                linewidth=1.6,
                color=colors[band],
                label=f"{band} excluded by artifact gate",
                zorder=4,
            )
    fit = result["intrinsic_scatter"]
    nu = np.geomspace(580.0, 1550.0, 400)
    model = fit["gamma_ref_mhz"] * (nu / fit["nu_ref_mhz"]) ** fit["alpha"]
    scatter = fit["intrinsic_log_scatter"]
    ax.plot(
        nu,
        model,
        color="black",
        lw=1.6,
        label=rf"cross-band fit: $\alpha={fit['alpha']:.2f}\pm{fit['alpha_err']:.2f}$",
    )
    ax.fill_between(
        nu,
        model * np.exp(-scatter),
        model * np.exp(scatter),
        color="black",
        alpha=0.12,
        label="intrinsic-scatter envelope",
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\nu$ (MHz)")
    ax.set_ylabel(r"$\gamma_1$ (MHz)")
    ax.grid(alpha=0.2, which="both")
    ax.legend(frameon=False, fontsize=8)
    for suffix in ("png", "pdf"):
        fig.savefig(output_stem.with_suffix(f".{suffix}"))
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chime-json", type=Path, required=True)
    parser.add_argument("--dsa-json", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = build_result(args.chime_json, args.dsa_json)
    (args.output_dir / "chromatica_cross_band_fit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    _write_csv(result["points"], args.output_dir / "chromatica_cross_band_points.csv")
    _plot(result, args.output_dir / "chromatica_cross_band_fit")
    primary = result["intrinsic_scatter"]
    formal = result["formal_no_extra_scatter"]
    print(
        f"alpha={primary['alpha']:.6f} +/- {primary['alpha_err']:.6f}; "
        f"gamma_1GHz={primary['gamma_ref_mhz']:.6f} MHz; "
        f"intrinsic_log_scatter={primary['intrinsic_log_scatter']:.6f}; "
        f"formal_chi2/dof={formal['chi_square']:.3f}/{formal['dof']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
