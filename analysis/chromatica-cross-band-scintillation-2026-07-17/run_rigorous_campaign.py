#!/usr/bin/env python3
"""Run the common, fail-closed Chromatica scintillation campaign."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import platform
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

import matplotlib
import numpy as np
import scipy

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scintillation.scint_analysis import analysis, figure_manifest, window_refit
from scintillation.scint_analysis import chime_artifact_guards as guards
from scintillation.scint_analysis import config as config_module
from scintillation.scint_analysis.acf_covariance import simulate_thin_screen_spectrum
from scintillation.scint_analysis.pipeline import ScintillationAnalysis
from scintillation.scint_analysis.rigorous_campaign import (
    bootstrap_acf_fit,
    combine_uncertainties,
    fit_acf_contract,
    generalized_lorentzian_sensitivity,
    qualify_gates,
)

SEED = 20260717
N_BOOTSTRAP = 80
N_INJECTION = 32
BLOCK_LENGTH = 8
DSA_CONFIG = Path("scintillation/configs/bursts/chromatica_dsa.yaml")
CHIME_CAMPAIGN = Path(
    "analysis/window-tuning-campaign-2026-07-17/results/chromatica_hi_campaign.json"
)
BAND_POLICIES = {
    "CHIME/FRB": {
        "fit_span_mhz": 5.0,
        "fit_span_variants_mhz": [1.0, 2.0, 3.0, 5.0],
        "harmonic_spacing_mhz": 0.390625,
        "harmonic_halfwidth_mhz": 0.05,
    },
    "DSA-110": {
        "fit_span_mhz": 25.0,
        "fit_span_variants_mhz": [8.0, 12.0, 18.0, 25.0],
        "harmonic_spacing_mhz": None,
        "harmonic_halfwidth_mhz": None,
    },
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, list | tuple):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return _jsonable(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _fresh_dsa_config(
    root: Path,
    *,
    burst_window: tuple[int, int],
    off_window: tuple[int, int],
    num_subbands: int,
    fixed_slices: list[list[int]] | None = None,
) -> dict[str, Any]:
    cfg = config_module.load_config(root / DSA_CONFIG)
    cfg = copy.deepcopy(cfg)
    cfg.pop("stored_fits", None)
    cfg["analysis"]["rfi_masking"]["manual_burst_window"] = list(burst_window)
    cfg["analysis"]["rfi_masking"]["manual_noise_window"] = list(off_window)
    acf_cfg = cfg["analysis"].setdefault("acf", {})
    acf_cfg["num_subbands"] = int(num_subbands)
    acf_cfg["use_snr_subbanding"] = True
    acf_cfg["first_fit_lag"] = 1
    if fixed_slices is None:
        acf_cfg.pop("subband_channel_slices", None)
    else:
        acf_cfg["subband_channel_slices"] = fixed_slices
    cfg["analysis"].setdefault("noise", {})["disable_template"] = True
    cfg["analysis"].setdefault("fit_2d", {})["enable"] = False
    pipe_opts = cfg.setdefault("pipeline_options", {})
    pipe_opts["force_recalc"] = True
    pipe_opts["save_intermediate_steps"] = False
    pipe_opts["halt_after_acf"] = True
    pipe_opts["diagnostic_plots"] = {"enable": False}
    pipe_opts["cache_directory"] = str(root / ".cache" / "rigorous-chromatica")
    return cfg


def _prepare_dsa(
    root: Path,
    *,
    burst_window: tuple[int, int],
    off_window: tuple[int, int],
    num_subbands: int,
    fixed_slices: list[list[int]] | None = None,
) -> dict[str, Any]:
    cfg = _fresh_dsa_config(
        root,
        burst_window=burst_window,
        off_window=off_window,
        num_subbands=num_subbands,
        fixed_slices=fixed_slices,
    )
    pipe = ScintillationAnalysis(cfg)
    pipe.run()
    if not pipe.acf_results or not pipe.acf_results.get("subband_acfs"):
        raise RuntimeError("DSA preparation produced no subband ACFs")
    return {
        "band": "DSA-110",
        "window": {"burst_lims": list(burst_window), "off_lims": list(off_window)},
        "config": cfg,
        "spectrum": pipe.masked_spectrum,
        "burst_lims": tuple(pipe.burst_lims),
        "off_lims": tuple(pipe.off_pulse_lims),
        "acf": pipe.acf_results,
    }


def _chime_weights(payload: dict[str, Any]) -> np.ndarray:
    record = payload["windows"]["weights"]
    values = np.asarray(record["values"], dtype=float)
    expected = int(record["t1"]) - int(record["t0"])
    if values.size != expected:
        raise ValueError("CHIME campaign weight span is inconsistent")
    return values


def _prepare_chime(
    *,
    burst_window: tuple[int, int],
    off_window: tuple[int, int],
    num_subbands: int,
    fixed_slices: list[list[int]] | None = None,
    time_weights: np.ndarray | None = None,
) -> dict[str, Any]:
    result = window_refit.refit(
        "chromatica_hi",
        burst_window,
        off_window,
        [],
        first_fit_lag=1,
        time_weights=time_weights,
        subband_channel_slices=fixed_slices,
        num_subbands=num_subbands,
        return_prepared=True,
    )
    common = result["common_acf"]
    acf = {
        **common,
        "subband_center_freqs_mhz": np.asarray(result["center_freqs"], dtype=float),
    }
    return {
        "band": "CHIME/FRB",
        "window": {
            "burst_lims": list(burst_window),
            "off_lims": list(off_window),
            "estimator": "matched-weight" if time_weights is not None else "boxcar",
        },
        "config": result["_prepared_config"],
        "spectrum": result["_prepared_spectrum"],
        "burst_lims": tuple(burst_window),
        "off_lims": tuple(off_window),
        "acf": acf,
    }


def _acf_value(acf: dict[str, Any], name: str, index: int) -> np.ndarray:
    return np.asarray(acf[name][index], dtype=float)


def _fit_one(
    prepared: dict[str, Any],
    index: int,
    *,
    fit_span_mhz: float,
    first_positive_lag: int = 1,
) -> dict[str, Any]:
    acf = prepared["acf"]
    policy = BAND_POLICIES[prepared["band"]]
    lags = _acf_value(acf, "subband_lags_mhz", index)
    available = float(np.nanmax(lags[lags > 0]))
    fit_span = min(float(fit_span_mhz), available)
    return fit_acf_contract(
        lags,
        _acf_value(acf, "subband_acfs", index),
        _acf_value(acf, "subband_acfs_err", index),
        channel_width_mhz=float(acf["subband_channel_widths_mhz"][index]),
        fit_range_mhz=fit_span,
        first_positive_lag=first_positive_lag,
        max_components=3,
        harmonic_spacing_mhz=policy["harmonic_spacing_mhz"],
        harmonic_halfwidth_mhz=policy["harmonic_halfwidth_mhz"],
    )


def _off_slice_starts(prepared: dict[str, Any], *, max_slices: int = 12) -> list[int]:
    burst_width = max(prepared["burst_lims"][1] - prepared["burst_lims"][0], 4)
    lo = prepared["off_lims"][0] + 2
    hi = prepared["off_lims"][1] - burst_width
    if hi <= lo:
        return []
    return list(range(lo, hi, burst_width + 4))[:max_slices]


def _artifact_controls(
    prepared: dict[str, Any], index: int, central: dict[str, Any]
) -> dict[str, Any]:
    if not central.get("fit_ok"):
        return {
            "off_pulse_null": guards.off_pulse_null_verdict(None, []),
            "low_lag_stability": guards.low_lag_stability_verdict(None, {}),
        }
    acf = prepared["acf"]
    channel_width = float(acf["subband_channel_widths_mhz"][index])
    channel_slice = tuple(int(value) for value in acf["subband_channel_slices"][index])
    fit_span = BAND_POLICIES[prepared["band"]]["fit_span_mhz"]
    spectrum = prepared["spectrum"]
    off_widths = []
    max_lag_bins = max(3, int(fit_span / channel_width) + 1)
    for start in _off_slice_starts(prepared):
        width = max(prepared["burst_lims"][1] - prepared["burst_lims"][0], 4)
        off_spectrum = spectrum.get_spectrum((start, start + width))[slice(*channel_slice)]
        acf_obj = analysis.calculate_acf(
            off_spectrum,
            channel_width,
            off_burst_spectrum_mean=None,
            max_lag_bins=max_lag_bins,
            first_fit_lag=1,
        )
        if acf_obj is None:
            continue
        fit = fit_acf_contract(
            acf_obj.lags,
            acf_obj.acf,
            acf_obj.err,
            channel_width_mhz=channel_width,
            fit_range_mhz=min(fit_span, float(np.max(acf_obj.lags))),
            first_positive_lag=1,
            max_components=3,
            harmonic_spacing_mhz=BAND_POLICIES[prepared["band"]][
                "harmonic_spacing_mhz"
            ],
            harmonic_halfwidth_mhz=BAND_POLICIES[prepared["band"]][
                "harmonic_halfwidth_mhz"
            ],
        )
        if fit.get("fit_ok"):
            off_widths.append(float(fit["components"]["bandwidth"]["gamma_mhz"]))

    on_width = float(central["components"]["bandwidth"]["gamma_mhz"])
    excision_widths = {}
    for first_lag in (2, 3, 4):
        fit = _fit_one(
            prepared,
            index,
            fit_span_mhz=fit_span,
            first_positive_lag=first_lag,
        )
        excision_widths[first_lag - 1] = (
            float(fit["components"]["bandwidth"]["gamma_mhz"])
            if fit.get("fit_ok")
            else None
        )
    return {
        "off_pulse_slice_starts": _off_slice_starts(prepared),
        "off_pulse_widths_mhz": off_widths,
        "off_pulse_null": guards.off_pulse_null_verdict(on_width, off_widths),
        "low_lag_stability": guards.low_lag_stability_verdict(
            on_width, excision_widths
        ),
    }


def _normalization_record(prepared: dict[str, Any], index: int) -> dict[str, Any]:
    acf = prepared["acf"]
    channel_slice = tuple(int(value) for value in acf["subband_channel_slices"][index])
    spectrum = prepared["spectrum"]
    on = spectrum.get_spectrum(prepared["burst_lims"])[slice(*channel_slice)]
    off = spectrum.get_spectrum(prepared["off_lims"])[slice(*channel_slice)]
    mean_on = float(np.ma.mean(on))
    mean_off = float(np.ma.mean(off))
    source_mean = mean_on - mean_off
    passed = bool(
        np.isfinite(mean_on)
        and np.isfinite(mean_off)
        and np.isfinite(source_mean)
        and source_mean > 0
        and on.count() >= 20
        and off.count() >= 20
    )
    return {
        "pass": passed,
        "mean_on": mean_on,
        "mean_off": mean_off,
        "source_mean": source_mean,
        "on_unmasked_channels": int(on.count()),
        "off_unmasked_channels": int(off.count()),
        "acf_denominator": "(mean_on - mean_off)^2",
    }


def _matched_injection(
    prepared: dict[str, Any],
    index: int,
    central: dict[str, Any],
    *,
    seed: int,
    coverage_half_width_mhz: float | None,
) -> dict[str, Any]:
    normalization = _normalization_record(prepared, index)
    starts = _off_slice_starts(prepared)
    if not central.get("fit_ok") or not normalization["pass"] or len(starts) < 3:
        return {
            "pass": False,
            "reason": "central fit, normalization, or three off-pulse slices unavailable",
            "n_off_slices": len(starts),
            "seed": seed,
        }
    acf = prepared["acf"]
    channel_width = float(acf["subband_channel_widths_mhz"][index])
    channel_slice = tuple(int(value) for value in acf["subband_channel_slices"][index])
    n_channels = channel_slice[1] - channel_slice[0]
    truth_gamma = float(central["components"]["bandwidth"]["gamma_mhz"])
    truth_m = min(float(central["components"]["m_narrow"]["value"]), 1.0)
    fit_span = BAND_POLICIES[prepared["band"]]["fit_span_mhz"]
    max_lag_bins = max(3, int(fit_span / channel_width) + 1)
    burst_width = max(prepared["burst_lims"][1] - prepared["burst_lims"][0], 4)
    rng = np.random.default_rng(seed)
    recovered_gamma = []
    recovered_m = []
    interval_covers = []
    spectrum = prepared["spectrum"]
    for trial in range(N_INJECTION):
        start = starts[trial % len(starts)]
        noise = spectrum.get_spectrum((start, start + burst_width))[slice(*channel_slice)]
        gain = simulate_thin_screen_spectrum(
            rng,
            n_channels,
            truth_gamma / channel_width,
            mod_index=truth_m,
            snr=np.inf,
        )
        injected = np.ma.MaskedArray(
            noise.data + normalization["source_mean"] * gain,
            mask=np.ma.getmaskarray(noise),
        )
        acf_obj = analysis.calculate_acf(
            injected,
            channel_width,
            off_burst_spectrum_mean=normalization["mean_off"],
            max_lag_bins=max_lag_bins,
            first_fit_lag=1,
        )
        if acf_obj is None:
            continue
        fit = fit_acf_contract(
            acf_obj.lags,
            acf_obj.acf,
            acf_obj.err,
            channel_width_mhz=channel_width,
            fit_range_mhz=min(fit_span, float(np.max(acf_obj.lags))),
            first_positive_lag=1,
            max_components=3,
            harmonic_spacing_mhz=BAND_POLICIES[prepared["band"]][
                "harmonic_spacing_mhz"
            ],
            harmonic_halfwidth_mhz=BAND_POLICIES[prepared["band"]][
                "harmonic_halfwidth_mhz"
            ],
        )
        if not fit.get("fit_ok"):
            continue
        width = fit["components"]["bandwidth"]
        recovered_gamma.append(float(width["gamma_mhz"]))
        recovered_m.append(float(fit["components"]["m_narrow"]["value"]))
        interval_covers.append(
            bool(
                coverage_half_width_mhz is not None
                and np.isfinite(coverage_half_width_mhz)
                and abs(float(width["gamma_mhz"]) - truth_gamma)
                <= coverage_half_width_mhz
            )
        )
    success_fraction = len(recovered_gamma) / N_INJECTION
    median_gamma = float(np.median(recovered_gamma)) if recovered_gamma else None
    median_m = float(np.median(recovered_m)) if recovered_m else None
    gamma_bias = abs(median_gamma / truth_gamma - 1.0) if median_gamma else None
    m_bias = abs(median_m / truth_m - 1.0) if median_m else None
    coverage = float(np.mean(interval_covers)) if interval_covers else 0.0
    passed = bool(
        success_fraction >= 0.8
        and gamma_bias is not None
        and gamma_bias <= 0.15
        and m_bias is not None
        and m_bias <= 0.15
        and coverage >= 0.6
    )
    return {
        "pass": passed,
        "seed": seed,
        "n_trials": N_INJECTION,
        "n_success": len(recovered_gamma),
        "success_fraction": success_fraction,
        "n_off_slices": len(starts),
        "truth_gamma_mhz": truth_gamma,
        "truth_m": truth_m,
        "median_gamma_mhz": median_gamma,
        "median_m": median_m,
        "median_gamma_relative_bias": gamma_bias,
        "median_m_relative_bias": m_bias,
        "coverage_half_width_mhz": coverage_half_width_mhz,
        "recovered_gamma_q16_q50_q84_mhz": (
            [float(value) for value in np.quantile(recovered_gamma, [0.16, 0.5, 0.84])]
            if recovered_gamma
            else None
        ),
        "recovered_m_q16_q50_q84": (
            [float(value) for value in np.quantile(recovered_m, [0.16, 0.5, 0.84])]
            if recovered_m
            else None
        ),
        "final_interval_coverage": coverage,
        "criteria": {
            "min_success_fraction": 0.8,
            "max_median_relative_bias": 0.15,
            "min_final_interval_coverage": 0.6,
        },
    }


def _variant_summary(
    central: dict[str, Any], variants: list[dict[str, Any]]
) -> dict[str, Any]:
    good = [row for row in variants if row.get("fit", {}).get("fit_ok")]
    if not central.get("fit_ok") or not good:
        return {"stable": False, "reason": "central or variant fits unavailable"}
    central_n = int(central["n_preferred"])
    central_gamma = float(central["components"]["bandwidth"]["gamma_mhz"])
    gammas = [float(row["fit"]["components"]["bandwidth"]["gamma_mhz"]) for row in good]
    model_counts = Counter(int(row["fit"]["n_preferred"]) for row in good)
    same_model_fraction = model_counts[central_n] / len(good)
    max_shift = max(abs(value / central_gamma - 1.0) for value in gammas)
    success_fraction = len(good) / len(variants) if variants else 0.0
    gamma_stable = bool(success_fraction >= 0.8 and max_shift <= 0.35)
    model_stable = bool(success_fraction >= 0.8 and same_model_fraction >= 0.6)

    def modulation_half_range(name: str) -> float | None:
        values = []
        for row in good:
            record = row["fit"]["components"][name]
            if record.get("value") is not None:
                values.append(float(record["value"]))
        return 0.5 * (max(values) - min(values)) if len(values) > 1 else None

    return {
        "stable": gamma_stable,
        "gamma_stable": gamma_stable,
        "model_stable": model_stable,
        "n_variants": len(variants),
        "n_success": len(good),
        "success_fraction": success_fraction,
        "central_n_preferred": central_n,
        "model_counts": {str(key): value for key, value in sorted(model_counts.items())},
        "same_model_fraction": same_model_fraction,
        "maximum_fractional_gamma_shift": max_shift,
        "gamma_systematic_half_range_mhz": 0.5 * (max(gammas) - min(gammas)),
        "m_narrow_systematic_half_range": modulation_half_range("m_narrow"),
        "m_broad_systematic_half_range": modulation_half_range("m_broad"),
        "m_total_systematic_half_range": modulation_half_range("m_total"),
        "criteria": {
            "min_success_fraction": 0.8,
            "min_same_model_fraction": 0.6,
            "max_fractional_gamma_shift": 0.35,
        },
    }


def _compact_fit(fit: dict[str, Any]) -> dict[str, Any]:
    """Drop repeated ACF arrays from variant records while retaining decisions."""
    if not fit.get("fit_ok"):
        return {
            "fit_ok": False,
            "reason": fit.get("reason"),
            "n_fit_points": fit.get("n_fit_points"),
        }
    return {
        key: fit[key]
        for key in (
            "fit_ok",
            "n_preferred",
            "criterion",
            "delta_bic",
            "f_test_p",
            "selected_bic",
            "selected_redchi",
            "n_fit_points",
            "components",
        )
    }


def _window_variants(root: Path, band: str, central: dict[str, Any]) -> list[dict[str, Any]]:
    fixed = [list(row) for row in central["acf"]["subband_channel_slices"]]
    prepared = []
    if band == "DSA-110":
        definitions = [
            ((1247, 1272), (0, 1181), "on_wide"),
            ((1250, 1269), (0, 1181), "on_core"),
            ((1249, 1270), (100, 1181), "off_late"),
        ]
        for burst, off, label in definitions:
            item = _prepare_dsa(
                root,
                burst_window=burst,
                off_window=off,
                num_subbands=4,
                fixed_slices=fixed,
            )
            item["variant_label"] = label
            prepared.append(item)
    else:
        payload = json.loads((root / CHIME_CAMPAIGN).read_text())
        definitions = []
        for row in payload["variants"]:
            windows = row["windows"]
            key = (tuple(windows["burst_lims"]), tuple(windows["off_lims"]))
            if key not in [(item[0], item[1]) for item in definitions]:
                definitions.append((key[0], key[1], windows.get("label", "boxcar")))
        for number, (burst, off, label) in enumerate(definitions):
            item = _prepare_chime(
                burst_window=burst,
                off_window=off,
                num_subbands=4,
                fixed_slices=fixed,
                time_weights=None,
            )
            item["variant_label"] = f"{label}_{number}"
            prepared.append(item)
    return prepared


def _partition_diagnostic(root: Path, band: str, central_payload: dict[str, Any]) -> dict[str, Any]:
    records = []
    for count in (2, 3, 4):
        if count == 4:
            prepared = central_payload
        elif band == "DSA-110":
            prepared = _prepare_dsa(
                root,
                burst_window=tuple(central_payload["burst_lims"]),
                off_window=tuple(central_payload["off_lims"]),
                num_subbands=count,
            )
        else:
            payload = json.loads((root / CHIME_CAMPAIGN).read_text())
            prepared = _prepare_chime(
                burst_window=tuple(central_payload["burst_lims"]),
                off_window=tuple(central_payload["off_lims"]),
                num_subbands=count,
                time_weights=_chime_weights(payload),
            )
        rows = []
        centers = np.asarray(prepared["acf"]["subband_center_freqs_mhz"], dtype=float)
        for index in range(len(prepared["acf"]["subband_acfs"])):
            fit = _fit_one(
                prepared,
                index,
                fit_span_mhz=BAND_POLICIES[band]["fit_span_mhz"],
            )
            rows.append(
                {
                    "index": index,
                    "center_frequency_mhz": float(centers[index]),
                    "fit_ok": bool(fit.get("fit_ok")),
                    "gamma_mhz": (
                        float(fit["components"]["bandwidth"]["gamma_mhz"])
                        if fit.get("fit_ok")
                        else None
                    ),
                    "n_preferred": fit.get("n_preferred"),
                }
            )
        records.append({"requested_num_subbands": count, "subbands": rows})
    return {
        "role": "campaign-level diagnostic; boundaries are not mapped one-to-one",
        "partitions": records,
    }


def _uncertainty_record(
    central: dict[str, Any], bootstrap: dict[str, Any], variants: dict[str, Any]
) -> dict[str, Any] | None:
    if not central.get("fit_ok") or not bootstrap.get("gamma_mhz"):
        return None
    bandwidth = central["components"]["bandwidth"]
    systematic = variants.get("gamma_systematic_half_range_mhz")
    if systematic is None:
        return None
    return combine_uncertainties(
        covariance_sigma=float(bandwidth["covariance_sigma_mhz"]),
        bootstrap_q16=float(bootstrap["gamma_mhz"]["q16"]),
        bootstrap_q84=float(bootstrap["gamma_mhz"]["q84"]),
        systematic_half_range=float(systematic),
    )


def _run_band(root: Path, band: str) -> dict[str, Any]:
    if band == "DSA-110":
        central_prepared = _prepare_dsa(
            root,
            burst_window=(1249, 1270),
            off_window=(0, 1181),
            num_subbands=4,
        )
    else:
        payload = json.loads((root / CHIME_CAMPAIGN).read_text())
        central_prepared = _prepare_chime(
            burst_window=tuple(payload["windows"]["burst_lims"]),
            off_window=tuple(payload["windows"]["off_lims"]),
            num_subbands=4,
            time_weights=_chime_weights(payload),
        )
    window_prepared = _window_variants(root, band, central_prepared)
    policy = BAND_POLICIES[band]
    centers = np.asarray(
        central_prepared["acf"]["subband_center_freqs_mhz"], dtype=float
    )
    subbands = []
    for index in range(len(central_prepared["acf"]["subband_acfs"])):
        central = _fit_one(
            central_prepared,
            index,
            fit_span_mhz=policy["fit_span_mhz"],
        )
        fit_variants = []
        for span in policy["fit_span_variants_mhz"]:
            for first_lag in (1, 2, 3):
                if span == policy["fit_span_mhz"] and first_lag == 1:
                    continue
                fit_variants.append(
                    {
                        "kind": "fit_policy",
                        "fit_span_mhz": span,
                        "first_positive_lag": first_lag,
                        "fit": _compact_fit(
                            _fit_one(
                                central_prepared,
                                index,
                                fit_span_mhz=span,
                                first_positive_lag=first_lag,
                            )
                        ),
                    }
                )
        for prepared in window_prepared:
            fit_variants.append(
                {
                    "kind": "window",
                    "label": prepared["variant_label"],
                    "window": prepared["window"],
                    "fit": _compact_fit(
                        _fit_one(
                            prepared,
                            index,
                            fit_span_mhz=policy["fit_span_mhz"],
                        )
                    ),
                }
            )
        variants = _variant_summary(central, fit_variants)
        acf = central_prepared["acf"]
        bootstrap = bootstrap_acf_fit(
            _acf_value(acf, "subband_lags_mhz", index),
            _acf_value(acf, "subband_acfs", index),
            _acf_value(acf, "subband_acfs_err", index),
            channel_width_mhz=float(acf["subband_channel_widths_mhz"][index]),
            fit_range_mhz=policy["fit_span_mhz"],
            first_positive_lag=1,
            max_components=3,
            n_bootstrap=N_BOOTSTRAP,
            block_length=BLOCK_LENGTH,
            seed=SEED + index + (100 if band == "DSA-110" else 0),
            harmonic_spacing_mhz=policy["harmonic_spacing_mhz"],
            harmonic_halfwidth_mhz=policy["harmonic_halfwidth_mhz"],
        )
        controls = _artifact_controls(central_prepared, index, central)
        normalization = _normalization_record(central_prepared, index)
        alternative = generalized_lorentzian_sensitivity(central)
        uncertainty = _uncertainty_record(central, bootstrap, variants)
        coverage_half_width = uncertainty["total_sigma"] if uncertainty else None
        injection = _matched_injection(
            central_prepared,
            index,
            central,
            seed=SEED + 1000 + index + (100 if band == "DSA-110" else 0),
            coverage_half_width_mhz=coverage_half_width,
        )
        gates = {
            "normalization": normalization["pass"],
            "fit_quality": bool(
                central.get("fit_ok")
                and central["components"]["bandwidth"]["eligible"]
            ),
            "off_pulse_null": controls["off_pulse_null"].get("null_pass"),
            "low_lag_stability": controls["low_lag_stability"].get("stable"),
            "bootstrap_stability": bootstrap.get("stable"),
            "variant_stability": variants.get("stable"),
            "alternative_shape": alternative.get("pass"),
            "matched_injection": injection.get("pass"),
        }
        qualification = qualify_gates(gates)
        accepted = bool(qualification["qualified"] and uncertainty is not None)
        if central.get("fit_ok"):
            central["components"]["bandwidth"]["admitted"] = accepted
            central["components"]["bandwidth"]["total_sigma_mhz"] = (
                uncertainty["total_sigma"] if uncertainty else None
            )
            for name, bootstrap_name, systematic_name in (
                ("m_narrow", "m_narrow", "m_narrow_systematic_half_range"),
                ("m_broad", "m_broad", "m_broad_systematic_half_range"),
                ("m_total", "m_total", "m_total_systematic_half_range"),
            ):
                record = central["components"][name]
                boot = bootstrap.get(bootstrap_name)
                systematic = variants.get(systematic_name)
                admitted = bool(
                    accepted
                    and record.get("eligible")
                    and boot is not None
                    and systematic is not None
                    and bootstrap.get("model_stable")
                    and variants.get("model_stable")
                )
                record["admitted"] = admitted
                record["total_sigma"] = None
                if admitted:
                    record["total_sigma"] = combine_uncertainties(
                        covariance_sigma=float(record["covariance_sigma"]),
                        bootstrap_q16=float(boot["q16"]),
                        bootstrap_q84=float(boot["q84"]),
                        systematic_half_range=float(systematic),
                    )["total_sigma"]
        subbands.append(
            {
                "index": index,
                "center_frequency_mhz": float(centers[index]),
                "channel_width_mhz": float(acf["subband_channel_widths_mhz"][index]),
                "num_channels": int(acf["subband_num_channels"][index]),
                "channel_slice": [int(value) for value in acf["subband_channel_slices"][index]],
                "central_fit": central,
                "normalization": normalization,
                "artifact_controls": controls,
                "bootstrap": bootstrap,
                "alternative_shape": alternative,
                "matched_injection": injection,
                "variant_stability": variants,
                "fit_variants": fit_variants,
                "bandwidth_uncertainty": uncertainty,
                "qualification": qualification,
                "accepted_for_cross_band": accepted,
            }
        )
    return {
        "band": band,
        "central_window": central_prepared["window"],
        "fit_policy": policy,
        "subbands": subbands,
        "subband_partition_diagnostic": _partition_diagnostic(
            root, band, central_prepared
        ),
    }


def _plot_band(record: dict[str, Any], output: Path) -> list[Path]:
    subbands = record["subbands"]
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.2), constrained_layout=True)
    for ax, subband in zip(axes.flat, subbands, strict=False):
        fit = subband["central_fit"]
        if fit.get("fit_ok"):
            x = np.asarray(fit["fit_lags_mhz"])
            y = np.asarray(fit["fit_acf"])
            error = np.asarray(fit["fit_err"])
            ax.errorbar(x, y, yerr=error, fmt=".", ms=2, lw=0.4, alpha=0.55)
            ax.plot(x, fit["model_acf"], color="black", lw=1.2)
            width = fit["components"]["bandwidth"]
            state = "admitted" if width.get("admitted") else "excluded"
            sigma = width.get("total_sigma_mhz")
            label = rf"$\gamma_1={width['gamma_mhz']:.3f}$ MHz"
            if sigma is not None:
                label += rf" $\pm {sigma:.3f}$"
            ax.text(0.97, 0.94, f"{label}\n{state}", transform=ax.transAxes,
                    ha="right", va="top", fontsize=8)
        ax.axhline(0.0, color="0.5", lw=0.5)
        ax.set_xlabel(r"$\Delta\nu$ (MHz)")
        ax.set_ylabel("normalized ACF")
        ax.text(0.03, 0.94, f"{subband['center_frequency_mhz']:.1f} MHz",
                transform=ax.transAxes, ha="left", va="top", fontsize=8)
    stem = output / ("chime_rigorous_acf_fits" if record["band"] == "CHIME/FRB"
                     else "dsa_rigorous_acf_fits")
    paths = []
    for suffix in ("png", "pdf"):
        path = stem.with_suffix(f".{suffix}")
        fig.savefig(path, dpi=240, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def _plot_modulation(records: list[dict[str, Any]], output: Path) -> list[Path]:
    fig, ax = plt.subplots(figsize=(6.2, 4.3), constrained_layout=True)
    colors = {"CHIME/FRB": "#4c78a8", "DSA-110": "#f58518"}
    markers = {"m_narrow": "o", "m_broad": "s", "m_total": "^"}
    for band_record in records:
        band = band_record["band"]
        for name, marker in markers.items():
            admitted_x, admitted_y, admitted_err = [], [], []
            excluded_x, excluded_y = [], []
            for subband in band_record["subbands"]:
                fit = subband["central_fit"]
                if not fit.get("fit_ok"):
                    continue
                item = fit["components"][name]
                if item.get("value") is None:
                    continue
                if item.get("admitted"):
                    admitted_x.append(subband["center_frequency_mhz"])
                    admitted_y.append(item["value"])
                    admitted_err.append(item["total_sigma"])
                else:
                    excluded_x.append(subband["center_frequency_mhz"])
                    excluded_y.append(item["value"])
            label = f"{band} {name.replace('_', ' ')}"
            if admitted_x:
                ax.errorbar(admitted_x, admitted_y, yerr=admitted_err, fmt=marker,
                            color=colors[band], capsize=2, label=label)
            if excluded_x:
                ax.scatter(excluded_x, excluded_y, marker="x", color=colors[band],
                           alpha=0.65, label=f"{label} excluded")
    ax.axhline(1.0, color="0.4", ls="--", lw=0.8, label=r"$m=1$")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\nu$ (MHz)")
    ax.set_ylabel("modulation index")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, fontsize=7, ncol=2)
    stem = output / "chromatica_modulation_qualification"
    paths = []
    for suffix in ("png", "pdf"):
        path = stem.with_suffix(f".{suffix}")
        fig.savefig(path, dpi=240, bbox_inches="tight")
        paths.append(path)
    plt.close(fig)
    return paths


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    chime = _run_band(root, "CHIME/FRB")
    dsa = _run_band(root, "DSA-110")
    raw_dsa = Path(dsa["subbands"][0]["normalization"].get("raw_path", ""))
    dsa_input = Path(
        os.path.expandvars(config_module.load_config(root / DSA_CONFIG)["input_data_path"])
    )
    chime_input = Path(
        os.path.expandvars(window_refit._base_config("chromatica_hi")["input_data_path"])
    )
    payload = {
        "schema": "flits.rigorous-scintillation-campaign/v1",
        "analysis": "Chromatica common CHIME/FRB and DSA-110 scintillation campaign",
        "status": "measurement_campaign_complete",
        "seed": SEED,
        "n_bootstrap": N_BOOTSTRAP,
        "n_matched_injection": N_INJECTION,
        "block_length": BLOCK_LENGTH,
        "bands": {"CHIME/FRB": chime, "DSA-110": dsa},
        "provenance": {
            "command": (
                "NUMBA_DISABLE_JIT=1 uv run --frozen python "
                "analysis/chromatica-cross-band-scintillation-2026-07-17/"
                "run_rigorous_campaign.py --root $PWD --output-dir "
                "analysis/chromatica-cross-band-scintillation-2026-07-17/rigorous"
            ),
            "git_head": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=root, text=True
            ).strip(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "dsa_input_path": str(dsa_input.resolve()),
            "dsa_input_sha256": _sha256(dsa_input),
            "chime_input_path": str(chime_input.resolve()),
            "chime_input_sha256": _sha256(chime_input),
            "runner_sha256": _sha256(Path(__file__)),
        },
    }
    del raw_dsa
    result_path = output / "chromatica_rigorous_scintillation.json"
    result_path.write_text(json.dumps(_jsonable(payload), indent=2, sort_keys=True) + "\n")
    figures = []
    for band_record in (chime, dsa):
        figures.extend(_plot_band(band_record, output))
    figures.extend(_plot_modulation([chime, dsa], output))
    for path in figures:
        if path.suffix == ".png":
            figure_manifest.register_figure(
                output,
                path.name,
                "Each fitted curve follows the central ACF peak; admitted/excluded labels "
                "and modulation markers agree with the rigorous campaign JSON.",
                campaign="Chromatica rigorous cross-band scintillation",
            )
    summary = {
        band: {
            "accepted": sum(
                row["accepted_for_cross_band"] for row in record["subbands"]
            ),
            "total": len(record["subbands"]),
        }
        for band, record in (("CHIME/FRB", chime), ("DSA-110", dsa))
    }
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
