"""Factor-aware spectral ACF fitting primitives.

This module owns three contracts that were previously spread across several
analysis paths:

* sub-band boundaries are defined once, in channel space, with their physical
  frequency bounds recorded;
* Lorentzian amplitudes are parameterized as ``m_i**2``;
* additive multiple-scale fits and multiplicative multiple-screen models have
  distinct names and distinct total-modulation formulae.

The additive model is a phenomenological decomposition.  It must not be
interpreted as the physical ACF of independent multiplicative screens, whose
cross terms are retained by :func:`multiplicative_lorentzian_acf`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass
from typing import Literal

import numpy as np
from lmfit import Model
from lmfit.models import ConstantModel
from scipy.stats import f as f_dist

SUPPORTED_CHIME_UPCHANNEL_FACTORS = (16, 32, 64, 128, 256, 512)
CHIME_COARSE_CHANNEL_WIDTH_MHZ = 0.390625

SubbandMode = Literal["equal_channels", "equal_snr", "fixed"]
LorentzianParameterization = Literal["phenomenological_sum", "multiplicative_screens"]


@dataclass(frozen=True)
class Subband:
    """One contiguous frequency-channel interval, using Python slice bounds."""

    start: int
    stop: int
    center_frequency_mhz: float
    bandwidth_mhz: float
    allocation_weight: float

    @property
    def channel_count(self) -> int:
        return self.stop - self.start


@dataclass(frozen=True)
class SubbandPlan:
    """Auditable sub-band plan for one factor-tagged spectrum."""

    mode: SubbandMode
    upchannel_factor: int | None
    channel_width_mhz: float
    subbands: tuple[Subband, ...]
    weight_definition: str

    def to_dict(self) -> dict:
        return {
            "mode": self.mode,
            "upchannel_factor": self.upchannel_factor,
            "channel_width_mhz": self.channel_width_mhz,
            "weight_definition": self.weight_definition,
            "subbands": [asdict(item) | {"channel_count": item.channel_count} for item in self.subbands],
        }


def validate_chime_factor_grid(upchannel_factor: int, channel_width_mhz: float) -> None:
    """Require a supported factor and its nominal CHIME fine-channel spacing."""

    factor = int(upchannel_factor)
    if factor not in SUPPORTED_CHIME_UPCHANNEL_FACTORS:
        raise ValueError(
            f"unsupported CHIME upchannel factor {factor}; "
            f"expected one of {SUPPORTED_CHIME_UPCHANNEL_FACTORS}"
        )
    expected = CHIME_COARSE_CHANNEL_WIDTH_MHZ / factor
    if not np.isclose(float(channel_width_mhz), expected, rtol=1.0e-6, atol=1.0e-12):
        raise ValueError(
            f"channel width {channel_width_mhz:.12g} MHz does not match U={factor} "
            f"nominal width {expected:.12g} MHz"
        )


def _validate_fixed_slices(
    fixed_slices: Sequence[Sequence[int]], n_channels: int
) -> list[tuple[int, int]]:
    slices = [(int(start), int(stop)) for start, stop in fixed_slices]
    if not slices:
        raise ValueError("fixed_slices must contain at least one interval")
    cursor = 0
    for start, stop in slices:
        if start != cursor or not (start < stop <= n_channels):
            raise ValueError(
                "fixed sub-band slices must be contiguous, non-empty, and cover "
                f"channels from zero; invalid interval {(start, stop)} after {cursor}"
            )
        cursor = stop
    if cursor != n_channels:
        raise ValueError(
            f"fixed sub-band slices stop at {cursor}, not the channel count {n_channels}"
        )
    return slices


def _weighted_partition(weights: np.ndarray, n_subbands: int) -> list[tuple[int, int]]:
    """Partition positive per-channel weights into contiguous near-equal totals."""

    total = float(np.sum(weights))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("equal-S/N sub-bands require positive finite total S/N² weight")
    cumulative = np.cumsum(weights)
    boundaries = [0]
    for index in range(1, n_subbands):
        target = total * index / n_subbands
        stop = int(np.searchsorted(cumulative, target, side="left") + 1)
        minimum = boundaries[-1] + 1
        maximum = weights.size - (n_subbands - index)
        boundaries.append(min(max(stop, minimum), maximum))
    boundaries.append(weights.size)
    return list(zip(boundaries[:-1], boundaries[1:], strict=True))


def spectrum_difference_noise_rms(
    off_pulse_power: np.ndarray | np.ma.MaskedArray,
    on_pulse_power: np.ndarray | np.ma.MaskedArray,
    *,
    time_weights: Sequence[float] | None = None,
) -> np.ma.MaskedArray:
    """Per-channel uncertainty of ``mean(on) - mean(off)``.

    The calculation uses the measured off-pulse standard deviation and the
    exact valid-sample weights in each channel. It assumes independent time
    samples; a campaign with measured temporal covariance must supply a
    covariance-aware replacement rather than relabel this estimate.
    """

    off = np.ma.masked_invalid(np.ma.asarray(off_pulse_power, dtype=float))
    on = np.ma.masked_invalid(np.ma.asarray(on_pulse_power, dtype=float))
    if off.ndim != 2 or on.ndim != 2 or off.shape[0] != on.shape[0]:
        raise ValueError("on- and off-pulse power must be channel-by-time arrays")
    if time_weights is None:
        base_weights = np.ones(on.shape[1], dtype=float)
    else:
        base_weights = np.asarray(time_weights, dtype=float)
        if base_weights.shape != (on.shape[1],):
            raise ValueError("time_weights must match the on-pulse time dimension")
        if np.any(~np.isfinite(base_weights)) or np.any(base_weights < 0):
            raise ValueError("time_weights must be finite and non-negative")
    valid_on = ~np.ma.getmaskarray(on)
    weights = valid_on * base_weights[None, :]
    sum_weights = np.sum(weights, axis=1)
    sum_squared_weights = np.sum(weights**2, axis=1)
    on_variance_factor = np.divide(
        sum_squared_weights,
        sum_weights**2,
        out=np.full(on.shape[0], np.nan),
        where=sum_weights > 0,
    )
    off_count = np.ma.count(off, axis=1)
    off_variance_factor = np.divide(
        1.0,
        off_count,
        out=np.full(on.shape[0], np.nan),
        where=off_count > 0,
    )
    off_std = np.ma.std(off, axis=1, ddof=1)
    noise = off_std * np.sqrt(on_variance_factor + off_variance_factor)
    invalid = (
        (off_count < 2)
        | (sum_weights <= 0)
        | ~np.isfinite(np.ma.filled(noise, np.nan))
        | (np.ma.filled(noise, 0.0) <= 0)
    )
    return np.ma.MaskedArray(np.ma.filled(noise, np.nan), mask=invalid)


def build_subband_plan(
    frequencies_mhz: Sequence[float],
    n_subbands: int,
    *,
    mode: SubbandMode,
    signal: Sequence[float] | np.ma.MaskedArray | None = None,
    noise_rms: Sequence[float] | np.ma.MaskedArray | None = None,
    fixed_slices: Sequence[Sequence[int]] | None = None,
    upchannel_factor: int | None = None,
) -> SubbandPlan:
    """Build contiguous equal-channel, equal-S/N, or fixed sub-bands.

    Equal-S/N uses the additive information weight
    ``q_i = (max(signal_i, 0) / noise_rms_i)**2``.  For independent channels,
    the sum of ``q_i`` is the squared signal-to-noise ratio.  Masked or invalid
    channels contribute zero.  This is intentionally not the historical
    equal-total-signal approximation.
    """

    frequencies = np.asarray(frequencies_mhz, dtype=float)
    if frequencies.ndim != 1 or frequencies.size < 2:
        raise ValueError("frequencies_mhz must be a one-dimensional grid with at least two channels")
    if not np.all(np.isfinite(frequencies)):
        raise ValueError("frequencies_mhz contains non-finite values")
    spacing = np.abs(np.diff(frequencies))
    channel_width = float(np.median(spacing))
    if not np.allclose(spacing, channel_width, rtol=1.0e-3, atol=1.0e-12):
        raise ValueError("sub-band planning requires a uniform frequency grid")
    if upchannel_factor is not None:
        validate_chime_factor_grid(upchannel_factor, channel_width)

    n_channels = frequencies.size
    n_subbands = int(n_subbands)
    if not 1 <= n_subbands <= n_channels:
        raise ValueError("n_subbands must be between one and the channel count")

    if fixed_slices is not None:
        if mode != "fixed":
            raise ValueError("fixed_slices requires mode='fixed'")
        slices = _validate_fixed_slices(fixed_slices, n_channels)
        weights = np.ones(n_channels, dtype=float)
        weight_definition = "fixed owner-reviewed channel intervals"
    elif mode == "equal_channels":
        boundaries = np.linspace(0, n_channels, n_subbands + 1, dtype=int)
        slices = list(zip(boundaries[:-1], boundaries[1:], strict=True))
        weights = np.ones(n_channels, dtype=float)
        weight_definition = "channel count"
    elif mode == "equal_snr":
        if signal is None or noise_rms is None:
            raise ValueError("equal_snr mode requires per-channel signal and noise_rms")
        signal_ma = np.ma.masked_invalid(np.ma.asarray(signal, dtype=float))
        noise_ma = np.ma.masked_invalid(np.ma.asarray(noise_rms, dtype=float))
        if signal_ma.shape != frequencies.shape or noise_ma.shape != frequencies.shape:
            raise ValueError("signal and noise_rms must match the frequency grid")
        invalid = (
            np.ma.getmaskarray(signal_ma)
            | np.ma.getmaskarray(noise_ma)
            | (np.ma.filled(noise_ma, 0.0) <= 0)
        )
        weights = np.square(
            np.clip(np.ma.filled(signal_ma, 0.0), 0.0, None)
            / np.where(invalid, 1.0, np.ma.filled(noise_ma, 1.0))
        )
        weights[invalid] = 0.0
        slices = _weighted_partition(weights, n_subbands)
        weight_definition = "sum of independent-channel S/N squared"
    else:
        raise ValueError(f"unsupported sub-band mode {mode!r}")

    subbands = []
    for start, stop in slices:
        local_freq = frequencies[start:stop]
        subbands.append(
            Subband(
                start=int(start),
                stop=int(stop),
                center_frequency_mhz=float(np.mean(local_freq)),
                bandwidth_mhz=float((stop - start) * channel_width),
                allocation_weight=float(np.sum(weights[start:stop])),
            )
        )
    return SubbandPlan(
        mode=mode,
        upchannel_factor=int(upchannel_factor) if upchannel_factor is not None else None,
        channel_width_mhz=channel_width,
        subbands=tuple(subbands),
        weight_definition=weight_definition,
    )


def lorentzian_component(x, gamma_mhz: float, modulation_index: float):
    """One Lorentzian ACF component with zero-lag amplitude ``m²``."""

    lag = np.asarray(x, dtype=float)
    return modulation_index**2 / (1.0 + (lag / gamma_mhz) ** 2)


def lorentzian_from_amplitude(lag_mhz, amplitude: float, gamma_mhz: float):
    """One Lorentzian parameterized by its zero-lag amplitude."""

    lag = np.asarray(lag_mhz, dtype=float)
    return amplitude / (1.0 + (lag / gamma_mhz) ** 2)


def summed_lorentzian_acf(
    lag_mhz,
    gamma_mhz: Sequence[float],
    modulation_indices: Sequence[float],
    *,
    baseline: float = 0.0,
):
    """Phenomenological sum of Lorentzian scales plus one constant."""

    lag = np.asarray(lag_mhz, dtype=float)
    gammas = np.asarray(gamma_mhz, dtype=float)
    mods = np.asarray(modulation_indices, dtype=float)
    if gammas.shape != mods.shape or gammas.ndim != 1 or gammas.size == 0:
        raise ValueError("gamma_mhz and modulation_indices must be non-empty matched vectors")
    out = np.full(lag.shape, float(baseline), dtype=float)
    for gamma, modulation in zip(gammas, mods, strict=True):
        out += lorentzian_component(lag, gamma, modulation)
    return out


def multiplicative_lorentzian_acf(
    lag_mhz,
    gamma_mhz: Sequence[float],
    modulation_indices: Sequence[float],
    *,
    baseline: float = 0.0,
):
    """Independent multiplicative-screen ACF, retaining every cross term.

    ``ACF + 1 = product_i(1 + m_i² rho_i)``.  For two screens this expands to
    the two Lorentzians plus the required product term.
    """

    lag = np.asarray(lag_mhz, dtype=float)
    gammas = np.asarray(gamma_mhz, dtype=float)
    mods = np.asarray(modulation_indices, dtype=float)
    if gammas.shape != mods.shape or gammas.ndim != 1 or gammas.size == 0:
        raise ValueError("gamma_mhz and modulation_indices must be non-empty matched vectors")
    product = np.ones(lag.shape, dtype=float)
    for gamma, modulation in zip(gammas, mods, strict=True):
        product *= 1.0 + lorentzian_component(lag, gamma, modulation)
    return product - 1.0 + float(baseline)


def total_modulation_index(
    modulation_indices: Sequence[float],
    *,
    parameterization: LorentzianParameterization,
) -> float:
    """Return modulation implied by the model's zero-lag excess."""

    mods = np.asarray(modulation_indices, dtype=float)
    if mods.ndim != 1 or mods.size == 0 or np.any(~np.isfinite(mods)) or np.any(mods < 0):
        raise ValueError("modulation_indices must be a non-empty finite non-negative vector")
    if parameterization == "phenomenological_sum":
        variance_fraction = float(np.sum(mods**2))
    elif parameterization == "multiplicative_screens":
        variance_fraction = float(np.prod(1.0 + mods**2) - 1.0)
    else:
        raise ValueError(f"unsupported Lorentzian parameterization {parameterization!r}")
    return float(np.sqrt(max(variance_fraction, 0.0)))


def total_modulation_index_uncertainty(
    modulation_indices: Sequence[float],
    covariance: np.ndarray,
    *,
    parameterization: LorentzianParameterization,
) -> float:
    """Propagate the full component covariance into the total modulation."""

    mods = np.asarray(modulation_indices, dtype=float)
    cov = np.asarray(covariance, dtype=float)
    if cov.shape != (mods.size, mods.size) or np.any(~np.isfinite(cov)):
        return float("nan")
    total = total_modulation_index(mods, parameterization=parameterization)
    if total == 0:
        return 0.0 if np.allclose(cov, 0.0) else float("nan")
    if parameterization == "phenomenological_sum":
        gradient = mods / total
    elif parameterization == "multiplicative_screens":
        factors = 1.0 + mods**2
        gradient = np.array(
            [
                modulation * np.prod(np.delete(factors, index)) / total
                for index, modulation in enumerate(mods)
            ]
        )
    else:
        raise ValueError(f"unsupported Lorentzian parameterization {parameterization!r}")
    variance = float(gradient @ cov @ gradient)
    return float(np.sqrt(max(variance, 0.0))) if np.isfinite(variance) else float("nan")


def noise_corrected_modulation_index(
    mean_signal: float,
    observed_variance: float,
    noise_variance: float,
    *,
    mean_signal_err: float | None = None,
    observed_variance_err: float | None = None,
    noise_variance_err: float | None = None,
) -> dict:
    """Noise-debiased direct modulation estimate.

    Returns no point estimate when the intrinsic variance
    ``observed_variance - noise_variance`` is non-positive.
    """

    mean_signal = float(mean_signal)
    observed_variance = float(observed_variance)
    noise_variance = float(noise_variance)
    if not (
        np.isfinite(mean_signal)
        and mean_signal > 0
        and np.isfinite(observed_variance)
        and observed_variance >= 0
        and np.isfinite(noise_variance)
        and noise_variance >= 0
    ):
        raise ValueError("mean and variances must be finite; mean positive; variances non-negative")
    intrinsic_variance = observed_variance - noise_variance
    modulation = (
        float(np.sqrt(intrinsic_variance) / mean_signal)
        if intrinsic_variance > 0
        else None
    )
    error_terms = (mean_signal_err, observed_variance_err, noise_variance_err)
    if any(item is not None for item in error_terms):
        if not all(item is not None and np.isfinite(item) and item >= 0 for item in error_terms):
            raise ValueError("all three non-negative uncertainty inputs are required together")
        if intrinsic_variance > 0:
            variance_err = np.hypot(observed_variance_err, noise_variance_err)
            dm_d_variance = 1.0 / (2.0 * mean_signal * np.sqrt(intrinsic_variance))
            dm_d_mean = -np.sqrt(intrinsic_variance) / mean_signal**2
            modulation_err = float(
                np.hypot(dm_d_variance * variance_err, dm_d_mean * mean_signal_err)
            )
        else:
            modulation_err = None
    else:
        modulation_err = None
    return {
        "modulation_index": modulation,
        "modulation_index_err": modulation_err,
        "intrinsic_variance": float(intrinsic_variance),
        "nonnegative_intrinsic_variance": float(max(intrinsic_variance, 0.0)),
        "detected": bool(intrinsic_variance > 0),
        "definition": "sqrt(max(observed variance - noise variance, 0)) / mean signal",
    }


def _n_lorentzian_model(n_components: int):
    model = ConstantModel(prefix="c_")
    for index in range(n_components):
        model += Model(lorentzian_component, prefix=f"l{index}_")
    return model


def _stderr(parameter) -> float:
    return float(parameter.stderr) if parameter.stderr is not None else float("nan")


def _component_record(result, index: int, gamma_lower: float, gamma_upper: float) -> dict:
    gamma_parameter = result.params[f"l{index}_gamma_mhz"]
    modulation_parameter = result.params[f"l{index}_modulation_index"]
    gamma = abs(float(gamma_parameter.value))
    gamma_err = _stderr(gamma_parameter)
    if gamma <= 1.01 * gamma_lower:
        status = "unresolved_by_channelization"
    elif gamma >= 0.99 * gamma_upper:
        status = "lower_limit_from_fit_span"
    elif not np.isfinite(gamma_err):
        status = "uncertainty_unavailable"
    else:
        status = "measured"
    return {
        "dnu_mhz": gamma,
        "m": abs(float(modulation_parameter.value)),
        "dnu_err": gamma_err,
        "m_err": _stderr(modulation_parameter),
        "width_status": status,
        "measurement_admissible": status == "measured",
    }


def fit_lorentzian_components(
    lags_mhz,
    acf,
    *,
    max_components: int = 3,
    acf_err=None,
    delta_bic_strong: float = 6.0,
    p_threshold: float = 0.05,
    channel_width_mhz: float | None = None,
    upchannel_factor: int | None = None,
) -> dict:
    """Fit and select 1..N phenomenological Lorentzian components.

    Only positive lags are fitted because a mirrored ACF does not provide
    independent samples.  An added component requires both strong Bayesian
    information criterion improvement and the legacy nested F-test check.
    The latter is retained for compatibility and labeled approximate because
    the null amplitude lies on a parameter boundary.
    """

    lags = np.asarray(lags_mhz, dtype=float)
    values = np.asarray(acf, dtype=float)
    errors = np.asarray(acf_err, dtype=float) if acf_err is not None else None
    if lags.shape != values.shape or lags.ndim != 1:
        raise ValueError("lags_mhz and acf must be matched one-dimensional arrays")
    if errors is not None and errors.shape != values.shape:
        raise ValueError("acf_err must match the ACF")
    max_components = int(max_components)
    if max_components < 1:
        raise ValueError("max_components must be positive")
    if upchannel_factor is not None:
        if channel_width_mhz is None:
            raise ValueError("factor-tagged fits require channel_width_mhz")
        validate_chime_factor_grid(upchannel_factor, channel_width_mhz)

    keep = np.isfinite(lags) & np.isfinite(values) & (lags > 0)
    if errors is not None:
        keep &= np.isfinite(errors) & (errors > 0)
    lags, values = lags[keep], values[keep]
    if errors is not None:
        errors = errors[keep]
    if lags.size < max(5, 2 * max_components + 2):
        raise ValueError("too few independent positive-lag ACF samples")
    order = np.argsort(lags)
    lags, values = lags[order], values[order]
    if errors is not None:
        errors = errors[order]

    span = float(np.max(lags))
    fine = float(np.min(np.diff(np.unique(lags))))
    peak = float(max(np.max(values) - np.median(values[-max(3, values.size // 5) :]), 1.0e-3))
    weights = None if errors is None else 1.0 / errors
    fits = []
    for count in range(1, max_components + 1):
        model = _n_lorentzian_model(count)
        params = model.make_params()
        params["c_c"].set(value=float(np.median(values[-max(3, values.size // 5) :])))
        seeds = (
            np.geomspace(max(2.0 * fine, span / 100.0), 0.5 * span, count)
            if count > 1
            else np.array([max(2.0 * fine, 0.2 * span)])
        )
        for index, gamma in enumerate(seeds):
            params[f"l{index}_gamma_mhz"].set(value=float(gamma), min=0.25 * fine, max=span)
            params[f"l{index}_modulation_index"].set(
                value=float(np.sqrt(peak / count)), min=0.0
            )
        try:
            result = model.fit(values, params, x=lags, weights=weights)
        except Exception:
            fits.append({"n": count, "success": False, "bic": np.inf, "chi2": np.inf})
            continue
        gamma_lower = 0.25 * fine
        gamma_upper = span
        components = sorted(
            (
                _component_record(result, index, gamma_lower, gamma_upper)
                for index in range(count)
            ),
            key=lambda item: item["dnu_mhz"],
            reverse=True,
        )
        component_mods = [item["m"] for item in components]
        modulation_names = [
            f"l{index}_modulation_index" for index in range(count)
        ]
        if result.covar is not None and all(
            name in result.var_names for name in modulation_names
        ):
            indices = [result.var_names.index(name) for name in modulation_names]
            modulation_covariance = np.asarray(result.covar)[np.ix_(indices, indices)]
        else:
            modulation_covariance = np.full((count, count), np.nan)
        fits.append(
            {
                "n": count,
                "success": bool(result.success),
                "bic": float(result.bic),
                "aic": float(result.aic),
                "chi2": float(result.chisqr),
                "redchi": float(result.redchi),
                "n_params": int(result.nvarys),
                "ndata": int(result.ndata),
                "constant": float(result.params["c_c"].value),
                "constant_err": _stderr(result.params["c_c"]),
                "components": components,
                "modulation_parameterization": "phenomenological_sum",
                "zero_lag_excess": float(np.sum(np.square(component_mods))),
                "m_total_equivalent": total_modulation_index(
                    component_mods, parameterization="phenomenological_sum"
                ),
                "m_total_equivalent_err": total_modulation_index_uncertainty(
                    component_mods,
                    modulation_covariance,
                    parameterization="phenomenological_sum",
                ),
                "m_total_is_physical_screen_model": count == 1,
            }
        )

    delta_bic = {}
    f_test = {}
    preferred = 1
    for count in range(2, max_components + 1):
        previous, current = fits[count - 2], fits[count - 1]
        if not (previous.get("success") and current.get("success")):
            break
        improvement = previous["bic"] - current["bic"]
        delta_bic[count] = float(improvement)
        dof_num = current["n_params"] - previous["n_params"]
        dof_den = current["ndata"] - current["n_params"]
        if (
            dof_num > 0
            and dof_den > 0
            and current["chi2"] > 0
            and previous["chi2"] > current["chi2"]
        ):
            statistic = (
                (previous["chi2"] - current["chi2"]) / dof_num
            ) / (current["chi2"] / dof_den)
            probability = float(f_dist.sf(statistic, dof_num, dof_den))
        else:
            probability = 1.0
        f_test[count] = probability
        if (
            preferred == count - 1
            and improvement > delta_bic_strong
            and probability < p_threshold
        ):
            preferred = count
        else:
            break

    return {
        "n_preferred": preferred,
        "fits": fits,
        "delta_bic": delta_bic,
        "f_test": f_test,
        "criterion": (
            f"BIC delta>{delta_bic_strong:g} and approximate nested F-test "
            f"p<{p_threshold:g}; both required"
        ),
        "selection_scope": (
            "diagnostic_only: per-lag errors do not encode correlated ACF-lag covariance; "
            "use acf_evidence for science-facing model comparison"
        ),
        "fit_domain": "positive lags only; zero lag excluded",
        "lag_unit": "MHz",
        "upchannel_factor": int(upchannel_factor) if upchannel_factor is not None else None,
        "channel_width_mhz": (
            float(channel_width_mhz) if channel_width_mhz is not None else None
        ),
    }
